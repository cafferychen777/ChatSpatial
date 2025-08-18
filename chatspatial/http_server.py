"""
HTTP Server wrapper for ChatSpatial MCP
Implements Streamable HTTP Transport with SSE (Server-Sent Events)
"""

import asyncio
import json
import logging
import uuid
from typing import Dict, Any, Optional
from pathlib import Path
import argparse

from fastapi import FastAPI, HTTPException, Request, Response
from fastapi.responses import StreamingResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import uvicorn

from .server import mcp, data_store, adapter

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Create FastAPI app
app = FastAPI(
    title="ChatSpatial MCP Server",
    description="Model Context Protocol server for spatial transcriptomics analysis",
    version="1.0.0"
)

# Configure CORS - restrict to localhost by default
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:*", "http://127.0.0.1:*"],  # Only localhost
    allow_credentials=True,
    allow_methods=["GET", "POST", "OPTIONS"],
    allow_headers=["*"],
    expose_headers=["X-Session-Id"]
)

# Session management
sessions: Dict[str, Dict[str, Any]] = {}


class MCPRequest(BaseModel):
    """MCP request format"""
    jsonrpc: str = "2.0"
    method: str
    params: Optional[Dict[str, Any]] = None
    id: Optional[str] = None


class MCPResponse(BaseModel):
    """MCP response format"""
    jsonrpc: str = "2.0"
    result: Optional[Any] = None
    error: Optional[Dict[str, Any]] = None
    id: Optional[str] = None


@app.get("/")
async def root():
    """Root endpoint with server information"""
    return {
        "name": "ChatSpatial MCP Server",
        "version": "1.0.0",
        "protocol": "MCP",
        "transports": ["stdio", "sse", "http"],
        "endpoints": {
            "rpc": "/rpc",
            "sse": "/sse",
            "health": "/health",
            "sessions": "/sessions"
        }
    }


@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return {
        "status": "healthy",
        "sessions": len(sessions),
        "datasets": len(data_store)
    }


@app.post("/sessions")
async def create_session():
    """Create a new session"""
    session_id = str(uuid.uuid4())
    sessions[session_id] = {
        "id": session_id,
        "created_at": asyncio.get_event_loop().time(),
        "last_activity": asyncio.get_event_loop().time(),
        "data_store": {}  # Session-specific data store
    }
    
    logger.info(f"Created new session: {session_id}")
    
    return JSONResponse(
        content={"session_id": session_id},
        headers={"X-Session-Id": session_id}
    )


@app.get("/sessions/{session_id}")
async def get_session(session_id: str):
    """Get session information"""
    if session_id not in sessions:
        raise HTTPException(status_code=404, detail="Session not found")
    
    session = sessions[session_id]
    return {
        "id": session["id"],
        "created_at": session["created_at"],
        "last_activity": session["last_activity"],
        "datasets": list(session.get("data_store", {}).keys())
    }


@app.delete("/sessions/{session_id}")
async def delete_session(session_id: str):
    """Delete a session"""
    if session_id not in sessions:
        raise HTTPException(status_code=404, detail="Session not found")
    
    del sessions[session_id]
    logger.info(f"Deleted session: {session_id}")
    
    return {"message": "Session deleted"}


@app.post("/rpc")
async def handle_rpc(request: MCPRequest, session_id: Optional[str] = None):
    """Handle RPC requests"""
    try:
        # Update session activity
        if session_id and session_id in sessions:
            sessions[session_id]["last_activity"] = asyncio.get_event_loop().time()
        
        # Route the request to appropriate handler
        method = request.method
        params = request.params or {}
        
        # Handle different MCP methods
        if method == "initialize":
            result = await handle_initialize(params)
        elif method == "tools/list":
            # Tools are registered directly with the mcp object
            result = await mcp.list_tools()
        elif method == "tools/call":
            result = await handle_tool_call(params, session_id)
        elif method == "resources/list":
            result = await adapter.handle_resource_list()
        elif method == "resources/read":
            result = await adapter.handle_resource_read(params.get("uri", ""))
        elif method == "prompts/list":
            result = await adapter.handle_prompt_list()
        elif method == "prompts/get":
            prompt_name = params.get("name", "")
            prompt = await adapter.prompt_manager.get_prompt(prompt_name)
            if not prompt:
                raise ValueError(f"Prompt '{prompt_name}' not found.")
            result = {
                "name": prompt.name,
                "description": prompt.description,
                "arguments": prompt.arguments
            }
        else:
            raise ValueError(f"Unknown method: {method}")
        
        return MCPResponse(
            result=result,
            id=request.id
        )
        
    except Exception as e:
        logger.error(f"RPC error: {str(e)}", exc_info=True)
        
        # Format error response
        error_response = {
            "code": -32603,
            "message": str(e),
            "data": {"method": request.method}
        }
        
        return MCPResponse(
            error=error_response,
            id=request.id
        )


@app.get("/sse")
async def handle_sse(session_id: Optional[str] = None):
    """Handle Server-Sent Events for streaming responses"""
    async def event_generator():
        """Generate SSE events"""
        try:
            # Send initial connection event
            yield f"data: {json.dumps({'type': 'connection', 'status': 'connected'})}\n\n"
            
            # Keep connection alive
            while True:
                # Send heartbeat every 30 seconds
                await asyncio.sleep(30)
                yield f"data: {json.dumps({'type': 'heartbeat', 'timestamp': asyncio.get_event_loop().time()})}\n\n"
                
        except asyncio.CancelledError:
            logger.info("SSE connection closed")
            raise
    
    return StreamingResponse(
        event_generator(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
            "X-Accel-Buffering": "no"  # Disable nginx buffering
        }
    )


async def handle_initialize(params: Dict[str, Any]) -> Dict[str, Any]:
    """Handle initialization request"""
    # Per MCP spec, protocolVersion should be date-based.
    # Capabilities should be objects, not booleans.
    client_info = params.get("clientInfo", {})
    logger.info(f"Initialize request from {client_info.get('name', 'unknown client')}")
    
    return {
        "protocolVersion": "2025-06-18",
        "serverInfo": {
            "name": "ChatSpatial",
            "version": "1.0.0"
        },
        "capabilities": {
            "tools": {"listChanged": False},
            "resources": {},
            "prompts": {},
            "logging": {}
        }
    }


async def handle_tool_call(params: Dict[str, Any], session_id: Optional[str]) -> Any:
    """Handle tool call with session-specific data store"""
    tool_name = params.get("name", "")
    tool_params = params.get("arguments", {})
    
    # Get the tool handler
    if tool_name not in mcp._tool_handlers:
        raise ValueError(f"Unknown tool: {tool_name}")
    
    handler = mcp._tool_handlers[tool_name]
    
    # Use session-specific data store if available
    if session_id and session_id in sessions:
        # Temporarily replace global data_store with session-specific one
        global data_store
        original_store = data_store
        data_store = sessions[session_id]["data_store"]
        
        try:
            # Call the tool
            result = await handler.call(**tool_params)
            return result
        finally:
            # Restore original data store
            data_store = original_store
    else:
        # Use global data store
        result = await handler.call(**tool_params)
        return result


# Security middleware
@app.middleware("http")
async def security_middleware(request: Request, call_next):
    """Security middleware for origin validation"""
    # Check Origin header for non-GET requests
    if request.method != "GET":
        origin = request.headers.get("origin", "")
        if origin and not any(origin.startswith(allowed) for allowed in ["http://localhost:", "http://127.0.0.1:"]):
            return JSONResponse(
                status_code=403,
                content={"error": "Forbidden: Invalid origin"}
            )
    
    # Add security headers
    response = await call_next(request)
    response.headers["X-Content-Type-Options"] = "nosniff"
    response.headers["X-Frame-Options"] = "DENY"
    response.headers["X-XSS-Protection"] = "1; mode=block"
    
    return response


# Rate limiting (simple in-memory implementation)
request_counts: Dict[str, list] = {}

@app.middleware("http")
async def rate_limit_middleware(request: Request, call_next):
    """Simple rate limiting middleware"""
    # Get client IP
    client_ip = request.client.host if request.client else "unknown"
    
    # Initialize request count for IP
    if client_ip not in request_counts:
        request_counts[client_ip] = []
    
    # Remove old requests (older than 1 minute)
    current_time = asyncio.get_event_loop().time()
    request_counts[client_ip] = [
        t for t in request_counts[client_ip] 
        if current_time - t < 60
    ]
    
    # Check rate limit (100 requests per minute)
    if len(request_counts[client_ip]) >= 100:
        return JSONResponse(
            status_code=429,
            content={"error": "Too many requests"}
        )
    
    # Record this request
    request_counts[client_ip].append(current_time)
    
    return await call_next(request)


def create_cli():
    """Create command-line interface"""
    parser = argparse.ArgumentParser(
        description="ChatSpatial MCP HTTP Server"
    )
    parser.add_argument(
        "--host",
        default="127.0.0.1",
        help="Host to bind to (default: 127.0.0.1 for security)"
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8000,
        help="Port to bind to (default: 8000)"
    )
    parser.add_argument(
        "--reload",
        action="store_true",
        help="Enable auto-reload for development"
    )
    parser.add_argument(
        "--allow-external",
        action="store_true",
        help="Allow external connections (not recommended)"
    )
    
    return parser


def main():
    """Run the HTTP server"""
    parser = create_cli()
    args = parser.parse_args()
    
    # Security warning for external connections
    if args.allow_external and args.host != "127.0.0.1":
        logger.warning(
            "WARNING: Allowing external connections. "
            "This is not recommended for security reasons."
        )
    
    # Configure and run server
    logger.info(f"Starting ChatSpatial HTTP server on {args.host}:{args.port}")
    
    uvicorn.run(
        "chatspatial.http_server:app",
        host=args.host,
        port=args.port,
        reload=args.reload,
        log_level="info"
    )


if __name__ == "__main__":
    main()