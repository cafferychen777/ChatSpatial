#!/usr/bin/env python3
"""
Example HTTP client for ChatSpatial MCP server
Demonstrates how to use the HTTP transport
"""

import asyncio
import json
import aiohttp
from typing import Dict, Any, Optional


class ChatSpatialHTTPClient:
    """HTTP client for ChatSpatial MCP server"""
    
    def __init__(self, base_url: str = "http://localhost:8000"):
        self.base_url = base_url
        self.session_id: Optional[str] = None
        
    async def create_session(self) -> str:
        """Create a new session"""
        async with aiohttp.ClientSession() as session:
            async with session.post(f"{self.base_url}/sessions") as resp:
                data = await resp.json()
                self.session_id = data["session_id"]
                return self.session_id
    
    async def call_rpc(self, method: str, params: Optional[Dict[str, Any]] = None) -> Any:
        """Call an RPC method"""
        request = {
            "jsonrpc": "2.0",
            "method": method,
            "params": params or {},
            "id": "1"
        }
        
        headers = {}
        if self.session_id:
            headers["X-Session-Id"] = self.session_id
        
        async with aiohttp.ClientSession() as session:
            async with session.post(
                f"{self.base_url}/rpc",
                json=request,
                headers=headers
            ) as resp:
                response = await resp.json()
                
                if "error" in response:
                    raise Exception(f"RPC Error: {response['error']}")
                
                return response.get("result")
    
    async def call_tool(self, tool_name: str, **kwargs) -> Any:
        """Call a tool"""
        return await self.call_rpc("tools/call", {
            "name": tool_name,
            "arguments": kwargs
        })
    
    async def list_tools(self) -> list:
        """List available tools"""
        return await self.call_rpc("tools/list")
    
    async def list_resources(self) -> list:
        """List available resources"""
        return await self.call_rpc("resources/list")
    
    async def read_resource(self, uri: str) -> str:
        """Read a resource"""
        return await self.call_rpc("resources/read", {"uri": uri})
    
    async def list_prompts(self) -> list:
        """List available prompts"""
        return await self.call_rpc("prompts/list")
    
    async def get_prompt(self, name: str, arguments: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """Get a prompt"""
        return await self.call_rpc("prompts/get", {
            "name": name,
            "arguments": arguments
        })


async def demo_basic_usage():
    """Demonstrate basic HTTP client usage"""
    print("=== ChatSpatial HTTP Client Demo ===\n")
    
    # Create client
    client = ChatSpatialHTTPClient()
    
    # Create session
    print("1. Creating session...")
    session_id = await client.create_session()
    print(f"   Session ID: {session_id}\n")
    
    # List tools
    print("2. Listing available tools...")
    tools = await client.list_tools()
    print(f"   Found {len(tools)} tools")
    for tool in tools[:3]:
        print(f"   - {tool['name']}: {tool.get('annotations', {}).get('title', tool['name'])}")
    print("   ...\n")
    
    # List resources
    print("3. Listing resources...")
    resources = await client.list_resources()
    print(f"   Found {len(resources)} resources")
    for resource in resources[:3]:
        print(f"   - {resource['uri']}: {resource['name']}")
    if len(resources) > 3:
        print("   ...\n")
    
    # List prompts
    print("4. Listing prompts...")
    prompts = await client.list_prompts()
    print(f"   Found {len(prompts)} prompts")
    for prompt in prompts[:3]:
        print(f"   - {prompt['name']}: {prompt['description']}")
    print("   ...\n")


async def demo_spatial_analysis():
    """Demonstrate spatial analysis workflow"""
    print("\n=== Spatial Analysis Workflow Demo ===\n")
    
    client = ChatSpatialHTTPClient()
    session_id = await client.create_session()
    print(f"Session created: {session_id}\n")
    
    try:
        # Load test data
        print("1. Loading spatial data...")
        # This would normally load real data
        # result = await client.call_tool(
        #     "load_data",
        #     data_path="/path/to/data.h5ad",
        #     data_type="10x_visium"
        # )
        print("   (Skipping actual data load in demo)\n")
        
        # Get a prompt for spatial analysis
        print("2. Getting spatial analysis prompt...")
        prompt = await client.get_prompt("analyze-spatial-expression", {
            "genes": ["CD3D", "CD3E", "CD8A"],
            "method": "moran"
        })
        print(f"   Prompt: {prompt['description']}")
        print(f"   Suggested action: {prompt['messages'][0]['content']}\n")
        
        # List resources again (would show loaded data)
        print("3. Checking available resources...")
        resources = await client.list_resources()
        for resource in resources:
            if "dataset" in resource['uri']:
                print(f"   - {resource['uri']}: {resource['name']}")
        
    except Exception as e:
        print(f"   Error: {e}")


async def demo_sse_streaming():
    """Demonstrate SSE streaming connection"""
    print("\n=== SSE Streaming Demo ===\n")
    
    print("Connecting to SSE endpoint...")
    
    async with aiohttp.ClientSession() as session:
        async with session.get("http://localhost:8000/sse") as resp:
            print(f"Connected! Status: {resp.status}")
            print("Receiving events (press Ctrl+C to stop):\n")
            
            count = 0
            async for line in resp.content:
                if line.startswith(b'data: '):
                    data = json.loads(line[6:])
                    print(f"Event: {data}")
                    
                    count += 1
                    if count >= 3:  # Stop after 3 events for demo
                        print("\nDemo complete!")
                        break


async def main():
    """Run all demos"""
    # Basic usage
    await demo_basic_usage()
    
    # Spatial analysis workflow
    await demo_spatial_analysis()
    
    # SSE streaming (commented out as it's long-running)
    # await demo_sse_streaming()
    
    print("\n=== Demo Complete ===")


if __name__ == "__main__":
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        print("\nDemo interrupted by user")