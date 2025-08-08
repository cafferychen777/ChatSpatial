# Security Policy

## Supported Versions

We provide security updates for the following versions of ChatSpatial:

| Version | Supported          |
| ------- | ------------------ |
| 0.3.x   | :white_check_mark: |
| 0.2.x   | :x:                |
| < 0.2   | :x:                |

## Reporting a Vulnerability

We take the security of ChatSpatial seriously. If you discover a security vulnerability, please follow these steps:

### Private Disclosure Process

1. **Do not** open a public GitHub issue for security vulnerabilities
2. Send an email to the maintainers with:
   - A detailed description of the vulnerability
   - Steps to reproduce the issue
   - Potential impact assessment
   - Any suggested fixes or mitigations

### What to Expect

- **Acknowledgment**: We will acknowledge receipt of your vulnerability report within 48 hours
- **Initial Assessment**: We will provide an initial assessment within 5 business days
- **Updates**: We will keep you informed of our progress throughout the process
- **Resolution**: We aim to resolve critical security issues within 30 days

### Scope

This security policy applies to:

- The main ChatSpatial codebase
- Official Docker images
- Dependencies and third-party integrations
- Model Context Protocol (MCP) server implementations

### Out of Scope

- Issues in third-party dependencies (please report directly to those projects)
- General software bugs that don't have security implications
- Issues requiring physical access to systems

## Security Best Practices

When using ChatSpatial:

1. **Environment Variables**: Store sensitive data (API keys, database credentials) in environment variables, not in code
2. **Network Security**: Use HTTPS for all communications when deploying in production
3. **Data Privacy**: Be mindful of sensitive biological data and comply with relevant regulations
4. **Dependencies**: Regularly update dependencies to get security patches
5. **Access Control**: Implement appropriate authentication and authorization for production deployments

## Security Features

ChatSpatial includes several security features:

- Input validation and sanitization
- Secure handling of file uploads and data processing
- MCP protocol compliance with security best practices
- Error handling that doesn't expose sensitive information

## Contributing to Security

If you're contributing code to ChatSpatial, please:

- Follow secure coding practices
- Run security linters and tests
- Be cautious with external dependencies
- Document any security-related changes

Thank you for helping to keep ChatSpatial secure!