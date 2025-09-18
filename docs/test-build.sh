#!/bin/bash

# Simple script to test Jekyll build locally
# This script checks if the Jekyll site builds successfully

echo "Testing Jekyll build for ChatSpatial documentation..."

# Check if we're in the right directory
if [ ! -f "_config.yml" ]; then
    echo "Error: _config.yml not found. Are you in the docs directory?"
    exit 1
fi

# Check if Gemfile exists
if [ ! -f "Gemfile" ]; then
    echo "Error: Gemfile not found."
    exit 1
fi

# Try to build the site (dry run)
echo "Attempting to build Jekyll site..."

# For systems with newer Ruby, try bundle install first
if command -v bundle &> /dev/null; then
    echo "Bundle command found, trying to install dependencies..."
    bundle install --quiet 2>/dev/null || echo "Bundle install failed, but continuing..."
fi

# Try Jekyll build with error checking
if command -v jekyll &> /dev/null; then
    echo "Jekyll found, attempting build..."
    jekyll build --dry-run 2>&1 | head -20
else
    echo "Jekyll not found globally, trying with bundle exec..."
    if command -v bundle &> /dev/null; then
        bundle exec jekyll build --dry-run 2>&1 | head -20
    else
        echo "Neither jekyll nor bundle found. You may need to install Ruby/Jekyll."
        echo "For macOS: brew install ruby"
        echo "Then: gem install bundler jekyll"
    fi
fi

echo ""
echo "Build test completed. Check output above for any errors."
echo "If successful, you can run 'bundle exec jekyll serve' to preview locally."