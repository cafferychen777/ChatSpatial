# GitHub Pages Configuration Guide

This document outlines the GitHub Pages configuration for ChatSpatial documentation.

## Repository Settings

To configure GitHub Pages for this repository:

1. Go to **Settings** → **Pages** in your GitHub repository
2. Set **Source** to "GitHub Actions"
3. The documentation will be automatically deployed when changes are pushed to the `main` branch

## Custom Domain (Optional)

If you want to use a custom domain:

1. Add a `CNAME` file to the `docs/` directory with your domain
2. Configure DNS settings for your domain
3. Update the `url` field in `docs/_config.yml`

## Branch Configuration

- **Main branch**: `main` - Production deployments
- **Feature branches**: Preview builds are generated for PRs affecting documentation

## Workflow Files

- `.github/workflows/pages.yml` - Main deployment workflow
- `.github/workflows/docs-preview.yml` - PR preview builds

## Troubleshooting

### Common Issues

1. **Build failures**: Check the Actions tab for detailed error logs
2. **Broken links**: Ensure all internal links use relative paths
3. **Missing assets**: Verify all images and CSS files are committed
4. **Theme issues**: Check that Just the Docs theme version is compatible

### Debug Steps

1. Run `docs/test-build.sh` locally to test the build
2. Check Jekyll configuration with `bundle exec jekyll doctor`
3. Review GitHub Actions logs for detailed error messages

## Local Development

```bash
cd docs
bundle install
bundle exec jekyll serve --livereload
```

Visit http://localhost:4000 to preview the site locally.

## Deployment Status

The documentation is available at: https://cafferychen777.github.io/ChatSpatial/

Check the deployment status in the **Environments** section of the repository.