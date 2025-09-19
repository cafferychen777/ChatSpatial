# ChatSpatial Documentation Theme Customization

This document describes the custom styling and branding applied to the ChatSpatial documentation site built with Just the Docs theme.

## 🎨 Design System Overview

ChatSpatial documentation features a modern, professional design system inspired by spatial transcriptomics visualization and bioinformatics tools. The theme emphasizes:

- **Accessibility**: WCAG 2.1 AA compliance
- **Mobile-first**: Responsive design for all devices
- **Brand consistency**: Cohesive visual identity
- **User experience**: Enhanced readability and navigation

## 📁 File Structure

```
docs/
├── _sass/                      # Custom SCSS modules
│   ├── color_schemes/
│   │   └── chatspatial.scss   # Custom color scheme for Just the Docs
│   ├── custom.scss            # Main custom styles and CSS variables
│   ├── responsive.scss        # Mobile-first responsive design
│   └── components.scss        # Reusable UI components
├── _layouts/
│   └── default.html           # Custom layout with JS integration
├── assets/
│   ├── css/
│   │   └── style.scss         # Main stylesheet entry point
│   ├── js/
│   │   └── custom.js          # Interactive functionality
│   └── images/
│       ├── logo.svg           # Main ChatSpatial logo
│       └── logo-icon.svg      # Icon-only logo variant
├── _config.yml                # Updated Jekyll configuration
└── style-guide.md             # Design system documentation
```

## 🎯 Brand Identity

### Color Palette

Our color scheme reflects spatial biology and modern bioinformatics:

- **Primary Blue**: `#0ea5e9` → `#0369a1` (Spatial data visualization)
- **Secondary Purple**: `#d946ef` → `#a21caf` (Cell type highlighting)
- **Accent Green**: `#22c55e` → `#15803d` (Gene expression gradients)
- **Neutral Grays**: Optimized for readability and accessibility

### Logo Design

The ChatSpatial logo combines:
- Spatial grid dots representing tissue coordinates
- Connection lines showing spatial relationships  
- AI chat bubble indicating conversational interface
- Modern typography with brand colors

### Typography

- **Body Text**: Inter font family for excellent readability
- **Code**: JetBrains Mono for technical content
- **Responsive scaling**: Mobile-first font size adjustments

## 🧩 Component System

### Buttons
- **Primary**: Gradient backgrounds with hover animations
- **Outline**: Clean borders with smooth transitions
- **Sizes**: Small, regular, and large variants
- **States**: Hover, focus, and disabled styles

### Cards
- **Shadow system**: Layered shadows for depth
- **Hover effects**: Subtle lift animations
- **Grid layouts**: Responsive feature showcases

### Alerts & Callouts
- **Contextual styling**: Info, success, warning, error states
- **Just the Docs integration**: Enhanced callout styling
- **Icon support**: Visual indicators for different message types

### Code Blocks
- **Language labels**: Clear identification of code type
- **Copy buttons**: One-click code copying
- **Syntax highlighting**: Enhanced color scheme
- **Responsive**: Horizontal scrolling on mobile

## 📱 Responsive Design

### Breakpoints
- **Mobile**: < 500px (touch-friendly, single column)
- **Tablet**: 500px - 750px (adaptive grid, touch optimized)
- **Desktop**: 750px - 1064px (multi-column, hover states)
- **Large**: > 1064px (full layout, maximum width)

### Mobile Optimizations
- Touch targets ≥ 44px for accessibility
- Simplified navigation patterns
- Optimized typography scaling
- Reduced animation for better performance

## ✨ Interactive Features

### JavaScript Enhancements
- **Tab system**: Organized content switching
- **Accordion components**: Collapsible content sections
- **Copy functionality**: Code block copying
- **Tooltip system**: Contextual help
- **Smooth scrolling**: Enhanced anchor navigation
- **Theme toggle**: Light/dark mode switching (planned)

### Animation System
- **Fade-in effects**: Intersection Observer-based
- **Progress bars**: Animated loading indicators
- **Hover states**: Smooth transitions
- **Reduced motion**: Respects user preferences

## 🔧 Configuration

### Jekyll Configuration (_config.yml)
```yaml
# Logo and branding
logo: "/ChatSpatial/assets/images/logo.svg"
color_scheme: "chatspatial"

# Enhanced callouts
callouts_level: loud
callouts:
  highlight:
    color: yellow
  important:
    color: blue
  # ... additional callout configurations
```

### Custom CSS Variables
```scss
:root {
  // Brand colors
  --cs-primary-500: #0ea5e9;
  --cs-secondary-500: #d946ef;
  --cs-accent-500: #22c55e;
  
  // Typography
  --cs-font-sans: 'Inter', sans-serif;
  --cs-font-mono: 'JetBrains Mono', monospace;
  
  // Shadows and effects
  --cs-shadow-sm: 0 1px 2px 0 rgb(0 0 0 / 0.05);
  --cs-shadow-md: 0 4px 6px -1px rgb(0 0 0 / 0.1);
  // ... additional variables
}
```

## 🚀 Usage Examples

### Feature Grid
```html
<div class="features-grid">
  <div class="feature-card">
    <div class="feature-icon">
      <!-- SVG icon -->
    </div>
    <h3>Feature Title</h3>
    <p>Feature description</p>
  </div>
</div>
```

### Enhanced Buttons
```markdown
[Primary Button]{: .btn .btn-primary}
[Outline Button]{: .btn .btn-outline}
[Success Button]{: .btn .btn-success}
```

### Custom Alerts
```html
<div class="alert alert-info">
  <div class="alert-content">
    <h4>Information</h4>
    <p>Alert message content</p>
  </div>
</div>
```

## 🎨 Style Guide

Visit `/style-guide.html` to see all components in action, including:
- Complete color palette
- Typography specimens  
- Button variations
- Card layouts
- Alert styles
- Progress bars
- Form elements
- Responsive demonstrations

## 🛠️ Development

### Building the Site
```bash
# Install dependencies
bundle install

# Serve locally
bundle exec jekyll serve

# Build for production
bundle exec jekyll build
```

### Adding New Components
1. Define styles in `_sass/components.scss`
2. Add responsive rules in `_sass/responsive.scss`  
3. Update color variables in `_sass/color_schemes/chatspatial.scss`
4. Document in `style-guide.md`

### Customizing Colors
Edit `_sass/color_schemes/chatspatial.scss` to modify:
- Brand color values
- Just the Docs color mappings
- Callout styling
- Component color assignments

## 📈 Performance Considerations

- **CSS optimization**: Modular SCSS architecture
- **Image formats**: SVG logos for crisp scaling
- **Font loading**: System font fallbacks
- **Animation**: GPU-accelerated transforms
- **Mobile**: Touch-optimized interactions

## ♿ Accessibility Features

- **Color contrast**: WCAG AA compliant ratios
- **Focus indicators**: Visible keyboard navigation
- **Screen readers**: Semantic HTML structure
- **Touch targets**: Minimum 44px size
- **Motion**: Reduced motion preference support
- **Alt text**: Descriptive image alternatives

## 🔮 Future Enhancements

- [ ] Dark mode theme toggle
- [ ] Advanced search styling
- [ ] Interactive tutorials
- [ ] Animation library integration
- [ ] Performance monitoring
- [ ] A/B testing framework

---

**Note**: This custom theme builds upon Just the Docs while maintaining compatibility with Jekyll and GitHub Pages. All customizations follow best practices for maintainability and performance.