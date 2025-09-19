// Wait for DOM to be loaded
document.addEventListener('DOMContentLoaded', function() {
  // Check if mermaid is available
  if (typeof mermaid !== 'undefined') {
    mermaid.initialize({
      startOnLoad: true,
      theme: 'default',
      themeVariables: {
        primaryColor: '#e1f5fe',
        primaryTextColor: '#000',
        primaryBorderColor: '#1976d2',
        lineColor: '#1976d2',
        secondaryColor: '#f3e5f5',
        tertiaryColor: '#fff3e0'
      }
    });
    
    // Initialize mermaid diagrams
    mermaid.init(undefined, document.querySelectorAll('.language-mermaid, .mermaid'));
    
    // Re-render mermaid diagrams when page content changes
    const observer = new MutationObserver(() => {
      mermaid.init(undefined, document.querySelectorAll('.language-mermaid, .mermaid'));
    });
    
    observer.observe(document.body, {
      childList: true,
      subtree: true
    });
  }
});
