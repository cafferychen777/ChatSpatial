document$.subscribe(() => {
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
  
  // Re-render mermaid diagrams when page content changes
  const observer = new MutationObserver(() => {
    mermaid.init(undefined, document.querySelectorAll('.mermaid'));
  });
  
  observer.observe(document.body, {
    childList: true,
    subtree: true
  });
});
