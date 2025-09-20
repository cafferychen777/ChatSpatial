"""
Color utilities for visualization.
"""

import matplotlib.cm as cm
import matplotlib.colors as mcolors
import pandas as pd


def _ensure_categorical_colors(adata, key):
    """
    Ensure categorical data has proper color mapping for visualization.

    This function checks if a categorical feature has associated colors in adata.uns,
    and generates them if missing. This is essential for proper visualization of
    categorical data like cell types, clusters, etc.

    Args:
        adata: AnnData object
        key: Column name in adata.obs that contains categorical data

    Returns:
        bool: True if colors were ensured (either existed or were created), False if not categorical
    """
    # Check if the key exists and is categorical
    if key not in adata.obs.columns:
        return False

    if not pd.api.types.is_categorical_dtype(adata.obs[key]):
        return False

    # Check if colors already exist
    colors_key = f"{key}_colors"
    if colors_key in adata.uns:
        # Colors already exist, ensure they match the number of categories
        categories = adata.obs[key].cat.categories
        if len(adata.uns[colors_key]) >= len(categories):
            return True

    # Generate colors for categories
    categories = adata.obs[key].cat.categories
    n_categories = len(categories)

    # Choose colormap based on number of categories
    if n_categories <= 10:
        # Use tab10 for small number of categories
        colors = [mcolors.to_hex(cm.tab10(i)) for i in range(n_categories)]
    elif n_categories <= 20:
        # Use tab20 for medium number of categories
        colors = [mcolors.to_hex(cm.tab20(i / 20)) for i in range(n_categories)]
    else:
        # Use a combination of colormaps for large number of categories
        # Cycle through tab20, tab20b, and tab20c
        colors = []
        cmaps = [cm.tab20, cm.tab20b, cm.tab20c]
        for i in range(n_categories):
            cmap_idx = i // 20
            color_idx = (i % 20) / 20
            if cmap_idx < len(cmaps):
                colors.append(mcolors.to_hex(cmaps[cmap_idx](color_idx)))
            else:
                # Fallback to hsv colormap for very large numbers
                colors.append(mcolors.to_hex(cm.hsv(i / n_categories)))

    # Store colors in uns
    adata.uns[colors_key] = colors[:n_categories]

    return True
