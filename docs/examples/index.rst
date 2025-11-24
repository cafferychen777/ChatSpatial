Examples and Case Studies
=========================

Real-world examples demonstrating ChatSpatial's capabilities with actual spatial transcriptomics datasets.

.. grid:: 1 2 2 3
    :gutter: 3

    .. grid-item-card:: 📊 Datasets
        :link: datasets/README
        :link-type: doc

        Example datasets and data sources

    .. grid-item-card:: 📝 Guides
        :link: guides/README
        :link-type: doc

        Practical guides for specific use cases

    .. grid-item-card:: 🔬 Workflows
        :link: workflows/README
        :link-type: doc

        Complete analysis pipelines

Featured Examples
-----------------

Mouse Brain Visium Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Complete end-to-end analysis of mouse brain spatial transcriptomics data:

- Data loading and quality control
- Spatial domain identification with SpaGCN
- Cell type annotation using reference data
- Visualization of results

:doc:`View Tutorial → <workflows/case_studies/mouse_brain_visium>`

Multi-Sample Integration
~~~~~~~~~~~~~~~~~~~~~~~~~

Learn how to integrate and compare multiple spatial samples:

- Batch effect correction with Harmony
- Cross-sample comparative analysis
- Integrated visualization

Cell Communication Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Identify and visualize cell-cell interactions:

- Ligand-receptor pair analysis with LIANA
- Spatial communication patterns
- Network visualization

Using These Examples
--------------------

Each example includes:

1. **Dataset Description** - What data is being analyzed
2. **Scientific Question** - Research goal of the analysis
3. **Complete Workflow** - Step-by-step analysis commands
4. **Result Interpretation** - How to understand the output
5. **Visualizations** - Publication-ready figures

Example Data
------------

Example datasets can be accessed via:

- Built-in ChatSpatial datasets
- Public repositories (10X Genomics, SpatialDB)
- Your own data (following the same patterns)

See :doc:`datasets/README` for details.

Contributing Examples
---------------------

Have a great analysis workflow to share? We'd love to include it!

See our `Contributing Guide <https://github.com/cafferychen777/ChatSpatial/blob/main/CONTRIBUTING.md>`_ for details.

Additional Resources
--------------------

- :doc:`../tutorials/index` - Step-by-step learning guides
- :doc:`../reference/api/index` - Detailed tool documentation
- :doc:`../reference/quick-reference/README` - Command cheat sheets

.. toctree::
   :maxdepth: 2
   :hidden:

   datasets/README
   guides/README
   workflows/README
