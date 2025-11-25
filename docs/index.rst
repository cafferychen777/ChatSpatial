ChatSpatial Documentation
==========================

**Chat with your spatial transcriptomics data. No coding required.**

ChatSpatial is a Model Context Protocol (MCP) server that integrates 60+ spatial analysis methods into Claude. Analyze your data through natural language conversation.

----

Getting Started
------------------

Choose your path:

.. grid:: 2
    :gutter: 3

    .. grid-item-card:: New Users
        :link: quickstart
        :link-type: doc

        **5-minute setup** - Get ChatSpatial running and analyze your first dataset

        Perfect for: Researchers, biologists, anyone new to ChatSpatial

    .. grid-item-card:: Advanced Users
        :link: advanced/methods-reference
        :link-type: doc

        **Detailed documentation** → All tools, parameters, and configuration

        Perfect for: Power users, developers, troubleshooting

----

Quick Links
-----------

**Essential:**

- :doc:`quickstart` - Get running in 5 minutes
- :doc:`examples` - Real analysis workflows
- :doc:`advanced/methods-reference` - All 60+ tools

**Support:**

- :doc:`advanced/troubleshooting` - Problem solving
- :doc:`advanced/faq` - Common questions
- `GitHub Issues <https://github.com/cafferychen777/ChatSpatial/issues>`_ - Report bugs

----

What Can ChatSpatial Do?
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Analysis Type
     - Natural Language Example
   * - **Load Data**
     - "Load my Visium dataset"
   * - **Spatial Domains**
     - "Identify tissue regions"
   * - **Cell Types**
     - "Annotate cell types using reference"
   * - **Deconvolution**
     - "Deconvolve spots with Cell2location"
   * - **Communication**
     - "Analyze cell-cell interactions"
   * - **Trajectories**
     - "Find developmental paths"
   * - **Enrichment**
     - "Run pathway analysis"
   * - **Visualization**
     - "Create spatial heatmap"

**60+ methods across 12 analysis categories** - all through conversation!

----

Architecture Overview
---------------------

.. image:: _static/architecture.png
   :width: 100%
   :align: center
   :alt: ChatSpatial Architecture

**ChatSpatial solves the data format barrier problem** through schema-enforced MCP execution, providing a unified interface to 60+ spatial transcriptomics methods across Python and R ecosystems.

----

Example Conversation
--------------------

.. code-block:: text

    User: "Load /path/to/visium_data.h5ad and identify spatial domains"

    ChatSpatial: Loaded 3,456 spots, 18,078 genes
                Identified 7 spatial domains using SpaGCN
                Generated visualization

    User: "Find marker genes for domain 3 and show me what cell type it is"

    ChatSpatial: Found 23 significant markers (adj. p < 0.05)
                Top markers: GFAP, S100B, AQP4
                Domain 3 shows astrocyte signature
                Created expression heatmap

**That's ChatSpatial. Natural conversation to scientific results.**

----

Learning Paths
--------------

**Path 1: First-Time User (30 minutes)**

1. :doc:`quickstart` - Setup and first analysis
2. :doc:`Basic Analysis Examples <examples>` - Learn fundamentals
3. Try your own data!

**Path 2: Intermediate User (1 hour)**

1. :doc:`Cell Type Analysis <examples>` - Annotation methods
2. :doc:`Deconvolution <examples>` - Estimate compositions
3. :doc:`Complete Workflows <examples>` - Multi-step analysis

**Path 3: Advanced User (2+ hours)**

1. :doc:`Advanced Analysis <examples>` - Trajectories, CNV, etc.
2. :doc:`advanced/methods-reference` - Deep dive into all tools
3. Combine multiple methods for publication-quality analysis

----

Key Features
------------

- **Natural Language Interface**: Analyze spatial data using conversational commands with Claude
- **60+ Methods**: Cell type annotation, spatial domains, trajectories, communication, and more
- **Multiple Formats**: H5AD, 10x Visium, Slide-seq, MERFISH, seqFISH
- **No Coding Required**: Just chat naturally with Claude
- **Publication Ready**: High-quality visualizations and comprehensive analysis

----

Community & Support
-------------------

- **GitHub Issues**: `Report bugs and request features <https://github.com/cafferychen777/ChatSpatial/issues>`_
- **Discussions**: `Community Q&A <https://github.com/cafferychen777/ChatSpatial/discussions>`_
- **Documentation**: Browse the guides below

----

.. toctree::
   :maxdepth: 1
   :caption: Core Documentation
   :hidden:

   quickstart
   examples

.. toctree::
   :maxdepth: 1
   :caption: Advanced Documentation
   :hidden:

   advanced/methods-reference
   advanced/installation
   advanced/configuration
   advanced/troubleshooting
   advanced/faq

----

**Ready to start?** See :doc:`quickstart`

Made for the spatial transcriptomics community
