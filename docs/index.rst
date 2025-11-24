ChatSpatial Documentation
==========================

**Agentic Workflow Orchestration for Spatial Transcriptomics Analysis**

.. grid:: 2
    :gutter: 3

    .. grid-item-card:: 🚀 Get Started
        :link: getting-started/installation
        :link-type: doc

        Install ChatSpatial and start analyzing spatial data with natural language commands.

    .. grid-item-card:: 📚 View on GitHub
        :link: https://github.com/cafferychen777/ChatSpatial

        Check out the source code and contribute to the project.

----

Key Features
------------

- **Natural Language Interface**: Analyze spatial data using conversational commands with Claude and other LLMs
- **Comprehensive Analysis Tools**: Cell type annotation, spatial domains, trajectory analysis, and advanced spatial statistics
- **Multiple Data Formats**: Support for H5AD, 10X Visium, Slide-seq, MERFISH, and other spatial transcriptomics formats

What is ChatSpatial?
--------------------

ChatSpatial is an agentic workflow orchestration platform that eliminates the "implementation tax" in spatial transcriptomics research. Built on the Model Context Protocol (MCP), it integrates 60 state-of-the-art methods from fragmented Python and R ecosystems into a unified conversational interface. Unlike autonomous agents, ChatSpatial prioritizes human-steered discovery through schema-enforced tool-calling that ensures reproducibility and analytical fidelity.

Advanced Capabilities
~~~~~~~~~~~~~~~~~~~~~

- **Advanced Algorithms**: Integration with state-of-the-art spatial analysis methods
- **Visualization**: Built-in plotting and visualization capabilities

Supported Analysis Types
~~~~~~~~~~~~~~~~~~~~~~~~~

- **Cell Type Annotation**: Marker-based, reference-based (Tangram, scANVI), and ML approaches
- **Spatial Domain Identification**: SpaGCN, STAGATE, Leiden/Louvain clustering
- **Cell Communication**: LIANA, CellPhoneDB, and spatial interaction analysis
- **Trajectory Analysis**: RNA velocity, pseudotime, and developmental trajectories
- **Spatial Statistics**: Moran's I, spatial autocorrelation, and neighborhood analysis
- **Deconvolution**: Cell2location, RCTD, Stereoscope for cell type proportions

Getting Started
---------------

Prerequisites
~~~~~~~~~~~~~

- Python 3.10 or higher (required for MCP)
- Claude Desktop or compatible MCP client
- Git for installation

Quick Installation
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Clone the repository
    git clone https://github.com/cafferychen777/ChatSpatial.git
    cd ChatSpatial

    # Create virtual environment (strongly recommended)
    python3 -m venv chatspatial_env
    source chatspatial_env/bin/activate  # macOS/Linux

    # Install dependencies (recommended: full installation)
    pip install -e ".[full]"

    # Run the MCP server
    python -m chatspatial

For detailed installation instructions, see the :doc:`Installation Guide <getting-started/installation>`.

Your First Analysis
~~~~~~~~~~~~~~~~~~~

Once installed, you can start analyzing spatial data immediately:

.. code-block:: text

    Load the mouse brain Visium dataset and show me a spatial plot of the top variable genes

.. code-block:: text

    Identify spatial domains in my data using spagcn method and visualize the results

.. code-block:: text

    Perform cell type annotation using marker genes and show the spatial distribution

Documentation Structure
-----------------------

.. grid:: 2
    :gutter: 3

    .. grid-item-card:: 🎯 Getting Started
        :link: getting-started/index
        :link-type: doc

        Installation, quick start, and basic configuration

    .. grid-item-card:: 📖 Tutorials
        :link: tutorials/index
        :link-type: doc

        Step-by-step guides for core analysis workflows

    .. grid-item-card:: 📘 Reference
        :link: reference/index
        :link-type: doc

        API documentation and technical reference

    .. grid-item-card:: 💡 Examples
        :link: examples/index
        :link-type: doc

        Real-world analysis examples and case studies

Community and Support
---------------------

- **GitHub Issues**: `Report bugs and request features <https://github.com/cafferychen777/ChatSpatial/issues>`_
- **Discussions**: `Community Q&A and discussions <https://github.com/cafferychen777/ChatSpatial/discussions>`_
- **Documentation**: This site and inline help

License
-------

ChatSpatial is distributed under the MIT License. See `LICENSE <https://github.com/cafferychen777/ChatSpatial/blob/main/LICENSE>`_ for more information.

----

.. toctree::
   :maxdepth: 2
   :caption: Contents
   :hidden:

   getting-started/index
   tutorials/index
   reference/index
   examples/index
   resources/README

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
