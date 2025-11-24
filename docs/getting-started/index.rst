Getting Started
===============

Welcome to ChatSpatial! This section will help you get up and running quickly.

.. grid:: 1 2 2 3
    :gutter: 3

    .. grid-item-card:: 📦 Installation
        :link: installation
        :link-type: doc

        Install ChatSpatial on your system

    .. grid-item-card:: 🚀 Quick Start
        :link: quick-start
        :link-type: doc

        Get started with your first analysis

    .. grid-item-card:: ⚙️ Configuration
        :link: configuration
        :link-type: doc

        Configure ChatSpatial for your needs

Installation Options
--------------------

ChatSpatial can be installed in several ways:

- **Full Installation** (Recommended): Includes all analysis methods and features
- **Core Installation**: Basic functionality only
- **Development Installation**: For contributors

Prerequisites
-------------

Before installing ChatSpatial, ensure you have:

- Python 3.10 or higher
- Git (for cloning the repository)
- Virtual environment tool (venv, conda, etc.)

Quick Installation
------------------

The fastest way to get started:

.. code-block:: bash

    git clone https://github.com/cafferychen777/ChatSpatial.git
    cd ChatSpatial
    python3 -m venv chatspatial_env
    source chatspatial_env/bin/activate
    pip install -e ".[full]"

Next Steps
----------

After installation:

1. Follow the :doc:`quick-start` guide for your first analysis
2. Configure ChatSpatial with :doc:`configuration`
3. Explore the :doc:`../tutorials/index` for detailed workflows

.. toctree::
   :maxdepth: 2
   :hidden:

   installation
   quick-start
   configuration
