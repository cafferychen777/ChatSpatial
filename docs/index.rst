ChatSpatial
===========

**MCP server for spatial transcriptomics analysis via natural language.**

Analyze your spatial data by chatting with Claude. No coding required.

.. only:: en

   .. raw:: html

      <div style="margin: 2rem 0; padding: 1.5rem; background: linear-gradient(135deg, #f0f9ff 0%, #e0f2fe 100%); border-radius: 12px; border-left: 4px solid #1565C0;">
        <p style="margin: 0; font-size: 1.1rem; color: #1e293b;">
          <strong>60+ methods</strong> across 15 analysis categories.
          Supports 10x Visium, Xenium, Slide-seq, MERFISH, and more.
        </p>
      </div>

.. only:: zh

   .. raw:: html

      <div style="margin: 2rem 0; padding: 1.5rem; background: linear-gradient(135deg, #f0f9ff 0%, #e0f2fe 100%); border-radius: 12px; border-left: 4px solid #1565C0;">
        <p style="margin: 0; font-size: 1.1rem; color: #1e293b;">
          覆盖 <strong>15 类分析的 60+ 方法</strong>。
          支持 10x Visium、Xenium、Slide-seq、MERFISH 等。
        </p>
      </div>

----

Getting Started
---------------

.. grid:: 3
    :gutter: 3

    .. grid-item-card:: Installation
        :link: installation
        :link-type: doc
        :class-card: sd-rounded-3

        **Set up the environment**

        Install ChatSpatial and prepare your Python environment.

    .. grid-item-card:: Quick Start
        :link: quickstart
        :link-type: doc
        :class-card: sd-rounded-3

        **Your first analysis**

        Load data and start analyzing with natural language.

    .. grid-item-card:: Docker / GHCR
        :link: docker
        :link-type: doc
        :class-card: sd-rounded-3

        **Run without local dependency resolution**

        Pull the GHCR image, mount data, and configure MCP.

    .. grid-item-card:: Concepts
        :link: concepts
        :link-type: doc
        :class-card: sd-rounded-3

        **Understand the methods**

        When to use which method and why.

    .. grid-item-card:: Examples
        :link: examples
        :link-type: doc
        :class-card: sd-rounded-3

        **See what's possible**

        Natural language commands for every analysis type.

----

Reference
---------

.. grid:: 2
    :gutter: 2

    .. grid-item-card:: Methods Reference
        :link: advanced/methods-reference
        :link-type: doc

        All 20 MCP tools with parameters

    .. grid-item-card:: Configuration
        :link: advanced/configuration
        :link-type: doc

        MCP client configuration

    .. grid-item-card:: Data Management
        :link: advanced/data-management
        :link-type: doc

        External dataset registry and fetch workflow

    .. grid-item-card:: Contributing
        :link: contributing
        :link-type: doc

        Contributor entry point for docs and code changes

----

Support
-------

- :doc:`advanced/troubleshooting` — Common issues and solutions
- :doc:`advanced/faq` — Frequently asked questions
- `GitHub Issues <https://github.com/cafferychen777/ChatSpatial/issues>`_ — Report bugs

----

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Getting Started

   installation
   docker
   quickstart
   concepts
   examples
   contributing

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Reference

   advanced/methods-reference
   advanced/configuration
   advanced/data-management

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Support

   advanced/troubleshooting
   advanced/faq
