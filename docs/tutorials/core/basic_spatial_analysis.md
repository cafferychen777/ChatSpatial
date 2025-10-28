---
layout: default
title: Basic Spatial Analysis
parent: Core Tutorials
grand_parent: Tutorials
nav_order: 1
description: "Your first spatial transcriptomics analysis with ChatSpatial"
---

# Your First Spatial Analysis with ChatSpatial

This tutorial shows you how to analyze spatial transcriptomics data using **natural language** in Claude Desktop. No coding required - just ask ChatSpatial what you want to do!

## What You'll Learn

By the end of this tutorial, you'll know how to:

1. **Ask ChatSpatial to load your data** using simple questions
2. **Request preprocessing** with natural language
3. **Discover spatial domains** by asking about tissue structure
4. **Identify cell types** through conversational requests
5. **Explore cell communication** with intuitive questions
6. **Get beautiful visualizations** automatically

## Before You Start

### What You Need
- ✅ **Claude Desktop** installed and running
- ✅ **ChatSpatial configured** (see [Getting Started](../../getting-started/))
- ✅ **Your spatial data** (we'll use a mouse brain example)
- ✅ **Basic curiosity** about your tissue!

### No Coding Experience? Perfect!
This tutorial is designed for biologists, not programmers. You'll interact with ChatSpatial exactly like chatting with a knowledgeable colleague who happens to be an expert in spatial transcriptomics analysis.

## Step 1: Loading Your Data

### 💬 Start the Conversation

Open Claude Desktop and start chatting with ChatSpatial. Here's how a typical conversation begins:

---

**👤 You:** "Hi! I have a mouse brain Visium dataset that I'd like to analyze. Can you help me load it?"

**🤖 Claude:** "I'd be happy to help you analyze your mouse brain Visium data! Let me load it for you. What's the file path to your dataset?"

**👤 You:** "The file is at `/Users/myname/data/mouse_brain_visium.h5ad`"

**🤖 Claude:** "Perfect! Let me load your mouse brain Visium dataset now."

*[ChatSpatial automatically loads the data]*

**🤖 Claude:** "Great! I've successfully loaded your dataset. Here's what I found:
- **2,698 spots** across the tissue
- **31,053 genes** measured
- **Spatial coordinates** are included
- **Platform**: 10x Visium
- **Total UMI counts**: 4,235,892

Your data looks good and ready for analysis! What would you like to explore first?"

---

### 🎯 What Just Happened?

ChatSpatial automatically:
- ✅ Detected the file format (H5AD)
- ✅ Loaded the expression matrix
- ✅ Found the spatial coordinates
- ✅ Calculated basic statistics
- ✅ Validated data quality

### 💡 Pro Tips for Data Loading

**Instead of technical commands, just ask naturally:**

- ❌ Don't say: "Execute load_data() with parameters..."
- ✅ Do say: "Can you load my Visium data from this file?"

> ⚠️ **CRITICAL**: Always use **absolute paths** (starting with `/`) when specifying file locations. ChatSpatial's MCP server cannot find files with relative paths like `./data.h5ad` or `data.h5ad`.

**Path Examples:**
- ✅ **Correct**: `/Users/myname/Downloads/card_spatial.h5ad`
- ✅ **Correct**: `/home/user/projects/spatial_data.h5ad`
- ❌ **Wrong**: `data.h5ad` or `./folder/data.h5ad`

**Common questions you can ask:**
- "What's in my dataset?"
- "How many cells/spots do I have?"
- "Is my data good quality?"
- "Can you show me a summary?"

## Step 2: Preparing Your Data for Analysis

### 💬 Continuing the Conversation

Now that your data is loaded, let's prepare it for analysis:

---

**👤 You:** "Great! Now I'd like to preprocess this data. Can you clean it up and prepare it for analysis?"

**🤖 Claude:** "Absolutely! I'll preprocess your mouse brain data to ensure high-quality analysis. This includes quality control, normalization, and finding the most informative genes. Let me do that now."

*[ChatSpatial runs comprehensive preprocessing]*

**🤖 Claude:** "Preprocessing complete! Here's what I did and found:

**Quality Control:**
- Filtered out low-quality spots (kept 2,651 of 2,698 spots)
- Removed rarely expressed genes (kept 18,078 of 31,053 genes)
- Identified 2,000 highly variable genes for analysis

**Data Processing:**
- Normalized gene expression across spots
- Applied log transformation
- Computed principal components (50 PCs)
- Built neighbor graphs for clustering

**Initial Clustering:**
- Found 8 preliminary clusters using Leiden algorithm
- Clusters show good spatial organization

Would you like to see a visualization of these clusters, or shall we move on to identifying spatial domains?"

---

### 🎯 What ChatSpatial Did Automatically

**Quality Control:**
- ✅ Removed spots with too few or too many genes
- ✅ Filtered out rarely expressed genes
- ✅ Calculated quality metrics (mitochondrial genes, etc.)

**Normalization:**
- ✅ Normalized total counts per spot
- ✅ Log-transformed the data
- ✅ Identified highly variable genes

**Dimensionality Reduction:**
- ✅ Performed PCA analysis
- ✅ Built neighbor graphs
- ✅ Created UMAP embedding

### 💡 Natural Ways to Ask for Preprocessing

**Instead of technical parameters, try:**
- "Can you clean up my data?"
- "Please prepare my data for analysis"
- "I want to normalize and filter my dataset"
- "Can you find the most variable genes?"

**Follow-up questions you might ask:**
- "How many spots were filtered out?"
- "Can you show me the quality metrics?"
- "What do the clusters look like?"
- "Can I see a UMAP plot?"

## Step 3: Discovering Spatial Domains

### 💬 Finding Tissue Architecture

Let's continue our conversation to discover the spatial organization of the brain tissue:

---

**👤 You:** "This looks good! Now I'm curious about the spatial organization of this brain tissue. Can you identify different spatial domains or regions?"

**🤖 Claude:** "Excellent question! The spatial organization is one of the most fascinating aspects of brain tissue. Let me identify distinct spatial domains using SpaGCN, which is particularly good at finding spatially coherent regions in brain tissue."

*[ChatSpatial analyzes spatial domains]*

**🤖 Claude:** "Fantastic! I've identified **7 distinct spatial domains** in your mouse brain tissue:

**Spatial Domain Analysis Results:**
- **Domain 0** (423 spots): Cortical layer regions - likely outer cortex
- **Domain 1** (381 spots): White matter regions - myelinated areas
- **Domain 2** (356 spots): Hippocampal formation - memory center
- **Domain 3** (298 spots): Deep cortical layers - inner cortex
- **Domain 4** (267 spots): Ventricular zones - stem cell regions
- **Domain 5** (234 spots): Corpus callosum - connecting fibers
- **Domain 6** (192 spots): Boundary regions - transition zones

**Quality Metrics:**
- Spatial coherence: 0.84 (excellent spatial organization)
- Silhouette score: 0.58 (good cluster separation)

Would you like me to create a beautiful spatial map showing these domains?"

**👤 You:** "Yes, please show me the spatial map!"

**🤖 Claude:** "Here's your spatial domain map! 🧠✨"

*[Beautiful spatial visualization appears showing the brain tissue with different colored regions]*

**🤖 Claude:** "This map shows the distinct anatomical regions of your mouse brain. Each color represents a different spatial domain. You can see how the domains follow the natural architecture of brain tissue - the layered cortex, the central white matter, and the hippocampal structures. Pretty amazing how the algorithm captured the biological organization!"

---

### 🎯 What You Just Discovered

**Biological Insights:**
- ✅ **7 anatomically meaningful regions** identified
- ✅ **Cortical layering** clearly visible
- ✅ **White matter** vs **gray matter** separation
- ✅ **Hippocampal structures** distinguished
- ✅ **Natural tissue boundaries** respected

### 💡 Natural Ways to Explore Spatial Domains

**Questions you can ask:**
- "What spatial regions are in my tissue?"
- "Can you find different anatomical areas?"
- "How is my tissue organized spatially?"
- "Show me the spatial domains"
- "Are there distinct regions in this sample?"

**Follow-up explorations:**
- "Which genes define each domain?"
- "Can you compare domains statistically?"
- "What's unique about domain 3?"
- "Show me domain boundaries more clearly"

See also: [Getting Started](../../getting-started/)

