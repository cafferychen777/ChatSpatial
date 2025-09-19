---
layout: default
title: Trajectory Analysis
parent: Advanced Tutorials
grand_parent: Tutorials
nav_order: 4
description: "RNA velocity and trajectory analysis in spatial data"
---

# RNA Velocity and Trajectory Analysis with ChatSpatial

Learn how to uncover cellular dynamics, developmental trajectories, and temporal processes in your spatial transcriptomics data through natural conversation with ChatSpatial.

## What You Will Discover

By the end of this tutorial, you will know how to:

1. **Ask about cell dynamics** in your tissue using natural language
2. **Compute RNA velocity** to see directional changes in gene expression
3. **Infer developmental trajectories** and cell fate transitions
4. **Explore spatial-temporal patterns** across tissue regions
5. **Interpret biological significance** of velocity and trajectory results
6. **Choose the right method** for your specific research question

## Before You Start

### What You Need
- ✅ **Spatial data loaded** in ChatSpatial
- ✅ **Basic preprocessing completed** (see [Basic Analysis](../core/basic_spatial_analysis.md))
- ✅ **Interest in** cellular dynamics and development

### What Are RNA Velocity & Trajectories?

**RNA Velocity** reveals the **direction and speed** of cellular changes by measuring the balance between newly transcribed (unspliced) and mature (spliced) RNA molecules. Think of it as seeing which direction cells are "moving" in gene expression space.

**Trajectory Analysis** reconstructs **developmental paths** and cell state transitions, helping you understand how cells differentiate, mature, or respond to stimuli over time.

## Biological Background

### Understanding RNA Velocity

RNA velocity is based on a simple biological principle:

1. **Genes are transcribed** → unspliced pre-mRNA is created
2. **Pre-mRNA is processed** → spliced mature mRNA is produced  
3. **mRNA is degraded** → both forms eventually disappear

The **ratio of unspliced to spliced RNA** tells us:
- **High unspliced/spliced**: Gene is being actively turned ON
- **Low unspliced/spliced**: Gene is being turned OFF
- **Balanced ratio**: Gene expression is in steady state

### Why Spatial RNA Velocity Matters

Traditional RNA velocity analyzes cells in isolation, but **spatial velocity** reveals:
- **Where** developmental processes occur in tissue
- **How** cellular transitions spread spatially
- **Whether** neighboring cells follow similar trajectories
- **Which** tissue regions drive developmental processes

## Step 1: Starting with RNA Velocity Analysis

### 💬 The Conversation Begins

---

**👤 You:** "I have a developing brain tissue dataset and I am curious about cellular dynamics. Can you analyze RNA velocity to see how cells are changing over time?"

**🤖 ChatSpatial:** "Good question! RNA velocity analysis is good for understanding developmental dynamics in brain tissue. Let me check your data and compute RNA velocity for you.

First, let me examine what type of velocity data you have available..."

*[ChatSpatial analyzes data structure]*

**🤖 ChatSpatial:** "I found that your dataset contains both **spliced and unspliced** RNA counts, which is suitable for direct velocity computation. Here are your options:

**Available Velocity Methods:**
1. **scVelo (recommended)** - Fast and reliable, suitable for most analyses
2. **VELOVI** - Deep learning approach for velocity computation, more accurate but slower
3. **SIRV** - If you have paired scRNA-seq data

Which velocity method would you prefer, or should I start with scVelo?"

**👤 You:** "Let us start with scVelo. I want to understand the basics first."

**🤖 ChatSpatial:** "Good choice! scVelo is good for getting started. Let me compute RNA velocity using the stochastic mode, which works well for most spatial datasets."

*[ChatSpatial computes RNA velocity]*

**🤖 ChatSpatial:** "RNA velocity computation is complete for your brain tissue. Here is what I found:

**RNA Velocity Results:**
- **2,651 spots** analyzed
- **Velocity computed** for 18,078 genes  
- **Mode**: Stochastic (captures cell state transitions)
- **Average velocity magnitude**: 0.34 (moderate dynamics)
- **Top dynamic genes**: Neurod1, Dlx2, Gsx2, Pax6, Tbr1

**Key Insights:**
- Strong velocity patterns in **ventricular zones** (high neurogenesis)
- Directional flow from **stem cell regions** toward **cortical layers**
- Active differentiation in **intermediate zones**
- Stable expression in **mature cortical areas**

Would you like to see a velocity field visualization or explore specific regions?"

---

### 🎯 What Just Happened?

ChatSpatial automatically:
- ✅ Detected spliced/unspliced RNA data layers
- ✅ Chose appropriate scVelo parameters for spatial data
- ✅ Computed velocity vectors for all genes and spots
- ✅ Built transition probability matrices
- ✅ Identified regions with high cellular dynamics
- ✅ Highlighted biologically relevant patterns

## Step 2: Exploring Velocity Patterns

### 💬 Diving Deeper into Dynamics

---

**👤 You:** "This is interesting! Can you show me where the most dynamic changes are happening in my tissue?"

**🤖 ChatSpatial:** "Absolutely! Let me create visualizations showing cellular dynamics across your brain tissue."

*[Velocity visualization appears]*

**🤖 ChatSpatial:** "Here is your RNA velocity map! 🧠⚡

**What you are seeing:**
- **Arrows** show velocity direction (where cells are heading)
- **Arrow length** indicates velocity magnitude (speed of change)
- **Colors** represent velocity confidence (red = high, blue = low)
- **Background** shows your spatial tissue structure

**Key Biological Patterns:**
1. **Ventricular Zone**: Strong outward velocities (cells leaving stem state)
2. **Subventricular Zone**: Directional flow toward cortex
3. **Intermediate Zone**: Cells transitioning to mature neurons
4. **Cortical Plate**: Low velocity (mature, stable cells)

This matches known brain development patterns. Your tissue shows active neurogenesis with cells migrating from germinal zones to their final cortical destinations."

**👤 You:** "Can you tell me which specific genes are driving these velocity patterns?"

**🤖 ChatSpatial:** "Good question! Let me identify the top velocity genes and their biological significance:

**Top Velocity-Driving Genes:**

**Neurogenesis & Differentiation:**
- **Neurod1** (Neural differentiation factor) - High velocity in transition zones
- **Pax6** (Neural progenitor marker) - Decreasing in differentiating areas
- **Tbr1** (Cortical neuron marker) - Increasing toward cortical plate

**Regional Specification:**
- **Dlx2** (Interneuron specification) - Active in ganglionic eminence
- **Gsx2** (Progenitor maintenance) - High in ventricular zones
- **Emx2** (Cortical patterning) - Gradient across cortical areas

**Migration & Maturation:**
- **Dcx** (Neuronal migration) - High in migratory paths
- **Map2** (Neuronal maturation) - Increasing in cortical layers
- **Gfap** (Glial differentiation) - Active in astrocyte-rich regions

These genes tell the story of brain development: stem cells (Pax6) differentiating into neurons (Neurod1), migrating to cortex (Dcx), and maturing (Map2)!"

---

### 💡 Natural Ways to Explore Velocity

**Questions you can ask:**
- "Where are the most dynamic changes happening?"
- "Which genes have the highest velocity?"
- "Can you show me velocity patterns in specific regions?"
- "How confident is the velocity computation?"
- "What biological processes do these patterns represent?"

**Follow-up explorations:**
- "Compare velocity between different conditions"
- "Show me genes that are speeding up vs slowing down"
- "Which areas have the most directional flow?"

## Step 3: Trajectory Inference and Pseudotime

### 💬 Understanding Developmental Paths

---

**👤 You:** "This velocity analysis is informative! Now I want to understand the actual developmental trajectories. Can you trace how cells progress through different states?"

**🤖 ChatSpatial:** "Good next step! Trajectory analysis will show you the developmental paths that cells follow. Let me infer pseudotime and cellular trajectories using the velocity information we just computed.

For trajectory inference, I can use:
- **CellRank** (recommended with velocity data) - Uses velocity to guide trajectory inference
- **Palantir** - Good general-purpose trajectory method
- **DPT** - Simple diffusion-based approach

Let me start with CellRank since we have velocity data."

*[ChatSpatial performs trajectory analysis]*

**🤖 ChatSpatial:** "I have reconstructed the developmental trajectories in your brain tissue using CellRank, which combines RNA velocity with spatial information.

**Trajectory Analysis Results:**

**Developmental Timeline:**
- **Pseudotime range**: 0.0 (earliest) to 1.0 (most mature)
- **7 distinct cell states** identified along trajectories
- **3 major developmental branches** discovered

**Cell State Progression:**
1. **Neural stem cells** (pseudotime 0.0-0.2) → Ventricular zones
2. **Neural progenitors** (pseudotime 0.2-0.4) → Subventricular zones  
3. **Migrating neurons** (pseudotime 0.4-0.6) → Intermediate zones
4. **Cortical neurons** (pseudotime 0.6-0.8) → Cortical layers
5. **Mature neurons** (pseudotime 0.8-1.0) → Upper cortical layers

**Spatial-Temporal Insights:**
- Trajectories flow from **inner (ventricular) to outer (cortical)** regions
- **Multiple parallel paths** suggest different neuronal subtypes
- **Branch points** occur at key developmental decision points
- **Spatial organization** mirrors temporal progression

This captures the inside-out pattern of mammalian cortical development."

**👤 You:** "This is very informative! Can you show me these trajectories on the actual tissue?"

**🤖 ChatSpatial:** "Absolutely! Let me create a spatial trajectory map."

*[Spatial-temporal visualization appears]*

**🤖 ChatSpatial:** "Here is your developmental trajectory map! 🌈🧠

**What you are seeing:**
- **Color gradient** represents pseudotime (purple = early, yellow = late)
- **Trajectory paths** show developmental flows
- **Branch points** indicate cell fate decisions
- **Arrows** show predicted future cell states

**Biological Insights:**
- The **radial pattern** matches known neural development
- **Multiple streams** suggest parallel differentiation programs  
- **Branch points** occur at anatomically correct locations
- **Timing matches** known developmental stages

Your tissue captures a snapshot of brain development in action."

---

## Step 4: Detailed Velocity Analysis with VELOVI

### 💬 Deep Learning Approach to Velocity Computation

---

**👤 You:** "This traditional approach is good! I am curious about the newer deep learning methods. Can you try VELOVI for better velocity computation?"

**🤖 ChatSpatial:** "Good idea! VELOVI uses deep learning to compute RNA velocity more accurately by learning complex velocity patterns. It is particularly good for noisy data and can capture subtler velocity dynamics. 

**Important**: VELOVI is specifically for **velocity computation** - after computing velocities with VELOVI, we will still use trajectory methods like CellRank or Palantir to infer the actual developmental paths. Let me run VELOVI for velocity computation."

*[ChatSpatial runs VELOVI deep learning model]*

**🤖 ChatSpatial:** "VELOVI velocity computation is complete. Here is how the improved velocities compare to traditional scVelo:

**VELOVI Velocity Results:**

**Technical Performance:**
- **Model training**: Converged after 800 epochs
- **Latent dimensions**: 10 (captured main variation)  
- **Velocity uncertainty**: Computed (shows confidence)
- **GPU acceleration**: Used for faster training

**Improved Velocity Quality:**
- **More refined velocity vectors** - smoother, less noisy
- **Better uncertainty estimates** - shows which regions are reliable
- **Subtle pattern discovery** - found additional velocity transitions
- **Improved gene rankings** - more biologically relevant

**Key Velocity Improvements from scVelo:**
1. **Smoother velocity fields** - less jagged vectors
2. **Better handling of noise** - more robust to technical artifacts
3. **Uncertainty quantification** - tells you where to trust velocity estimates
4. **Novel transition detection** - discovered intermediate velocity states

**When to Use VELOVI for Velocity:**
- Complex datasets with multiple cell types
- Noisy or challenging data
- Need for velocity uncertainty estimates
- Publication-quality velocity analyses

Now let me use these improved VELOVI velocities with CellRank for trajectory inference!"

*[ChatSpatial runs CellRank with VELOVI velocities]*

**🤖 ChatSpatial:** "I have now computed trajectories using the improved VELOVI velocities. The trajectory inference shows cleaner developmental paths due to the more accurate velocity estimates from VELOVI."

**👤 You:** "That is helpful! What about the uncertainty - which regions should I trust most?"

**🤖 ChatSpatial:** "Good question! The uncertainty information is valuable for interpretation:

**High Confidence Regions (trust these results):**
- **Ventricular zones**: Clear stem cell dynamics
- **Cortical plate**: Well-defined mature states
- **Major transition zones**: Strong velocity signals

**Medium Confidence Regions (interpret carefully):**
- **Boundary areas**: Mixed cell populations
- **Low-cell-density spots**: Limited data
- **Tissue edges**: Potential artifacts

**Lower Confidence Regions (use caution):**
- **Sparse areas**: Few cells for reliable inference
- **Highly heterogeneous spots**: Mixed signals
- **Technical noise regions**: Processing artifacts

**Practical Recommendations:**
1. **Focus biological interpretation** on high-confidence regions
2. **Validate uncertain areas** with additional experiments  
3. **Use confidence maps** to guide follow-up studies
4. **Combine with spatial context** for better interpretation

This uncertainty information helps you prioritize which findings to pursue experimentally!"

---

## Step 5: Method Comparison and Selection

### 💬 Choosing the Right Approach

---

**👤 You:** "I have seen scVelo and VELOVI now. How do I know which method is best for my research question?"

**🤖 ChatSpatial:** "Good question! Let me break down when to use each method:

## Method Comparison Guide

### **Velocity Computation Methods**

#### **scVelo** 🚀
**Best for:**
- First-time RNA velocity analysis
- Quick exploratory studies  
- Well-behaved, clean datasets
- Standard developmental processes
- Teaching and demonstrations

**Advantages:**
- Fast computation (minutes)
- Well-established and validated
- Extensive documentation
- Works with most datasets
- Good for beginners

**Limitations:**
- Can be noisy with complex data
- Limited uncertainty quantification
- Basic handling of technical artifacts

#### **VELOVI** 🧠
**Best for:**
- Complex, multi-lineage systems
- Noisy or challenging datasets
- Publication-quality velocity computation
- Need uncertainty estimates
- Detailed research projects

**Advantages:**
- More accurate velocity estimates
- Built-in uncertainty quantification
- Better noise handling
- Can discover subtle patterns
- Deep learning powered

**Limitations:**
- Slower (requires training)
- More complex to interpret
- Higher computational requirements
- Newer (less established)

#### **SIRV** 🔬
**Best for:**
- Spatial data WITHOUT velocity layers
- Have paired scRNA-seq reference
- Cross-dataset integration
- Validation studies

**Advantages:**
- Works without spliced/unspliced data
- Leverages reference datasets
- Good for data integration

**Limitations:**
- Requires reference scRNA-seq data
- More complex setup
- Dependent on reference quality

### **Trajectory Inference Methods**

#### **CellRank** 🎯
**Best for:**
- When velocity data is available
- Complex transition dynamics
- Absorbing states identification
- Publication-quality trajectory analysis

**Advantages:**
- Uses velocity to guide trajectories
- Identifies terminal states automatically
- Handles complex branching patterns
- Well-validated method

**Limitations:**
- Requires velocity computation first
- More computationally intensive
- Complex parameter tuning

#### **Palantir** 🌟
**Best for:**
- General-purpose trajectory inference
- Multi-branched developmental processes
- When velocity data is limited
- Robust trajectory computation

**Advantages:**
- Works with or without velocity
- Handles complex branching well
- Good for exploratory analysis
- Stable across different datasets

**Limitations:**
- May not utilize velocity optimally
- Requires manual root selection
- Less integration with spatial info

#### **DPT** 📏
**Best for:**
- Simple trajectory inference
- Quick exploratory analysis
- Linear or simple branching processes
- Fallback when other methods fail

**Advantages:**
- Fast and simple
- Well-established method
- Minimal parameter tuning
- Good for simple systems

**Limitations:**
- Does not use velocity information
- Limited for complex trajectories
- Less complex than modern methods

## **My Recommendations:**

### **Typical Workflow:**
1. **Velocity Computation**: scVelo (quick) OR VELOVI (accurate)
2. **Trajectory Inference**: CellRank (with velocity) OR Palantir (general) OR DPT (simple)

### **Specific Recommendations:**

**For Beginners:** 
- **Velocity**: Start with **scVelo** to learn concepts
- **Trajectory**: Use **CellRank** to utilize velocity information

**For Publication:** 
- **Velocity**: Use **VELOVI** for accurate velocity computation
- **Trajectory**: Use **CellRank** for detailed trajectory inference

**For Exploration:** 
- **Velocity**: Try **both scVelo and VELOVI** and compare
- **Trajectory**: Compare **CellRank and Palantir** results

**For Integration:** 
- **Velocity**: Consider **SIRV** if you have references
- **Trajectory**: Use **Palantir** for robust inference

**Your dataset** (brain development) would benefit from:
1. **scVelo first** → **CellRank** - quick velocity + trajectory insights
2. **VELOVI** → **CellRank** - refined velocity + trajectory analysis
3. **Compare both workflows** - validate findings"

---

## Step 6: Interpreting Biological Significance

### 💬 Making Sense of Your Results

---

**👤 You:** "I have all these velocity and trajectory results now. How do I interpret them biologically? What do they mean for my research?"

**🤖 ChatSpatial:** "Good question! Let me help you interpret these results in biological context:

## 🧬 Biological Interpretation Guide

### **RNA Velocity Patterns Mean:**

**High Velocity Regions:**
- **Active cellular transitions** happening
- **Developmental hotspots** or **response zones**
- **Unstable cell states** (transitioning)
- **Key regulatory events** occurring

**Low Velocity Regions:**
- **Stable, mature cell states**
- **Homeostatic maintenance** modes
- **Terminal differentiation** achieved
- **Quiescent cell populations**

**Directional Flow Patterns:**
- **Developmental progressions** (stem → mature)
- **Migration pathways** (source → destination)
- **Response cascades** (stimulus → response)
- **Tissue organization** principles

### **Your Brain Development Results Show:**

**Biological Processes:**
1. **Neurogenesis**: Active stem cell differentiation
2. **Migration**: Radial movement toward cortex  
3. **Specification**: Cell fate determination
4. **Maturation**: Functional neuron development
5. **Organization**: Layer-specific positioning

**Spatial Organization:**
- **Inside-out development**: Earlier neurons in deeper layers
- **Radial glia scaffolding**: Guiding migration paths
- **Regionalization**: Different areas develop differently
- **Temporal waves**: Coordinated developmental timing

### **Research Implications:**

**What Your Data Reveals:**
1. **Normal development patterns** - your tissue is developing correctly
2. **Timing of processes** - when different events occur
3. **Spatial coordination** - how location affects development
4. **Key regulatory genes** - which genes control transitions

**Follow-up Research Questions:**
- Which **signaling pathways** regulate these transitions?
- How do **perturbations** (disease/drugs) affect these patterns?
- What **cell-cell interactions** guide development?
- Which **transcription factors** drive state changes?

**Experimental Validation:**
1. **In situ hybridization** - validate top velocity genes
2. **Immunostaining** - confirm protein expression patterns
3. **Perturbation studies** - test predicted regulators
4. **Single-cell validation** - verify computational predictions"

**👤 You:** "This is helpful! Can you suggest specific experiments to validate these findings?"

**🤖 ChatSpatial:** "Absolutely! Here are targeted experiments to validate your velocity and trajectory findings:

## 🔬 Validation Experiment Guide

### **Immediate Validation (1-2 weeks):**

**1. In Situ Hybridization**
- Target: **Top 5 velocity genes** (Neurod1, Pax6, Tbr1, Dlx2, Dcx)
- Prediction: Expression gradients matching trajectory paths
- Outcome: Validate spatial expression patterns

**2. Immunofluorescence**
- Markers: **Pax6** (stem), **Tbr1** (neurons), **Map2** (mature)
- Prediction: Spatial gradient from ventricular → cortical
- Outcome: Confirm protein-level trajectory progression

### **Extended Validation (1-2 months):**

**3. Time-course Analysis**
- Method: Multiple developmental timepoints
- Prediction: Earlier timepoints show more ventricular velocity
- Outcome: Validate temporal dynamics

**4. Perturbation Studies**
- Target: **Key velocity regulators** (Pax6, Neurod1)
- Method: Knockdown/overexpression
- Prediction: Altered trajectory patterns
- Outcome: Test functional importance

### **Complex Validation (3-6 months):**

**5. Live Imaging**
- Method: Time-lapse of developing tissue
- Prediction: Cell movements match trajectory predictions
- Outcome: Direct observation of dynamics

**6. Single-cell Validation**
- Method: scRNA-seq time series
- Prediction: Individual cells follow predicted paths
- Outcome: Cell-level trajectory validation

**Priority Ranking:**
1. **In situ + IF** (essential, quick)
2. **Perturbation** (functional validation)  
3. **Time-course** (temporal validation)
4. **Live imaging** (direct validation)

Would you like detailed protocols for any of these experiments?"

---

## Common Questions and Troubleshooting

### 💬 "My data does not have spliced/unspliced layers. Can I still do velocity analysis?"

**🤖 ChatSpatial:** "Yes! You have several options:

**Option 1: SIRV Method**
- Use a reference scRNA-seq dataset with velocity layers
- ChatSpatial can transfer velocity from reference to your spatial data
- Works well for similar tissue types

**Option 2: Trajectory Without Velocity**
- Focus on trajectory analysis using expression patterns only  
- Methods like Palantir and DPT work without velocity
- Still provides pseudotime and developmental insights

**Option 3: Expression Dynamics**
- Analyze temporal patterns in gene expression
- Use spatial patterns to infer developmental processes
- Less direct but still informative

Just ask: 'Can you analyze trajectories without velocity data?'"

### 💬 "The velocity patterns look noisy. How can I improve them?"

**🤖 ChatSpatial:** "Noisy velocity is common! Try these approaches:

**Parameter Optimization:**
- Ask: 'Can you try velocity with different parameters?'
- Increase filtering stringency
- Adjust neighborhood parameters

**Method Switching:**
- Try VELOVI: 'Can you run VELOVI for smoother results?'
- Use different velocity modes: 'Try deterministic mode'

**Data Quality:**
- Check: 'What is the quality of my spliced/unspliced data?'
- Consider preprocessing improvements

**Spatial Smoothing:**
- Ask: 'Can you spatially smooth the velocity results?'
- Incorporate spatial information more strongly"

### 💬 "How do I know if my trajectories are biologically meaningful?"

**🤖 ChatSpatial:** "Good question! Look for these validation signs:

**Biological Consistency:**
- Trajectories follow known developmental patterns
- Gene expression changes match literature
- Spatial organization makes anatomical sense
- Timing aligns with known processes

**Technical Quality:**
- High velocity confidence scores
- Consistent results across methods
- Stable across parameter changes
- Reproducible in replicate samples

**Experimental Support:**
- Marker genes match predicted states
- Perturbation studies confirm predictions
- Independent datasets show similar patterns

Ask: 'How confident should I be in these trajectory results?'"

## Parameter Selection Guide

### 🎯 Choosing the Right Parameters

**RNA Velocity Parameters:**

```text
Mode Selection:
- "stochastic": Good for most analyses (recommended)
- "deterministic": For clean, consistent data  
- "dynamical": For complex temporal dynamics

Number of PCs:
- 30-50: Standard for most datasets
- Higher: For complex, multi-lineage data
- Lower: For simple systems or small datasets
```

**Trajectory Parameters:**

```text
Method Selection:
- "cellrank": Best when velocity is available
- "palantir": Good general-purpose method
- "dpt": Simple fallback method

Spatial Weight:
- 0.3-0.5: Balanced spatial-expression integration
- Higher: More spatial influence  
- Lower: More expression-based
```

**VELOVI Parameters (for velocity computation only):**

```text
Training Epochs:
- 500-1000: Standard training for velocity computation
- More: For complex velocity dynamics
- Less: For quick velocity exploration

Latent Dimensions:
- 10-20: Standard for most velocity data
- Higher: For very complex velocity systems
- Lower: For simple velocity dynamics

Note: After VELOVI computes velocities, use trajectory methods:
- CellRank (recommended with VELOVI velocities)
- Palantir (general trajectory inference)
- DPT (simple pseudotime)
```

## Next Steps and Additional Applications

### Recommended Analysis Workflow

### 🔄 Complete Velocity + Trajectory Pipeline

**Step 1: Choose Velocity Method**
- **Quick exploration**: "Compute RNA velocity with scVelo"
- **High accuracy**: "Compute RNA velocity with VELOVI"
- **No velocity layers**: "Use SIRV with reference data"

**Step 2: Trajectory Inference** 
- **With velocity**: "Infer trajectories with CellRank using the velocity data"
- **General purpose**: "Infer trajectories with Palantir"
- **Simple analysis**: "Compute pseudotime with DPT"

**Step 3: Integration & Validation**
- **Compare methods**: "Compare CellRank vs Palantir trajectories"
- **Spatial validation**: "Show trajectory patterns on tissue"
- **Biological validation**: "Identify marker genes along trajectories"

### 🚀 Continue Your Analysis

**Ready to explore more?** Try these additional analyses:

- **"Can you compare velocity patterns between conditions?"**
- **"Show me how trajectories change with drug treatment"**
- **"Which pathways regulate these developmental transitions?"**
- **"Can you predict future cell states based on current velocity?"**
- **"How do spatial constraints affect trajectory progression?"**

### 🔬 Integration with Other Analyses

**Combine trajectory analysis with:**
- [Cell Communication Analysis](../analysis/cell_communication_analysis.md) - How signaling guides trajectories
- [Cell Type Annotation](../analysis/cell_type_annotation.md) - Identity along trajectories  
- [Visualization Tutorial](../core/visualization_tutorial.md) - Trajectory visualization plots

### 📚 Additional Topics

**Explore further approaches:**
- **Multi-modal trajectory analysis** (RNA + protein + chromatin)
- **Cross-species trajectory comparison** 
- **Disease trajectory analysis** (normal vs pathological)
- **Drug response trajectories** (treatment effects over time)

## Summary

RNA velocity and trajectory analysis reveal the **dynamic story** hidden in your spatial transcriptomics data. Through natural conversation with ChatSpatial, you can:

- **Uncover cellular dynamics** without complex programming
- **Discover developmental pathways** through intuitive questions
- **Compare different methods** to find the best approach
- **Interpret biological significance** with expert guidance
- **Design validation experiments** based on computational predictions

**Remember**: Every trajectory tells a story about cellular journeys through development, disease, or response. ChatSpatial helps you read that story and understand its biological significance.

**Key Takeaway**: The combination of RNA velocity (showing direction) and trajectory analysis (showing paths) provides insights into how cellular processes unfold in space and time. Your tissue is not just a static snapshot - it is a dynamic system in motion, and these analyses help you see that motion clearly.

Start your trajectory exploration today: *"Can you analyze the developmental dynamics in my tissue?"*