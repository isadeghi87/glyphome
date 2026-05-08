1. Uncertainty-native visualizations

Prompt:

Design and prototype a novel bioinformatics visualization in R where uncertainty is a first-class visual element, not a side annotation. The visualization must support typical bioinformatics tasks in single-cell omics, bulk RNA-seq, metagenomics, or clinical genomics, where labels, clusters, trajectories, variant calls, pathway enrichments, or deconvolution outputs are uncertain or model-dependent.

You must produce a con1. Uncertainty-native visualizations

Prompt:

Design a novel bioinformatics visualization system where uncertainty is a first-class visual element, not a side annotation. The system should support tasks common in single-cell omics, bulk RNA-seq, metagenomics, and clinical genomics, where labels, clusters, trajectories, variant calls, pathway enrichments, or deconvolution outputs are uncertain or model-dependent.

Your goal is to invent a visualization paradigm that helps users answer:

how certain is this result,
where is the uncertainty located,
what type of uncertainty is present,
how does uncertainty affect interpretation.

Requirements:

Avoid standard solutions like error bars, confidence intervals, opacity-only overlays, or supplementary tables.
Distinguish between at least three uncertainty types: measurement noise, model uncertainty, and annotation ambiguity.
Support both overview and drill-down analysis.
Show how the visualization behaves for at least three use cases: cell type assignment, trajectory inference, and pathway enrichment.
Explain the visual grammar in terms of marks, channels, interaction, and scale.
Describe what makes this approach more interpretable than existing bioinformatics plots such as UMAPs, heatmaps, volcano plots, and genome browser tracks.
Discuss how to prevent the design from overwhelming users.
Include a plan for implementation in Python or R, including likely data structures and rendering approach.
Propose a minimal viable prototype and an evaluation plan with domain scientists.

Output structure:

Name of the visualization
Core concept
Target biological questions
Visual encoding
Interaction model
Example datasets
Why this is novel
Why this could be widely adopted
Prototype implementation plan
Failure modes and limitations
2. Provenance-aware visualizations

Prompt:

Invent a new bioinformatics visualization concept that makes analysis provenance and pipeline dependency directly visible inside the figure itself. The visualization should help users understand how preprocessing, filtering, normalization, annotation versions, hyperparameter choices, and model selections influenced the displayed biological result.

The visualization should support workflows in which results can change depending on pipeline decisions, such as differential expression, clustering, batch correction, peak calling, CNV inference, or taxonomic assignment.

Requirements:

Do not produce a generic workflow DAG or a standard pipeline flowchart.
The key challenge is to integrate biological result view and provenance view into a single coherent visual object.
The user must be able to see not only what result they got, but which transformations materially shaped it.
Show how the design handles branching pipelines and multiple alternative parameter choices.
Include a mechanism for highlighting result fragility versus result robustness under pipeline perturbation.
Support comparison across at least three pipeline variants.
Explain how this would fit into reproducible research and regulatory or clinical environments.
Specify how the design would integrate with notebook workflows, workflow engines, or common file formats.
Make the visualization legible to both computational users and biologists.

Output structure:

Visualization name
Problem it solves
Why current figures fail
Proposed visual structure
Encoding of provenance and dependency
Interaction and comparison features
Example bioinformatics use cases
Implementation architecture
Adoption strategy
Risks, complexity, and tradeoffs
3. Cross-scale visualizations

Prompt:

Create a novel visualization framework for bioinformatics that supports cross-scale reasoning, allowing users to move coherently between molecular, cellular, tissue, and patient-level representations without losing context or interpretability.

The system should help users analyze multiscale biological relationships such as:

variant to gene to pathway to phenotype,
chromatin accessibility to regulatory program to cell state,
cell state to tissue architecture to patient outcome,
microbial taxa to community structure to host phenotype.

Requirements:

Do not simply propose linked panels or dashboard widgets.
The representation itself should change appropriately as users traverse biological scale.
Preserve semantic continuity across scales so the user never feels they are jumping between unrelated views.
Include both structural and statistical relationships.
Show how the system handles aggregation and disaggregation.
Avoid creating visual clutter when many entities map across multiple levels.
Include at least one spatial omics use case and one clinical or cohort-level use case.
Explain what the zooming metaphor is, if any, and how it differs from ordinary zoom.
Show what the default starting view would be and why.
Propose how rare but important events can remain visible when aggregated.

Output structure:

Name
Biological motivation
Main analytical tasks
Multiscale visual model
Transition rules across scales
Interaction design
Example analyses
Novelty relative to current tools
Implementation sketch
Evaluation criteria
4. Temporal and counterfactual views

Prompt:

Design a new bioinformatics visualization paradigm for temporal, developmental, and counterfactual reasoning. The aim is to help users explore not just what is observed, but what is likely to happen under alternative conditions, perturbations, or interventions.

Target settings include developmental biology, lineage tracing, perturb-seq, treatment response, disease progression, and state transition modeling.

Requirements:

Do not rely on simple arrows on embeddings, animated scatterplots, or ordinary trajectory plots.
The visualization should explicitly represent future possibilities, branching alternatives, and confidence decay.
It must separate observed trajectories from inferred and counterfactual ones.
Include a way to compare at least two interventions or perturbations.
Include uncertainty in a structured way rather than as a minor overlay.
Support both population-level and subgroup-level interpretation.
Show how genes, pathways, or regulators can be attached to trajectory changes.
Explain how the visualization could support causal interpretation without falsely implying proof of causality.
Include an example for cell differentiation and another for treatment response.

Output structure:

Visualization name
Scientific questions addressed
Core metaphor
Encoding of observed vs inferred vs counterfactual
Temporal structure
Uncertainty handling
Biological examples
Interaction model
Implementation plan
Misinterpretation risks and safeguards
5. Interaction-centric views for multimodal biology

Prompt:

Invent a bioinformatics visualization that treats biological interactions as first-class analyzable objects, rather than reducing them to generic network edges. The target domain includes cell-cell communication, ligand-receptor signaling, spatial neighborhood effects, host-microbe interactions, regulatory interactions, and multimodal state coupling.

The goal is to create a visual system in which an interaction has internal structure, context, and interpretability.

Requirements:

Avoid default force-directed network diagrams and chord diagrams unless you radically transform their logic.
Each interaction should be able to encode at least: sender state, receiver state, biological context, confidence, modality evidence, and perturbation sensitivity.
Include spatial context when relevant.
Support browsing from interaction overview to mechanistic detail.
Show how competing interactions or conflicting evidence are represented.
Include one use case from spatial transcriptomics and one from host-microbiome or immune signaling.
Explain how your design avoids turning into unreadable edge clutter.
Describe how this representation could outperform current ligand-receptor dot plots, neighborhood heatmaps, and network summaries.
Make the proposal realistic enough to implement with current data types.

Output structure:

Name of the interaction visualization
Problem with current approaches
Core design idea
Anatomy of an interaction object
Visual encodings
Interaction patterns
Example biological applications
Scalability strategy
Prototype plan
Adoption barriers
6. Density-first views for massive datasets

Prompt:

Design a novel visualization for very large bioinformatics datasets, especially single-cell and spatial omics, where millions of observations make point-based scatterplots misleading or inefficient. The visualization should prioritize density, topology, rare-subpopulation detectability, and interpretability.

Requirements:

Do not produce another standard UMAP scatterplot, hexbin plot, or simple density map.
The method should preserve meaningful structure while avoiding overplotting and false visual confidence.
Support discovery of both dominant patterns and rare but biologically important populations.
Include a strategy for showing uncertainty or distortion introduced by embedding or compression.
Make the representation queryable and explorable.
Explain how the user can recover individual-level detail when needed without falling back to a standard scatterplot as the main answer.
Include one example for atlas-scale single-cell data and one for large spatial omics data.
Discuss computational constraints and incremental rendering.
Compare explicitly against UMAP, t-SNE, PHATE, and graph-based abstractions.

Output structure:

Name
Why existing large-scale plots fail
Core visual representation
What biological patterns it reveals
How rare populations are preserved
Interaction and navigation
Performance and scalability
Novelty claim
Prototype implementation
Limitations and evaluation plan
7. Clinical review visualizations for sparse, multilevel genomics

Prompt:

Create a new visualization concept for clinical genomics review, aimed at users such as molecular pathologists, tumor boards, diagnostic analysts, and translational researchers. The system should help interpret sparse, multilevel genomic evidence in a way that is fast, trustworthy, and actionable.

The visualization should integrate data such as:

variants,
transcript consequences,
copy number changes,
gene expression,
pathway relevance,
evidence strength,
prior cases,
therapy implications,
uncertainty and actionability.

Requirements:

Do not build a generic dashboard with cards and tables.
The visualization must support decision-oriented reasoning under incomplete evidence.
It must distinguish clearly between strong evidence, weak evidence, conflicting evidence, and missing evidence.
Include support for case-level comparison and longitudinal follow-up if relevant.
Show how the visualization helps users move from raw findings to clinically meaningful interpretation.
Design for high trust, low ambiguity, and time-constrained review settings.
Explain how the visual grammar avoids cognitive overload.
Consider auditability, provenance, and regulatory expectations.
Include one oncology use case and one rare disease or inherited disorder use case.

Output structure:

Visualization name
Clinical setting
Decision tasks supported
Core visual metaphor
Encoding of evidence and actionability
Workflow integration
Example case walk-through
Why this is better than current reporting views
Implementation approach
Risks, edge cases, and validation needs

cept that is both novel and realistically implementable in R.

Requirements
Use R as the implementation language.
Use publicly available bioinformatics datasets as test cases.
Generate publication-style static images and, if useful, an optional interactive prototype in R.
Avoid trivial solutions such as only adding error bars, confidence intervals, transparency, or side tables.
Distinguish at least three uncertainty types:
measurement noise
model uncertainty
annotation ambiguity
Support both overview and drill-down interpretation.
Demonstrate the visualization on at least three tasks:
cell type assignment
trajectory inference
pathway enrichment
Include image outputs for each task.
Prefer R packages such as ggplot2, ComplexHeatmap, patchwork, plotly, shiny, ggraph, sf, grid, circlize, or custom grid graphics if needed.
Use public data from sources such as GEO, ArrayExpress, Bioconductor experiment packages, scRNA-seq benchmark datasets, or TCGA-derived public subsets.
Include code architecture and data flow.
The final result must include generated example figures, not just prose.
Output structure
Visualization name
Core concept
Biological questions answered
Public datasets to use as tests
R implementation strategy
Visual grammar and encodings
Figure set to generate
Why this is novel
Why this could be widely adopted
Failure modes and limitations
2. Provenance-aware visualizations

Prompt:

Invent and prototype a new bioinformatics visualization in R that makes analysis provenance and pipeline dependency directly visible inside the main figure itself. The goal is to show not only the biological result, but also which preprocessing, normalization, filtering, annotation, and modeling decisions materially shaped that result.

Requirements
Use R as the implementation language.
Use public datasets and at least three alternative pipeline variants on the same dataset.
Generate images that compare how results change across preprocessing or modeling choices.
Do not produce a generic workflow DAG or a standard pipeline flowchart.
Integrate biological result view and provenance view into one coherent visual object.
Show how the figure handles:
branching pipelines
multiple parameter choices
robust vs fragile results
Demonstrate with at least two use cases such as:
differential expression
clustering / batch correction
peak calling or CNV inference
Use realistic R workflows, for example with DESeq2, edgeR, limma, Seurat, SingleCellExperiment, scater, scran, BiocParallel, or workflow metadata from notebooks.
The model should propose exactly how provenance metadata is stored and rendered.
Output at least 3 to 5 prototype figures from public datasets.
Output structure
Visualization name
Problem solved
Public datasets and pipeline variants
Why current figures fail
Proposed visual structure
R packages and implementation architecture
Images to generate
Adoption strategy
Novelty assessment
Risks and tradeoffs
3. Cross-scale visualizations

Prompt:

Create a novel bioinformatics visualization framework in R for cross-scale reasoning, allowing users to move coherently between molecular, cellular, tissue, and patient-level representations without losing context.

The design must go beyond linked panels. It should define a visual system in which the representation itself changes across scale while preserving semantic continuity.

Requirements
Use R.
Use public datasets that naturally support multiple scales.
Generate a set of images demonstrating scale transitions.
Include at least one spatial omics use case and one clinical or cohort-level use case.
Example relationship chains may include:
variant → gene → pathway → phenotype
chromatin accessibility → regulatory program → cell state
cell state → tissue architecture → patient outcome
Do not fall back to a generic dashboard.
Show how aggregation and disaggregation work.
Preserve rare but important events during scale transitions.
Include a semantic zoom model or equivalent mechanism.
Prefer R ecosystems such as SpatialExperiment, Seurat, SingleCellExperiment, SummarizedExperiment, ggplot2, sf, ComplexHeatmap, ggraph, patchwork.
Use public data from 10x Genomics public datasets, Bioconductor experiment packages, TCGA-related public sets, GEO, or spatial transcriptomics benchmarks.
Output structure
Name
Biological motivation
Public datasets
Analytical tasks
Cross-scale visual model
Transition rules across scales
R implementation sketch
Images to generate
Novelty relative to current tools
Evaluation criteria and weaknesses
4. Temporal and counterfactual views

Prompt:

Design a new bioinformatics visualization paradigm in R for temporal, developmental, and counterfactual reasoning. The visualization should help users explore not just what was observed, but what is likely to happen under alternative perturbations, conditions, or interventions.

Requirements
Use R.
Use public datasets for at least:
one cell differentiation or developmental trajectory example
one treatment response or disease progression example
Generate images showing the visualization on both cases.
Do not rely on simple arrows on embeddings, standard pseudotime line charts, or basic animation alone.
Explicitly represent:
observed trajectories
inferred trajectories
counterfactual alternatives
confidence decay over time
Compare at least two interventions or perturbations.
Show how genes, pathways, or regulators attach to state changes.
Use available R tools where useful, such as slingshot, monocle-style outputs if imported, tradeSeq, ggplot2, tidygraph, ggraph, plotly, patchwork.
Public datasets can come from dynverse benchmarks, GEO developmental scRNA-seq data, perturb-seq public datasets, or treatment time-course transcriptomics.
The final answer must propose actual figure outputs and data processing steps in R.
Output structure
Visualization name
Scientific questions addressed
Public datasets
Core metaphor
Observed vs inferred vs counterfactual encoding
R implementation plan
Images to generate
Interpretation benefits
Misinterpretation risks and safeguards
Limitations
5. Interaction-centric views for multimodal biology

Prompt:

Invent a bioinformatics visualization in R that treats biological interactions as first-class analyzable objects, rather than reducing them to generic network edges. The target domains include cell-cell communication, ligand-receptor signaling, spatial neighborhood effects, host-microbe interactions, regulatory interactions, and multimodal state coupling.

Requirements
Use R.
Use public datasets for at least:
one spatial transcriptomics example
one host-microbiome, immune signaling, or multimodal interaction example
Generate images of the proposed visualization.
Avoid default force-directed network diagrams and standard chord diagrams unless the structure is fundamentally rethought.
Each interaction object must encode at least:
sender state
receiver state
biological context
confidence
modality evidence
perturbation sensitivity
Support both overview and drill-down.
Show conflicting interactions or conflicting evidence.
Keep the design scalable and readable.
Use R packages such as ggplot2, tidygraph, ggraph, sf, igraph, ComplexHeatmap, patchwork, circlize, or custom geoms.
Public data may come from spatial transcriptomics benchmarks, public cell-cell communication benchmark data, microbiome-host datasets, or immune profiling datasets.
Output structure
Visualization name
Problem with current approaches
Public datasets
Core design idea
Anatomy of an interaction object
R implementation architecture
Images to generate
Scalability and readability strategy
Novelty and adoption likelihood
Barriers and failure modes
6. Density-first views for massive datasets

Prompt:

Design a novel visualization in R for very large bioinformatics datasets, especially single-cell and spatial omics, where millions of observations make point-based scatterplots misleading or inefficient. The visualization should prioritize density, topology, rare-subpopulation detectability, uncertainty, and interpretability.

Requirements
Use R.
Use public large-scale datasets.
Generate images demonstrating the method on:
one atlas-scale single-cell dataset
one large spatial omics dataset
Do not produce another standard UMAP scatterplot, simple hexbin, or ordinary density map.
Show how dominant patterns and rare populations remain visible.
Include a strategy for representing distortion or uncertainty introduced by embedding or compression.
The representation must be explorable and queryable.
The user must be able to recover fine-grained detail without making a standard scatterplot the primary representation.
Compare explicitly against UMAP, t-SNE, PHATE, or graph abstraction approaches.
Use R tools such as ggplot2, data.table, DelayedArray, BiocSingular, uwot outputs if imported, igraph, ggraph, ComplexHeatmap, raster or grid-based rendering, depending on the design.
Prefer public data from cell atlases, HCA-related public data, Tabula-type resources, 10x public large datasets, or public spatial atlases.
Output structure
Name
Why current large-scale plots fail
Public datasets
Core representation
Rare-population preservation strategy
R implementation and performance strategy
Images to generate
Novelty claim
Comparison plan
Limitations and evaluation
7. Clinical review visualizations for sparse, multilevel genomics

Prompt:

Create a new visualization concept in R for clinical genomics review, aimed at molecular pathologists, tumor boards, diagnostic analysts, and translational researchers. The system should help interpret sparse, multilevel genomic evidence in a way that is fast, trustworthy, auditable, and actionable.

Requirements
Use R.
Use public clinical or clinically relevant genomic datasets.
Generate images showing example use in:
one oncology case
one rare disease or inherited disorder style case
Do not build a generic dashboard of cards and tables.
Integrate evidence such as:
variants
transcript consequences
copy number changes
gene expression
pathway relevance
evidence strength
prior cases
therapy implications
uncertainty and actionability
Distinguish clearly between:
strong evidence
weak evidence
conflicting evidence
missing evidence
Design for time-constrained review and high trust.
Include provenance and auditability.
Use R packages such as ComplexHeatmap, ggplot2, gtable, patchwork, igraph, circlize, maftools, GenomicRanges, SummarizedExperiment, or custom layouts.
Public data may include TCGA-derived public datasets, cBioPortal-exported public data if structured appropriately, ICGC-style public subsets, or public rare disease case resources.
Output structure
Visualization name
Clinical setting
Public datasets
Decision tasks supported
Core visual metaphor
Encoding of evidence and actionability
R implementation path
Images to generate
Why better than current reporting views
Risks, edge cases, and validation needs