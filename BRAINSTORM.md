# Visualization package — review and roadmap to publication

## 1. State of the package (May 2026)

Seven self-contained R packages, each producing 3 multi-panel figures.
Code style is consistent (named palettes, deterministic synthetic data,
patchwork composition, `<pkg>_save_demo()` entry point). Quality is good
but several panels still rely on simplifications that should be tightened
before submission to a bioinformatics journal.

### Per-package audit

| Package | Strength | Weakness to fix |
|---|---|---|
| uncertainr | Tri-channel uncertainty halos are genuinely novel | Decomposed UMAP panels share the same point set — add **cross-channel correlation map** and **per-cell uncertainty fingerprint** |
| provtrace | Consensus heatmap × method agreement is informative | Add a **fragility score panel** — per gene, distribution of ranks across pipelines |
| scalefusion | Four resolution rows is a clean idea | A4 (compartment) currently re-uses the same encoding as A3 — replace with a **scale-disagreement diagnostic** (which scale changes the answer) |
| counterpath | Counterfactual fans are clear | Add a **gene-attached counterfactual** panel: which transcripts diverge under the alternative treatment |
| intercapsule | 6-D capsule glyph is the most novel object in the suite | Crowded with 10+ pairs — add **density-aware bundling** and a **capsule cross-section** legend |
| densweave | Density terrain + valley call-out is good | Add a **rare-event significance ring** computed from density quantile permutations |
| clinladder | Multi-modal evidence ladder is decision-oriented | Add a **case-vs-cohort** comparison rung (where does this patient sit relative to similar cases) |

### Gaps not covered by any current package

These are the *novel additions* implemented in this revision:

- **Survival / prognosis** has no first-class visual element anywhere.
- **Per-patient cohort narrative** across modalities and time is missing.
- **Multi-objective treatment decisions** (efficacy × toxicity × confidence) are not visualized.

## 2. Three new packages added

### 2.1 `survivora` — prognostic landscape cartography
Treats *hazard* as a topographic surface over a low-dimensional embedding.
- **Panel A**: Embedding with hazard-ratio terrain (smoothed local Cox hazard) as the
  background and bootstrap **stability contours** showing where the prognostic call is robust.
- **Panel B**: Locality-resolved Kaplan–Meier small multiples (one per landscape region).
- **Panel C**: Per-cluster hazard-ratio forest with bootstrap CI and prognostic
  reproducibility score.

What is novel: survival is rendered as a *spatial channel of the embedding*, not a side
plot. Stability is shown as a contour, not an error bar. This is not a UMAP, not a forest
plot, not a KM curve in isolation — it is their fusion.

### 2.2 `cohortprism` — spectral patient prisms
Each patient is rendered as a vertical prism stacking continuous molecular bands, discrete
mutation ticks, and a temporal clinical-event lane.
- **Panel A**: Cohort prism wall with patients ordered by similarity dendrogram, response
  ribbon along the top.
- **Panel B**: Single-patient drill-down prism with annotated events.
- **Panel C**: Cohort-vs-case projection — where this patient sits in continuous
  multi-modal space.

What is novel: replaces OncoPrint's discrete grid with **continuous spectral bands** and
**temporal clinical events** in the same compact glyph. Designed for tumor board / molecular
review where both molecular state and clinical history must be inspected together.

### 2.3 `riskweave` — multi-objective decision frontier
Treats treatment selection as a Pareto problem.
- **Panel A**: Decision frontier in efficacy × toxicity space with confidence rings;
  Pareto-dominant options highlighted.
- **Panel B**: Modality-contribution rays — for each candidate treatment, why does the
  recommendation hold (genomic / clinical / imaging contributions).
- **Panel C**: Patient drift — how the recommendation shifts when a single piece of
  evidence is perturbed (sensitivity to missing data).

What is novel: clinical decision support normally uses tables and probabilities; here the
decision *space itself* is the figure, with sensitivity made visual.

## 3. Manuscript framing

Suggested target: *Bioinformatics* / *Genome Biology* / *Briefings in Bioinformatics*
methods or applications track.

Recommended figure ladder:
- Figure 1 — paradigm overview (one composite panel per package family).
- Figure 2 — uncertainty-native suite (`uncertainr`, `densweave`).
- Figure 3 — provenance and scale (`provtrace`, `scalefusion`).
- Figure 4 — interaction and dynamics (`intercapsule`, `counterpath`).
- Figure 5 — clinical decision suite (`clinladder`, `survivora`, `cohortprism`,
  `riskweave`).
- Supplement — failure modes, ablation, per-package code architecture.

Required before submission:
- Replace synthetic data with at least one public dataset per package
  (TCGA-LIHC for survivora, MSK-IMPACT subset for cohortprism, public IO trial
  data for riskweave).
- User study: 5–8 domain experts rate interpretability vs. matched conventional
  plots (UMAP, OncoPrint, decision tree).
- Runtime benchmarks for atlas-scale inputs.
- CRAN / Bioconductor packaging audit.
