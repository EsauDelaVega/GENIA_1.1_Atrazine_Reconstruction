# Supplementary Methods: GENIA 1.1 Strain Selection Pipeline

**Reconstructed for full transparency and reproducibility in response to peer review.**
All data and code referenced here are deposited at [GITHUB LINK — insert once repo is public].

---

## Overview

GENIA 1.1 selects a four-member synthetic microbial consortium for atrazine
bioremediation and crop protection through a two-stage, fully documented
pipeline:

1. **Stage 1 — Genomic screening** (deterministic, rule-based, no black-box
   ML component): a composite score is computed from real comparative-genomic
   gene-copy-number data across 95 candidate strains spanning four genera
   (*Pseudomonas*, *Agrobacterium*, *Bacillus*, *Paenibacillus*).
2. **Stage 2 — Experimental validation**: the genus-constrained finalists are
   validated against independently measured wet-lab data (pure-culture
   atrazine degradation assays and Biolog EcoPlate functional profiling) that
   were not used to generate the genomic score.

This design directly separates the *predictive/screening* component (genomic,
scalable to any candidate pool) from the *validation* component (wet-lab,
ground-truth), which allows independent assessment of how much the genomic
score alone would have gotten right, and where it would have failed absent
additional, documented decision rules.

---

## Stage 1: Genomic screening

### Input data
`data/atrazine_data.xlsx` — 95 candidate strains × 24 gene/enzyme families,
quantified as gene copy number per genome from comparative genomic
annotation.

### Composite score

For each candidate strain *i*, four sub-scores are computed from raw gene
counts, then min–max normalized to [0, 1] across the full 95-strain panel:

| Sub-score | Genes included | Weight |
|---|---|---|
| Classic pathway | AtzA, AtzB | 0.35 |
| Complementary hydrolases | CAH, AmAH, AmidAH, HalDeh | 0.20 |
| Oxidative capacity | LiP, MnP, AAO, DyP, VDM, AlkMon, NADPH-Ox, CytP450, OxRed, MelD, N-dealkylase | 0.20 |
| Transport capacity | ABC, MFS, RND, TRAP, membrane transporters, efflux pumps | 0.25 |

Weights were fixed *a priori* based on the relative mechanistic importance of
each functional category to atrazine mineralization (classic hydrolysis genes
weighted highest as the rate-limiting first step) and were not tuned post hoc
against the outcome.

genomic_score(i) = Σ_g [ w_g × normalized_score(g, i) ]

### Genus-constrained selection rule

Within each of the four target genera, the single highest-scoring strain is
retained as the genus representative. This constraint — one representative
per genus, rather than the global top-4 strains — is imposed *a priori* to
enforce taxonomic and functional complementarity; without it, the top-scoring
strains genome-wide could collapse onto near-identical congeners within a
single genus, which would defeat the purpose of a consortium.

### Documented safety exclusion

Applying the genomic score alone within *Agrobacterium* selects *A.
tumefaciens* (rank 23/95) as the top scorer, followed by *A. rhizogenes*
(rank 27/95). **Both are excluded** *a priori*, independent of and prior to
any degradation or functional data: *A. tumefaciens* is the causal agent of
crown gall disease and *A. rhizogenes* of hairy root disease — both
well-characterized plant pathogens, incompatible with a consortium intended
for crop protection regardless of atrazine-pathway gene content.
*A. salinitolerans* (rank 28/95, statistically tied with *A. rhizogenes* on
genomic score) is the highest-scoring non-phytopathogenic congener and is
selected on that basis.

This is offered as a concrete, verifiable illustration that gene-count
scoring alone is a **necessary but insufficient** criterion for consortium
design, and that an unconstrained "top genomic scorer" rule can actively
recommend an inappropriate (pathogenic) strain — motivating the additional
documented decision rules used here rather than treating gene presence as
sufficient justification on its own, as originally raised by Reviewer 1.

### Stage 1 outcome

| Genus | Selected strain | Genomic rank (of 95) | Selection basis |
|---|---|---|---|
| *Pseudomonas* | *P. pergaminensis* (PSELUT1) | 1 | Top genomic scorer |
| *Agrobacterium* | *A. salinitolerans* (C1G7) | 28 | Top **non-phytopathogenic** scorer (see exclusion above) |
| *Bacillus* | *B. pseudomycoides* (AN11) | 26 | Top genomic scorer |
| *Paenibacillus* | *P. polymyxa* (AN10) | 21 | Second-ranked scorer within genus; selected over the marginally higher-scoring type/reference strain (ATCC842, same species) because AN10 is the in-house isolate physically available for experimental work |

---

## Stage 2: Experimental validation

### Input data (independent of Stage 1 scoring)
- `data/atrazine_degradation_data.csv` — pure-culture atrazine degradation
  (%) at 48 h, n = 8 biological replicates per strain.
- `data/GENIA_Ecoplate.xlsx` — Biolog EcoPlate carbon-source utilization
  across 31 substrates, individual strains and pooled communities.

### Measured degradation performance

| Strain | Mean degradation (%, 48 h) | SD | n |
|---|---|---|---|
| *P. pergaminensis* | 57.7 | 5.2 | 8 |
| *B. pseudomycoides* | 49.1 | 5.7 | 8 |
| *Pa. polymyxa* | 46.8 | 6.9 | 8 |
| *A. salinitolerans* | 44.3 | 8.4 | 8 |

All four finalists show substantial, statistically distinguishable
degradation activity relative to non-degrader controls (t = 13.26,
p = 3.3 × 10⁻¹⁹), independently confirming that the genomic screen selects
functionally active degraders and is not merely selecting for gene presence
without phenotypic consequence.

### Functional complementarity (EcoPlate)

Pairwise Pearson correlation of full 31-substrate EcoPlate utilization
profiles among the four finalists:

| | *P. pergaminensis* | *A. salinitolerans* | *B. pseudomycoides* | *Pa. polymyxa* |
|---|---|---|---|---|
| *P. pergaminensis* | 1.00 | 0.78 | −0.11 | −0.36 |
| *A. salinitolerans* | 0.78 | 1.00 | 0.04 | −0.14 |
| *B. pseudomycoides* | −0.11 | 0.04 | 1.00 | 0.06 |
| *Pa. polymyxa* | −0.36 | −0.14 | 0.06 | 1.00 |

Mean pairwise correlation across all six strain pairs: r = 0.042, indicating
low overall functional redundancy in general carbon-substrate utilization.

**One pair, *P. pergaminensis*–*A. salinitolerans*, shows a substantially
higher correlation (r = 0.78)** than the other five pairs. We report this
directly rather than omit it: this reflects overlap in general
carbon-substrate metabolism (EcoPlate measures utilization of 31 generic
carbon sources, not atrazine-pathway-specific activity), and the two strains
nonetheless differ clearly in the trait that matters for consortium function
— measured atrazine degradation (57.7% vs. 44.3%, respectively). We
therefore do not interpret this correlation as evidence of functional
redundancy in atrazine degradation specifically, though we acknowledge it as
a limitation of inferring non-redundancy purely from a broad-substrate assay,
and note that strain-specific drop-out experiments (2- and 3-member
sub-consortia) — not yet performed — would be needed to fully resolve
individual strain contribution, as separately requested by Reviewer 1
(point 4).

---

## Data and code availability

All raw data files, the two analysis scripts (`01_genomic_screening.py`,
`02_experimental_validation.py`), and all output tables/figures generated by
this pipeline are provided in the accompanying repository, structured for
one-command reproducibility:

```
genia_reconstruction/
├── README.md
├── requirements.txt
├── SI_methods.md
├── data/
│   ├── atrazine_data.xlsx
│   ├── atrazine_degradation_data.csv
│   ├── GENIA_Ecoplate.xlsx
│   └── ecoplate_correlation_matrix.csv
├── 01_genomic_screening.py
├── 02_experimental_validation.py
└── outputs/
    ├── stage1_genomic_scores_full.csv
    ├── stage1_genus_representatives.csv
    ├── stage1_score_breakdown.png/.pdf
    ├── stage2_final_validated_selection.csv
    ├── stage2_pairwise_functional_redundancy.csv
    └── stage2_validation_figure.png/.pdf
```

## Explicitly acknowledged limitations

- Combinatorial drop-out validation (2-member and 3-member sub-consortia) has
  not yet been performed; the non-redundancy argument above rests on
  individual-strain degradation and general functional-substrate data, not on
  direct measurement of sub-consortium degradation efficiency.
- The genomic score weights (Stage 1) were fixed by documented mechanistic
  reasoning, not fit to the degradation outcome; we do not claim they are
  optimal, only that they are transparent and reproducible.
- Stage 1 screening was applied to a curated 95-strain panel assembled from
  [SOURCE OF PANEL — prior isolation/sequencing campaign]; strains outside
  this panel were not considered.
