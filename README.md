# GENIA 1.1 — Reproducible strain-selection pipeline (atrazine)

**Companion repository to [GENIA_Framework_Biodegradation](https://github.com/EsauDelaVega/GENIA_Framework_Biodegradation)** (the general GENIA v2.0 ML framework: GAT, Node2Vec, Random Forest for multi-pollutant SynCom design).

This repository documents, independently and transparently, the specific
**four-strain atrazine consortium selection** (GENIA 1.1) used in
"Assembling the Fantastic Four" (submitted to *Environmental Science &
Technology*), in direct response to peer-review requests for reproducibility
and methodological transparency (manuscript es-2026-07925e).

Where the main GENIA framework repository documents the general-purpose GNN/
Node2Vec/Random Forest architecture applicable to arbitrary pollutants and
candidate panels, this repository documents the **specific, deterministic,
rule-based screening logic** actually used to arrive at the atrazine
four-member consortium reported in the manuscript, plus its experimental
validation against real wet-lab data. It intentionally does not use a
black-box model at this stage — every score, weight, and exclusion rule is
stated explicitly and is traceable to a data file in `data/`.

## What this is

A two-stage, fully deterministic and documented selection process:

1. **`01_genomic_screening.py`** — scores 95 candidate strains from real
   comparative-genomic gene-copy-number data (KOfam/eggNOG-based annotation),
   applies a genus-diversity constraint, and documents an explicit
   biosafety exclusion rule.
2. **`02_experimental_validation.py`** — validates the four selected strains
   against independently measured wet-lab data (pure-culture degradation
   assays, Biolog EcoPlate functional profiling).

See `SI_methods.md` for the full methodology write-up, including exact
scoring formula, weights, and an explicit discussion of limitations.

## Reproduce

```bash
pip install -r requirements.txt
python 01_genomic_screening.py
python 02_experimental_validation.py
```

Outputs (tables + figures) are written to `outputs/`.

## Data provenance

All files in `data/` are real measurements:
- `atrazine_data.xlsx` — comparative genomic gene-copy-number annotation, 95
  candidate strains.
- `atrazine_degradation_data.csv` — pure-culture atrazine degradation assays,
  n = 8 replicates/strain.
- `GENIA_Ecoplate.xlsx`, `ecoplate_correlation_matrix.csv` — Biolog EcoPlate
  functional carbon-substrate utilization profiling.

All scores in this pipeline trace directly to a raw measurement in `data/`;
see `SI_methods.md` for the exact provenance of each input file.

## License / citation

[Add license and citation info once linked to the published manuscript.]

## Relationship to the main GENIA framework

| | GENIA_Framework_Biodegradation | This repository (GENIA 1.1 atrazine) |
|---|---|---|
| Scope | General-purpose ML framework (any pollutant, any candidate panel) | Specific atrazine 4-strain consortium selection |
| Model | GAT + Node2Vec + Random Forest ensemble | Deterministic weighted rule-based scoring |
| Purpose | Reusable methodology | Reproducibility/transparency artifact for peer review |
| Status | Active development, v2.0 | Frozen snapshot tied to manuscript es-2026-07925e |
