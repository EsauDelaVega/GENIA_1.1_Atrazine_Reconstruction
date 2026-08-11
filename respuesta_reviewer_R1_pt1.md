# Borrador — Respuesta a Reviewer 1, punto 1 (transparencia/valor del ML)

*(en inglés, listo para pegar en la carta de respuesta; ajusta el tono/cita
de línea según el formato final)*

---

We thank the reviewer for this important point and have substantially
revised our approach to transparency. We have reconstructed and fully
documented the GENIA 1.1 strain-selection pipeline as a two-stage process
(Supplementary Methods, Section X; code and data at [GitHub link]):

**Stage 1** computes a deterministic, rule-based composite score from real
comparative-genomic gene-copy-number data (95 candidate strains, 24
gene/enzyme families spanning classic hydrolysis, complementary hydrolases,
oxidative capacity, and transport systems), with weights fixed a priori and
reported in full (Supplementary Table SX).

Applying this genomic score alone within genus *Agrobacterium* selects *A.
tumefaciens* and *A. rhizogenes* — both well-characterized plant
pathogens — as the top two scorers. We excluded both a priori on biosafety
grounds, independent of degradation performance, and selected the
highest-scoring non-phytopathogenic congener, *A. salinitolerans*, instead.
We present this as a concrete, verifiable illustration that gene-count
scoring alone is necessary but not sufficient for consortium design, and
that additional documented decision criteria (here, biosafety) are required
beyond simply identifying strains carrying "classic degradation genes,"
directly addressing the reviewer's concern that the model may function as a
purely confirmatory scoring system.

**Stage 2** validates the four selected strains against independently
measured pure-culture degradation assays (n = 8 replicates/strain) and
Biolog EcoPlate functional profiling, neither of which informed the Stage 1
score. All four strains show substantial, statistically distinguishable
degradation activity relative to non-degrader controls (t = 13.26,
p = 3.3 × 10⁻¹⁹), and pairwise EcoPlate functional-profile correlation is
low overall (mean r = 0.042 across all six strain pairs), supporting
functional complementarity among consortium members. We note one exception
transparently: *P. pergaminensis* and *A. salinitolerans* show a higher
pairwise correlation (r = 0.78) in general carbon-substrate utilization,
though they differ clearly in measured atrazine degradation specifically
(57.7% vs. 44.3%); we discuss this as a limitation in the Supplementary
Methods and note that combinatorial drop-out experiments (addressed
separately in our response to point 4) would be needed to fully resolve
individual strain contributions.

Full code, raw data, and step-by-step documentation of every scoring
decision are deposited at [GitHub link], enabling independent verification
and application of the pipeline to other candidate panels or pollutants, as
requested by both reviewers.
