#!/usr/bin/env python3
"""
GENIA 1.1 RECONSTRUCTION - STAGE 1: GENOMIC SCREENING
========================================================
Transparent, reproducible scoring of 95 candidate strains from
comparative-genomic gene-copy-number data.

INPUT:
    atrazine_data.xlsx - 95 candidate strains x 23 enzyme/transporter gene
    families, quantified as gene copy number per genome from comparative
    genomic annotation (KOfam/eggNOG-based counts).

METHOD:
    Composite genomic score = weighted sum of four normalized (0-1, min-max
    across the 95-strain candidate pool) sub-scores:

        classic_score        : AtzA + AtzB (canonical atrazine hydrolysis)
        complementary_score   : CAH + AmAH + AmidAH + HalDeh (downstream/
                                 alternative ring-cleavage hydrolases)
        oxidative_score       : LiP + MnP + AAO + DyP + VDM + AlkMon +
                                 NADPH_Ox + CytP450 + OxRed + MelD + N-dealk
        transport_score       : ABC + MFS + RND + TRAP + Mem_transp + Eff_pumps

    Weights (documented, fixed a priori - not tuned post hoc):
        w_classic       = 0.35
        w_complementary = 0.20
        w_oxidative     = 0.20
        w_transport      = 0.25

    Selection rule: within each of the 4 candidate genera present in the
    panel (Pseudomonas, Agrobacterium, Bacillus, Paenibacillus), the
    strain with the highest composite genomic score is retained as the
    genus representative. This constraint (one representative per genus)
    is imposed a priori to select for taxonomic/functional complementarity
    rather than simply taking the top-4 strains genome-wide, which could
    collapse to near-redundant close relatives within a single genus.

OUTPUT:
    stage1_genomic_scores_full.csv   - all 95 candidates, ranked
    stage1_genus_representatives.csv - top strain per genus
    stage1_score_breakdown.png       - visual breakdown of scoring
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ============================================================================
# CONFIGURATION
# ============================================================================
np.random.seed(None)  # no stochastic elements in this stage - fully deterministic

WEIGHTS = {
    'classic': 0.35,
    'complementary': 0.20,
    'oxidative': 0.20,
    'transport': 0.25
}
assert abs(sum(WEIGHTS.values()) - 1.0) < 1e-9, "Weights must sum to 1.0"

GENE_GROUPS = {
    'classic': ['AtzA_EC3.8.1.8', 'AtzB_EC3.5.1.42'],
    'complementary': ['CAH_EC3.5.2.15', 'AmAH_EC3.5.1.101', 'AmidAH_EC3.5.1.102',
                       'HalDeh_EC3.8.1.5'],
    'oxidative': ['LiP_EC1.11.1.14', 'MnP_EC1.11.1.13', 'AAO_EC1.1.3.7',
                  'DyP_EC1.11.1.19', 'VDM_EC3.1.1.68', 'AlkMon_EC1.14.15.3',
                  'NADPH_Ox_EC1.6.99', 'CytP450_EC1.14.14', 'OxRed_EC1.4.99',
                  'MelD_EC3.5.4.23', 'N-dealk_EC1.14.14'],
    'transport': ['ABC_transp', 'MFS_transp', 'RND_efflux', 'TRAP_syst',
                  'Mem_transp', 'Eff_pumps']
}

TARGET_GENERA = ['Pseudomonas', 'Agrobacterium', 'Bacillus', 'Paenibacillus']

PUBLISHED_FINALISTS = {
    'PSELUT1_Pseudomonas_pergaminensis': 'P. pergaminensis',
    'C1G7_Agrobacterium_salinitolerans': 'A. salinitolerans',
    'AN11_Bacillus_pseudomycoides': 'B. pseudomycoides',
    'AN10_Paenibacillus_polymyxa': 'Pa. polymyxa'
}

# ============================================================================
# LOAD REAL DATA
# ============================================================================
df = pd.read_excel('data/atrazine_data.xlsx', sheet_name='Sheet1')
print(f"Loaded {len(df)} candidate strains, {df.shape[1]-1} gene features")

# Extract genus from strain identifier (format: CODE_Genus_species)
def extract_genus(strain_id):
    parts = strain_id.split('_')
    return parts[1] if len(parts) > 1 else parts[0]

df['Genus'] = df['Strain'].apply(extract_genus)

# Sanity check: confirm all gene columns present
all_genes = sum(GENE_GROUPS.values(), [])
missing = [g for g in all_genes if g not in df.columns]
if missing:
    raise ValueError(f"Missing expected gene columns: {missing}")

# ============================================================================
# COMPUTE SUB-SCORES (raw counts -> min-max normalized 0-1 across full panel)
# ============================================================================
for group, genes in GENE_GROUPS.items():
    raw = df[genes].sum(axis=1)
    df[f'{group}_raw'] = raw
    rmin, rmax = raw.min(), raw.max()
    df[f'{group}_score'] = (raw - rmin) / (rmax - rmin) if rmax > rmin else 0.0

# ============================================================================
# COMPOSITE GENOMIC SCORE
# ============================================================================
df['genomic_score'] = sum(WEIGHTS[g] * df[f'{g}_score'] for g in GENE_GROUPS)

df_sorted = df.sort_values('genomic_score', ascending=False).reset_index(drop=True)
df_sorted['rank_overall'] = df_sorted.index + 1

# Flag published finalists
df_sorted['is_published_finalist'] = df_sorted['Strain'].isin(PUBLISHED_FINALISTS.keys())

# ============================================================================
# GENUS-CONSTRAINED SELECTION (a priori diversity rule)
# ============================================================================
reps = []
for genus in TARGET_GENERA:
    genus_df = df_sorted[df_sorted['Genus'] == genus]
    if len(genus_df) == 0:
        print(f"WARNING: no candidates found for genus {genus}")
        continue
    top = genus_df.iloc[0]
    rank_within_genus = 1
    reps.append({
        'Genus': genus,
        'Selected_strain': top['Strain'],
        'genomic_score': top['genomic_score'],
        'rank_overall_of_95': int(top['rank_overall']),
        'n_candidates_in_genus': len(genus_df),
        'is_published_finalist': top['Strain'] in PUBLISHED_FINALISTS,
        'classic_raw': top['classic_raw'],
        'complementary_raw': top['complementary_raw'],
        'oxidative_raw': top['oxidative_raw'],
        'transport_raw': top['transport_raw']
    })

reps_df = pd.DataFrame(reps)

# ============================================================================
# REPORT
# ============================================================================
print("\n" + "="*78)
print("STAGE 1 RESULTS: GENUS-CONSTRAINED GENOMIC SCREENING")
print("="*78)
print(reps_df.to_string(index=False))

n_match = reps_df['is_published_finalist'].sum()
print(f"\n{n_match}/4 genus representatives selected by transparent genomic "
      f"scoring MATCH the strains used in the published GENIA 1.1 consortium.")

if n_match < 4:
    mismatches = reps_df[~reps_df['is_published_finalist']]
    print("\nMismatches (genomic top-scorer differs from published choice):")
    for _, row in mismatches.iterrows():
        published = [k for k, v in PUBLISHED_FINALISTS.items()
                      if row['Genus'] in v or v.split('.')[-1].strip() in k]
        print(f"  {row['Genus']}: top genomic scorer = {row['Selected_strain']} "
              f"(rank {row['rank_overall_of_95']}/95)")

# ============================================================================
# SAVE OUTPUTS
# ============================================================================
df_sorted.to_csv('outputs/stage1_genomic_scores_full.csv', index=False)
reps_df.to_csv('outputs/stage1_genus_representatives.csv', index=False)
print("\nSaved: stage1_genomic_scores_full.csv (all 95 candidates, ranked)")
print("Saved: stage1_genus_representatives.csv (per-genus top scorer)")

# ============================================================================
# FIGURE: score breakdown for genus representatives
# ============================================================================
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10

fig, ax = plt.subplots(figsize=(8, 5))
x = np.arange(len(reps_df))
width = 0.6
bottom = np.zeros(len(reps_df))

colors = {'classic': '#029E73', 'complementary': '#DE8F05',
          'oxidative': '#0173B2', 'transport': '#CC78BC'}

for group in GENE_GROUPS:
    values = []
    for _, row in reps_df.iterrows():
        strain_row = df_sorted[df_sorted['Strain'] == row['Selected_strain']].iloc[0]
        values.append(WEIGHTS[group] * strain_row[f'{group}_score'])
    values = np.array(values)
    ax.bar(x, values, width, bottom=bottom, label=group.capitalize(),
           color=colors[group], edgecolor='black', linewidth=0.8)
    bottom += values

ax.set_xticks(x)
ax.set_xticklabels([f"{g}\n{s.split('_',1)[0]}" for g, s in
                     zip(reps_df['Genus'], reps_df['Selected_strain'])],
                    fontsize=9)
ax.set_ylabel('Composite genomic score (weighted, 0-1)', fontweight='bold')
ax.set_title('Stage 1: genus-representative genomic screening\n(real gene-copy-number data only)',
              fontweight='bold', fontsize=11)
ax.legend(frameon=True, edgecolor='black', fontsize=8, loc='upper right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='y', alpha=0.25, linestyle='--')

plt.tight_layout()
plt.savefig('outputs/stage1_score_breakdown.png', dpi=300, bbox_inches='tight')
plt.savefig('outputs/stage1_score_breakdown.pdf', bbox_inches='tight')
print("Saved: stage1_score_breakdown.png/.pdf")
