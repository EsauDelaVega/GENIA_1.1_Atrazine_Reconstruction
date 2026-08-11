#!/usr/bin/env python3
"""
GENIA 1.1 RECONSTRUCTION - STAGE 2: EXPERIMENTAL VALIDATION
========================================================
Validates the four finalist strains using measured wet-lab data:

    atrazine_degradation_data.csv  - measured atrazine degradation (%),
                                      n=8 replicates/strain, individual
                                      pure-culture assays
    GENIA_Ecoplate.xlsx             - measured Biolog EcoPlate substrate
                                      utilization (31 carbon sources),
                                      individual strains + pooled
                                      communities (GENIA_1, GENIA_2)
    ecoplate_correlation_matrix.csv - pairwise Pearson correlation of
                                      EcoPlate metabolic profiles between
                                      strains (functional redundancy proxy)

This stage does NOT re-derive the selection - it documents, with real
data, how the four genus representatives selected in Stage 1 (after
applying the a priori safety/non-pathogenicity exclusion of A. tumefaciens,
see note below) perform experimentally and how functionally
non-redundant they are relative to one another.

SAFETY EXCLUSION CRITERION (documented, not data-driven):
    Two of the top-scoring genomic candidates within genus Agrobacterium
    were excluded a priori on biosafety grounds, before any experimental
    (degradation/EcoPlate) data were consulted:
        - A. tumefaciens (rank 23/95 overall) - causal agent of crown
          gall disease
        - A. rhizogenes  (rank 27/95 overall) - causal agent of hairy
          root disease
    Both are well-characterized plant pathogens and therefore unsuitable
    for a consortium intended for crop protection applications,
    irrespective of atrazine-pathway gene content. A. salinitolerans
    (rank 28/95, effectively tied with A. rhizogenes) is the highest-
    scoring non-phytopathogenic congener in the candidate panel and was
    selected on that basis. This exclusion rule is independent of, and
    precedes, any degradation-performance data.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10

FINALISTS = {
    'P_pergaminensis': 'P. pergaminensis',
    'A_salinitolerans': 'A. salinitolerans',
    'B_pseudomycoides': 'B. pseudomycoides',
    'Pa_polymyxa': 'Pa. polymyxa'
}

# ============================================================================
# 1. REAL MEASURED DEGRADATION
# ============================================================================
deg = pd.read_csv('data/atrazine_degradation_data.csv')
deg_finalists = deg[deg['Strain'].isin(FINALISTS.keys())]

deg_summary = deg_finalists.groupby('Strain')['Degradation_pct'].agg(
    ['mean', 'std', 'count']).reset_index()
deg_summary.columns = ['Strain', 'degradation_mean', 'degradation_std', 'n']
print("="*78)
print("MEASURED ATRAZINE DEGRADATION (pure culture, individual strains)")
print("="*78)
print(deg_summary.to_string(index=False))

# ============================================================================
# 2. REAL ECOPLATE FUNCTIONAL DIVERSITY
# ============================================================================
eco = pd.read_excel('data/GENIA_Ecoplate.xlsx')
eco_finalist_cols = [c for c in eco.columns if c in FINALISTS]

eco_summary = []
for strain in eco_finalist_cols:
    vals = eco[strain].values
    active = (vals > 0.25).sum()  # substrate considered "utilized" above blank threshold
    mean_activity = vals.mean()
    # Shannon diversity of substrate utilization profile
    p = vals / vals.sum()
    p = p[p > 0]
    shannon = -np.sum(p * np.log(p))
    eco_summary.append({
        'Strain': strain,
        'mean_substrate_activity': mean_activity,
        'n_substrates_active': active,
        'shannon_diversity': shannon
    })
eco_summary = pd.DataFrame(eco_summary)
print("\n" + "="*78)
print("MEASURED ECOPLATE FUNCTIONAL PROFILE (31 carbon sources)")
print("="*78)
print(eco_summary.to_string(index=False))

# ============================================================================
# 3. REAL PAIRWISE FUNCTIONAL REDUNDANCY (EcoPlate correlation matrix)
# ============================================================================
corr = pd.read_csv('data/ecoplate_correlation_matrix.csv', index_col=0)
finalist_corr = corr.loc[eco_finalist_cols, eco_finalist_cols]
print("\n" + "="*78)
print("PAIRWISE FUNCTIONAL CORRELATION AMONG FINALISTS (EcoPlate profiles)")
print("="*78)
print(finalist_corr.round(3).to_string())

pairs = []
cols = finalist_corr.columns.tolist()
for i in range(len(cols)):
    for j in range(i+1, len(cols)):
        pairs.append({
            'strain_1': cols[i], 'strain_2': cols[j],
            'pearson_r': finalist_corr.iloc[i, j]
        })
pairs_df = pd.DataFrame(pairs).sort_values('pearson_r', ascending=False)
mean_redundancy = pairs_df['pearson_r'].mean()
print(f"\nMean pairwise functional correlation among the 4 finalists: "
      f"r = {mean_redundancy:.3f}")
print("(Lower r = more functionally complementary / non-redundant substrate use)")

# ============================================================================
# 4. COMBINE INTO FINAL DOCUMENTED SELECTION TABLE
# ============================================================================
genomic = pd.read_csv('outputs/stage1_genomic_scores_full.csv')

# Direct, explicit strain-code lookup (species code -> full genomic isolate id)
GENOMIC_ID = {
    'P_pergaminensis': 'PSELUT1_Pseudomonas_pergaminensis',
    'A_salinitolerans': 'C1G7_Agrobacterium_salinitolerans',
    'B_pseudomycoides': 'AN11_Bacillus_pseudomycoides',
    'Pa_polymyxa': 'AN10_Paenibacillus_polymyxa'
}

final_table = []
for code, label in FINALISTS.items():
    g_row = genomic[genomic['Strain'] == GENOMIC_ID[code]]
    deg_row = deg_summary[deg_summary['Strain'] == code]
    eco_row = eco_summary[eco_summary['Strain'] == code]

    final_table.append({
        'Strain': label,
        'genomic_isolate_id': GENOMIC_ID[code],
        'genomic_score_stage1': g_row['genomic_score'].values[0] if len(g_row) else np.nan,
        'genomic_rank_of_95': int(g_row['rank_overall'].values[0]) if len(g_row) else np.nan,
        'measured_degradation_pct_mean': deg_row['degradation_mean'].values[0] if len(deg_row) else np.nan,
        'measured_degradation_pct_std': deg_row['degradation_std'].values[0] if len(deg_row) else np.nan,
        'n_replicates': int(deg_row['n'].values[0]) if len(deg_row) else np.nan,
        'ecoplate_substrates_active_of_31': eco_row['n_substrates_active'].values[0] if len(eco_row) else np.nan,
        'ecoplate_shannon_diversity': eco_row['shannon_diversity'].values[0] if len(eco_row) else np.nan,
    })

final_df = pd.DataFrame(final_table)
final_df.to_csv('outputs/stage2_final_validated_selection.csv', index=False)
pairs_df.to_csv('outputs/stage2_pairwise_functional_redundancy.csv', index=False)

print("\n" + "="*78)
print("SAVED: stage2_final_validated_selection.csv")
print("SAVED: stage2_pairwise_functional_redundancy.csv")
print("="*78)

# ============================================================================
# FIGURE: measured degradation + functional non-redundancy, side by side
# ============================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

# Panel A: measured degradation (real data, jittered points)
colors_list = ['#7EAED9', '#F4A582', '#92D1B3', '#B8A9D9']
np.random.seed(0)
for i, (code, label) in enumerate(FINALISTS.items()):
    vals = deg_finalists[deg_finalists['Strain'] == code]['Degradation_pct'].values
    mean, sem = vals.mean(), stats.sem(vals)
    ax1.bar(i, mean, width=0.6, color=colors_list[i], edgecolor='black',
             linewidth=1.2, alpha=0.85, zorder=2)
    ax1.errorbar(i, mean, yerr=sem, fmt='none', color='black', capsize=5,
                 elinewidth=1.5, zorder=3)
    jitter = np.random.normal(i, 0.06, len(vals))
    ax1.scatter(jitter, vals, color='white', edgecolors='black', s=35,
               linewidth=1, zorder=4)

ax1.set_xticks(range(4))
ax1.set_xticklabels(list(FINALISTS.values()), rotation=30, ha='right', fontsize=9,
                    fontstyle='italic')
ax1.set_ylabel('Measured atrazine\ndegradation (%, 48h)', fontweight='bold')
ax1.set_title('a  Real degradation data\n(pure culture, n=8/strain)', loc='left', fontweight='bold')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.grid(axis='y', alpha=0.25, linestyle='--')

# Panel B: pairwise functional correlation heatmap
im = ax2.imshow(finalist_corr.values, cmap='RdBu_r', vmin=-1, vmax=1)
ax2.set_xticks(range(4))
ax2.set_yticks(range(4))
short_labels = [FINALISTS[c] for c in eco_finalist_cols]
ax2.set_xticklabels(short_labels, rotation=30, ha='right', fontsize=8, fontstyle='italic')
ax2.set_yticklabels(short_labels, fontsize=8, fontstyle='italic')
for i in range(4):
    for j in range(4):
        ax2.text(j, i, f'{finalist_corr.values[i,j]:.2f}', ha='center', va='center',
                 fontsize=8, color='white' if abs(finalist_corr.values[i,j])>0.5 else 'black')
ax2.set_title('b  Functional redundancy\n(EcoPlate profile correlation)', loc='left', fontweight='bold')
plt.colorbar(im, ax=ax2, fraction=0.046, pad=0.04, label='Pearson r')

plt.tight_layout()
plt.savefig('outputs/stage2_validation_figure.png', dpi=300, bbox_inches='tight')
plt.savefig('outputs/stage2_validation_figure.pdf', bbox_inches='tight')
print("SAVED: stage2_validation_figure.png/.pdf")
