import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

# Load data
immune = pd.read_csv('data/metadata/immune_markers.tsv', sep='\t', index_col=0)
ft = pd.read_csv('data/exported/table/feature-table.tsv', sep='\t', skiprows=1, index_col=0)
tax = pd.read_csv('data/exported/taxonomy/taxonomy.tsv', sep='\t', index_col=0)

# Get genus level
def get_genus(taxon):
    parts = taxon.split(';')
    for p in parts:
        if 'g__' in p and not p.strip().endswith('g__'):
            return p.strip().replace('g__', '').strip()
    return 'Unassigned'

tax['genus'] = tax['Taxon'].apply(get_genus)
ft_genus = ft.copy()
ft_genus.index = ft_genus.index.map(lambda x: tax.loc[x, 'genus'] if x in tax.index else 'Unassigned')
ft_genus = ft_genus.groupby(ft_genus.index).sum()
ft_rel = ft_genus.div(ft_genus.sum(axis=0), axis=1) * 100

# Top 15 genera
top15 = ft_rel.mean(axis=1).nlargest(15).index
ft_top = ft_rel.loc[top15].T

# Align samples
common = immune.index.intersection(ft_top.index)
immune_aligned = immune.loc[common, ['sIgA', 'TNF_alpha', 'IL6']]
ft_aligned = ft_top.loc[common]

# Compute Spearman correlations
corr_matrix = pd.DataFrame(index=top15, columns=['sIgA', 'TNF_alpha', 'IL6'])
pval_matrix = pd.DataFrame(index=top15, columns=['sIgA', 'TNF_alpha', 'IL6'])

for genus in top15:
    for marker in ['sIgA', 'TNF_alpha', 'IL6']:
        r, p = stats.spearmanr(ft_aligned[genus], immune_aligned[marker])
        corr_matrix.loc[genus, marker] = r
        pval_matrix.loc[genus, marker] = p

corr_matrix = corr_matrix.astype(float)
pval_matrix = pval_matrix.astype(float)

# Plot heatmap
fig, ax = plt.subplots(figsize=(7, 9))

sns.heatmap(corr_matrix,
            cmap='RdBu_r',
            center=0,
            vmin=-1, vmax=1,
            annot=True, fmt='.2f',
            linewidths=0.5,
            ax=ax,
            cbar_kws={'label': 'Spearman r'})

# Mark significant correlations
for i in range(len(top15)):
    for j, marker in enumerate(['sIgA', 'TNF_alpha', 'IL6']):
        if pval_matrix.iloc[i, j] < 0.05:
            ax.text(j + 0.5, i + 0.85, '*', ha='center',
                    fontsize=14, color='black', fontweight='bold')

ax.set_xticklabels(['sIgA', 'TNF-α', 'IL-6'], fontsize=11)
ax.set_title('Microbiome — immune marker correlations\n(Spearman r, * p < 0.05)',
             fontweight='bold', fontsize=12)
ax.set_ylabel('')

plt.tight_layout()
plt.savefig('results/figures/figure5_immune_correlation.png',
            dpi=300, bbox_inches='tight')
plt.close()
print("Figure 5 saved")