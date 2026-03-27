import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

os.makedirs('results/figures', exist_ok=True)

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

meta = pd.read_csv('data/metadata/metadata.tsv', sep='\t', index_col=0)
colors = {'case': '#E8593C', 'control': '#1D9E75'}

# Figure 1: Shannon Diversity
shannon = pd.read_csv('data/exported/shannon/alpha-diversity.tsv', sep='\t', index_col=0)
shannon = shannon.join(meta)

fig, ax = plt.subplots(figsize=(6, 5))
for group, data in shannon.groupby('case_control'):
    ax.scatter([group]*len(data), data['shannon_entropy'],
               color=colors[group], alpha=0.7, s=80, zorder=3)
for group, data in shannon.groupby('case_control'):
    mean = data['shannon_entropy'].mean()
    x = 0 if group == 'case' else 1
    ax.hlines(mean, x-0.2, x+0.2, colors=colors[group], linewidth=2.5)
ax.set_xlabel('Group')
ax.set_ylabel('Shannon diversity index')
ax.set_title('Alpha diversity: Celiac disease vs controls', fontweight='bold')
ax.set_xticks([0, 1])
ax.set_xticklabels(['Celiac (case)', 'Control'])
plt.tight_layout()
plt.savefig('results/figures/figure1_shannon_diversity.png', dpi=300, bbox_inches='tight')
plt.close()
print("Figure 1 saved")

# Figure 2: Bray-Curtis PCoA
with open('data/exported/bray_curtis_pcoa/ordination.txt') as f:
    lines = f.readlines()

prop_idx = [i for i, l in enumerate(lines) if l.startswith('Proportion explained')][0]
props = list(map(float, lines[prop_idx + 1].strip().split('	')))

site_start = [i for i, l in enumerate(lines) if l.startswith('Site')][0] + 1
site_lines = []
for line in lines[site_start:]:
    if line.strip() == '' or line.startswith('Biplot') or line.startswith('Site constraints'):
        break
    site_lines.append(line.strip().split('\t'))

pcoa_df = pd.DataFrame(site_lines)
pcoa_df.columns = ['sample'] + [f'PC{i}' for i in range(1, len(pcoa_df.columns))]
pcoa_df = pcoa_df.set_index('sample').astype(float)
pcoa_df = pcoa_df.join(meta)

fig, ax = plt.subplots(figsize=(7, 6))
for group, data in pcoa_df.groupby('case_control'):
    ax.scatter(data['PC1'], data['PC2'],
               color=colors[group], label=group.capitalize(),
               s=100, alpha=0.8, edgecolors='white', linewidth=0.5)
ax.set_xlabel(f'PC1 ({props[0]*100:.1f}% variance explained)')
ax.set_ylabel(f'PC2 ({props[1]*100:.1f}% variance explained)')
ax.set_title('Beta diversity: Bray-Curtis PCoA', fontweight='bold')
ax.legend(title='Group', frameon=False)
ax.axhline(0, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)
ax.axvline(0, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('results/figures/figure2_bray_curtis_pcoa.png', dpi=300, bbox_inches='tight')
plt.close()
print("Figure 2 saved")

# Figure 3: Taxonomy barplot
ft = pd.read_csv('data/exported/table/feature-table.tsv', sep='\t', skiprows=1, index_col=0)
tax = pd.read_csv('data/exported/taxonomy/taxonomy.tsv', sep='\t', index_col=0)

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

top10 = ft_rel.mean(axis=1).nlargest(10).index
ft_top = ft_rel.loc[top10].copy()
ft_top.loc['Other'] = 100 - ft_top.sum()

sample_order = meta.sort_values('case_control').index
ft_top = ft_top[sample_order]

fig, ax = plt.subplots(figsize=(12, 6))
palette = sns.color_palette('tab20', len(ft_top))
bottom = np.zeros(len(ft_top.columns))
for i, (genus, row) in enumerate(ft_top.iterrows()):
    ax.bar(range(len(ft_top.columns)), row.values,
           bottom=bottom, color=palette[i], label=genus, width=0.8)
    bottom += row.values

ax.set_xticks(range(len(ft_top.columns)))
ax.set_xticklabels([f"{s}\n({meta.loc[s,'case_control']})" for s in ft_top.columns],
                   rotation=45, ha='right', fontsize=8)
ax.set_ylabel('Relative abundance (%)')
ax.set_title('Taxonomic composition: top 10 genera', fontweight='bold')
ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', frameon=False, fontsize=9)

n_cases = (meta['case_control'] == 'case').sum()
ax.axvline(n_cases - 0.5, color='black', linewidth=1.5, linestyle='--', alpha=0.7)
ax.text(n_cases/2 - 0.5, 105, 'Celiac cases', ha='center', fontsize=10, fontweight='bold', color='#E8593C')
ax.text(n_cases + (len(ft_top.columns)-n_cases)/2 - 0.5, 105, 'Controls', ha='center', fontsize=10, fontweight='bold', color='#1D9E75')

plt.tight_layout()
plt.savefig('results/figures/figure3_taxonomy_barplot.png', dpi=300, bbox_inches='tight')
plt.close()
print("Figure 3 saved")

print("\nAll figures saved to results/figures/")