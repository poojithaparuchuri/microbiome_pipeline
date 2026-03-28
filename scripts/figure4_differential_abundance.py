import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

lfc = pd.read_csv('data/exported/ancombc_genus/lfc_slice.csv', index_col=0)
se = pd.read_csv('data/exported/ancombc_genus/se_slice.csv', index_col=0)
qval = pd.read_csv('data/exported/ancombc_genus/q_val_slice.csv', index_col=0)

col = [c for c in lfc.columns if 'control' in c.lower()][0]

sig = qval[col] < 0.05
lfc_sig = lfc.loc[sig, col]
se_sig = se.loc[sig, col]

def clean_name(name):
    if name.startswith('Unassigned'):
        return 'Unassigned'
    parts = name.split(';')
    for p in reversed(parts):
        p = p.strip()
        if 'g__' in p and not p.endswith('g__'):
            return p.replace('g__', '').strip()
        if 'f__' in p and not p.endswith('f__'):
            return p.replace('f__', '').strip() + ' (family)'
        if 'd__Archaea' in p:
            return 'Archaea (domain)'
    return name.split(';')[0].strip()

name_map = {i: clean_name(i) for i in lfc_sig.index}
lfc_sig.index = [name_map[i] for i in lfc_sig.index]
se_sig = se_sig.loc[sig]
se_sig.index = [clean_name(i) for i in se_sig.index]

order = lfc_sig.sort_values()
se_ordered = se_sig.reindex(order.index).fillna(0)
colors = ['#1D9E75' if v > 0 else '#E8593C' for v in order.values]

fig, ax = plt.subplots(figsize=(8, max(4, len(order) * 0.6 + 2)))

ax.barh(range(len(order)), order.values,
        color=colors, alpha=0.85, height=0.6)

ax.errorbar(order.values, range(len(order)),
            xerr=se_ordered.values * 1.96,
            fmt='none', color='black', capsize=3, linewidth=1)

ax.set_yticks(range(len(order)))
ax.set_yticklabels(order.index, fontsize=10)
ax.axvline(0, color='black', linewidth=0.8)
ax.set_xlabel('Log fold change (control vs case)', fontsize=11)
ax.set_title('Differentially abundant genera\nCeliac cases vs controls (ANCOM-BC, q < 0.05)',
             fontweight='bold', fontsize=12)

legend = [Patch(color='#1D9E75', label='Enriched in controls'),
          Patch(color='#E8593C', label='Enriched in cases')]
ax.legend(handles=legend, frameon=False, loc='lower right')

ax.text(0.02, 0.02, f'n = {len(order)} significant genera',
        transform=ax.transAxes, fontsize=9, color='gray')

plt.tight_layout()
plt.savefig('results/figures/figure4_differential_abundance.png',
            dpi=300, bbox_inches='tight')
plt.close()
print("Figure 4 saved")