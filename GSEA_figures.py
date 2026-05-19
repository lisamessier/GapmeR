import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np

# Load data
df = pd.read_csv('gsea_sat2_vs_mock_significant.csv')
df = df[df['NES'] < 0].copy()
df = df[df['p.adjust'] < 0.05]
df['-log10padj'] = -np.log10(df['p.adjust'])

# Cancer-relevant gene sets
selected_ids = [
    'HALLMARK_MYC_TARGETS_V1',
    'HALLMARK_E2F_TARGETS',
    'HALLMARK_G2M_CHECKPOINT',
    'REACTOME_M_PHASE',
    'REACTOME_DNA_REPLICATION',
    'REACTOME_CELL_CYCLE_CHECKPOINTS',
    'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
    'HALLMARK_MTORC1_SIGNALING',
    'HALLMARK_DNA_REPAIR',
    'REACTOME_TRANSLATION',
    'WONG_EMBRYONIC_STEM_CELL_CORE',
    'HALLMARK_MYC_TARGETS_V2',
    'REACTOME_TRANSCRIPTIONAL_REGULATION_BY_TP53',
    'REACTOME_CELLULAR_SENESCENCE',
    'WP_HALLMARK_OF_CANCER_SUSTAINING_PROLIFERATIVE_SIGNALING',
]

plot_df = df[df['Description'].isin(selected_ids)].copy()
plot_df = plot_df.sort_values('NES', ascending=True)

# Clean labels
def clean(s):
    s = s.replace('HALLMARK_', '').replace('REACTOME_', '').replace('WP_', '').replace('WONG_', '')
    return s.replace('_', ' ').title()

plot_df['label'] = plot_df['Description'].apply(clean)

# Settings
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans'],
    'axes.spines.right': False,
    'axes.spines.top': False,
})

cmap = mcolors.LinearSegmentedColormap.from_list('blue_pink', ['#0D1B8E', '#FF1493'])
norm = mcolors.Normalize(vmin=plot_df['-log10padj'].min(), vmax=plot_df['-log10padj'].max())

sizes = 40 + 280 * (plot_df['setSize'] - plot_df['setSize'].min()) / \
        (plot_df['setSize'].max() - plot_df['setSize'].min() + 1)

# Make figure
fig, ax = plt.subplots(figsize=(9, 7), constrained_layout=True)

sc = ax.scatter(plot_df['NES'], range(len(plot_df)),
                c=plot_df['-log10padj'], cmap=cmap, norm=norm,
                s=sizes, alpha=0.9, zorder=3, edgecolors='white', linewidths=0.4)

ax.set_yticks(range(len(plot_df)))
ax.set_yticklabels(plot_df['label'], fontsize=9)
ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=10)
ax.grid(axis='x', linestyle='--', linewidth=0.4, alpha=0.5, zorder=0)
ax.tick_params(labelsize=9)

# Colorbar
cb = plt.colorbar(sc, ax=ax, shrink=0.4, pad=0.02, aspect=15)
cb.set_label('–log₁₀(adj. p-value)', fontsize=8)
cb.ax.tick_params(labelsize=7)

# Intersection size legend
for sz, lbl in [(40, 'small'), (160, 'med'), (320, 'large')]:
    ax.scatter([], [], s=sz, color='grey', alpha=0.6, label=lbl)
ax.legend(title='Set size', title_fontsize=7, fontsize=7,
          loc='lower right', framealpha=0.6)

# Panel label
ax.text(0.5, 1.02, 'HSAT2 KD', transform=ax.transAxes,
        fontsize=12, fontweight='bold', va='bottom', ha='center')

# Save
fig.savefig('gsea_bubbleplot.jpg', bbox_inches='tight', dpi=300)
fig.savefig('gsea_bubbleplot.svg', bbox_inches='tight', dpi=300)
print('Saved gsea_bubbleplot.jpg and gsea_bubbleplot.svg')