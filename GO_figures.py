import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np

# Load data
hsat2 = pd.read_excel('2026_hsat2_downregulated_pathways_log2fc_0.4.xlsx')
hsat2['gene_ratio'] = hsat2['intersection_size'] / hsat2['query_size']

# Remove 'cytoplasm' before selecting top terms b/c it warps the plot
hsat2_filtered = hsat2[hsat2['term_name'] != 'cytoplasm']

# Settings
SOURCE_LABELS = {
    'GO:BP': 'GO:BP', 'GO:CC': 'GO:CC', 'GO:MF': 'GO:MF',
    'REAC': 'Reactome', 'KEGG': 'KEGG', 'WP': 'WikiPathways',
    'TF': 'TF', 'MIRNA': 'miRNA', 'CORUM': 'CORUM',
}

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans'],
    'axes.spines.right': False,
    'axes.spines.top': False,
})

cmap_hsat2 = mcolors.LinearSegmentedColormap.from_list('blue_pink', ['#0D1B8E', '#FF1493'])

# Bubble plot function
def make_bubbleplot(ax, df, cmap, top_n=15):
    df = df.nsmallest(top_n, 'adjusted_p_value').copy()
    df = df.sort_values('negative_log10_of_adjusted_p_value', ascending=True)
    df['lbl'] = df['term_name'].str[:55]

    norm = mcolors.Normalize(vmin=df['adjusted_p_value'].min(),
                              vmax=df['adjusted_p_value'].max())

    sizes = 40 + 280 * (df['intersection_size'] - df['intersection_size'].min()) / \
            (df['intersection_size'].max() - df['intersection_size'].min() + 1)

    sc = ax.scatter(df['gene_ratio'], range(len(df)),
                    c=df['adjusted_p_value'], cmap=cmap, norm=norm,
                    s=sizes, alpha=0.9, zorder=3, edgecolors='white', linewidths=0.4)

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df['lbl'], fontsize=8)
    ax.set_xlabel('Gene ratio', fontsize=9)
    ax.grid(axis='x', linestyle='--', linewidth=0.4, alpha=0.5, zorder=0)
    ax.tick_params(axis='both', labelsize=8)

    # Source database labels on right axis
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks(range(len(df)))
    ax2.set_yticklabels(
        [SOURCE_LABELS.get(s, s) for s in df['source']],
        fontsize=6.5, color='#555555'
    )
    ax2.spines['right'].set_visible(True)
    ax2.spines['top'].set_visible(False)
    ax2.tick_params(right=False)

    return sc, norm

# Build figure
fig, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)

sc, norm = make_bubbleplot(ax, hsat2_filtered, cmap_hsat2)

# Colorbar
sm = cm.ScalarMappable(norm=norm, cmap=cmap_hsat2)
sm.set_array([])
cb = plt.colorbar(sm, ax=ax, shrink=0.4, pad=0.18, aspect=15)
cb.set_label('Adjusted p-value', fontsize=8)
cb.ax.tick_params(labelsize=7)

# Intersection size legend
for sz, lbl in [(40, 'low'), (160, 'mid'), (320, 'high')]:
    ax.scatter([], [], s=sz, color='grey', alpha=0.6, label=lbl)
ax.legend(title='Intersection\nsize', title_fontsize=7, fontsize=7,
          loc='upper right', framealpha=0.6)

# Panel label
ax.text(0.5, 1.02, 'HSAT2 KD', transform=ax.transAxes,
        fontsize=12, fontweight='bold', va='bottom', ha='center')

# Save
fig.savefig('hsat2_bubbleplot.svg', bbox_inches='tight', dpi=300)
fig.savefig('hsat2_bubbleplot.jpg', bbox_inches='tight', dpi=300)
print('Saved hsat2_bubbleplot.svg and hsat2_bubbleplot.jpg')