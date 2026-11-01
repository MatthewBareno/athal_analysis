"""Generate publication figure replicating the Pharmacogenetic Biogeographic Grouping System layout:
Panel (a): Pairwise Principal Component scatter matrix (PC1 to PC4) forming lower triangle.
Panel (b): Multi-class OvR ROC-AUC vs. Principal Components (1 to 10) per ancestry group.
Bottom Right: Group legend.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

# 1001 Genomes ancestry groups and configuration
GROUP_CONFIG = {
    'relict': {'label': 'Relict (Ancestral)', 'color': '#dc2626'},
    'spain': {'label': 'Spain / Iberia', 'color': '#ea580c'},
    'western_europe': {'label': 'Western Europe', 'color': '#4338ca'},
    'germany': {'label': 'Germany', 'color': '#047857'},
    'central_europe': {'label': 'Central Europe', 'color': '#16a34a'},
    'italy_balkan_caucasus': {'label': 'Italy / Balkans / Caucasus', 'color': '#9333ea'},
    'south_sweden': {'label': 'South Sweden', 'color': '#2563eb'},
    'north_sweden': {'label': 'North Sweden', 'color': '#06b6d4'},
    'asia': {'label': 'Asia', 'color': '#db2777'},
    'admixed': {'label': 'Admixed', 'color': '#64748b'}
}

ORDERED_GROUPS = [
    'relict', 'spain', 'western_europe', 'germany', 'central_europe',
    'italy_balkan_caucasus', 'south_sweden', 'north_sweden', 'asia', 'admixed'
]


def compute_ovr_auc(df, max_pcs=10):
    ks = list(range(1, max_pcs + 1))
    auc_data = {g: [] for g in ORDERED_GROUPS}
    
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    for k in ks:
        pc_cols = [f'PC{i}' for i in range(1, k + 1)]
        X = df[pc_cols].values
        
        for g in ORDERED_GROUPS:
            y = (df['admixture_group'] == g).astype(int).values
            y_probs = np.zeros(len(y))
            for train_idx, test_idx in cv.split(X, y):
                clf = LogisticRegression(class_weight='balanced', max_iter=1000, random_state=42)
                clf.fit(X[train_idx], y[train_idx])
                y_probs[test_idx] = clf.predict_proba(X[test_idx])[:, 1]
            auc = roc_auc_score(y, y_probs)
            auc_data[g].append(auc)
            
    return pd.DataFrame(auc_data, index=ks)


def create_figure(df, df_auc):
    os.makedirs('results/popgen', exist_ok=True)
    
    # Match the publication dimensions from the reference paper
    fig = plt.figure(figsize=(15.0, 8.8), dpi=300)
    
    # Outer layout: Left side (Panel a), Right side (Panel b + Legend)
    outer = gridspec.GridSpec(1, 2, width_ratios=[1.15, 1.0], wspace=0.24, left=0.06, right=0.96, top=0.93, bottom=0.08)
    
    # Left Grid: 3 rows, 3 cols for PC lower triangle
    gs_left = gridspec.GridSpecFromSubplotSpec(3, 3, subplot_spec=outer[0], hspace=0.28, wspace=0.28)
    
    # Right Grid: Top = AUC plot (Panel b), Bottom = Legend
    gs_right = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[1], height_ratios=[1.35, 0.65], hspace=0.32)
    
    # --- PANEL (a): PC Lower Triangular Matrix ---
    pc_pairs = [
        (0, 0, 'PC1', 'PC2'),
        (1, 0, 'PC1', 'PC3'),
        (1, 1, 'PC2', 'PC3'),
        (2, 0, 'PC1', 'PC4'),
        (2, 1, 'PC2', 'PC4'),
        (2, 2, 'PC3', 'PC4')
    ]
    
    axes_pc = {}
    for r, c, xpc, ypc in pc_pairs:
        ax = fig.add_subplot(gs_left[r, c])
        axes_pc[(r, c)] = ax
        
        # Plot points by group
        for grp in ORDERED_GROUPS:
            sub = df[df['admixture_group'] == grp]
            if sub.empty:
                continue
            ax.scatter(
                sub[xpc], sub[ypc],
                color=GROUP_CONFIG[grp]['color'],
                alpha=0.78,
                s=24,
                linewidth=0.35,
                edgecolors='#0f172a'
            )
            
        ax.grid(True, linestyle='-', alpha=0.45, color='#e2e8f0', linewidth=0.7)
        ax.set_axisbelow(True)
        ax.tick_params(axis='both', which='major', labelsize=7.0, pad=1.5, length=3, color='#94a3b8')
        
        # Border styling: light gray borders
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        for spine in ['left', 'bottom']:
            ax.spines[spine].set_color('#cbd5e1')
            ax.spines[spine].set_linewidth(0.8)
            
    # Add Row Labels (PC2, PC3, PC4) to the left of the leftmost plots
    axes_pc[(0, 0)].set_ylabel('PC2', fontsize=11, fontweight='bold', labelpad=8)
    axes_pc[(1, 0)].set_ylabel('PC3', fontsize=11, fontweight='bold', labelpad=8)
    axes_pc[(2, 0)].set_ylabel('PC4', fontsize=11, fontweight='bold', labelpad=8)
    
    # Add Column Labels (PC1, PC2, PC3) at the bottom of the bottom row plots
    axes_pc[(2, 0)].set_xlabel('PC1', fontsize=11, fontweight='bold', labelpad=8)
    axes_pc[(2, 1)].set_xlabel('PC2', fontsize=11, fontweight='bold', labelpad=8)
    axes_pc[(2, 2)].set_xlabel('PC3', fontsize=11, fontweight='bold', labelpad=8)
    
    # Add (a) label above top-left subplot
    axes_pc[(0, 0)].text(-0.35, 1.22, '(a)', transform=axes_pc[(0, 0)].transAxes,
                          fontsize=14, fontweight='bold', va='top', ha='left')
    
    # --- PANEL (b): Multi-PC AUC Curve ---
    ax_b = fig.add_subplot(gs_right[0])
    
    ks = df_auc.index.values
    for grp in ORDERED_GROUPS:
        cfg = GROUP_CONFIG[grp]
        ax_b.plot(
            ks, df_auc[grp].values,
            color=cfg['color'],
            marker='o',
            markersize=7.0,
            linewidth=2.0,
            alpha=0.92,
            label=f"{cfg['label']}"
        )
        
    ax_b.set_xlabel('Principal Components', fontsize=11, fontweight='bold', labelpad=7)
    ax_b.set_ylabel('AUC', fontsize=11, fontweight='bold', labelpad=7)
    ax_b.set_xticks(ks)
    ax_b.set_xlim(0.65, 10.35)
    ax_b.set_ylim(0.38, 1.04)
    ax_b.set_yticks([0.4, 0.6, 0.8, 1.0])
    ax_b.grid(True, linestyle='-', alpha=0.45, color='#e2e8f0', linewidth=0.7)
    ax_b.set_axisbelow(True)
    
    # Clean spines
    ax_b.spines['top'].set_visible(False)
    ax_b.spines['right'].set_visible(False)
    ax_b.spines['left'].set_color('#cbd5e1')
    ax_b.spines['bottom'].set_color('#cbd5e1')
    ax_b.tick_params(axis='both', which='major', labelsize=9.5, length=4, color='#94a3b8')
    
    # Add (b) label above AUC subplot
    ax_b.text(-0.10, 1.15, '(b)', transform=ax_b.transAxes,
              fontsize=14, fontweight='bold', va='top', ha='left')
              
    # --- PANEL (b) Bottom: Legend (Clean, unbordered matching reference) ---
    ax_leg = fig.add_subplot(gs_right[1])
    ax_leg.axis('off')
    
    handles, labels = ax_b.get_legend_handles_labels()
    
    ax_leg.legend(
        handles, labels,
        loc='center left',
        bbox_to_anchor=(0.12, 0.5),
        ncol=1,
        frameon=False,
        fontsize=9.5,
        handletextpad=0.8,
        labelspacing=0.4,
        markerscale=1.2
    )
    
    out_png = 'results/popgen/pca_auc_biogeographic_system.png'
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"Generated exact replica figure: {out_png}")
    return out_png


def main():
    print("Loading data...")
    df = pd.read_csv('data/processed/athal_1003_pca_kmeans.csv')
    print("Computing One-vs-Rest AUC across PCs 1..10...")
    df_auc = compute_ovr_auc(df, max_pcs=10)
    df_auc.to_csv('results/popgen/pca_ovr_auc_results.csv')
    print("Generating publication plot matching requested layout...")
    create_figure(df, df_auc)
    print("Done!")


if __name__ == '__main__':
    main()
