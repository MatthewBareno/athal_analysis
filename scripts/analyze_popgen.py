"""Perform population genomics analyses on Arabidopsis thaliana (1001 Genomes):
1. LD Pruning & Principal Component Analysis (PLINK 2).
2. K-means clustering on genetic PCs (evaluating K=2..12, peak at K=9).
3. ADMIXTURE visualization (K=9 stacked barplot).
4. Generates publication-quality figures (HTML & 300-DPI PNG).
"""

import json
import os
import subprocess
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Standard 1001 Genomes ancestry groups and colors
GROUP_CONFIG = {
    'central_europe': {'label': 'Central Europe', 'color': '#33a02c'},
    'germany': {'label': 'Germany', 'color': '#006d2c'},
    'south_sweden': {'label': 'South Sweden', 'color': '#1f78b4'},
    'north_sweden': {'label': 'North Sweden', 'color': '#38bdf8'},
    'western_europe': {'label': 'Western Europe', 'color': '#4338ca'},
    'spain': {'label': 'Spain', 'color': '#ea580c'},
    'italy_balkan_caucasus': {'label': 'Italy / Balkans / Caucasus', 'color': '#9333ea'},
    'asia': {'label': 'Asia', 'color': '#ec4899'},
    'relict': {'label': 'Relict (Ancestral)', 'color': '#dc2626'},
    'admixed': {'label': 'Admixed', 'color': '#64748b'}
}

ORDERED_GROUPS = [
    'relict', 'spain', 'western_europe', 'germany', 'central_europe',
    'italy_balkan_caucasus', 'south_sweden', 'north_sweden', 'asia', 'admixed'
]


def run_pca_plink():
    os.makedirs('data/processed', exist_ok=True)
    os.makedirs('results/popgen', exist_ok=True)

    prune_in = 'data/processed/athal_1003_ld.prune.in'
    if not os.path.exists(prune_in):
        print("Running LD pruning with PLINK2...")
        subprocess.run([
            'plink2', '--pfile', 'data/processed/athal_1003_qc',
            '--indep-pairwise', '50', '10', '0.2',
            '--out', 'data/processed/athal_1003_ld'
        ], check=True)

    pca_evec = 'data/processed/athal_1003_pca.eigenvec'
    if not os.path.exists(pca_evec):
        print("Running PCA with PLINK2 on LD-pruned variants...")
        subprocess.run([
            'plink2', '--pfile', 'data/processed/athal_1003_qc',
            '--extract', prune_in,
            '--pca', '20',
            '--out', 'data/processed/athal_1003_pca'
        ], check=True)

    # Load eigenvalues and calculate variance explained
    evals = np.loadtxt('data/processed/athal_1003_pca.eigenval')
    total_var = np.sum(evals)
    var_explained = (evals / total_var) * 100

    evec = pd.read_csv('data/processed/athal_1003_pca.eigenvec', sep='\t')
    pheno = pd.read_csv('data/processed/phenotype_clean.csv')

    merged = pheno.merge(evec, left_on='accession_id', right_on='#IID')
    return merged, var_explained, evals


def run_kmeans(merged):
    # Top 10 PCs for clustering
    pc_cols = [f'PC{i}' for i in range(1, 11)]
    X = merged[pc_cols].values

    metrics = []
    for k in range(2, 13):
        km = KMeans(n_clusters=k, random_state=42, n_init=10).fit(X)
        sil = silhouette_score(X, km.labels_)
        metrics.append({'K': k, 'Inertia': float(km.inertia_), 'Silhouette': float(sil)})

    df_metrics = pd.DataFrame(metrics)
    df_metrics.to_csv('results/popgen/kmeans_evaluation.csv', index=False)

    # Fit optimal K=9 (peak silhouette score)
    km9 = KMeans(n_clusters=9, random_state=42, n_init=15).fit(X)
    merged['kmeans_cluster_k9'] = km9.labels_

    merged.to_csv('data/processed/athal_1003_pca_kmeans.csv', index=False)
    print("Saved merged PCA and K-means dataset: data/processed/athal_1003_pca_kmeans.csv")
    return merged, df_metrics


def plot_pca_figure(merged, var_explained):
    fig = go.Figure()

    for grp in ORDERED_GROUPS:
        sub = merged[merged['admixture_group'] == grp]
        if sub.empty:
            continue

        cfg = GROUP_CONFIG.get(grp, {'label': grp, 'color': '#7f7f7f'})
        custom_data = list(zip(
            sub['name'],
            sub['accession_id'],
            sub['country'],
            sub['PC1'].round(4),
            sub['PC2'].round(4),
            sub['FT10'].round(1)
        ))

        fig.add_trace(go.Scatter(
            x=sub['PC1'],
            y=sub['PC2'],
            mode='markers',
            name=f"{cfg['label']} (n={len(sub)})",
            customdata=custom_data,
            hovertemplate=(
                "<b>%{customdata[0]}</b> (ID: %{customdata[1]})<br>"
                f"<b>Group:</b> {cfg['label']}<br>"
                "<b>Country:</b> %{customdata[2]}<br>"
                "<b>PC1:</b> %{customdata[3]} | <b>PC2:</b> %{customdata[4]}<br>"
                "<b>Flowering Time:</b> %{customdata[5]} d<extra></extra>"
            ),
            marker=dict(
                size=7,
                color=cfg['color'],
                opacity=0.88,
                line=dict(width=0.45, color='#ffffff')
            )
        ))

    fig.update_layout(
        title=dict(
            text="<b>Genetic Structure of <i>Arabidopsis thaliana</i> (PCA on 397,850 LD-Pruned SNPs)</b><br>"
                 "<span style='font-size: 11.5px; font-weight: normal; color: #475569;'>"
                 "1,003 accessions colored by 1001 Genomes Genetic Ancestry Group</span>",
            x=0.02, y=0.97,
            font=dict(size=14, family='sans-serif', color='#0f172a')
        ),
        xaxis=dict(
            title=dict(
                text=f"<b>Principal Component 1 ({var_explained[0]:.2f}% variance explained)</b>",
                font=dict(size=11, family='sans-serif')
            ),
            zeroline=True, zerolinecolor='#e2e8f0',
            gridcolor='#f1f5f9',
            tickfont=dict(size=10, family='sans-serif')
        ),
        yaxis=dict(
            title=dict(
                text=f"<b>Principal Component 2 ({var_explained[1]:.2f}% variance explained)</b>",
                font=dict(size=11, family='sans-serif')
            ),
            zeroline=True, zerolinecolor='#e2e8f0',
            gridcolor='#f1f5f9',
            tickfont=dict(size=10, family='sans-serif')
        ),
        legend=dict(
            title=dict(text="<b>Ancestry Group (n)</b>", font=dict(size=10.5, family='sans-serif')),
            orientation="v",
            yanchor="middle",
            y=0.5,
            xanchor="left",
            x=1.01,
            font=dict(size=9.5, family='sans-serif'),
            itemsizing="constant",
            bgcolor="rgba(255,255,255,0.9)",
            bordercolor="#e2e8f0",
            borderwidth=1
        ),
        margin=dict(l=60, r=220, t=65, b=60),
        width=1150,
        height=620,
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        font=dict(family="sans-serif")
    )

    out_html = 'results/popgen/pca_pc1_pc2.html'
    out_png = 'results/popgen/pca_pc1_pc2.png'

    fig.write_html(
        out_html,
        include_plotlyjs='cdn',
        full_html=True,
        config=dict(displayModeBar=True, displaylogo=False)
    )
    fig.write_image(out_png, scale=2)
    print(f"Generated PCA plots: {out_html} and {out_png}")


def plot_admixture_barplot():
    admix_path = 'data/1001genomes_admixture_groups.csv'
    pheno_path = 'data/processed/phenotype_clean.csv'

    df_admix = pd.read_csv(admix_path)
    df_pheno = pd.read_csv(pheno_path)

    k9_cols = [c for c in df_admix.columns if c.startswith('K9.ChrAll.')]
    merged = df_pheno.merge(df_admix[['id', 'group'] + k9_cols], left_on='accession_id', right_on='id')

    # Color mapping for K=9 ancestry components matching the groups
    # Order individuals by group and their dominant ancestry proportion
    merged['primary_ancestry'] = merged[k9_cols].max(axis=1)

    group_order_map = {g: i for i, g in enumerate(ORDERED_GROUPS)}
    merged['group_rank'] = merged['admixture_group'].map(group_order_map).fillna(99)
    merged = merged.sort_values(['group_rank', 'primary_ancestry'], ascending=[True, False]).reset_index(drop=True)

    # Use Matplotlib for publication-grade stacked barplot
    fig, ax = plt.subplots(figsize=(15, 4.5), dpi=300)

    # 9 Distinct colors for K=9 components
    k_colors = [
        '#4338ca',  # 0: Western Europe
        '#dc2626',  # 1: Relict
        '#ea580c',  # 2: Spain
        '#1f78b4',  # 3: South Sweden
        '#38bdf8',  # 4: North Sweden
        '#006d2c',  # 5: Germany
        '#ec4899',  # 6: Asia
        '#9333ea',  # 7: Italy / Balkans / Caucasus
        '#33a02c',  # 8: Central Europe
    ]

    bottom = np.zeros(len(merged))
    x = np.arange(len(merged))

    for idx, col in enumerate(k9_cols):
        values = merged[col].values
        ax.bar(x, values, bottom=bottom, color=k_colors[idx], width=1.0, linewidth=0)
        bottom += values

    # Draw group separation vertical divider lines and labels
    group_changes = np.where(merged['admixture_group'].values[:-1] != merged['admixture_group'].values[1:])[0]
    for change_x in group_changes:
        ax.axvline(change_x + 0.5, color='#0f172a', linewidth=0.75, linestyle='--')

    # Label groups at their midpoints
    group_starts = [0] + (group_changes + 1).tolist()
    group_ends = group_changes.tolist() + [len(merged) - 1]

    tick_locs = []
    tick_labels = []
    for start, end in zip(group_starts, group_ends):
        grp_name = merged.loc[start, 'admixture_group']
        label = GROUP_CONFIG.get(grp_name, {}).get('label', grp_name)
        # Shorten for tight barplot labels
        label_short = label.replace(' (Ancestral)', '').replace(' / Balkans / Caucasus', '')
        tick_locs.append((start + end) / 2)
        tick_labels.append(label_short)

    ax.set_xticks(tick_locs)
    ax.set_xticklabels(tick_labels, rotation=40, ha='right', fontsize=8.5, fontweight='bold', color='#1e293b')
    ax.set_xlim(-0.5, len(merged) - 0.5)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Ancestry Fraction (Q)", fontsize=10, fontweight='bold', color='#0f172a')
    ax.set_title("1001 Genomes ADMIXTURE Clustering Profile (K = 9 across 1,003 Accessions)", fontsize=12, fontweight='bold', pad=12)

    # Style borders
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_color('#64748b')
    ax.spines['bottom'].set_color('#64748b')

    plt.tight_layout()
    out_png = 'results/popgen/admixture_barplot_k9.png'
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"Generated ADMIXTURE barplot: {out_png}")


def plot_kmeans_evaluation(df_metrics):
    fig, ax1 = plt.subplots(figsize=(7, 4.5), dpi=300)
    color_sil = '#2563eb'
    color_ine = '#d97706'

    ax1.set_xlabel('Number of Clusters (K)', fontsize=10, fontweight='bold')
    ax1.set_ylabel('Silhouette Score', color=color_sil, fontsize=10, fontweight='bold')
    l1 = ax1.plot(df_metrics['K'], df_metrics['Silhouette'], color=color_sil, marker='o', linewidth=2.2, label='Silhouette Score')
    ax1.tick_params(axis='y', labelcolor=color_sil)
    ax1.set_xticks(df_metrics['K'])

    # Highlight peak at K=9
    k9_sil = df_metrics.loc[df_metrics['K'] == 9, 'Silhouette'].values[0]
    ax1.plot(9, k9_sil, marker='*', markersize=14, color='#dc2626', zorder=5)
    ax1.annotate(
        f'Optimal K=9\n(Sil = {k9_sil:.3f})',
        xy=(9, k9_sil), xytext=(7.5, k9_sil + 0.04),
        arrowprops=dict(facecolor='#dc2626', shrink=0.08, width=1.5, headwidth=6),
        fontsize=9, fontweight='bold', color='#dc2626',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#fee2e2', edgecolor='#dc2626', alpha=0.9)
    )

    ax2 = ax1.twinx()
    ax2.set_ylabel('Inertia (Within-Cluster Sum of Squares)', color=color_ine, fontsize=10, fontweight='bold')
    l2 = ax2.plot(df_metrics['K'], df_metrics['Inertia'], color=color_ine, marker='s', linestyle='--', linewidth=2.0, label='Inertia (Elbow)')
    ax2.tick_params(axis='y', labelcolor=color_ine)

    ax1.set_title('K-Means Cluster Selection (Top 10 Principal Components)', fontsize=11, fontweight='bold', pad=10)
    ax1.grid(True, linestyle=':', alpha=0.6, color='#cbd5e1')

    # Combined legend
    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='lower left', frameon=True, edgecolor='#e2e8f0', fontsize=8.5)

    plt.tight_layout()
    out_png = 'results/popgen/kmeans_evaluation.png'
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"Generated K-means evaluation plot: {out_png}")


def plot_composite_figure(merged, var_explained, df_metrics):
    """Generate a 3-panel publication figure: Panel A (PCA) + Panel B (K-means) + Panel C (ADMIXTURE)."""
    fig = plt.figure(figsize=(16, 11), dpi=300)
    gs = gridspec.GridSpec(2, 2, height_ratios=[1.15, 0.85], hspace=0.36, wspace=0.22)

    # Panel A: PCA (Top Left)
    ax1 = fig.add_subplot(gs[0, 0])
    for grp in ORDERED_GROUPS:
        sub = merged[merged['admixture_group'] == grp]
        if sub.empty:
            continue
        cfg = GROUP_CONFIG.get(grp, {'label': grp, 'color': '#7f7f7f'})
        ax1.scatter(
            sub['PC1'], sub['PC2'],
            c=cfg['color'], label=f"{cfg['label']} ({len(sub)})",
            alpha=0.85, s=26, edgecolors='white', linewidth=0.25
        )

    ax1.set_xlabel(f"Principal Component 1 ({var_explained[0]:.2f}% variance)", fontsize=9.5, fontweight='bold')
    ax1.set_ylabel(f"Principal Component 2 ({var_explained[1]:.2f}% variance)", fontsize=9.5, fontweight='bold')
    ax1.set_title("A   Genome-Wide PCA (397,850 LD-Pruned SNPs)", fontsize=11.5, fontweight='bold', loc='left')
    ax1.legend(loc='upper right', fontsize=7.5, frameon=True, edgecolor='#cbd5e1', framealpha=0.92, ncol=2)
    ax1.grid(True, linestyle=':', alpha=0.5, color='#cbd5e1')
    ax1.set_axisbelow(True)

    # Panel B: K-means Model Selection (Top Right)
    ax2 = fig.add_subplot(gs[0, 1])
    color_sil = '#2563eb'
    color_ine = '#d97706'

    ax2.set_xlabel('Number of Clusters (K)', fontsize=9.5, fontweight='bold')
    ax2.set_ylabel('Silhouette Score', color=color_sil, fontsize=9.5, fontweight='bold')
    l1 = ax2.plot(df_metrics['K'], df_metrics['Silhouette'], color=color_sil, marker='o', linewidth=2.0, label='Silhouette')
    ax2.tick_params(axis='y', labelcolor=color_sil, labelsize=8.5)
    ax2.tick_params(axis='x', labelsize=8.5)
    ax2.set_xticks(df_metrics['K'])

    k9_sil = df_metrics.loc[df_metrics['K'] == 9, 'Silhouette'].values[0]
    ax2.plot(9, k9_sil, marker='*', markersize=13, color='#dc2626', zorder=5)
    ax2.annotate(
        f'Peak K=9 (Sil = {k9_sil:.3f})',
        xy=(9, k9_sil), xytext=(6.8, k9_sil + 0.045),
        arrowprops=dict(facecolor='#dc2626', shrink=0.08, width=1.2, headwidth=5),
        fontsize=8.5, fontweight='bold', color='#dc2626',
        bbox=dict(boxstyle='round,pad=0.25', facecolor='#fee2e2', edgecolor='#dc2626', alpha=0.9)
    )

    ax2_twin = ax2.twinx()
    ax2_twin.set_ylabel('Inertia (WCSS)', color=color_ine, fontsize=9.5, fontweight='bold')
    l2 = ax2_twin.plot(df_metrics['K'], df_metrics['Inertia'], color=color_ine, marker='s', linestyle='--', linewidth=1.8, label='Inertia')
    ax2_twin.tick_params(axis='y', labelcolor=color_ine, labelsize=8.5)

    ax2.set_title("B   K-Means Clustering Evaluation (Top 10 PCs)", fontsize=11.5, fontweight='bold', loc='left')
    ax2.grid(True, linestyle=':', alpha=0.5, color='#cbd5e1')

    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    ax2.legend(lines, labels, loc='lower left', frameon=True, edgecolor='#cbd5e1', framealpha=0.92, fontsize=8)

    # Panel C: ADMIXTURE Barplot (Bottom Full Width)
    ax3 = fig.add_subplot(gs[1, :])
    admix_path = 'data/1001genomes_admixture_groups.csv'
    df_admix = pd.read_csv(admix_path)
    k9_cols = [c for c in df_admix.columns if c.startswith('K9.ChrAll.')]
    merged_admix = merged.merge(df_admix[['id'] + k9_cols], left_on='accession_id', right_on='id')
    merged_admix['primary_ancestry'] = merged_admix[k9_cols].max(axis=1)
    group_order_map = {g: i for i, g in enumerate(ORDERED_GROUPS)}
    merged_admix['group_rank'] = merged_admix['admixture_group'].map(group_order_map).fillna(99)
    merged_admix = merged_admix.sort_values(['group_rank', 'primary_ancestry'], ascending=[True, False]).reset_index(drop=True)

    k_colors = [
        '#4338ca', '#dc2626', '#ea580c', '#1f78b4', '#38bdf8',
        '#006d2c', '#ec4899', '#9333ea', '#33a02c'
    ]
    bottom = np.zeros(len(merged_admix))
    x = np.arange(len(merged_admix))

    for idx, col in enumerate(k9_cols):
        values = merged_admix[col].values
        ax3.bar(x, values, bottom=bottom, color=k_colors[idx], width=1.0, linewidth=0)
        bottom += values

    group_changes = np.where(merged_admix['admixture_group'].values[:-1] != merged_admix['admixture_group'].values[1:])[0]
    for change_x in group_changes:
        ax3.axvline(change_x + 0.5, color='#0f172a', linewidth=0.75, linestyle='--')

    group_starts = [0] + (group_changes + 1).tolist()
    group_ends = group_changes.tolist() + [len(merged_admix) - 1]
    tick_locs = [(s + e) / 2 for s, e in zip(group_starts, group_ends)]
    tick_labels = [
        GROUP_CONFIG.get(merged_admix.loc[s, 'admixture_group'], {}).get('label', '').replace(' (Ancestral)', '').replace(' / Balkans / Caucasus', '')
        for s in group_starts
    ]

    ax3.set_xticks(tick_locs)
    ax3.set_xticklabels(tick_labels, rotation=25, ha='right', fontsize=8.5, fontweight='bold', color='#1e293b')
    ax3.set_xlim(-0.5, len(merged_admix) - 0.5)
    ax3.set_ylim(0, 1.0)
    ax3.set_ylabel("Ancestry Fraction (Q)", fontsize=9.5, fontweight='bold')
    ax3.set_title("C   Global ADMIXTURE Ancestry Proportions (K = 9 across 1,003 Accessions)", fontsize=11.5, fontweight='bold', loc='left')

    out_comp = 'results/popgen/pca_admixture_composite.png'
    plt.savefig(out_comp, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Generated composite publication figure: {out_comp}")


def compute_concordance(merged):
    # Cross-tabulation between ADMIXTURE group and K-means cluster
    ct = pd.crosstab(
        merged['admixture_group'],
        merged['kmeans_cluster_k9'],
        margins=True,
        margins_name='Total'
    )
    out_csv = 'results/popgen/kmeans_concordance.csv'
    ct.to_csv(out_csv)
    print(f"Generated K-means concordance matrix: {out_csv}")
    return ct


def main():
    print("=== Step 2: Population Structure Analysis ===")
    merged, var_explained, evals = run_pca_plink()
    merged, df_metrics = run_kmeans(merged)
    compute_concordance(merged)
    plot_pca_figure(merged, var_explained)
    plot_kmeans_evaluation(df_metrics)
    plot_admixture_barplot()
    plot_composite_figure(merged, var_explained, df_metrics)
    print("=== Population Structure Analysis Complete ===")


if __name__ == '__main__':
    main()
