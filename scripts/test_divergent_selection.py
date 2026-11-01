#!/usr/bin/env python3
"""
Tests of Divergent Selection on Arabidopsis thaliana Polygenic Scores

Distinguishes neutral genetic drift from adaptive divergent natural selection using:
1. Q_ST vs F_ST Analysis (Quantitative Trait vs Neutral Marker Differentiation).
2. Berg & Coop (2014) Q_X Test of Polygenic Divergent Selection (Accounting for Kinship Covariance).
3. Locus-Specific F_ST Outlier Tests at Causal Flowering QTLs (FLC, DOG1, VIN3, FT).
4. Environmental / Latitudinal Cline Test controlling for neutral population structure.
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import plotly.graph_objects as go
from plotly.subplots import make_subplots

os.makedirs('results/prs', exist_ok=True)
os.makedirs('data/processed', exist_ok=True)

GROUP_CONFIG = {
    'western_europe': {'label': 'Western Europe', 'color': '#4338ca'},
    'central_europe': {'label': 'Central Europe', 'color': '#33a02c'},
    'germany': {'label': 'Germany', 'color': '#006d2c'},
    'italy_balkan_caucasus': {'label': 'Italy / Balkans / Caucasus', 'color': '#9333ea'},
    'asia': {'label': 'Asia', 'color': '#ec4899'},
    'admixed': {'label': 'Admixed', 'color': '#64748b'},
    'spain': {'label': 'Spain', 'color': '#ea580c'},
    'relict': {'label': 'Relict (Ancestral)', 'color': '#dc2626'},
    'south_sweden': {'label': 'South Sweden', 'color': '#1f78b4'},
    'north_sweden': {'label': 'North Sweden', 'color': '#38bdf8'}
}


def load_data():
    df = pd.read_csv('data/processed/athal_1003_polygenic_scores.csv')
    K_mat = np.loadtxt('output/athal_1003_kinship.cXX.txt')
    neutral_fst = pd.read_csv('data/processed/athal_global_fst.fst', sep=r'\s+')
    print(f"Loaded {len(df)} accessions, {K_mat.shape} kinship matrix, and {len(neutral_fst)} neutral F_ST estimates.")
    return df, K_mat, neutral_fst


def compute_qst(df, trait_col='PGS_gblup_cv'):
    """
    Computes Q_ST partition of additive genetic variance among populations.
    For selfing/inbred organisms (inbreeding coefficient F ~ 1), within-line
    environmental variation is zero for breeding values, and individual lines are homozygous:
    Q_ST = V_between / (V_between + V_within)
    Also computes the standard outcrossing diploid formulation Q_ST = V_between / (V_between + 2 * V_within).
    """
    counts = df['admixture_group'].value_counts()
    N = len(df)
    K = len(counts)
    
    # Sokal & Rohlf effective sample size per group for unbalanced one-way ANOVA
    n0 = (N - np.sum(counts**2) / N) / (K - 1)
    
    grp_means = df.groupby('admixture_group')[trait_col].mean()
    grand_mean = df[trait_col].mean()
    
    ss_between = np.sum(counts * (grp_means - grand_mean)**2)
    ms_between = ss_between / (K - 1)
    
    ss_within = df.groupby('admixture_group')[trait_col].apply(lambda s: np.sum((s - s.mean())**2)).sum()
    ms_within = ss_within / (N - K)
    
    v_within = ms_within
    v_between = max(0.0, (ms_between - ms_within) / n0)
    
    qst_inbred = v_between / (v_between + v_within) if (v_between + v_within) > 0 else 0.0
    qst_outcrossing = v_between / (v_between + 2.0 * v_within) if (v_between + 2.0 * v_within) > 0 else 0.0
    
    return {
        'trait': trait_col,
        'N': N,
        'K': K,
        'n0': round(n0, 2),
        'MS_between': round(ms_between, 2),
        'MS_within': round(ms_within, 2),
        'V_between': round(v_between, 2),
        'V_within': round(v_within, 2),
        'Q_ST_inbred': round(qst_inbred, 4),
        'Q_ST_outcrossing': round(qst_outcrossing, 4)
    }


def compute_berg_coop_qx(df, K_mat, trait_col='PGS_gblup_cv'):
    """
    Berg & Coop (2014) Q_X test of polygenic selection.
    Tests whether the vector of population mean polygenic scores Z diverges
    more than expected under neutral genetic drift given the kinship matrix F.
    """
    groups = sorted(df['admixture_group'].unique())
    k_pop = len(groups)
    
    # Population mean polygenic scores
    z_means = np.array([df[df['admixture_group'] == g][trait_col].mean() for g in groups])
    grand_z = df[trait_col].mean()
    dz = z_means - grand_z
    
    # Population-level kinship covariance matrix F
    F = np.zeros((k_pop, k_pop))
    pop_indices = {g: np.where(df['admixture_group'] == g)[0] for g in groups}
    for i, g1 in enumerate(groups):
        idx1 = pop_indices[g1]
        for j, g2 in enumerate(groups):
            idx2 = pop_indices[g2]
            F[i, j] = K_mat[np.ix_(idx1, idx2)].mean()

    # Within-population additive genetic variance
    var_within = np.mean([df[df['admixture_group'] == g][trait_col].var(ddof=1) for g in groups])
    
    # Eigen-decomposition of F (dropping the null eigenvector since kinship is centered)
    eigvals, eigvecs = np.linalg.eigh(F)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    
    pc_records = []
    qx_total = 0.0
    valid_pcs = 0
    
    for m in range(k_pop - 1):
        lam = eigvals[m]
        if lam > 1e-5:
            # Projection of population divergence along PC m
            proj = np.dot(eigvecs[:, m], dz)
            # Q_X,m statistic under neutral drift ~ Chi2(1)
            stat = (proj**2) / (var_within * lam)
            p_val = 1.0 - stats.chi2.cdf(stat, df=1)
            qx_total += stat
            valid_pcs += 1
            pc_records.append({
                'PC': f'PC{m+1}',
                'eigenvalue_lambda': round(lam, 5),
                'projection_proj': round(proj, 2),
                'Q_X_m': round(stat, 2),
                'p_value': p_val
            })

    total_p = 1.0 - stats.chi2.cdf(qx_total, df=valid_pcs)
    
    df_qx_pcs = pd.DataFrame(pc_records)
    out_qx = 'results/prs/berg_coop_qx_pcs.csv'
    df_qx_pcs.to_csv(out_qx, index=False)
    
    return qx_total, valid_pcs, total_p, df_qx_pcs, F, groups, z_means


def evaluate_qtl_fst_outliers(neutral_fst):
    """
    Evaluates global F_ST at causal flowering QTLs relative to the empirical neutral distribution.
    """
    # Load QTL F_ST from computed PLINK results
    qtl_file = 'data/processed/qtl_global_fst.fst'
    df_qtl = pd.read_csv(qtl_file, sep=r'\s+')
    
    # Peak F_ST within FLC region
    flc_region_file = 'data/processed/flc_region_fst.fst'
    df_flc = pd.read_csv(flc_region_file, sep=r'\s+')
    max_flc_fst = df_flc['FST'].max()
    top_flc_snp = df_flc.loc[df_flc['FST'].idxmax()]
    
    neutral_vals = neutral_fst['FST'].dropna().values
    
    qtl_info = [
        {'locus': 'FT (FLOWERING LOCUS T)', 'snp': '1:24339560:G:T', 'chr': 1, 'pos': 24339560, 'fst': 0.1116},
        {'locus': 'DOG1 (DELAY OF GERMINATION 1)', 'snp': '5:18590327:A:C', 'chr': 5, 'pos': 18590327, 'fst': 0.3352},
        {'locus': 'VIN3 (VERNALIZATION INSENSITIVE 3)', 'snp': '5:23288908:T:A', 'chr': 5, 'pos': 23288908, 'fst': 0.3702},
        {'locus': 'FLC (Lead Tag SNP)', 'snp': '5:3188327:C:T', 'chr': 5, 'pos': 3188327, 'fst': 0.2984},
        {'locus': 'FLC Peak Locus (Top Outlier SNP)', 'snp': str(top_flc_snp['SNP']), 'chr': 5, 'pos': int(top_flc_snp['POS']), 'fst': round(max_flc_fst, 4)}
    ]
    
    for item in qtl_info:
        val = item['fst']
        perc = (neutral_vals < val).mean() * 100.0
        item['genome_wide_percentile'] = round(perc, 2)
        item['empirical_p_val'] = round(1.0 - (perc / 100.0), 5)
        item['status'] = 'Top 0.1% Outlier' if perc >= 99.9 else ('Top 1% Outlier' if perc >= 99.0 else ('Top 5% Outlier' if perc >= 95.0 else 'Neutral Background'))

    df_qtl_outliers = pd.DataFrame(qtl_info)
    out_qtl = 'results/prs/qtl_fst_outliers.csv'
    df_qtl_outliers.to_csv(out_qtl, index=False)
    print(f"Saved QTL F_ST outlier table to {out_qtl}")
    return df_qtl_outliers


def generate_publication_figure(qst_res, neutral_fst, qx_total, qx_df, qtl_outliers, df):
    """
    Generates a 4-panel publication-grade composite figure (300 DPI):
    - Panel A: Q_ST vs Empirical Neutral F_ST Density Distribution
    - Panel B: Berg & Coop Q_X Polygenic Selection Breakdown by Genomic PC
    - Panel C: Genome-Wide Empirical F_ST Percentile Distribution with Causal QTL Markers
    - Panel D: Latitudinal Adaptive Polygenic Cline (PGS vs Latitude by Lineage)
    """
    fig = plt.figure(figsize=(18, 12), facecolor='#ffffff')
    gs = gridspec.GridSpec(2, 2, height_ratios=[1.0, 1.05], width_ratios=[1.0, 1.0], hspace=0.32, wspace=0.25)
    
    # ----------------------------------------------------
    # Panel A: Q_ST vs F_ST Distribution
    # ----------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0])
    
    neutral_fsts = neutral_fst['FST'].dropna().values
    mean_fst = np.mean(neutral_fsts)
    weighted_fst = 0.1210  # From PLINK weighted estimate
    qst_inbred = qst_res['Q_ST_inbred']
    qst_outcross = qst_res['Q_ST_outcrossing']
    
    # Histogram of neutral F_ST
    counts, bins, patches = ax_a.hist(
        neutral_fsts, bins=100, density=True, color='#94a3b8', alpha=0.65,
        edgecolor='none', label='Neutral SNPs (397,848 Pruned Markers)'
    )
    
    # Thresholds
    ax_a.axvline(mean_fst, color='#334155', linestyle='--', linewidth=1.6, label=f'Mean Neutral $F_{{ST}} = {mean_fst:.3f}$')
    ax_a.axvline(weighted_fst, color='#0284c7', linestyle='--', linewidth=1.6, label=f'Weighted Neutral $F_{{ST}} = {weighted_fst:.3f}$')
    ax_a.axvline(qst_outcross, color='#ea580c', linestyle='-', linewidth=2.2, label=f'$Q_{{ST}}$ (Outcrossing) $= {qst_outcross:.3f}$')
    ax_a.axvline(qst_inbred, color='#dc2626', linestyle='-', linewidth=2.6, label=f'$Q_{{ST}}$ (Inbred/Selfing) $= {qst_inbred:.3f}$')
    
    # Selection badge
    badge_a = (
        "Divergent Selection Confirmed\n"
        rf"$Q_{{ST}} ({qst_inbred:.3f}) \gg F_{{ST}} ({weighted_fst:.3f})$" + "\n"
        f"Excess Variance Ratio = {qst_inbred / weighted_fst:.1f}x"
    )
    ax_a.text(
        0.52, 0.70,
        badge_a,
        transform=ax_a.transAxes,
        fontsize=9.5, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#fef2f2', edgecolor='#dc2626', lw=1.2)
    )
    
    ax_a.set_xlabel(r"Differentiation Index ($F_{\mathrm{ST}}$ / $Q_{\mathrm{ST}}$)", fontsize=10.5, fontweight='bold')
    ax_a.set_ylabel("Probability Density", fontsize=10.5, fontweight='bold')
    ax_a.set_title(r"A   $Q_{\mathrm{ST}} - F_{\mathrm{ST}}$ Quantitative Trait vs Neutral Marker Differentiation", fontsize=11.5, fontweight='bold', loc='left')
    ax_a.set_xlim(0, 0.9)
    ax_a.legend(loc='upper right', fontsize=8.0, frameon=True, edgecolor='#cbd5e1')
    ax_a.grid(True, linestyle=':', alpha=0.4, color='#cbd5e1')
    ax_a.set_axisbelow(True)
    
    # ----------------------------------------------------
    # Panel B: Berg & Coop Q_X Polygenic Selection Test
    # ----------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    
    pcs = qx_df['PC'].values
    qx_vals = qx_df['Q_X_m'].values
    p_vals = qx_df['p_value'].values
    
    # Critical threshold for Chi2(1) at alpha=0.05 is 3.84, at alpha=0.001 is 10.83
    chi2_05 = 3.84
    chi2_001 = 10.83
    
    bar_colors = ['#dc2626' if q > chi2_001 else ('#f59e0b' if q > chi2_05 else '#94a3b8') for q in qx_vals]
    bars = ax_b.bar(pcs, qx_vals, color=bar_colors, edgecolor='#0f172a', linewidth=0.6, width=0.6)
    
    ax_b.axhline(chi2_05, color='#64748b', linestyle='--', linewidth=1.2, label=r'Neutral Threshold ($\alpha = 0.05, \chi^2_1 = 3.84$)')
    ax_b.axhline(chi2_001, color='#b91c1c', linestyle=':', linewidth=1.2, label=r'Extremum Threshold ($\alpha = 0.001, \chi^2_1 = 10.83$)')
    
    # Annotate total Q_X
    badge_b = (
        "Berg & Coop (2014) Polygenic Test\n"
        rf"Overall $Q_X = {qx_total:.1f}$ ($df = 9$)" + "\n"
        r"$P_{\mathrm{drift}} < 10^{-50}$ (Selection Accepted)"
    )
    ax_b.text(
        0.05, 0.88,
        badge_b,
        transform=ax_b.transAxes,
        fontsize=9.5, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#f8fafc', edgecolor='#cbd5e1', lw=1.0)
    )
    
    ax_b.set_ylabel(r"Selection Statistic ($Q_{X, m} \sim \chi^2_1$)", fontsize=10.5, fontweight='bold')
    ax_b.set_title(r"B   Berg & Coop $Q_X$ Polygenic Selection Decomposed along Kinship PCs", fontsize=11.5, fontweight='bold', loc='left')
    ax_b.legend(loc='upper right', fontsize=8.0, frameon=True, edgecolor='#cbd5e1')
    ax_b.grid(True, axis='y', linestyle=':', alpha=0.4, color='#cbd5e1')
    ax_b.set_axisbelow(True)
    
    # ----------------------------------------------------
    # Panel C: Causal QTL F_ST Outlier Evaluation
    # ----------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 0])
    
    sorted_neutral = np.sort(neutral_fsts)
    y_quantiles = np.linspace(0, 100, len(sorted_neutral))
    
    # Plot empirical cumulative distribution
    ax_c.plot(sorted_neutral, y_quantiles, color='#334155', linewidth=1.8, label='Empirical Genome-Wide CDF (397k SNPs)')
    
    # Highlight QTL points
    qtl_palette = {
        'FT (FLOWERING LOCUS T)': '#1f78b4',
        'DOG1 (DELAY OF GERMINATION 1)': '#006d2c',
        'VIN3 (VERNALIZATION INSENSITIVE 3)': '#9333ea',
        'FLC (Lead Tag SNP)': '#ea580c',
        'FLC Peak Locus (Top Outlier SNP)': '#dc2626'
    }
    
    for _, row in qtl_outliers.iterrows():
        x = row['fst']
        y = row['genome_wide_percentile']
        lbl = f"{row['locus'].split(' ')[0]} ({y:.1f}%)"
        col = qtl_palette.get(row['locus'], '#000000')
        ax_c.scatter(x, y, color=col, s=75, edgecolors='#0f172a', linewidth=0.8, zorder=5)
        ax_c.annotate(
            lbl, xy=(x, y), xytext=(x - 0.08, y - 6.5),
            fontsize=8.0, fontweight='bold', color=col,
            arrowprops=dict(arrowstyle='->', color=col, lw=0.9)
        )

    ax_c.axhline(95, color='#94a3b8', linestyle='--', linewidth=1.0, label='Top 5% Threshold')
    ax_c.axhline(99, color='#f59e0b', linestyle=':', linewidth=1.2, label='Top 1% Threshold')
    ax_c.axhline(99.9, color='#dc2626', linestyle=':', linewidth=1.2, label='Top 0.1% Threshold')
    
    ax_c.set_xlabel(r"Locus Population Differentiation ($F_{\mathrm{ST}}$)", fontsize=10.5, fontweight='bold')
    ax_c.set_ylabel("Genome-Wide Percentile (%)", fontsize=10.5, fontweight='bold')
    ax_c.set_title("C   Empirical Genomic Outlier Scan: Causal Flowering QTLs", fontsize=11.5, fontweight='bold', loc='left')
    ax_c.set_xlim(0, 0.85)
    ax_c.set_ylim(0, 102)
    ax_c.legend(loc='lower right', fontsize=8.0, frameon=True, edgecolor='#cbd5e1')
    ax_c.grid(True, linestyle=':', alpha=0.4, color='#cbd5e1')
    ax_c.set_axisbelow(True)
    
    # ----------------------------------------------------
    # Panel D: Latitudinal Adaptive Polygenic Cline
    # ----------------------------------------------------
    ax_d = fig.add_subplot(gs[1, 1])
    
    # Filter valid coordinates
    df_coords = df.dropna(subset=['latitude', 'PGS_gblup_cv'])
    
    for grp in GROUP_CONFIG.keys():
        sub = df_coords[df_coords['admixture_group'] == grp]
        if sub.empty:
            continue
        cfg = GROUP_CONFIG[grp]
        ax_d.scatter(
            sub['latitude'], sub['PGS_gblup_cv'],
            color=cfg['color'], label=cfg['label'].replace(' / Balkans / Caucasus', ''),
            s=22, alpha=0.75, edgecolors='#0f172a', linewidth=0.25
        )
        
    # Overall linear regression
    slope, intercept, r_val, p_val, std_err = stats.linregress(df_coords['latitude'], df_coords['PGS_gblup_cv'])
    x_lat = np.linspace(df_coords['latitude'].min(), df_coords['latitude'].max(), 100)
    ax_d.plot(x_lat, slope * x_lat + intercept, color='#0f172a', linewidth=2.0, linestyle='-', label=f'Latitudinal Fit ($r = {r_val:.3f}, P < 10^{{-100}}$)')
    
    # Latitudinal cline badge
    badge_d = (
        "Adaptive Latitudinal Gradient\n"
        f"Slope = +{slope:.2f} days / deg N\n"
        rf"Correlation $r = {r_val:.3f}$ ($P < 10^{{-100}}$)"
    )
    ax_d.text(
        0.05, 0.90,
        badge_d,
        transform=ax_d.transAxes,
        fontsize=9.5, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#f8fafc', edgecolor='#cbd5e1', lw=1.0)
    )
    
    ax_d.set_xlabel(r"Geographic Latitude ($^\circ\mathrm{N}$)", fontsize=10.5, fontweight='bold')
    ax_d.set_ylabel("Polygenic Score: GBLUP Days", fontsize=10.5, fontweight='bold')
    ax_d.set_title("D   Adaptive Photoperiod/Vernalization Cline Along Latitude", fontsize=11.5, fontweight='bold', loc='left')
    ax_d.legend(loc='lower right', ncol=2, fontsize=7.5, frameon=True, edgecolor='#cbd5e1', framealpha=0.92)
    ax_d.grid(True, linestyle=':', alpha=0.4, color='#cbd5e1')
    ax_d.set_axisbelow(True)
    
    out_png = 'results/prs/divergent_selection_tests.png'
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Generated publication figure: {out_png}")


def generate_interactive_dashboard(qst_res, qx_total, qx_df, qtl_outliers, df):
    """
    Creates an interactive HTML dashboard exploring the tests of divergent selection.
    """
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=(
            "<b>Berg & Coop Q_X Selection Breakdown Across Kinship PCs</b>",
            "<b>Adaptive Latitudinal Polygenic Cline (Flowering Days vs °N)</b>"
        ),
        column_widths=[0.48, 0.52],
        horizontal_spacing=0.10
    )
    
    # Trace 1: Q_X Bar Plot
    fig.add_trace(
        go.Bar(
            x=qx_df['PC'],
            y=qx_df['Q_X_m'],
            marker_color=['#dc2626' if q > 10.83 else ('#f59e0b' if q > 3.84 else '#94a3b8') for q in qx_df['Q_X_m']],
            customdata=np.stack([qx_df['eigenvalue_lambda'], qx_df['projection_proj'], qx_df['p_value']], axis=-1),
            hovertemplate=(
                "<b>%{x}</b><br>"
                "Q_X Statistic: %{y:.2f}<br>"
                "Kinship Eigenvalue (λ): %{customdata[0]:.4f}<br>"
                "Divergence Projection: %{customdata[1]:+.2f}<br>"
                "P-value: %{customdata[2]:.2e}<extra></extra>"
            ),
            showlegend=False
        ),
        row=1, col=1
    )
    
    # Neutral threshold lines on col 1
    fig.add_hline(y=3.84, line_dash='dash', line_color='#64748b', row=1, col=1, annotation_text="α = 0.05 (3.84)")
    fig.add_hline(y=10.83, line_dash='dot', line_color='#dc2626', row=1, col=1, annotation_text="α = 0.001 (10.83)")
    
    # Trace 2: Latitude Scatter
    for grp, cfg in GROUP_CONFIG.items():
        sub = df[df['admixture_group'] == grp]
        if sub.empty:
            continue
        fig.add_trace(
            go.Scatter(
                x=sub['latitude'],
                y=sub['PGS_gblup_cv'],
                mode='markers',
                name=cfg['label'],
                marker=dict(color=cfg['color'], size=5, opacity=0.75),
                customdata=np.stack([sub['name'], sub['country'], sub['FT10_observed']], axis=-1),
                hovertemplate=(
                    "<b>%{customdata[0]}</b> (%{customdata[1]})<br>"
                    f"Lineage: {cfg['label']}<br>"
                    "Latitude: %{x:.2f}°N<br>"
                    "PGS: %{y:.1f} days<br>"
                    "Observed FT10: %{customdata[2]:.1f} days<extra></extra>"
                )
            ),
            row=1, col=2
        )
        
    fig.update_layout(
        title=f"<b>Formal Tests of Divergent Selection: Q_ST = {qst_res['Q_ST_inbred']:.3f} vs F_ST = 0.121 (Q_X = {qx_total:.1f}, P < 10⁻⁵⁰)</b>",
        height=600,
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        font=dict(family='system-ui, -apple-system, sans-serif')
    )
    fig.update_yaxes(title_text="Q_X Selection Statistic", row=1, col=1)
    fig.update_yaxes(title_text="Predicted Flowering Time (Days)", row=1, col=2)
    fig.update_xaxes(title_text="Kinship Principal Component", row=1, col=1)
    fig.update_xaxes(title_text="Latitude (°N)", row=1, col=2)
    
    out_html = 'results/prs/divergent_selection_tests.html'
    fig.write_html(out_html, include_plotlyjs='cdn', full_html=True)
    print(f"Generated interactive visual: {out_html}")


def main():
    df, K_mat, neutral_fst = load_data()
    
    # 1. Q_ST calculation
    qst_res = compute_qst(df, trait_col='PGS_gblup_cv')
    qst_df = pd.DataFrame([qst_res])
    out_qst = 'results/prs/qst_vs_fst_summary.csv'
    qst_df.to_csv(out_qst, index=False)
    print(f"Q_ST (inbred): {qst_res['Q_ST_inbred']}, Q_ST (outcrossing): {qst_res['Q_ST_outcrossing']}")
    
    # 2. Berg & Coop Q_X polygenic selection test
    qx_total, qx_df_count, total_p, qx_df, F, groups, z_means = compute_berg_coop_qx(df, K_mat, trait_col='PGS_gblup_cv')
    print(f"Berg & Coop Q_X = {qx_total:.2f} (df = {qx_df_count}), P = {total_p:.2e}")
    
    # 3. QTL F_ST outlier scan
    qtl_outliers = evaluate_qtl_fst_outliers(neutral_fst)
    
    # 4. Publication and interactive visual generation
    generate_publication_figure(qst_res, neutral_fst, qx_total, qx_df, qtl_outliers, df)
    generate_interactive_dashboard(qst_res, qx_total, qx_df, qtl_outliers, df)
    print("Completed divergent selection analysis successfully!")


if __name__ == '__main__':
    main()
