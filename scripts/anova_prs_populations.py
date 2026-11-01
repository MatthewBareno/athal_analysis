#!/usr/bin/env python3
"""
ANOVA and Corrected Multiple Comparisons for Polygenic Scores Across Macro-Regions

Performs:
1. One-way ANOVA and Welch's ANOVA testing average polygenic score (PGS) difference
   across the 10 1001 Genomes macro-regions / genetic lineages.
2. Pairwise multiple comparisons:
   - Tukey's Honestly Significant Difference (HSD) test
   - Welch's two-sample t-tests (allowing unequal variances)
   - Multiple testing corrections: Bonferroni, Holm, and Benjamini-Hochberg (FDR)
   - Standardized effect sizes: Cohen's d and Hedges' g
3. Detailed summary statistics per population.
4. Publication-grade multi-panel figure (PNG @ 300 DPI) and interactive HTML report.
"""

import os
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from itertools import combinations

# Setup directory
os.makedirs('results/prs', exist_ok=True)

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
    csv_path = 'data/processed/athal_1003_polygenic_scores.csv'
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} accessions from {csv_path}")
    return df


def compute_group_summaries(df, score_col='PGS_gblup_cv'):
    records = []
    for grp, sub in df.groupby('admixture_group'):
        vals = sub[score_col].dropna()
        n = len(vals)
        mean_val = vals.mean()
        std_val = vals.std(ddof=1)
        sem_val = std_val / np.sqrt(n) if n > 1 else 0
        ci95_low, ci95_high = stats.t.interval(0.95, df=n-1, loc=mean_val, scale=sem_val) if n > 1 else (mean_val, mean_val)
        median_val = vals.median()
        q25, q75 = vals.quantile(0.25), vals.quantile(0.75)
        iqr_val = q75 - q25

        # Also get observed FT10 mean and 4-QTL PGS mean
        ft10_mean = sub['FT10_observed'].mean()
        qtl_mean = sub['PGS_4qtl_cv'].mean()

        cfg = GROUP_CONFIG.get(grp, {'label': grp, 'color': '#7f7f7f'})
        records.append({
            'admixture_group': grp,
            'label': cfg['label'],
            'n': n,
            'mean_pgs': round(mean_val, 2),
            'std_pgs': round(std_val, 2),
            'sem_pgs': round(sem_val, 2),
            'ci95_low': round(ci95_low, 2),
            'ci95_high': round(ci95_high, 2),
            'median_pgs': round(median_val, 2),
            'q25_pgs': round(q25, 2),
            'q75_pgs': round(q75, 2),
            'iqr_pgs': round(iqr_val, 2),
            'mean_ft10_observed': round(ft10_mean, 2),
            'mean_pgs_4qtl': round(qtl_mean, 2)
        })

    summary_df = pd.DataFrame(records).sort_values('mean_pgs', ascending=True).reset_index(drop=True)
    out_summary = 'results/prs/pgs_population_group_summary.csv'
    summary_df.to_csv(out_summary, index=False)
    print(f"Saved population group summary to {out_summary}")
    return summary_df


def run_anova(df, score_col='PGS_gblup_cv'):
    groups = [sub[score_col].values for _, sub in df.groupby('admixture_group')]
    
    # 1. Standard Fisher One-way ANOVA
    f_stat, p_val = stats.f_oneway(*groups)
    
    # OLS formulation for ANOVA sum of squares and degrees of freedom
    model = ols(f'{score_col} ~ C(admixture_group)', data=df).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    
    ss_between = anova_table.loc['C(admixture_group)', 'sum_sq']
    df_between = anova_table.loc['C(admixture_group)', 'df']
    ms_between = ss_between / df_between
    
    ss_resid = anova_table.loc['Residual', 'sum_sq']
    df_resid = anova_table.loc['Residual', 'df']
    ms_resid = ss_resid / df_resid
    
    ss_total = ss_between + ss_resid
    df_total = df_between + df_resid
    
    # Effect sizes
    eta_squared = ss_between / ss_total
    omega_squared = (ss_between - df_between * ms_resid) / (ss_total + ms_resid)
    
    # 2. Levene's test for homogeneity of variances
    levene_stat, levene_pval = stats.levene(*groups, center='median')
    
    # 3. Welch's ANOVA (does not assume equal variance)
    # Using scipy.stats.alexandergovern
    ag_res = stats.alexandergovern(*groups)
    
    anova_results = {
        'metric': [
            'F_statistic', 'p_value',
            'df_between', 'df_within', 'df_total',
            'SS_between', 'SS_within', 'SS_total',
            'MS_between', 'MS_within',
            'eta_squared', 'omega_squared',
            'levene_statistic', 'levene_p_value',
            'alexander_govern_stat', 'alexander_govern_p_value'
        ],
        'value': [
            f_stat, p_val,
            df_between, df_resid, df_total,
            ss_between, ss_resid, ss_total,
            ms_between, ms_resid,
            eta_squared, omega_squared,
            levene_stat, levene_pval,
            ag_res.statistic, ag_res.pvalue
        ]
    }
    anova_res_df = pd.DataFrame(anova_results)
    out_anova = 'results/prs/pgs_population_anova_stats.csv'
    anova_res_df.to_csv(out_anova, index=False)
    print(f"One-way ANOVA: F({int(df_between)}, {int(df_resid)}) = {f_stat:.2f}, p = {p_val:.2e}")
    print(f"Eta-squared = {eta_squared:.4f} (explained variance: {eta_squared*100:.1f}%)")
    print(f"Omega-squared = {omega_squared:.4f}")
    print(f"Levene's test: stat = {levene_stat:.2f}, p = {levene_pval:.2e}")
    print(f"Alexander-Govern (Welch ANOVA): A = {ag_res.statistic:.2f}, p = {ag_res.pvalue:.2e}")
    return anova_res_df, f_stat, p_val, eta_squared


def run_pairwise_comparisons(df, ordered_groups, score_col='PGS_gblup_cv'):
    # 1. Tukey's HSD
    tukey = pairwise_tukeyhsd(endog=df[score_col], groups=df['admixture_group'], alpha=0.05)
    tukey_df = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])
    
    # Index tukey results by sorted tuple for fast lookup
    tukey_map = {}
    for _, row in tukey_df.iterrows():
        pair_key = tuple(sorted([row['group1'], row['group2']]))
        tukey_map[pair_key] = {
            'tukey_diff': row['meandiff'],
            'tukey_p_adj': row['p-adj'],
            'tukey_ci_low': row['lower'],
            'tukey_ci_high': row['upper'],
            'tukey_reject': row['reject']
        }

    # 2. Pairwise Welch t-tests and effect sizes for all 45 combinations
    records = []
    pairs = list(combinations(ordered_groups, 2))
    
    for g1, g2 in pairs:
        s1 = df[df['admixture_group'] == g1][score_col]
        s2 = df[df['admixture_group'] == g2][score_col]
        n1, n2 = len(s1), len(s2)
        m1, m2 = s1.mean(), s2.mean()
        v1, v2 = s1.var(ddof=1), s2.var(ddof=1)
        diff = m2 - m1  # Direction: g2 - g1
        
        # Welch's t-test (unequal variance)
        t_res = stats.ttest_ind(s2, s1, equal_var=False)
        t_stat = t_res.statistic
        p_raw = t_res.pvalue
        df_welch = t_res.df
        
        # Standardized effect sizes
        # Pooled standard deviation
        s_pooled = np.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
        cohen_d = diff / s_pooled if s_pooled > 0 else 0
        # Hedges' g correction factor for sample size
        hedges_j = 1 - (3 / (4 * (n1 + n2) - 9))
        hedges_g = cohen_d * hedges_j
        
        # Tukey lookup
        tk = tukey_map.get(tuple(sorted([g1, g2])), {})
        
        records.append({
            'group1': g1,
            'label1': GROUP_CONFIG.get(g1, {}).get('label', g1),
            'n1': n1,
            'mean1': round(m1, 2),
            'group2': g2,
            'label2': GROUP_CONFIG.get(g2, {}).get('label', g2),
            'n2': n2,
            'mean2': round(m2, 2),
            'diff_mean_days': round(diff, 2),
            'welch_t_stat': round(t_stat, 3),
            'welch_df': round(df_welch, 1),
            'p_raw': p_raw,
            'cohen_d': round(cohen_d, 3),
            'hedges_g': round(hedges_g, 3),
            'tukey_p_adj': tk.get('tukey_p_adj', np.nan),
            'tukey_ci_low': tk.get('tukey_ci_low', np.nan),
            'tukey_ci_high': tk.get('tukey_ci_high', np.nan),
            'tukey_reject': tk.get('tukey_reject', False)
        })

    pairwise_df = pd.DataFrame(records)
    
    # Multiple testing corrections on raw p-values
    # 1. Bonferroni
    rej_bonf, p_bonf, _, _ = multipletests(pairwise_df['p_raw'], alpha=0.05, method='bonferroni')
    # 2. Holm-Bonferroni step-down
    rej_holm, p_holm, _, _ = multipletests(pairwise_df['p_raw'], alpha=0.05, method='holm')
    # 3. Benjamini-Hochberg FDR
    rej_fdr, p_fdr, _, _ = multipletests(pairwise_df['p_raw'], alpha=0.05, method='fdr_bh')
    
    pairwise_df['p_bonferroni'] = p_bonf
    pairwise_df['reject_bonferroni'] = rej_bonf
    pairwise_df['p_holm'] = p_holm
    pairwise_df['reject_holm'] = rej_holm
    pairwise_df['p_fdr_bh'] = p_fdr
    pairwise_df['reject_fdr_bh'] = rej_fdr
    
    # Significance stars based on FDR
    def get_stars(p):
        if p < 0.001:
            return '***'
        elif p < 0.01:
            return '**'
        elif p < 0.05:
            return '*'
        else:
            return 'ns'
            
    pairwise_df['sig_fdr'] = pairwise_df['p_fdr_bh'].apply(get_stars)
    pairwise_df['sig_bonf'] = pairwise_df['p_bonferroni'].apply(get_stars)
    
    out_pw = 'results/prs/pgs_population_pairwise_comparisons.csv'
    pairwise_df.to_csv(out_pw, index=False)
    print(f"Saved pairwise comparisons to {out_pw}")
    print(f"Significant pairs out of 45: Bonferroni = {rej_bonf.sum()}, Holm = {rej_holm.sum()}, FDR BH = {rej_fdr.sum()}, Tukey HSD = {pairwise_df['tukey_reject'].sum()}")
    return pairwise_df


def generate_publication_figure(df, summary_df, pairwise_df, anova_stats, ordered_groups):
    """
    Creates a multi-panel publication-grade figure (PNG @ 300 DPI):
    - Panel A: Ordered Forest / Distribution Plot of polygenic scores across populations
    - Panel B: Pairwise Statistical Significance Matrix (-log10 FDR q-values)
    - Panel C: Volcano Plot of Effect Size (Delta Days) vs -log10(p_FDR)
    """
    fig = plt.figure(figsize=(18, 12), facecolor='#ffffff')
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1.0], width_ratios=[1.15, 1.0], hspace=0.32, wspace=0.25)
    
    f_stat = anova_stats.loc[anova_stats['metric'] == 'F_statistic', 'value'].values[0]
    p_val = anova_stats.loc[anova_stats['metric'] == 'p_value', 'value'].values[0]
    eta_sq = anova_stats.loc[anova_stats['metric'] == 'eta_squared', 'value'].values[0]
    
    # ----------------------------------------------------
    # Panel A: Population Polygenic Score Distributions
    # ----------------------------------------------------
    ax_a = fig.add_subplot(gs[0, :])
    
    # Boxplot + Mean Point CI
    positions = np.arange(len(ordered_groups))
    group_data = [df[df['admixture_group'] == grp]['PGS_gblup_cv'].values for grp in ordered_groups]
    colors = [GROUP_CONFIG.get(grp, {}).get('color', '#64748b') for grp in ordered_groups]
    labels_clean = [GROUP_CONFIG.get(grp, {}).get('label', grp).replace(' / Balkans / Caucasus', '') for grp in ordered_groups]
    
    # Add subtle jitter points
    np.random.seed(42)
    for i, (grp, data) in enumerate(zip(ordered_groups, group_data)):
        jitter = np.random.normal(0, 0.08, size=len(data))
        ax_a.scatter(
            i + jitter, data,
            color=colors[i], alpha=0.30, s=18, edgecolors='none', zorder=2
        )
        
    # Overlay elegant boxplots
    bp = ax_a.boxplot(
        group_data, positions=positions, widths=0.45,
        patch_artist=True, showfliers=False, zorder=3,
        boxprops=dict(facecolor='white', edgecolor='#334155', linewidth=1.4, alpha=0.85),
        medianprops=dict(color='#0f172a', linewidth=2.0),
        whiskerprops=dict(color='#475569', linewidth=1.2),
        capprops=dict(color='#475569', linewidth=1.2)
    )
    
    # Overlay Group Mean and 95% CI marker
    for i, grp in enumerate(ordered_groups):
        sub_sum = summary_df[summary_df['admixture_group'] == grp].iloc[0]
        m = sub_sum['mean_pgs']
        low, high = sub_sum['ci95_low'], sub_sum['ci95_high']
        ax_a.errorbar(
            i, m, yerr=[[m - low], [high - m]],
            fmt='D', color=colors[i], markeredgecolor='#0f172a', markeredgewidth=1.2,
            markersize=7.5, elinewidth=2.2, capsize=4.5, capthick=1.8, zorder=4
        )
        # Annotate mean value above
        ax_a.text(
            i, high + 2.2, f"{m:.1f}d\n(n={sub_sum['n']})",
            ha='center', va='bottom', fontsize=8.0, fontweight='bold', color='#1e293b'
        )

    ax_a.set_xticks(positions)
    ax_a.set_xticklabels(labels_clean, fontsize=9.5, fontweight='bold', color='#0f172a')
    ax_a.set_ylabel("Cross-Validated Polygenic Score (Days to Flowering)", fontsize=10.5, fontweight='bold', color='#0f172a')
    ax_a.set_title(
        f"A   Polygenic Score Differentiation Across Macro-Regions  "
        f"(One-Way ANOVA: F(9, 993) = {f_stat:.1f}, P = {p_val:.2e}, η² = {eta_sq*100:.1f}%)",
        fontsize=12, fontweight='bold', loc='left', pad=12
    )
    ax_a.grid(True, axis='y', linestyle='--', alpha=0.4, color='#cbd5e1')
    ax_a.set_axisbelow(True)
    ax_a.set_ylim(45, 145)
    
    # ----------------------------------------------------
    # Panel B: Pairwise Statistical Significance Matrix
    # ----------------------------------------------------
    ax_b = fig.add_subplot(gs[1, 0])
    
    n_grps = len(ordered_groups)
    p_mat = np.zeros((n_grps, n_grps))
    diff_mat = np.zeros((n_grps, n_grps))
    star_mat = np.empty((n_grps, n_grps), dtype=object)
    star_mat.fill('')
    
    # Fill matrix: lower triangle = -log10(p_FDR), upper triangle = mean difference (days)
    for _, row in pairwise_df.iterrows():
        i = ordered_groups.index(row['group1'])
        j = ordered_groups.index(row['group2'])
        p_val_adj = max(row['p_fdr_bh'], 1e-45)
        log_p = -np.log10(p_val_adj)
        diff = row['diff_mean_days']
        star = row['sig_fdr']
        
        # Store in lower triangle: row > col
        row_idx = max(i, j)
        col_idx = min(i, j)
        p_mat[row_idx, col_idx] = log_p
        star_mat[row_idx, col_idx] = star
        diff_mat[row_idx, col_idx] = abs(diff)

    # Mask upper triangle and diagonal for heatmap
    mask = np.triu(np.ones_like(p_mat, dtype=bool))
    masked_p = np.ma.masked_where(mask, p_mat)
    
    cmap = plt.cm.YlOrRd
    cmap.set_bad('white')
    im = ax_b.imshow(masked_p, cmap=cmap, vmin=0, vmax=35, aspect='auto')
    
    # Annotate cells
    for r in range(n_grps):
        for c in range(r):
            val = p_mat[r, c]
            st = star_mat[r, c]
            d = diff_mat[r, c]
            txt_color = 'white' if val > 18 else '#0f172a'
            ax_b.text(c, r, f"{st}\n{d:.1f}d", ha='center', va='center',
                      fontsize=7.5, fontweight='bold', color=txt_color)
            
    cbar = fig.colorbar(im, ax=ax_b, fraction=0.035, pad=0.03)
    cbar.set_label(r"$-\log_{10}(P_{\mathrm{FDR}})$", fontsize=9.5, fontweight='bold')
    
    # Formatting axes
    short_labels = [l.split(' ')[0] if l != 'Italy / Balkans / Caucasus' else 'Italy' for l in labels_clean]
    ax_b.set_xticks(range(n_grps))
    ax_b.set_xticklabels(labels_clean, rotation=45, ha='right', fontsize=8.5, fontweight='bold')
    ax_b.set_yticks(range(n_grps))
    ax_b.set_yticklabels(labels_clean, fontsize=8.5, fontweight='bold')
    ax_b.set_title(
        r"B   Pairwise Contrasts Matrix ($-\log_{10} P_{\mathrm{FDR}}$ and Mean Difference)",
        fontsize=11.5, fontweight='bold', loc='left', pad=10
    )
    
    # ----------------------------------------------------
    # Panel C: Volcano Plot of Mean Difference vs Significance
    # ----------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 1])
    
    log_p_all = -np.log10(np.clip(pairwise_df['p_fdr_bh'], 1e-45, 1.0))
    diff_all = pairwise_df['diff_mean_days']
    sig_colors = []
    for p in pairwise_df['p_fdr_bh']:
        if p < 0.001:
            sig_colors.append('#dc2626')  # Dark red
        elif p < 0.05:
            sig_colors.append('#f59e0b')  # Amber
        else:
            sig_colors.append('#94a3b8')  # Slate grey

    ax_c.scatter(
        diff_all, log_p_all,
        c=sig_colors, s=55, alpha=0.85, edgecolors='#1e293b', linewidth=0.5, zorder=3
    )
    
    # Significance threshold lines
    ax_c.axhline(-np.log10(0.05), color='#64748b', linestyle='--', linewidth=1.2, label='FDR q = 0.05')
    ax_c.axhline(-np.log10(0.001), color='#b91c1c', linestyle=':', linewidth=1.2, label='FDR q = 0.001')
    
    # Annotate top contrasting pairs
    top_contrasts = pairwise_df.sort_values('p_raw').head(5)
    for _, top_row in top_contrasts.iterrows():
        x = top_row['diff_mean_days']
        y = -np.log10(max(top_row['p_fdr_bh'], 1e-45))
        lbl = f"{top_row['label2'].split(' ')[0]} vs {top_row['label1'].split(' ')[0]}"
        ax_c.annotate(
            lbl, xy=(x, y), xytext=(x - 4.5, y - 2.5),
            fontsize=7.5, fontweight='bold', color='#1e293b',
            arrowprops=dict(arrowstyle='->', color='#475569', lw=0.8)
        )

    ax_c.set_xlabel(r"Difference in Mean Polygenic Score ($\Delta$ Days)", fontsize=10, fontweight='bold')
    ax_c.set_ylabel(r"$-\log_{10}(P_{\mathrm{FDR}})$", fontsize=10, fontweight='bold')
    ax_c.set_title(
        r"C   Volcano Plot of Population Contrasts (45 Pairs)",
        fontsize=11.5, fontweight='bold', loc='left', pad=10
    )
    ax_c.legend(loc='upper left', fontsize=8.5, frameon=True, edgecolor='#cbd5e1')
    ax_c.grid(True, linestyle=':', alpha=0.5, color='#cbd5e1')
    ax_c.set_axisbelow(True)
    
    # Save figure
    out_png = 'results/prs/pgs_population_contrasts.png'
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Generated publication figure: {out_png}")


def generate_interactive_dashboard(df, summary_df, pairwise_df, anova_stats, ordered_groups):
    """
    Creates an interactive HTML dashboard:
    1. Interactive Box / Violin plots by population.
    2. Interactive Pairwise Contrast Heatmap with hover details.
    3. Pairwise comparison table.
    """
    labels = [GROUP_CONFIG.get(grp, {}).get('label', grp) for grp in ordered_groups]
    colors = [GROUP_CONFIG.get(grp, {}).get('color', '#64748b') for grp in ordered_groups]
    
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=(
            "<b>Polygenic Score Distributions by Macro-Region</b>",
            "<b>Pairwise Difference Matrix (|Δ Days| & FDR Significance)</b>"
        ),
        column_widths=[0.52, 0.48],
        horizontal_spacing=0.10
    )
    
    # 1. Boxplots
    for i, grp in enumerate(ordered_groups):
        sub_vals = df[df['admixture_group'] == grp]['PGS_gblup_cv']
        cfg = GROUP_CONFIG.get(grp, {})
        fig.add_trace(
            go.Box(
                y=sub_vals,
                name=cfg.get('label', grp),
                marker_color=cfg.get('color', '#64748b'),
                boxpoints='all',
                jitter=0.25,
                pointpos=0,
                marker=dict(size=4, opacity=0.5),
                showlegend=False
            ),
            row=1, col=1
        )
        
    # 2. Heatmap of pairwise differences
    n_grps = len(ordered_groups)
    diff_matrix = np.full((n_grps, n_grps), np.nan)
    hover_matrix = np.empty((n_grps, n_grps), dtype=object)
    
    for _, row in pairwise_df.iterrows():
        i = ordered_groups.index(row['group1'])
        j = ordered_groups.index(row['group2'])
        d = abs(row['diff_mean_days'])
        diff_matrix[i, j] = d
        diff_matrix[j, i] = d
        
        hov_text = (
            f"<b>{row['label1']} vs {row['label2']}</b><br>"
            f"Mean Difference: {row['diff_mean_days']:+.2f} days<br>"
            f"Welch t: {row['welch_t_stat']:.2f} (df={row['welch_df']:.1f})<br>"
            f"FDR q-value: {row['p_fdr_bh']:.2e} ({row['sig_fdr']})<br>"
            f"Bonferroni p: {row['p_bonferroni']:.2e}<br>"
            f"Tukey p-adj: {row['tukey_p_adj']:.4f}<br>"
            f"Cohen's d: {row['cohen_d']:.2f}"
        )
        hover_matrix[i, j] = hov_text
        hover_matrix[j, i] = hov_text

    for i in range(n_grps):
        diff_matrix[i, i] = 0.0
        hover_matrix[i, i] = f"<b>{labels[i]}</b><br>Identical Population (Δ = 0)"

    fig.add_trace(
        go.Heatmap(
            z=diff_matrix,
            x=labels,
            y=labels,
            text=hover_matrix,
            hoverinfo='text',
            colorscale='Blues',
            colorbar=dict(title='|Δ Days|', len=0.8, x=1.02),
            showscale=True
        ),
        row=1, col=2
    )
    
    f_stat = anova_stats.loc[anova_stats['metric'] == 'F_statistic', 'value'].values[0]
    p_val = anova_stats.loc[anova_stats['metric'] == 'p_value', 'value'].values[0]
    eta_sq = anova_stats.loc[anova_stats['metric'] == 'eta_squared', 'value'].values[0]

    fig.update_layout(
        title=f"<b>Arabidopsis thaliana Polygenic Score Population Differentiation (ANOVA F = {f_stat:.1f}, P = {p_val:.2e}, η² = {eta_sq*100:.1f}%)</b>",
        height=620,
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        font=dict(family='system-ui, -apple-system, sans-serif')
    )
    
    fig.update_yaxes(title_text="Predicted Flowering Time (Days)", row=1, col=1)
    fig.update_xaxes(tickangle=35, row=1, col=1)
    fig.update_xaxes(tickangle=35, row=1, col=2)
    fig.update_yaxes(autorange='reversed', row=1, col=2)
    
    out_html = 'results/prs/pgs_population_contrasts.html'
    fig.write_html(out_html, include_plotlyjs='cdn', full_html=True)
    print(f"Generated interactive dashboard: {out_html}")


def main():
    df = load_data()
    summary_df = compute_group_summaries(df, score_col='PGS_gblup_cv')
    ordered_groups = summary_df['admixture_group'].tolist()
    
    anova_stats, f_stat, p_val, eta_sq = run_anova(df, score_col='PGS_gblup_cv')
    pairwise_df = run_pairwise_comparisons(df, ordered_groups, score_col='PGS_gblup_cv')
    
    generate_publication_figure(df, summary_df, pairwise_df, anova_stats, ordered_groups)
    generate_interactive_dashboard(df, summary_df, pairwise_df, anova_stats, ordered_groups)
    print("Done analysis!")


if __name__ == '__main__':
    main()
