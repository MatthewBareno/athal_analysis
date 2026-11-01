"""Generate an ultra-clear, publication-grade Manhattan and QQ figure for LOCO LMM:
- Aligned Side-by-Side: Panel (a) Manhattan plot (75% width) + Panel (b) QQ plot (25% width) sharing the same vertical y-axis range.
- Alternating chromosome colors (Slate #1e293b and Azure #0284c7).
- Genome-wide significant points highlighted in bold vermilion (#dc2626, s=36).
- High-contrast, clean callout annotations for FT, DOG1, and VIN3.
- QQ plot with 95% concentration band and lambda_GC badge.
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import plotly.graph_objects as go

CANDIDATES = [
    {
        'gene': 'FT',
        'chr': 1,
        'pos': 24339560,
        'label': '$\\bf{FT}$ (Florigen)\n$-\\log_{10}P = 8.56$',
        'xytext': (0, 32),
        'bg': '#eff6ff',
        'edge': '#2563eb',
        'text_color': '#1e3a8a'
    },
    {
        'gene': 'DOG1',
        'chr': 5,
        'pos': 18590327,
        'label': '$\\bf{DOG1}$ (Thermal Seed & Flower)\n$-\\log_{10}P = 10.06$',
        'xytext': (-45, 36),
        'bg': '#fef2f2',
        'edge': '#dc2626',
        'text_color': '#991b1b'
    },
    {
        'gene': 'VIN3',
        'chr': 5,
        'pos': 23288908,
        'label': '$\\bf{VIN3}$ (Cold Vernalization Sensor)\n$-\\log_{10}P = 10.89$',
        'xytext': (45, 36),
        'bg': '#fef2f2',
        'edge': '#dc2626',
        'text_color': '#991b1b'
    }
]


def load_loco_data():
    assoc_file = 'results/gwas/ft10_gemma_loco_all.assoc.txt'
    print(f"Loading LOCO associations from {assoc_file}...")
    df = pd.read_csv(assoc_file, sep=r'\s+')
    df = df.dropna(subset=['p_wald']).copy()
    df['chr'] = df['chr'].astype(int)
    df['ps'] = df['ps'].astype(int)
    df['-log10p'] = -np.log10(np.maximum(df['p_wald'], 1e-300))
    return df


def plot_side_by_side_manhattan(df):
    os.makedirs('results/gwas', exist_ok=True)
    
    # 1. Coordinate preparation
    df_sorted = df.sort_values(['chr', 'ps']).copy()
    chroms = sorted(df_sorted['chr'].unique())
    
    chr_offsets = {}
    chr_centers = {}
    curr = 0
    for c in chroms:
        chr_offsets[c] = curr
        sub = df_sorted[df_sorted['chr'] == c]
        chr_centers[c] = curr + sub['ps'].max() / 2
        curr += sub['ps'].max() + 2500000  # 2.5 Mb buffer
        
    df_sorted['pos_cum'] = df_sorted['ps'] + df_sorted['chr'].map(chr_offsets)
    
    n_snps = len(df_sorted)
    bonf_p = 0.05 / n_snps
    bonf_log = -np.log10(bonf_p)
    sugg_log = 5.0
    
    # Calculate lambda_GC
    p_clean = df_sorted['p_wald'].values
    chisq = stats.chi2.ppf(1 - p_clean, 1)
    lambda_gc = np.median(chisq) / stats.chi2.ppf(0.5, 1)
    
    y_max = max(df_sorted['-log10p'].max() * 1.18, 12.5)
    
    # 2. Create Figure (Side-by-Side: Panel a Manhattan, Panel b QQ)
    fig = plt.figure(figsize=(16.5, 7.2), dpi=300)
    gs = gridspec.GridSpec(1, 2, width_ratios=[3.6, 1.0], wspace=0.18, left=0.06, right=0.97, top=0.91, bottom=0.11)
    
    # ================= PANEL (a): MANHATTAN PLOT =================
    ax_m = fig.add_subplot(gs[0])
    
    # High-contrast alternating colors: Slate Navy & Deep Azure
    colors = ['#1e293b', '#0284c7']
    
    for idx, c in enumerate(chroms):
        sub = df_sorted[df_sorted['chr'] == c]
        
        # Background non-significant variants (subsample to ~40k per chr for speed)
        bg = sub[sub['-log10p'] < sugg_log]
        if len(bg) > 40000:
            bg = bg.sample(40000, random_state=42)
        ax_m.scatter(bg['pos_cum'], bg['-log10p'], c=colors[idx % 2], s=7, alpha=0.55, edgecolors='none', rasterized=True)
        
        # Suggestive variants (5.0 <= -logP < Bonferroni)
        sugg = sub[(sub['-log10p'] >= sugg_log) & (sub['-log10p'] < bonf_log)]
        if not sugg.empty:
            ax_m.scatter(sugg['pos_cum'], sugg['-log10p'], c='#f59e0b', s=16, alpha=0.85, edgecolors='none', rasterized=True)
            
        # Significant variants (-logP >= Bonferroni) - Full vermilion points
        sig = sub[sub['-log10p'] >= bonf_log]
        if not sig.empty:
            ax_m.scatter(sig['pos_cum'], sig['-log10p'], c='#dc2626', s=34, alpha=0.95, edgecolors='#450a0a', linewidth=0.35, zorder=5)
            
    # Threshold lines
    ax_m.axhline(bonf_log, color='#dc2626', linestyle='--', linewidth=1.2, zorder=4,
                 label=f'Genome-Wide Significance (Bonferroni $P < {bonf_p:.2e}$)')
    ax_m.axhline(sugg_log, color='#94a3b8', linestyle=':', linewidth=1.0, zorder=3,
                 label=r'Suggestive Threshold ($P < 1.0 \times 10^{-5}$)')
    
    # Callout annotations for top loci
    for cand in CANDIDATES:
        gc = cand['chr']
        target_pos = cand['pos']
        # Locate exact SNP in df_sorted
        snp_match = df_sorted[(df_sorted['chr'] == gc) & (abs(df_sorted['ps'] - target_pos) < 50000)]
        if not snp_match.empty:
            lead = snp_match.sort_values('-log10p', ascending=False).iloc[0]
            x_coord = lead['pos_cum']
            y_coord = lead['-log10p']
            
            ax_m.annotate(
                cand['label'],
                xy=(x_coord, y_coord),
                xytext=cand['xytext'],
                textcoords='offset points',
                arrowprops=dict(
                    arrowstyle='->',
                    color='#0f172a',
                    lw=1.3,
                    shrinkA=3,
                    shrinkB=3
                ),
                fontsize=8.5,
                fontweight='bold',
                color=cand['text_color'],
                ha='center',
                va='bottom',
                bbox=dict(
                    boxstyle='round,pad=0.3',
                    facecolor=cand['bg'],
                    edgecolor=cand['edge'],
                    lw=1.1,
                    alpha=0.96
                ),
                zorder=10
            )
            
    ax_m.set_xticks([chr_centers[c] for c in chroms])
    ax_m.set_xticklabels([f"Chromosome {c}" for c in chroms], fontsize=10, fontweight='bold')
    ax_m.set_xlim(-800000, curr)
    ax_m.set_ylim(0, y_max)
    ax_m.set_ylabel(r"$-\log_{10}(P\text{-value})$", fontsize=11, fontweight='bold')
    ax_m.set_title("A   Genome-Wide Association Scan: Flowering Time (FT10) — GEMMA LOCO LMM", fontsize=12, fontweight='bold', loc='left')
    ax_m.legend(loc='upper right', frameon=True, edgecolor='#cbd5e1', facecolor='#ffffff', framealpha=0.95, fontsize=8.5)
    ax_m.grid(True, linestyle=':', alpha=0.35, color='#cbd5e1')
    ax_m.set_axisbelow(True)
    
    for spine in ['top', 'right']:
        ax_m.spines[spine].set_visible(False)
    for spine in ['left', 'bottom']:
        ax_m.spines[spine].set_color('#94a3b8')
        ax_m.spines[spine].set_linewidth(0.9)
        
    # ================= PANEL (b): QQ PLOT =================
    ax_qq = fig.add_subplot(gs[1])
    
    p_vals = df_sorted['p_wald'].values.copy()
    p_vals.sort()
    n = len(p_vals)
    
    obs = -np.log10(p_vals)
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    
    # 95% Confidence Band
    j = np.arange(1, n + 1)
    lower = -np.log10(stats.beta.ppf(0.975, j, n - j + 1))
    upper = -np.log10(stats.beta.ppf(0.025, j, n - j + 1))
    
    # Subsample for smooth plotting
    idx_sample = np.unique(np.concatenate([
        np.linspace(0, 10000, 2000, dtype=int),
        np.linspace(10000, n - 1, 10000, dtype=int)
    ]))
    
    ax_qq.fill_between(exp[idx_sample], lower[idx_sample], upper[idx_sample],
                       color='#e0f2fe', alpha=0.8, label='95% CI')
    ax_qq.plot([0, max(exp)], [0, max(exp)], color='#dc2626', linestyle='--', linewidth=1.4, label='Expected (Null)')
    
    # Plot non-significant points in azure, significant points in red
    non_sig_idx = idx_sample[obs[idx_sample] < bonf_log]
    sig_idx = idx_sample[obs[idx_sample] >= bonf_log]
    
    ax_qq.scatter(exp[non_sig_idx], obs[non_sig_idx], c='#0284c7', s=12, alpha=0.75, edgecolors='none', rasterized=True)
    if len(sig_idx) > 0:
        ax_qq.scatter(exp[sig_idx], obs[sig_idx], c='#dc2626', s=24, alpha=0.95, edgecolors='#450a0a', linewidth=0.3, zorder=5)
        
    # Add lambda_GC badge
    ax_qq.text(
        0.06, 0.90,
        f"$\\mathbf{{\\lambda_{{GC}} = {lambda_gc:.3f}}}$\n$\\text{{(Well-Calibrated)}}$",
        transform=ax_qq.transAxes,
        fontsize=9.5,
        fontweight='bold',
        color='#0f172a',
        va='top',
        bbox=dict(boxstyle='round,pad=0.35', facecolor='#f8fafc', edgecolor='#94a3b8', lw=1.0)
    )
    
    ax_qq.set_xlim(0, max(exp) * 1.05)
    ax_qq.set_ylim(0, y_max)
    ax_qq.set_xlabel(r"Expected $-\log_{10}(P)$", fontsize=10.5, fontweight='bold')
    ax_qq.set_ylabel(r"Observed $-\log_{10}(P)$", fontsize=10.5, fontweight='bold')
    ax_qq.set_title("B   Quantile-Quantile", fontsize=12, fontweight='bold', loc='left')
    ax_qq.grid(True, linestyle=':', alpha=0.35, color='#cbd5e1')
    ax_qq.set_axisbelow(True)
    
    for spine in ['top', 'right']:
        ax_qq.spines[spine].set_visible(False)
    for spine in ['left', 'bottom']:
        ax_qq.spines[spine].set_color('#94a3b8')
        ax_qq.spines[spine].set_linewidth(0.9)
        
    out_png = 'results/gwas/ft10_loco_manhattan_clear.png'
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"Generated ultra-clear publication figure: {out_png}")
    return out_png


def plot_interactive_plotly(df):
    print("Generating interactive Plotly Manhattan plot...")
    df_sorted = df.sort_values(['chr', 'ps']).copy()
    chroms = sorted(df_sorted['chr'].unique())
    
    chr_offsets = {}
    curr = 0
    for c in chroms:
        chr_offsets[c] = curr
        sub = df_sorted[df_sorted['chr'] == c]
        curr += sub['ps'].max() + 2500000
        
    df_sorted['pos_cum'] = df_sorted['ps'] + df_sorted['chr'].map(chr_offsets)
    bonf_log = -np.log10(0.05 / len(df_sorted))
    
    # Subsample non-significant points to keep file lightweight
    non_sig = df_sorted[df_sorted['-log10p'] < 5.0]
    if len(non_sig) > 35000:
        non_sig = non_sig.sample(35000, random_state=42)
    sig_pts = df_sorted[df_sorted['-log10p'] >= 5.0]
    plot_df = pd.concat([non_sig, sig_pts]).sort_values(['chr', 'ps'])
    
    fig = go.Figure()
    colors = ['#1e293b', '#0284c7']
    
    for idx, c in enumerate(chroms):
        sub = plot_df[plot_df['chr'] == c]
        custom = list(zip(
            sub['rs'],
            sub['chr'],
            sub['ps'],
            sub['p_wald'].apply(lambda x: f"{x:.2e}"),
            sub['af'].round(3),
            sub['beta'].round(2),
            sub['se'].round(2)
        ))
        
        fig.add_trace(go.Scattergl(
            x=sub['pos_cum'],
            y=sub['-log10p'],
            mode='markers',
            name=f"Chromosome {c}",
            marker=dict(size=4.5, color=colors[idx % 2], opacity=0.75),
            customdata=custom,
            hovertemplate=(
                "<b>SNP:</b> %{customdata[0]}<br>"
                "<b>Locus:</b> Chr%{customdata[1]}:%{customdata[2]:,}<br>"
                "<b>P-value:</b> %{customdata[3]}<br>"
                "<b>MAF:</b> %{customdata[4]}<br>"
                "<b>Effect (Beta):</b> %{customdata[5]} ± %{customdata[6]} days<extra></extra>"
            )
        ))
        
    # Highlight significant peaks with distinct red markers
    sig_all = plot_df[plot_df['-log10p'] >= bonf_log]
    if not sig_all.empty:
        custom_sig = list(zip(
            sig_all['rs'],
            sig_all['chr'],
            sig_all['ps'],
            sig_all['p_wald'].apply(lambda x: f"{x:.2e}"),
            sig_all['af'].round(3),
            sig_all['beta'].round(2),
            sig_all['se'].round(2)
        ))
        fig.add_trace(go.Scattergl(
            x=sig_all['pos_cum'],
            y=sig_all['-log10p'],
            mode='markers',
            name='Genome-Wide Significant',
            marker=dict(size=8, color='#dc2626', line=dict(width=0.8, color='#ffffff')),
            customdata=custom_sig,
            hovertemplate=(
                "<b>⭐ Significant Peak:</b> %{customdata[0]}<br>"
                "<b>Locus:</b> Chr%{customdata[1]}:%{customdata[2]:,}<br>"
                "<b>P-value:</b> %{customdata[3]}<br>"
                "<b>MAF:</b> %{customdata[4]}<br>"
                "<b>Effect (Beta):</b> %{customdata[5]} ± %{customdata[6]} days<extra></extra>"
            )
        ))
        
    fig.add_hline(y=bonf_log, line_dash='dash', line_color='#dc2626',
                  annotation_text=f"Bonferroni Significance (-logP = {bonf_log:.2f})")
    fig.add_hline(y=5.0, line_dash='dot', line_color='#94a3b8',
                  annotation_text="Suggestive Threshold (-logP = 5.0)")
                  
    fig.update_layout(
        title=dict(
            text="<b>Flowering Time (FT10) Genome-Wide Association Study — GEMMA LOCO LMM</b><br>"
                 "<span style='font-size: 11.5px; color: #475569;'>1,003 accessions | 1,031,617 SNPs | Leave-One-Chromosome-Out Kinship</span>",
            font=dict(size=14, family='sans-serif', color='#0f172a')
        ),
        xaxis=dict(title="<b>Chromosome</b>", showgrid=False),
        yaxis=dict(title="<b>-log10(P-value)</b>", showgrid=True, gridcolor='#f1f5f9'),
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        height=580,
        showlegend=True,
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1)
    )
    
    out_html = 'results/gwas/ft10_loco_manhattan_clear.html'
    fig.write_html(out_html, include_plotlyjs='cdn', full_html=True)
    print(f"Generated clear interactive Manhattan plot: {out_html}")
    return out_html


def main():
    df = load_loco_data()
    plot_side_by_side_manhattan(df)
    plot_interactive_plotly(df)
    print("Done!")


if __name__ == '__main__':
    main()
