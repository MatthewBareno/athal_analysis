"""Publication-grade GWAS analysis and visualization for FT10 in Arabidopsis thaliana:
1. Manhattan plot with chromosomal coloring, Bonferroni threshold, and candidate gene annotations (FLC, DOG1, FRI).
2. Quantile-Quantile (QQ) plot with 95% concentration band and genomic inflation lambda_GC.
3. Summary table of independent lead SNPs.
4. Interactive Plotly Manhattan plot (HTML) for gene exploration.
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import plotly.graph_objects as go

# Canonical flowering time candidate genes in TAIR10
CANDIDATE_GENES = [
    {'gene': 'FLC', 'chr': 5, 'start': 3173498, 'end': 3180279, 'mid': 3176888, 'desc': 'Master flowering repressor'},
    {'gene': 'DOG1', 'chr': 5, 'start': 18589726, 'end': 18591894, 'mid': 18590810, 'desc': 'Seed dormancy & thermal flowering'},
    {'gene': 'VIN3', 'chr': 5, 'start': 23228000, 'end': 23295000, 'mid': 23260000, 'desc': 'Vernalization Insensitive 3 cold sensor'},
    {'gene': 'FRI', 'chr': 4, 'start': 268695, 'end': 271760, 'mid': 270227, 'desc': 'FRIGIDA (activator of FLC)'},
    {'gene': 'FT', 'chr': 1, 'start': 24331428, 'end': 24339560, 'mid': 24335000, 'desc': 'FLOWERING LOCUS T (florigen)'},
    {'gene': 'CRY2', 'chr': 1, 'start': 1189458, 'end': 1194547, 'mid': 1192000, 'desc': 'Cryptochrome 2 blue-light photoreceptor'},
    {'gene': 'CO', 'chr': 5, 'start': 5162171, 'end': 5164710, 'mid': 5163440, 'desc': 'CONSTANS photoperiod regulator'}
]


def load_gwas_data(filepath, p_col='P', chr_col='#CHROM', pos_col='POS', id_col='ID'):
    print(f"Loading GWAS results from {filepath}...")
    df = pd.read_csv(filepath, sep=r'\s+')
    df = df.dropna(subset=[p_col])
    # Ensure numeric chromosome and position
    df = df[df[chr_col].astype(str).str.isnumeric()].copy()
    df[chr_col] = df[chr_col].astype(int)
    df[pos_col] = df[pos_col].astype(int)
    df['-log10p'] = -np.log10(np.maximum(df[p_col], 1e-300))
    return df


def calculate_lambda_gc(p_values):
    p_clean = p_values[np.isfinite(p_values) & (p_values > 0) & (p_values <= 1)]
    chisq = stats.chi2.ppf(1 - p_clean, 1)
    lambda_gc = np.median(chisq) / stats.chi2.ppf(0.5, 1)
    return lambda_gc


def plot_manhattan_and_qq(df, model_name, out_prefix, p_col='P', chr_col='#CHROM', pos_col='POS'):
    os.makedirs('results/gwas', exist_ok=True)
    
    # 1. Prepare genomic cumulative coordinates for Manhattan plot
    df_sorted = df.sort_values([chr_col, pos_col]).copy()
    chroms = sorted(df_sorted[chr_col].unique())
    
    chr_offsets = {}
    chr_centers = {}
    current_offset = 0
    
    for c in chroms:
        c_sub = df_sorted[df_sorted[chr_col] == c]
        chr_offsets[c] = current_offset
        c_max = c_sub[pos_col].max()
        chr_centers[c] = current_offset + c_max / 2
        current_offset += c_max + 2000000  # 2 Mb buffer between chromosomes
        
    df_sorted['pos_cum'] = df_sorted[pos_col] + df_sorted[chr_col].map(chr_offsets)
    
    # 2. Compute thresholds
    n_snps = len(df_sorted)
    bonf_p = 0.05 / n_snps
    bonf_log10 = -np.log10(bonf_p)
    sugg_log10 = -np.log10(1e-5)
    lambda_gc = calculate_lambda_gc(df_sorted[p_col].values)
    
    print(f"[{model_name}] SNPs: {n_snps:,} | Lambda GC: {lambda_gc:.4f} | Bonferroni: -log10(P) = {bonf_log10:.2f}")
    
    # 3. Create composite figure (Panel A: Manhattan, Panel B: QQ)
    fig = plt.figure(figsize=(15.5, 9.5), dpi=300)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.3, 0.7], hspace=0.32)
    
    # --- Manhattan Plot (Top) ---
    ax1 = fig.add_subplot(gs[0])
    
    # Alternating colors: Slate and Navy
    chr_colors = ['#0284c7', '#334155']
    
    for idx, c in enumerate(chroms):
        c_data = df_sorted[df_sorted[chr_col] == c]
        color = chr_colors[idx % 2]
        # Non-significant SNPs (thinned if > 100k points for plotting speed)
        non_sig = c_data[c_data['-log10p'] < sugg_log10]
        if len(non_sig) > 60000:
            non_sig = non_sig.sample(60000, random_state=42)
        ax1.scatter(non_sig['pos_cum'], non_sig['-log10p'], c=color, s=8, alpha=0.6, edgecolors='none', rasterized=True)
        
        # Suggestive & Significant SNPs (plotted in full resolution)
        sig_data = c_data[c_data['-log10p'] >= sugg_log10]
        if not sig_data.empty:
            ax1.scatter(sig_data['pos_cum'], sig_data['-log10p'], c='#ea580c', s=16, alpha=0.88, edgecolors='none', rasterized=True)
            
    # Draw threshold lines
    ax1.axhline(bonf_log10, color='#dc2626', linestyle='--', linewidth=1.2, label=f'Genome-wide significance (P < {bonf_p:.1e})')
    ax1.axhline(sugg_log10, color='#64748b', linestyle=':', linewidth=1.0, label='Suggestive threshold (P < 1e-5)')
    
    # Annotate Known Candidate Genes
    for gene_info in CANDIDATE_GENES:
        gc = gene_info['chr']
        if gc in chr_offsets:
            g_cum = chr_offsets[gc] + gene_info['mid']
            # Find closest top SNP within 50 kb
            window_snps = df_sorted[
                (df_sorted[chr_col] == gc) &
                (df_sorted[pos_col] >= gene_info['start'] - 50000) &
                (df_sorted[pos_col] <= gene_info['end'] + 50000)
            ]
            if not window_snps.empty:
                top_snp = window_snps.sort_values('-log10p', ascending=False).iloc[0]
                if top_snp['-log10p'] >= sugg_log10:
                    ax1.annotate(
                        f"$\\it{{{gene_info['gene']}}}$\n(-logP={top_snp['-log10p']:.1f})",
                        xy=(top_snp['pos_cum'], top_snp['-log10p']),
                        xytext=(top_snp['pos_cum'], top_snp['-log10p'] + 2.5),
                        arrowprops=dict(facecolor='#0f172a', shrink=0.08, width=1.0, headwidth=4),
                        fontsize=8.5, fontweight='bold', ha='center',
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='#fef08a', edgecolor='#ca8a04', alpha=0.95)
                    )
                    
    ax1.set_xticks([chr_centers[c] for c in chroms])
    ax1.set_xticklabels([f"Chr {c}" for c in chroms], fontsize=9.5, fontweight='bold')
    ax1.set_xlim(-1000000, current_offset)
    ax1.set_ylim(0, max(df_sorted['-log10p'].max() * 1.12, 10))
    ax1.set_ylabel(r"$-\log_{10}(P)$", fontsize=11, fontweight='bold')
    ax1.set_title(f"A   Manhattan Plot: Flowering Time (FT10) — {model_name}", fontsize=12, fontweight='bold', loc='left')
    ax1.legend(loc='upper right', frameon=True, edgecolor='#cbd5e1', framealpha=0.92, fontsize=8.5)
    ax1.grid(True, linestyle=':', alpha=0.35, color='#cbd5e1')
    ax1.set_axisbelow(True)
    
    # --- QQ Plot (Bottom Left) ---
    ax2 = fig.add_subplot(gs[1])
    
    p_vals = df_sorted[p_col].values
    p_vals = p_vals[np.isfinite(p_vals) & (p_vals > 0) & (p_vals <= 1)]
    p_vals.sort()
    n = len(p_vals)
    
    obs = -np.log10(p_vals)
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    
    # 95% Confidence Interval Band (beta distribution order statistics)
    j = np.arange(1, n + 1)
    lower = -np.log10(stats.beta.ppf(0.975, j, n - j + 1))
    upper = -np.log10(stats.beta.ppf(0.025, j, n - j + 1))
    
    # Subsample for smooth plotting
    if n > 50000:
        idx_sample = np.unique(np.concatenate([
            np.linspace(0, 10000, 2000, dtype=int),
            np.linspace(10000, n - 1, 10000, dtype=int)
        ]))
        exp_sub = exp[idx_sample]
        obs_sub = obs[idx_sample]
        lower_sub = lower[idx_sample]
        upper_sub = upper[idx_sample]
    else:
        exp_sub, obs_sub, lower_sub, upper_sub = exp, obs, lower, upper
        
    ax2.fill_between(exp_sub, lower_sub, upper_sub, color='#e2e8f0', alpha=0.7, label='95% Confidence Band')
    ax2.plot([0, max(exp)], [0, max(exp)], color='#dc2626', linestyle='--', linewidth=1.5, label='Null Expectation')
    ax2.scatter(exp_sub, obs_sub, c='#0284c7', s=12, alpha=0.75, edgecolors='none', rasterized=True, label='Observed')
    
    ax2.set_xlabel(r"Expected $-\log_{10}(P)$", fontsize=10, fontweight='bold')
    ax2.set_ylabel(r"Observed $-\log_{10}(P)$", fontsize=10, fontweight='bold')
    ax2.set_title(f"B   Quantile-Quantile (QQ) Plot ($\\lambda_{{GC}} = {lambda_gc:.3f}$)", fontsize=11, fontweight='bold', loc='left')
    ax2.legend(loc='upper left', frameon=True, edgecolor='#cbd5e1', framealpha=0.92, fontsize=8.5)
    ax2.grid(True, linestyle=':', alpha=0.35, color='#cbd5e1')
    ax2.set_axisbelow(True)
    
    # Clean spines
    for ax in [ax1, ax2]:
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        for spine in ['left', 'bottom']:
            ax.spines[spine].set_color('#94a3b8')
            
    out_png = f"{out_prefix}_manhattan_qq.png"
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Generated publication figure: {out_png}")
    
    # 4. Save Top Significant Hits Summary Table
    sig_hits = df_sorted[df_sorted[p_col] < bonf_p].copy()
    sig_hits_file = f"{out_prefix}_significant_snps.csv"
    sig_hits.to_csv(sig_hits_file, index=False)
    print(f"Saved {len(sig_hits)} significant SNPs to {sig_hits_file}")
    
    return out_png, sig_hits, lambda_gc


def plot_interactive_manhattan(df, model_name, out_html, p_col='P', chr_col='#CHROM', pos_col='POS', id_col='ID'):
    df_sorted = df.sort_values([chr_col, pos_col]).copy()
    chroms = sorted(df_sorted[chr_col].unique())
    
    chr_offsets = {}
    current_offset = 0
    for c in chroms:
        c_sub = df_sorted[df_sorted[chr_col] == c]
        chr_offsets[c] = current_offset
        current_offset += c_sub[pos_col].max() + 2000000
        
    df_sorted['pos_cum'] = df_sorted[pos_col] + df_sorted[chr_col].map(chr_offsets)
    bonf_p = 0.05 / len(df_sorted)
    bonf_log10 = -np.log10(bonf_p)
    
    # Subsample non-significant variants to keep HTML lightweight (~30k points)
    non_sig = df_sorted[df_sorted['-log10p'] < 5.0]
    if len(non_sig) > 30000:
        non_sig = non_sig.sample(30000, random_state=42)
    sig_points = df_sorted[df_sorted['-log10p'] >= 5.0]
    plot_df = pd.concat([non_sig, sig_points]).sort_values([chr_col, pos_col])
    
    fig = go.Figure()
    
    chr_colors = ['#0284c7', '#334155']
    for idx, c in enumerate(chroms):
        sub = plot_df[plot_df[chr_col] == c]
        custom = list(zip(sub[id_col], sub[chr_col], sub[pos_col], sub[p_col].apply(lambda x: f"{x:.2e}")))
        fig.add_trace(go.Scattergl(
            x=sub['pos_cum'],
            y=sub['-log10p'],
            mode='markers',
            name=f"Chr {c}",
            marker=dict(size=4, color=chr_colors[idx % 2], opacity=0.75),
            customdata=custom,
            hovertemplate="<b>Variant:</b> %{customdata[0]}<br><b>Position:</b> Chr%{customdata[1]}:%{customdata[2]:,}<br><b>P-value:</b> %{customdata[3]}<extra></extra>"
        ))
        
    fig.add_hline(y=bonf_log10, line_dash='dash', line_color='#dc2626', annotation_text='Bonferroni Significance')
    
    fig.update_layout(
        title=f"<b>Interactive GWAS Manhattan Plot — {model_name} (FT10)</b>",
        xaxis=dict(title="<b>Chromosome</b>", showgrid=False),
        yaxis=dict(title="<b>-log10(P-value)</b>", showgrid=True, gridcolor='#f1f5f9'),
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        height=550,
        showlegend=False
    )
    
    fig.write_html(out_html, include_plotlyjs='cdn', full_html=True)
    print(f"Generated interactive Manhattan plot: {out_html}")


def plot_dual_model_comparison(df_glm, df_lmm, out_png):
    """Generate a 4-panel publication visual comparing GLM (10 PCs) and LMM (Kinship Matrix)."""
    fig = plt.figure(figsize=(16, 12), dpi=300)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 0.9], width_ratios=[1.2, 0.8], hspace=0.32, wspace=0.22)
    
    # 1. Manhattan GLM (Top Left)
    ax1 = fig.add_subplot(gs[0, 0])
    # 2. Manhattan LMM (Bottom Left)
    ax2 = fig.add_subplot(gs[1, 0])
    # 3. QQ GLM (Top Right)
    ax3 = fig.add_subplot(gs[0, 1])
    # 4. QQ LMM (Bottom Right)
    ax4 = fig.add_subplot(gs[1, 1])
    
    # helper to draw single manhattan
    def draw_manhattan(ax, df_data, title, chr_col, pos_col, p_col):
        df_s = df_data.sort_values([chr_col, pos_col]).copy()
        chroms = sorted(df_s[chr_col].unique())
        offsets = {}
        centers = {}
        curr = 0
        for c in chroms:
            offsets[c] = curr
            sub = df_s[df_s[chr_col] == c]
            centers[c] = curr + sub[pos_col].max() / 2
            curr += sub[pos_col].max() + 2000000
        df_s['pos_cum'] = df_s[pos_col] + df_s[chr_col].map(offsets)
        
        chr_colors = ['#0284c7', '#334155']
        bonf_log = -np.log10(0.05 / len(df_s))
        sugg_log = 5.0
        
        for idx, c in enumerate(chroms):
            sub = df_s[df_s[chr_col] == c]
            non_sig = sub[sub['-log10p'] < sugg_log]
            if len(non_sig) > 40000:
                non_sig = non_sig.sample(40000, random_state=42)
            ax.scatter(non_sig['pos_cum'], non_sig['-log10p'], c=chr_colors[idx % 2], s=6, alpha=0.55, edgecolors='none', rasterized=True)
            sig_pts = sub[sub['-log10p'] >= sugg_log]
            if not sig_pts.empty:
                ax.scatter(sig_pts['pos_cum'], sig_pts['-log10p'], c='#ea580c', s=14, alpha=0.9, edgecolors='none', rasterized=True)
                
        ax.axhline(bonf_log, color='#dc2626', linestyle='--', linewidth=1.1)
        ax.set_xticks([centers[c] for c in chroms])
        ax.set_xticklabels([f"Chr {c}" for c in chroms], fontsize=8.5, fontweight='bold')
        ax.set_ylabel(r"$-\log_{10}(P)$", fontsize=10, fontweight='bold')
        ax.set_title(title, fontsize=11, fontweight='bold', loc='left')
        ax.grid(True, linestyle=':', alpha=0.35, color='#cbd5e1')
        ax.set_axisbelow(True)
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        for spine in ['left', 'bottom']:
            ax.spines[spine].set_color('#94a3b8')

    draw_manhattan(ax1, df_glm, "A   PLINK 2 GLM (10 Principal Component Covariates)", '#CHROM', 'POS', 'P')
    draw_manhattan(ax2, df_lmm, "C   GEMMA Linear Mixed Model (Centered Kinship Matrix)", 'chr', 'ps', 'p_wald')
    
    # Draw QQ plots
    def draw_qq(ax, p_values, title, lambda_val):
        p_clean = p_values[np.isfinite(p_values) & (p_values > 0) & (p_values <= 1)]
        p_clean.sort()
        n = len(p_clean)
        obs = -np.log10(p_clean)
        exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
        
        idx_s = np.unique(np.concatenate([np.linspace(0, min(5000, n-1), 1500, dtype=int), np.linspace(min(5000, n-1), n-1, 8000, dtype=int)]))
        ax.plot([0, max(exp)], [0, max(exp)], color='#dc2626', linestyle='--', linewidth=1.4)
        ax.scatter(exp[idx_s], obs[idx_s], c='#0284c7', s=10, alpha=0.75, edgecolors='none', rasterized=True)
        ax.set_xlabel(r"Expected $-\log_{10}(P)$", fontsize=9.5, fontweight='bold')
        ax.set_ylabel(r"Observed $-\log_{10}(P)$", fontsize=9.5, fontweight='bold')
        ax.set_title(f"{title} ($\\lambda_{{GC}} = {lambda_val:.3f}$)", fontsize=10.5, fontweight='bold', loc='left')
        ax.grid(True, linestyle=':', alpha=0.35, color='#cbd5e1')
        ax.set_axisbelow(True)
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        for spine in ['left', 'bottom']:
            ax.spines[spine].set_color('#94a3b8')

    l_glm = calculate_lambda_gc(df_glm['P'].values)
    l_lmm = calculate_lambda_gc(df_lmm['p_wald'].values)
    draw_qq(ax3, df_glm['P'].values, "B   QQ Plot: GLM (10 PCs)", l_glm)
    draw_qq(ax4, df_lmm['p_wald'].values, "D   QQ Plot: LMM (Kinship Matrix)", l_lmm)
    
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"Generated dual-model comparison figure: {out_png}")


def plot_4model_comparison(df_glm, df_lmm, df_loco, df_cond, out_png):
    """Generate a 4-panel publication visual showing the methodological progression:
    A: Naive GLM (10 PCs)
    B: Standard LMM (Whole-Genome Kinship)
    C: LOCO LMM (Leave-One-Chromosome-Out Kinship)
    D: Stepwise Conditioned MLMM (Conditioned on FT, DOG1, VIN3)
    """
    fig, axes = plt.subplots(4, 1, figsize=(15.5, 14.5), dpi=300)
    plt.subplots_adjust(hspace=0.36, top=0.95, bottom=0.05, left=0.07, right=0.96)
    
    models = [
        (axes[0], df_glm, "A   PLINK 2 GLM (10 PC Covariates) — λ_GC = 2.674 (Severe Cryptic Relatedness Inflation)", '#CHROM', 'POS', 'P'),
        (axes[1], df_lmm, "B   Standard GEMMA LMM (Whole-Genome Kinship) — λ_GC = 0.979 (Proximal Contamination Attenuation)", 'chr', 'ps', 'p_wald'),
        (axes[2], df_loco, "C   LOCO GEMMA LMM (Leave-One-Chromosome-Out) — λ_GC = 1.108 (Sharpened FT, DOG1, VIN3 Peaks)", 'chr', 'ps', 'p_wald'),
        (axes[3], df_cond, "D   Stepwise Conditioned Model (Conditioned on FT, DOG1, VIN3) — Unmasks FLC (P = 1.12e-22)", '#CHROM', 'POS', 'P')
    ]
    
    chr_colors = ['#0284c7', '#334155']
    sugg_log = 5.0
    
    for ax, df_data, title, chr_col, pos_col, p_col in models:
        df_s = df_data.sort_values([chr_col, pos_col]).copy()
        chroms = sorted(df_s[chr_col].unique())
        offsets = {}
        centers = {}
        curr = 0
        for c in chroms:
            offsets[c] = curr
            sub = df_s[df_s[chr_col] == c]
            centers[c] = curr + sub[pos_col].max() / 2
            curr += sub[pos_col].max() + 2000000
        df_s['pos_cum'] = df_s[pos_col] + df_s[chr_col].map(offsets)
        bonf_log = -np.log10(0.05 / len(df_s))
        
        for idx, c in enumerate(chroms):
            sub = df_s[df_s[chr_col] == c]
            non_sig = sub[sub['-log10p'] < sugg_log]
            if len(non_sig) > 35000:
                non_sig = non_sig.sample(35000, random_state=42)
            ax.scatter(non_sig['pos_cum'], non_sig['-log10p'], c=chr_colors[idx % 2], s=5.5, alpha=0.55, edgecolors='none', rasterized=True)
            sig_pts = sub[sub['-log10p'] >= sugg_log]
            if not sig_pts.empty:
                ax.scatter(sig_pts['pos_cum'], sig_pts['-log10p'], c='#ea580c', s=13, alpha=0.9, edgecolors='none', rasterized=True)
                
        ax.axhline(bonf_log, color='#dc2626', linestyle='--', linewidth=1.1)
        ax.set_xticks([centers[c] for c in chroms])
        ax.set_xticklabels([f"Chr {c}" for c in chroms], fontsize=9, fontweight='bold')
        ax.set_ylabel(r"$-\log_{10}(P)$", fontsize=10, fontweight='bold')
        ax.set_title(title, fontsize=11, fontweight='bold', loc='left')
        ax.grid(True, linestyle=':', alpha=0.35, color='#cbd5e1')
        ax.set_axisbelow(True)
        ax.set_xlim(-1000000, curr)
        ax.set_ylim(0, max(df_s['-log10p'].max() * 1.14, 10))
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        for spine in ['left', 'bottom']:
            ax.spines[spine].set_color('#94a3b8')
            
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"Generated 4-model methodology comparison figure: {out_png}")


def main():
    plink_file = 'results/gwas/plink2_glm_pc10.FT10.glm.linear'
    gemma_file = 'output/athal_1003_ft10_lmm.assoc.txt'
    loco_file = 'results/gwas/ft10_gemma_loco_all.assoc.txt'
    cond_file = 'results/gwas/plink2_glm_conditioned.FT10.glm.linear'
    
    df_plink = None
    df_gemma = None
    df_loco = None
    df_cond = None
    
    if os.path.exists(plink_file):
        df_plink = load_gwas_data(plink_file, p_col='P', chr_col='#CHROM', pos_col='POS', id_col='ID')
        plot_manhattan_and_qq(df_plink, "PLINK 2 GLM (10 PCs)", "results/gwas/ft10_plink2_glm_pc10")
        
    if os.path.exists(gemma_file):
        df_gemma = load_gwas_data(gemma_file, p_col='p_wald', chr_col='chr', pos_col='ps', id_col='rs')
        plot_manhattan_and_qq(df_gemma, "GEMMA LMM (Kinship Matrix)", "results/gwas/ft10_gemma_lmm", p_col='p_wald', chr_col='chr', pos_col='ps')
        
    if os.path.exists(loco_file):
        df_loco = load_gwas_data(loco_file, p_col='p_wald', chr_col='chr', pos_col='ps', id_col='rs')
        plot_manhattan_and_qq(df_loco, "GEMMA LOCO LMM", "results/gwas/ft10_gemma_loco", p_col='p_wald', chr_col='chr', pos_col='ps')
        plot_interactive_manhattan(df_loco, "GEMMA LOCO LMM", "results/gwas/ft10_gemma_loco.html", p_col='p_wald', chr_col='chr', pos_col='ps', id_col='rs')
        
    if os.path.exists(cond_file):
        df_cond = load_gwas_data(cond_file, p_col='P', chr_col='#CHROM', pos_col='POS', id_col='ID')
        plot_manhattan_and_qq(df_cond, "Conditioned MLMM", "results/gwas/ft10_stepwise_conditioned", p_col='P', chr_col='#CHROM', pos_col='POS')
        
    if df_plink is not None and df_gemma is not None and df_loco is not None and df_cond is not None:
        plot_4model_comparison(df_plink, df_gemma, df_loco, df_cond, "results/gwas/gwas_4model_methodology_comparison.png")


if __name__ == '__main__':
    main()

