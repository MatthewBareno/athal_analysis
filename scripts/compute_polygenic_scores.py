"""Compute polygenic score assignments for each of the 1,003 Arabidopsis thaliana accessions:
1. Whole-Genome GBLUP (Kinship-based genomic estimated breeding values, 5-fold CV and full sample).
2. Major-effect 4-QTL additive model (FT, DOG1, VIN3, FLC).
3. Hybrid multi-component model (4 QTLs + Kinship polygenic residual).
4. Generates individual accession score table (CSV) and publication figures (PNG & HTML).
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import KFold
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_absolute_error
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import plotly.graph_objects as go
import plotly.express as px

GROUP_CONFIG = {
    'relict': {'label': 'Relict', 'color': '#dc2626'},
    'spain': {'label': 'Spain', 'color': '#ea580c'},
    'western_europe': {'label': 'Western Europe', 'color': '#4338ca'},
    'germany': {'label': 'Germany', 'color': '#047857'},
    'central_europe': {'label': 'Central Europe', 'color': '#16a34a'},
    'italy_balkan_caucasus': {'label': 'Italy / Balkans', 'color': '#9333ea'},
    'south_sweden': {'label': 'South Sweden', 'color': '#2563eb'},
    'north_sweden': {'label': 'North Sweden', 'color': '#06b6d4'},
    'asia': {'label': 'Asia', 'color': '#db2777'},
    'admixed': {'label': 'Admixed', 'color': '#64748b'}
}

ORDERED_GROUPS = [
    'relict', 'spain', 'western_europe', 'germany', 'central_europe',
    'italy_balkan_caucasus', 'south_sweden', 'north_sweden', 'asia', 'admixed'
]


def run_pgs_pipeline():
    os.makedirs('results/prs', exist_ok=True)
    os.makedirs('data/processed', exist_ok=True)

    # 1. Load phenotype & metadata
    pheno_file = 'data/processed/phenotype_clean.csv'
    df_pheno = pd.read_csv(pheno_file)
    y = df_pheno['FT10'].values
    n = len(y)

    # 2. Load kinship matrix
    kinship_file = 'output/athal_1003_kinship.cXX.txt'
    K = np.loadtxt(kinship_file)

    # From GEMMA REMLE null model: Vg = 607.317, Ve = 20.405
    vg = 607.317
    ve = 20.405
    lambda_reg = ve / vg

    # 3. Load 4-QTL dosages (FT, DOG1, VIN3, FLC)
    qtl_file = 'data/processed/lead_4qtl_dosages.raw'
    df_qtl = pd.read_csv(qtl_file, sep=r'\s+')
    qtl_cols = [c for c in df_qtl.columns if ':' in c]
    for c in qtl_cols:
        df_qtl[c] = df_qtl[c].fillna(df_qtl[c].mean())
    X_qtl = df_qtl[qtl_cols].values

    # 4. Fit 5-Fold Cross-Validation for Out-Of-Sample Evaluation
    kf = KFold(n_splits=5, shuffle=True, random_state=42)

    y_pred_gblup_cv = np.zeros(n)
    y_pred_qtl_cv = np.zeros(n)
    y_pred_hybrid_cv = np.zeros(n)

    for train_idx, test_idx in kf.split(K):
        # A. 4-QTL Linear Regression
        lr = LinearRegression()
        lr.fit(X_qtl[train_idx], y[train_idx])
        y_pred_qtl_cv[test_idx] = lr.predict(X_qtl[test_idx])

        # B. Whole-Genome GBLUP
        K_tr_tr = K[np.ix_(train_idx, train_idx)]
        K_te_tr = K[np.ix_(test_idx, train_idx)]
        y_tr = y[train_idx]
        mu = np.mean(y_tr)
        y_tr_c = y_tr - mu

        inv = np.linalg.solve(K_tr_tr + lambda_reg * np.eye(len(train_idx)), y_tr_c)
        y_pred_gblup_cv[test_idx] = mu + K_te_tr @ inv

        # C. Hybrid Model (4 QTLs fixed + GBLUP on residuals)
        res_tr = y_tr - lr.predict(X_qtl[train_idx])
        inv_res = np.linalg.solve(K_tr_tr + lambda_reg * np.eye(len(train_idx)), res_tr)
        y_pred_hybrid_cv[test_idx] = lr.predict(X_qtl[test_idx]) + K_te_tr @ inv_res

    # 5. Full-Sample Genomic Estimated Breeding Values (GEBV)
    mu_full = np.mean(y)
    inv_full = np.linalg.solve(K + lambda_reg * np.eye(n), y - mu_full)
    gebv_u = K @ inv_full
    y_pred_gblup_full = mu_full + gebv_u

    lr_full = LinearRegression().fit(X_qtl, y)
    y_pred_qtl_full = lr_full.predict(X_qtl)

    # 6. Build Individual Accession Table
    df_scores = df_pheno[['accession_id', 'name', 'country', 'admixture_group', 'latitude', 'longitude', 'elevation']].copy()
    df_scores['FT10_observed'] = y.round(2)
    df_scores['PGS_gblup_cv'] = y_pred_gblup_cv.round(2)
    df_scores['PGS_gblup_full'] = y_pred_gblup_full.round(2)
    df_scores['PGS_4qtl_cv'] = y_pred_qtl_cv.round(2)
    df_scores['PGS_hybrid_cv'] = y_pred_hybrid_cv.round(2)
    df_scores['gblup_residual'] = (y - y_pred_gblup_cv).round(2)

    out_csv = 'data/processed/athal_1003_polygenic_scores.csv'
    df_scores.to_csv(out_csv, index=False)
    print(f"Saved individual accession polygenic scores to {out_csv}")

    # 7. Compute Model Performance Metrics
    metrics = []
    models_dict = {
        '4-QTL Additive Score (FT, DOG1, VIN3, FLC)': y_pred_qtl_cv,
        'Whole-Genome GBLUP (Kinship Matrix)': y_pred_gblup_cv,
        'Hybrid Model (4 QTLs + Polygenic Residual)': y_pred_hybrid_cv
    }

    for m_name, y_pred in models_dict.items():
        r2 = r2_score(y, y_pred)
        corr = stats.pearsonr(y, y_pred)[0]
        mae = mean_absolute_error(y, y_pred)
        metrics.append({
            'Model': m_name,
            'CV_R2': round(r2, 4),
            'CV_Pearson_r': round(corr, 4),
            'CV_MAE_days': round(mae, 2)
        })

    df_metrics = pd.DataFrame(metrics)
    out_metrics = 'results/prs/pgs_model_metrics.csv'
    df_metrics.to_csv(out_metrics, index=False)
    print("\nModel Evaluation (5-Fold Cross-Validation):")
    print(df_metrics.to_string(index=False))

    # 8. Create Publication Figure (3 Panels)
    plot_pgs_figures(df_scores, df_metrics)
    plot_interactive_pgs(df_scores)


def plot_pgs_figures(df_scores, df_metrics):
    fig = plt.figure(figsize=(16, 6.8), dpi=300)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.15, 1.0], wspace=0.22, left=0.06, right=0.96, top=0.88, bottom=0.12)

    # Panel A: Observed vs Predicted Scatter Plot
    ax1 = fig.add_subplot(gs[0])
    y_obs = df_scores['FT10_observed'].values
    y_pred = df_scores['PGS_gblup_cv'].values

    for grp in ORDERED_GROUPS:
        sub = df_scores[df_scores['admixture_group'] == grp]
        if sub.empty:
            continue
        cfg = GROUP_CONFIG.get(grp, {'label': grp, 'color': '#7f7f7f'})
        ax1.scatter(
            sub['PGS_gblup_cv'], sub['FT10_observed'],
            color=cfg['color'], label=f"{cfg['label']}",
            s=26, alpha=0.82, edgecolors='#0f172a', linewidth=0.3
        )

    # Linear fit and identity line
    slope, intercept, r_val, p_val, std_err = stats.linregress(y_pred, y_obs)
    x_vals = np.linspace(min(y_pred) - 2, max(y_pred) + 2, 100)
    ax1.plot(x_vals, x_vals, color='#94a3b8', linestyle='--', linewidth=1.2, label='Identity Line ($y = x$)')
    ax1.plot(x_vals, slope * x_vals + intercept, color='#dc2626', linewidth=1.8, label=f'Linear Fit ($R^2 = {r_val**2:.3f}$)')

    # Annotation badge
    mae = mean_absolute_error(y_obs, y_pred)
    ax1.text(
        0.05, 0.93,
        f"$\\mathbf{{5\\text{{-Fold Cross-Validation}}}}\n"
        f"R^2 = \\mathbf{{{r_val**2:.3f}}}\n"
        f"\\text{{Pearson }} r = \\mathbf{{{r_val:.3f}}}\n"
        f"\\text{{MAE}} = \\mathbf{{{mae:.2f}\\text{{ days}}}}$",
        transform=ax1.transAxes,
        fontsize=9.5, fontweight='bold', va='top',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#f8fafc', edgecolor='#cbd5e1', lw=1.0)
    )

    ax1.set_xlabel('Predicted Flowering Time: GBLUP (Days)', fontsize=10.5, fontweight='bold', labelpad=6)
    ax1.set_ylabel('Observed Flowering Time: FT10 (Days)', fontsize=10.5, fontweight='bold', labelpad=6)
    ax1.set_title('A   Out-of-Sample Genomic Prediction (GBLUP Kinship Model)', fontsize=11.5, fontweight='bold', loc='left')
    ax1.grid(True, linestyle=':', alpha=0.4, color='#cbd5e1')
    ax1.set_axisbelow(True)
    ax1.legend(loc='lower right', ncol=2, fontsize=7.5, frameon=True, edgecolor='#cbd5e1', framealpha=0.92)

    for spine in ['top', 'right']:
        ax1.spines[spine].set_visible(False)
    for spine in ['left', 'bottom']:
        ax1.spines[spine].set_color('#94a3b8')

    # Panel B: Polygenic Burden across Geographic Lineages
    ax2 = fig.add_subplot(gs[1])
    
    # Order groups by median predicted PGS
    box_data = []
    box_labels = []
    box_colors = []
    
    for grp in ORDERED_GROUPS:
        sub = df_scores[df_scores['admixture_group'] == grp]
        if not sub.empty:
            box_data.append(sub['PGS_gblup_cv'].values)
            box_labels.append(GROUP_CONFIG[grp]['label'])
            box_colors.append(GROUP_CONFIG[grp]['color'])

    bplot = ax2.boxplot(
        box_data,
        patch_artist=True,
        tick_labels=box_labels,
        medianprops=dict(color='#0f172a', linewidth=1.5),
        whiskerprops=dict(color='#64748b', linewidth=1.1),
        capprops=dict(color='#64748b', linewidth=1.1),
        flierprops=dict(marker='o', markersize=3, alpha=0.5, markerfacecolor='#94a3b8')
    )

    for patch, color in zip(bplot['boxes'], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
        patch.set_edgecolor('#0f172a')
        patch.set_linewidth(0.8)

    ax2.set_xticklabels(box_labels, rotation=35, ha='right', fontsize=8.5, fontweight='bold')
    ax2.set_ylabel('Predicted Genetic Breeding Value (Days)', fontsize=10.5, fontweight='bold', labelpad=6)
    ax2.set_title('B   Polygenic Flowering Burden by Geographic Lineage', fontsize=11.5, fontweight='bold', loc='left')
    ax2.grid(True, linestyle=':', alpha=0.4, color='#cbd5e1')
    ax2.set_axisbelow(True)

    for spine in ['top', 'right']:
        ax2.spines[spine].set_visible(False)
    for spine in ['left', 'bottom']:
        ax2.spines[spine].set_color('#94a3b8')

    out_png = 'results/prs/polygenic_score_evaluation.png'
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"Generated polygenic evaluation figure: {out_png}")


def plot_interactive_pgs(df_scores):
    fig = px.scatter(
        df_scores,
        x='PGS_gblup_cv',
        y='FT10_observed',
        color='admixture_group',
        color_discrete_map={k: v['color'] for k, v in GROUP_CONFIG.items()},
        hover_name='name',
        hover_data=['accession_id', 'country', 'admixture_group', 'PGS_gblup_cv', 'FT10_observed', 'gblup_residual'],
        labels={
            'PGS_gblup_cv': 'Predicted Flowering Time (GBLUP Days)',
            'FT10_observed': 'Observed Flowering Time (FT10 Days)',
            'admixture_group': 'Lineage'
        },
        title="<b>Arabidopsis thaliana Polygenic Score Validation (5-Fold CV GBLUP)</b>"
    )

    # Add identity line
    min_val = min(df_scores['PGS_gblup_cv'].min(), df_scores['FT10_observed'].min())
    max_val = max(df_scores['PGS_gblup_cv'].max(), df_scores['FT10_observed'].max())
    fig.add_shape(
        type='line',
        x0=min_val, y0=min_val, x1=max_val, y1=max_val,
        line=dict(color='#94a3b8', dash='dash', width=1.5)
    )

    fig.update_layout(
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        height=600,
        font=dict(family='sans-serif')
    )

    out_html = 'results/prs/polygenic_score_evaluation.html'
    fig.write_html(out_html, include_plotlyjs='cdn', full_html=True)
    print(f"Generated interactive polygenic score visual: {out_html}")


if __name__ == '__main__':
    run_pgs_pipeline()
