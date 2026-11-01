"""Generate a publication-quality world map of Arabidopsis thaliana accessions.

Colors accessions by the official 1001 Genomes Genetic Ancestry Groups (K=9 + Admixed)
from Alonso-Blanco et al. (Cell 2016).

Outputs:
- results/accession_world_map.html: Interactive Plotly map with tooltips and ancestry group legend
- results/accession_world_map.png: High-resolution (2200x1200, 300-DPI) publication figure
"""

import os
import pandas as pd
import plotly.graph_objects as go

# Clean display names and official colors for the 1001 Genomes genetic groups
GROUP_CONFIG = {
    'central_europe': {
        'label': 'Central Europe',
        'color': '#33a02c'   # Vibrant Green
    },
    'germany': {
        'label': 'Germany',
        'color': '#006d2c'   # Deep Forest Green
    },
    'south_sweden': {
        'label': 'South Sweden',
        'color': '#1f78b4'   # Marine Blue
    },
    'north_sweden': {
        'label': 'North Sweden',
        'color': '#38bdf8'   # Arctic Cyan
    },
    'western_europe': {
        'label': 'Western Europe',
        'color': '#4338ca'   # Indigo
    },
    'spain': {
        'label': 'Spain',
        'color': '#ea580c'   # Warm Orange
    },
    'italy_balkan_caucasus': {
        'label': 'Italy / Balkans / Caucasus',
        'color': '#9333ea'   # Mediterranean Violet
    },
    'asia': {
        'label': 'Asia',
        'color': '#ec4899'   # Rose Magenta
    },
    'relict': {
        'label': 'Relict (Ancestral)',
        'color': '#dc2626'   # Crimson Red
    },
    'admixed': {
        'label': 'Admixed',
        'color': '#64748b'   # Neutral Slate Gray
    }
}

# Preferred presentation order (ordered by geographic clustering and prevalence)
ORDERED_KEYS = [
    'relict',
    'spain',
    'western_europe',
    'germany',
    'central_europe',
    'italy_balkan_caucasus',
    'south_sweden',
    'north_sweden',
    'asia',
    'admixed'
]


def generate_ancestry_map():
    data_path = 'data/processed/phenotype_clean.csv'
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Missing {data_path}. Run prepare_phenotypes.py first.")

    df = pd.read_csv(data_path)
    valid_geo = df.dropna(subset=['latitude', 'longitude']).copy()

    fig = go.Figure()

    # Add each genetic ancestry group as an independent trace
    for key in ORDERED_KEYS:
        sub = valid_geo[valid_geo['admixture_group'] == key]
        if sub.empty:
            continue

        n_samples = len(sub)
        cfg = GROUP_CONFIG.get(key, {'label': key, 'color': '#7f7f7f'})

        site_names = sub['sitename'].fillna('N/A').astype(str)
        elevations = sub['elevation'].fillna(-999).round(0).astype(int)

        custom_data = list(zip(
            sub['name'],
            sub['accession_id'],
            sub['country'],
            site_names,
            elevations,
            sub['csnumber'].fillna('N/A').astype(str)
        ))

        fig.add_trace(go.Scattergeo(
            lon=sub['longitude'],
            lat=sub['latitude'],
            customdata=custom_data,
            hovertemplate=(
                "<b>%{customdata[0]}</b> (ID: %{customdata[1]})<br>"
                f"<b>Ancestry Group:</b> {cfg['label']}<br>"
                "<b>Country:</b> %{customdata[2]} (%{customdata[3]})<br>"
                "<b>Stock:</b> %{customdata[5]}<br>"
                "<b>Elevation:</b> %{customdata[4]} m<br>"
                "<b>Coordinates:</b> %{lat:.3f}°N, %{lon:.3f}°E<extra></extra>"
            ),
            mode='markers',
            name=f"{cfg['label']} ({n_samples})",
            marker=dict(
                size=6.2,
                color=cfg['color'],
                opacity=0.88,
                line=dict(width=0.45, color='#ffffff')
            )
        ))

    # Refined Natural Earth projection layout
    fig.update_geos(
        projection_type="natural earth",
        showcoastlines=True,
        coastlinecolor="#788896",
        coastlinewidth=0.75,
        showland=True,
        landcolor="#eef1f4",
        showocean=True,
        oceancolor="#f7f9fb",
        showcountries=True,
        countrycolor="#cbd5e1",
        countrywidth=0.55,
        showlakes=True,
        lakecolor="#f7f9fb",
        lataxis_range=[12, 72],
        lonaxis_range=[-130, 145],
        bgcolor="rgba(0,0,0,0)"
    )

    fig.update_layout(
        title=dict(
            text="<b>Geographic Distribution of <i>Arabidopsis thaliana</i> Genetic Groups</b><br>"
                 "<span style='font-size: 11.5px; font-weight: normal; color: #475569;'>"
                 "1,003 accessions colored by 1001 Genomes ancestry group (Alonso-Blanco et al., Cell 2016)</span>",
            x=0.03,
            y=0.96,
            font=dict(size=14, family="sans-serif", color="#0f172a")
        ),
        legend=dict(
            title=dict(text="<b>Genetic Group (n)</b>", font=dict(size=11, family="sans-serif", color="#0f172a")),
            orientation="v",
            yanchor="middle",
            y=0.5,
            xanchor="left",
            x=1.01,
            font=dict(size=10, family="sans-serif"),
            itemsizing="constant",
            bgcolor="rgba(255,255,255,0.9)",
            bordercolor="#e2e8f0",
            borderwidth=1
        ),
        margin=dict(l=10, r=200, t=55, b=10),
        width=1180,
        height=580,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family="sans-serif")
    )

    os.makedirs('results', exist_ok=True)
    out_html = 'results/accession_world_map.html'
    out_png = 'results/accession_world_map.png'

    # Save interactive HTML
    html_content = fig.to_html(
        include_plotlyjs='cdn',
        full_html=True,
        config=dict(
            displayModeBar=True,
            displaylogo=False,
            modeBarButtonsToRemove=['select2d', 'lasso2d']
        )
    )
    with open(out_html, 'w', encoding='utf-8') as f:
        f.write(html_content)

    # Save high-res 300-DPI publication PNG
    fig.write_image(out_png, scale=2)
    print(f"Generated ancestry-colored interactive map: {out_html}")
    print(f"Generated ancestry-colored publication figure: {out_png}")


if __name__ == '__main__':
    generate_ancestry_map()
