"""Generate a publication-quality world map of Arabidopsis thaliana accessions.

Visualizes the sampling locations of the 1,003 accessions across the native
and introduced range without phenotype stratification.
"""

import os
import pandas as pd
import plotly.graph_objects as go


def generate_accession_map():
    data_path = 'data/processed/phenotype_clean.csv'
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Missing {data_path}. Run prepare_phenotypes.py first.")

    df = pd.read_csv(data_path)
    # Filter to accessions with valid coordinates
    valid_geo = df.dropna(subset=['latitude', 'longitude']).copy()
    valid_geo = valid_geo.sort_values('country')

    # Clean metadata for tooltips
    site_names = valid_geo['sitename'].fillna('N/A').astype(str)
    elevations = valid_geo['elevation'].fillna(-999).round(0).astype(int)

    custom_data = list(zip(
        valid_geo['name'],
        valid_geo['accession_id'],
        valid_geo['country'],
        site_names,
        elevations,
        valid_geo['csnumber'].fillna('N/A').astype(str)
    ))

    fig = go.Figure()

    # Scattergeo trace for accession collection sites
    fig.add_trace(go.Scattergeo(
        lon=valid_geo['longitude'],
        lat=valid_geo['latitude'],
        customdata=custom_data,
        hovertemplate=(
            "<b>%{customdata[0]}</b> (ID: %{customdata[1]})<br>"
            "<b>Country:</b> %{customdata[2]}<br>"
            "<b>Site:</b> %{customdata[3]}<br>"
            "<b>Stock:</b> %{customdata[5]}<br>"
            "<b>Coordinates:</b> %{lat:.3f}°N, %{lon:.3f}°E<extra></extra>"
        ),
        mode='markers',
        marker=dict(
            size=5.5,
            color='#1b4d3e',  # Professional deep botanical green
            opacity=0.82,
            line=dict(width=0.6, color='#ffffff')
        ),
        name="Accession"
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
            text="<b>Geographic Distribution of <i>Arabidopsis thaliana</i> Accessions</b><br>"
                 "<span style='font-size: 11.5px; font-weight: normal; color: #475569;'>"
                 f"Collection sites for 1,003 1001-Genomes accessions (999 geo-referenced sites)</span>",
            x=0.03,
            y=0.96,
            font=dict(size=14, family="sans-serif", color="#0f172a")
        ),
        showlegend=False,
        margin=dict(l=10, r=10, t=55, b=10),
        width=1100,
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
    print(f"Generated clean accession map: {out_html}")
    print(f"Generated publication figure: {out_png}")


if __name__ == '__main__':
    generate_accession_map()
