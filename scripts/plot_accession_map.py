"""Generate a professional world map of Arabidopsis thaliana 1001 Genomes accessions.

Produces:
- results/accession_world_map.html: Interactive Plotly map with tooltips
- results/accession_world_map.png: High-resolution publication-quality figure (2200x1160)
"""

import os
import pandas as pd
import plotly.graph_objects as go


def generate_accession_maps():
    data_path = 'data/processed/phenotype_clean.csv'
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Missing {data_path}. Run prepare_phenotypes.py first.")

    df = pd.read_csv(data_path)
    # Filter to accessions with valid geographic coordinates
    valid_geo = df.dropna(subset=['latitude', 'longitude']).copy()
    valid_geo = valid_geo.sort_values('FT10')

    # Compact customdata for hovertemplate
    site_names = valid_geo['sitename'].fillna('N/A').astype(str)
    elevations = valid_geo['elevation'].fillna(-999).round(0).astype(int)

    custom_data = list(zip(
        valid_geo['name'],
        valid_geo['accession_id'],
        valid_geo['country'],
        valid_geo['FT10'].round(1),
        elevations,
        site_names
    ))

    fig = go.Figure()

    # Scattergeo trace
    fig.add_trace(go.Scattergeo(
        lon=valid_geo['longitude'],
        lat=valid_geo['latitude'],
        customdata=custom_data,
        hovertemplate=(
            "<b>%{customdata[0]}</b> (#%{customdata[1]})<br>"
            "<b>Country / Site:</b> %{customdata[2]} (%{customdata[5]})<br>"
            "<b>FT10 Flowering Time:</b> %{customdata[3]} days<br>"
            "<b>Coordinates:</b> %{lat:.2f}°N, %{lon:.2f}°E<extra></extra>"
        ),
        mode='markers',
        marker=dict(
            size=5.5,
            color=valid_geo['FT10'],
            colorscale='Viridis',
            cmin=valid_geo['FT10'].quantile(0.02),
            cmax=valid_geo['FT10'].quantile(0.98),
            colorbar=dict(
                title=dict(
                    text="<b>FT10 (Days)</b>",
                    font=dict(size=11, family="sans-serif")
                ),
                thickness=13,
                len=0.72,
                x=0.98,
                y=0.5,
                outlinewidth=0,
                tickfont=dict(size=10, family="sans-serif")
            ),
            opacity=0.88,
            line=dict(width=0.45, color='rgba(0, 0, 0, 0.7)')
        )
    ))

    # Professional Natural Earth projection layout
    fig.update_geos(
        projection_type="natural earth",
        showcoastlines=True,
        coastlinecolor="rgba(125, 135, 145, 0.6)",
        coastlinewidth=0.7,
        showland=True,
        landcolor="#f4f6f8",
        showocean=True,
        oceancolor="#eaf2f8",
        showcountries=True,
        countrycolor="rgba(175, 185, 195, 0.5)",
        countrywidth=0.5,
        showlakes=True,
        lakecolor="#eaf2f8",
        lataxis_range=[10, 72],
        lonaxis_range=[-130, 145],
        bgcolor="rgba(0,0,0,0)"
    )

    fig.update_layout(
        title=dict(
            text="<b>Geographic Distribution of <i>Arabidopsis thaliana</i> Accessions</b><br>"
                 "<span style='font-size: 11.5px; font-weight: normal; color: #555;'>"
                 f"1,003 1001-Genomes accessions colored by flowering time at 10 °C (FT10)</span>",
            x=0.03,
            y=0.96,
            font=dict(size=14, family="sans-serif")
        ),
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

    # Save high-res PNG image
    fig.write_image(out_png, scale=2)
    print(f"Generated interactive map: {out_html}")
    print(f"Generated high-resolution figure: {out_png}")


if __name__ == '__main__':
    generate_accession_maps()
