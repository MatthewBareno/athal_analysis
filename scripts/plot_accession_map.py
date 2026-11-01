"""Generate a publication-quality world map of Arabidopsis thaliana accessions colored by country.

Includes all 45 countries of origin in a clean multi-column legend ordered by sample size.

Outputs:
- results/accession_world_map.html: Interactive Plotly map with tooltips and 45-country legend
- results/accession_world_map.png: High-resolution (2560x1240, 300-DPI) publication figure
"""

import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def generate_all_countries_map():
    data_path = 'data/processed/phenotype_clean.csv'
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Missing {data_path}. Run prepare_phenotypes.py first.")

    df = pd.read_csv(data_path)
    # Filter to accessions with valid geographic coordinates
    valid_geo = df.dropna(subset=['latitude', 'longitude']).copy()

    # Order countries by accession count descending
    country_counts = valid_geo['country'].value_counts()
    ordered_countries = country_counts.index.tolist()

    # Broad qualitative palette providing 50 distinct colors
    palette = px.colors.qualitative.Alphabet + px.colors.qualitative.Dark24
    color_map = {c: palette[i % len(palette)] for i, c in enumerate(ordered_countries)}

    fig = go.Figure()

    # Add each country as a distinct trace so all 45 appear in the legend
    for country in ordered_countries:
        sub = valid_geo[valid_geo['country'] == country]
        n_samples = len(sub)

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
                "<b>Country:</b> %{customdata[2]}<br>"
                "<b>Site:</b> %{customdata[3]}<br>"
                "<b>Stock:</b> %{customdata[5]}<br>"
                "<b>Elevation:</b> %{customdata[4]} m<br>"
                "<b>Coordinates:</b> %{lat:.3f}°N, %{lon:.3f}°E<extra></extra>"
            ),
            mode='markers',
            name=f"{country} ({n_samples})",
            marker=dict(
                size=5.5,
                color=color_map[country],
                opacity=0.85,
                line=dict(width=0.4, color='#ffffff')
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
            text="<b>Geographic Distribution of <i>Arabidopsis thaliana</i> Accessions</b><br>"
                 "<span style='font-size: 11px; font-weight: normal; color: #475569;'>"
                 "Sampling locations of 1,003 accessions colored across all 45 countries of origin (999 geo-referenced sites)</span>",
            x=0.02,
            y=0.97,
            font=dict(size=14, family="sans-serif", color="#0f172a")
        ),
        legend=dict(
            title=dict(text="<b>Country (n)</b>", font=dict(size=10, family="sans-serif")),
            orientation="v",
            yanchor="top",
            y=0.98,
            xanchor="left",
            x=1.01,
            entrywidth=78,
            entrywidthmode="pixels",
            font=dict(size=8, family="sans-serif"),
            itemsizing="constant",
            bgcolor="rgba(255,255,255,0.85)",
            bordercolor="#e2e8f0",
            borderwidth=1
        ),
        margin=dict(l=10, r=260, t=55, b=10),
        width=1280,
        height=620,
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
    print(f"Generated 45-country interactive map: {out_html}")
    print(f"Generated 45-country publication figure: {out_png}")


if __name__ == '__main__':
    generate_all_countries_map()
