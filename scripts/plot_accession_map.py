"""Generate a publication-quality world map of Arabidopsis thaliana accessions grouped by macro-region.

Classifies all 1,003 accessions into 10 geographic macro-regions corresponding
to the primary sampling zones and genetic ancestry clusters in 1001 Genomes studies.

Outputs:
- results/accession_world_map.html: Interactive Plotly map with macro-region filtering
- results/accession_world_map.png: High-resolution (2200x1160, 300-DPI) publication figure
"""

import os
import pandas as pd
import plotly.graph_objects as go

# Standard biogeographic macro-regions for 1001 Genomes accessions
MACRO_REGIONS = {
    # Scandinavia & Baltic
    'SWE': 'Scandinavia & Baltic', 'NOR': 'Scandinavia & Baltic', 'FIN': 'Scandinavia & Baltic',
    'DEN': 'Scandinavia & Baltic', 'EST': 'Scandinavia & Baltic', 'LTU': 'Scandinavia & Baltic',

    # Western Europe
    'UK': 'Western Europe', 'FRA': 'Western Europe', 'GER': 'Western Europe',
    'NED': 'Western Europe', 'BEL': 'Western Europe', 'IRL': 'Western Europe',

    # Iberian Peninsula
    'ESP': 'Iberian Peninsula', 'POR': 'Iberian Peninsula', 'CPV': 'Iberian Peninsula',

    # South Europe & Balkans
    'ITA': 'South Europe & Balkans', 'BUL': 'South Europe & Balkans', 'ROU': 'South Europe & Balkans',
    'SRB': 'South Europe & Balkans', 'CRO': 'South Europe & Balkans', 'GRC': 'South Europe & Balkans',

    # Central Europe
    'CZE': 'Central Europe', 'AUT': 'Central Europe', 'SUI': 'Central Europe',
    'POL': 'Central Europe', 'SVK': 'Central Europe',

    # Eastern Europe & Russia
    'RUS': 'Eastern Europe & Russia', 'UKR': 'Eastern Europe & Russia',

    # North America
    'USA': 'North America', 'CAN': 'North America',

    # Caucasus & Near East
    'AZE': 'Caucasus & Near East', 'GEO': 'Caucasus & Near East', 'ARM': 'Caucasus & Near East',
    'TUR': 'Caucasus & Near East', 'IRN': 'Caucasus & Near East', 'LBN': 'Caucasus & Near East',

    # Central & East Asia
    'KGZ': 'Central & East Asia', 'KAZ': 'Central & East Asia', 'UZB': 'Central & East Asia',
    'TJK': 'Central & East Asia', 'AFG': 'Central & East Asia', 'IND': 'Central & East Asia',
    'CHN': 'Central & East Asia', 'JPN': 'Central & East Asia',

    # North Africa
    'MAR': 'North Africa'
}

# Colorblind-safe, distinct publication palette
REGION_COLORS = {
    'Western Europe': '#2563eb',          # Royal Blue
    'Scandinavia & Baltic': '#0284c7',    # Nordic Cyan
    'Iberian Peninsula': '#ea580c',       # Warm Terracotta
    'South Europe & Balkans': '#9333ea',  # Mediterranean Violet
    'Central Europe': '#16a34a',          # Continental Green
    'Eastern Europe & Russia': '#0d9488', # Boreal Teal
    'North America': '#dc2626',           # Crimson
    'Caucasus & Near East': '#d97706',    # Caucasus Ochre
    'Central & East Asia': '#db2777',     # Asian Rose
    'North Africa': '#854d0e'             # Saharan Bronze
}


def generate_macro_region_map():
    data_path = 'data/processed/phenotype_clean.csv'
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Missing {data_path}. Run prepare_phenotypes.py first.")

    df = pd.read_csv(data_path)
    valid_geo = df.dropna(subset=['latitude', 'longitude']).copy()
    valid_geo['macro_region'] = valid_geo['country'].map(MACRO_REGIONS).fillna('Other')

    # Compute accession count per region and order by prevalence
    region_counts = valid_geo['macro_region'].value_counts()
    ordered_regions = region_counts.index.tolist()

    fig = go.Figure()

    # Add each macro-region as a separate interactive trace
    for region in ordered_regions:
        sub = valid_geo[valid_geo['macro_region'] == region]
        n_samples = len(sub)
        color = REGION_COLORS.get(region, '#64748b')

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
                f"<b>Region:</b> {region}<br>"
                "<b>Country:</b> %{customdata[2]} (%{customdata[3]})<br>"
                "<b>Stock:</b> %{customdata[5]}<br>"
                "<b>Elevation:</b> %{customdata[4]} m<br>"
                "<b>Coordinates:</b> %{lat:.3f}°N, %{lon:.3f}°E<extra></extra>"
            ),
            mode='markers',
            name=f"{region} (n={n_samples})",
            marker=dict(
                size=6.2,
                color=color,
                opacity=0.85,
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
            text="<b>Geographic Distribution of <i>Arabidopsis thaliana</i> Accessions</b><br>"
                 "<span style='font-size: 11.5px; font-weight: normal; color: #475569;'>"
                 f"Sampling sites for 1,003 accessions across 10 biogeographic macro-regions (999 geo-referenced sites)</span>",
            x=0.03,
            y=0.96,
            font=dict(size=14, family="sans-serif", color="#0f172a")
        ),
        legend=dict(
            title=dict(text="<b>Macro-Region</b>", font=dict(size=11, family="sans-serif")),
            orientation="v",
            yanchor="middle",
            y=0.5,
            xanchor="left",
            x=1.01,
            font=dict(size=9.5, family="sans-serif"),
            itemsizing="constant",
            bgcolor="rgba(255,255,255,0.85)",
            bordercolor="#e2e8f0",
            borderwidth=1
        ),
        margin=dict(l=10, r=180, t=55, b=10),
        width=1160,
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
    print(f"Generated macro-region interactive map: {out_html}")
    print(f"Generated macro-region publication figure: {out_png}")


if __name__ == '__main__':
    generate_macro_region_map()
