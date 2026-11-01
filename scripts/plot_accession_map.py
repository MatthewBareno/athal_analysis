"""Generate a publication-quality world map of Arabidopsis thaliana accessions colored by country.

Produces:
- results/accession_world_map.html: Interactive Plotly map with country legend & hover tooltips
- results/accession_world_map.png: High-resolution publication-quality figure (2200x1160)
"""

import os
import pandas as pd
import plotly.graph_objects as go


def generate_country_accession_map():
    data_path = 'data/processed/phenotype_clean.csv'
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Missing {data_path}. Run prepare_phenotypes.py first.")

    df = pd.read_csv(data_path)
    valid_geo = df.dropna(subset=['latitude', 'longitude']).copy()

    # Calculate accession counts per country
    country_counts = valid_geo['country'].value_counts()
    
    # Countries with >= 10 accessions get dedicated distinct color categories;
    # Countries with < 10 accessions are grouped as "Other" to keep the publication legend legible.
    major_countries = country_counts[country_counts >= 10].index.tolist()
    
    def get_group(country):
        if country in major_countries:
            return f"{country} (n={country_counts[country]})"
        else:
            return f"Other ({len(country_counts) - len(major_countries)} countries, n={country_counts[country_counts < 10].sum()})"

    valid_geo['country_group'] = valid_geo['country'].apply(get_group)

    # Order groups by frequency
    group_order = [f"{c} (n={country_counts[c]})" for c in major_countries]
    other_label = f"Other ({len(country_counts) - len(major_countries)} countries, n={country_counts[country_counts < 10].sum()})"
    group_order.append(other_label)

    # High-contrast, colorblind-friendly publication palette
    palette = [
        '#1f77b4',  # SWE - Steel Blue
        '#e6550d',  # ESP - Deep Orange
        '#2ca02c',  # GER - Forest Green
        '#756bb1',  # ITA - Purple
        '#d62728',  # UK - Crimson Red
        '#31a354',  # RUS - Emerald
        '#6baed6',  # FRA - Sky Blue
        '#e7ba52',  # USA - Warm Amber
        '#9ecae1',  # CZE - Powder Blue
        '#fdae6b',  # BUL - Salmon
        '#bcbddc',  # AZE - Lavender
        '#e377c2',  # GEO - Rose
        '#17becf',  # NED - Cyan
        '#8c564b',  # POR - Sienna Brown
        '#7f7f7f',  # Other - Medium Gray
    ]
    color_map = {grp: palette[i % len(palette)] for i, grp in enumerate(group_order)}

    fig = go.Figure()

    # Add each country group as an independent trace (allows toggling in HTML)
    for grp in group_order:
        sub = valid_geo[valid_geo['country_group'] == grp]
        if sub.empty:
            continue

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
                "<b>Coordinates:</b> %{lat:.3f}°N, %{lon:.3f}°E<extra></extra>"
            ),
            mode='markers',
            name=grp,
            marker=dict(
                size=6.0,
                color=color_map[grp],
                opacity=0.85,
                line=dict(width=0.45, color='#ffffff')
            )
        ))

    # Professional Natural Earth projection layout
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
                 f"Sampling sites for 1,003 accessions colored by country of origin (999 geo-referenced sites)</span>",
            x=0.03,
            y=0.96,
            font=dict(size=14, family="sans-serif", color="#0f172a")
        ),
        legend=dict(
            title=dict(text="<b>Country of Origin</b>", font=dict(size=11, family="sans-serif")),
            orientation="v",
            yanchor="middle",
            y=0.5,
            xanchor="left",
            x=1.01,
            font=dict(size=9.5, family="sans-serif"),
            itemsizing="constant",
            bgcolor="rgba(255,255,255,0.8)",
            bordercolor="#e2e8f0",
            borderwidth=1
        ),
        margin=dict(l=10, r=160, t=55, b=10),
        width=1150,
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
    print(f"Generated country-colored interactive map: {out_html}")
    print(f"Generated country-colored publication figure: {out_png}")


if __name__ == '__main__':
    generate_country_accession_map()
