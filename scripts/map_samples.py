import pandas as pd
import plotly.express as px

df = pd.read_csv('data/accession_list.csv')
fig = px.scatter_geo(df, lat='latitude', lon='longitude', hover_name='name')
fig.write_html('results/sample_map.html')
