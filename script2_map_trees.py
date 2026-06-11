import pandas as pd
import folium
from folium.plugins import MarkerCluster
import requests

PARKS = {
    'Central': {
        'lat_min': 40.764, 'lon_min': -73.982,
        'lat_max': 40.800, 'lon_max': -73.949,
        'center': [40.782, -73.965],
    },
    'Prospect': {
        'lat_min': 40.656, 'lon_min': -73.978,
        'lat_max': 40.683, 'lon_max': -73.955,
        'center': [40.666, -73.969],
    },
    'Clove_Lakes': {
        'lat_min': 40.620, 'lon_min': -74.120,
        'lat_max': 40.637, 'lon_max': -74.102,
        'center': [40.628, -74.111],
    },
}

SPECIES_COLORS = [
    '#e6194b', '#3cb44b', '#4363d8', '#f58231', '#911eb4',
    '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990',
]
OTHER_COLOR = '#aaaaaa'

URL = 'https://data.cityofnewyork.us/resource/uvpi-gqnh.json'


def fetch_trees(park_name, bounds):
    where = (
        f"latitude > '{bounds['lat_min']}' AND latitude < '{bounds['lat_max']}' "
        f"AND longitude > '{bounds['lon_min']}' AND longitude < '{bounds['lon_max']}'"
    )
    params = {
        '$where': where,
        '$limit': 50000,
        '$select': 'latitude,longitude,spc_common,health,tree_dbh,status',
    }
    resp = requests.get(URL, params=params, timeout=30)
    resp.raise_for_status()
    df = pd.DataFrame(resp.json())
    df['park'] = park_name
    print(f"  {park_name}: {len(df)} trees")
    return df


def build_color_map(df):
    top = df['spc_common'].value_counts().head(len(SPECIES_COLORS)).index.tolist()
    return {species: SPECIES_COLORS[i] for i, species in enumerate(top)}


def build_map(df, color_map):
    center = [df['latitude'].astype(float).mean(), df['longitude'].astype(float).mean()]
    m = folium.Map(location=center, zoom_start=12, tiles='CartoDB positron')

    for park_name, bounds in PARKS.items():
        cluster = MarkerCluster(name=park_name).add_to(m)
        park_df = df[df['park'] == park_name]

        for _, row in park_df.iterrows():
            try:
                lat, lon = float(row['latitude']), float(row['longitude'])
            except (ValueError, KeyError):
                continue

            species = row.get('spc_common', 'unknown') or 'unknown'
            color = color_map.get(species, OTHER_COLOR)
            health = row.get('health', 'N/A') or 'N/A'
            dbh = row.get('tree_dbh', 'N/A') or 'N/A'

            folium.CircleMarker(
                location=[lat, lon],
                radius=5,
                color=color,
                fill=True,
                fill_color=color,
                fill_opacity=0.7,
                popup=folium.Popup(
                    f"<b>{species}</b><br>Park: {park_name}<br>Health: {health}<br>DBH: {dbh} in",
                    max_width=200,
                ),
                tooltip=species,
            ).add_to(cluster)

    # legend
    legend_html = '<div style="position:fixed;bottom:30px;left:30px;z-index:1000;background:white;padding:10px;border-radius:6px;font-size:12px;box-shadow:2px 2px 6px rgba(0,0,0,0.3)">'
    legend_html += '<b>Top species</b><br>'
    for species, color in color_map.items():
        legend_html += f'<span style="background:{color};display:inline-block;width:12px;height:12px;margin-right:4px;border-radius:50%"></span>{species}<br>'
    legend_html += f'<span style="background:{OTHER_COLOR};display:inline-block;width:12px;height:12px;margin-right:4px;border-radius:50%"></span>other<br>'
    legend_html += '</div>'
    m.get_root().html.add_child(folium.Element(legend_html))

    folium.LayerControl().add_to(m)
    return m


if __name__ == '__main__':
    print('Fetching tree data...')
    frames = [fetch_trees(name, bounds) for name, bounds in PARKS.items()]
    df = pd.concat(frames, ignore_index=True)
    print(f'Total: {len(df)} trees')

    color_map = build_color_map(df)
    m = build_map(df, color_map)
    m.save('trees_map.html')
    print('Saved to trees_map.html')
