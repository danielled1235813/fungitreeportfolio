"""
Fetch individual tree records from three sources for Central Park,
Prospect Park, and Clove Lakes. Outputs a unified CSV and Folium map.

Sources
-------
NYC Open Data hn5i-inap  : Prospect Park + Clove Lakes (Central Park not available)
OpenStreetMap Overpass   : all three parks  (natural=tree nodes)
iNaturalist              : all three parks  (research-grade plant observations)
"""

import time
import requests
import pandas as pd
import folium
from folium.plugins import MarkerCluster

PARKS = {
    'Central': {
        'lat_min': 40.764, 'lon_min': -73.982,
        'lat_max': 40.800, 'lon_max': -73.949,
        'center': [40.782, -73.965],
        'inat_place_id': 49955,
        'nyc_opendata': False,
    },
    'Prospect': {
        'lat_min': 40.656, 'lon_min': -73.978,
        'lat_max': 40.683, 'lon_max': -73.955,
        'center': [40.666, -73.969],
        'inat_place_id': 55174,
        'nyc_opendata': True,
    },
    'Clove_Lakes': {
        'lat_min': 40.620, 'lon_min': -74.120,
        'lat_max': 40.637, 'lon_max': -74.102,
        'center': [40.628, -74.111],
        'inat_place_id': 125420,
        'nyc_opendata': True,
    },
}

NYC_URL = 'https://data.cityofnewyork.us/resource/hn5i-inap.json'
OVERPASS_URL = 'https://overpass-api.de/api/interpreter'
INAT_URL = 'https://api.inaturalist.org/v1/observations'

PROVENANCE_COLORS = {
    'NYC Open Data - Forestry Tree Points': '#2196F3',
    'OpenStreetMap': '#4CAF50',
    'iNaturalist': '#FF9800',
}


# ── fetchers ──────────────────────────────────────────────────────────────────

def _parse_genusspecies(val):
    """'Acer nigrum - black maple' → ('Acer nigrum', 'black maple')"""
    if not val:
        return '', ''
    parts = str(val).split(' - ', 1)
    return parts[0].strip(), (parts[1].strip() if len(parts) > 1 else '')


def fetch_nyc_opendata(park_name, bounds):
    # build query string manually — requests encodes '$' which Socrata rejects
    # note: NW lat/lon first, then SE lat/lon per Socrata within_box spec
    records, offset, limit = [], 0, 5000
    while True:
        where = (
            f"within_box(location,{bounds['lat_max']},{bounds['lon_min']},"
            f"{bounds['lat_min']},{bounds['lon_max']})"
            " AND tpstructure!='Retired' AND tpstructure!='Stump'"
        )
        qs = (
            f"$where={where}&$limit={limit}&$offset={offset}"
            "&$select=genusspecies,dbh,tpcondition,tpstructure,location,createddate"
        )
        try:
            resp = requests.get(f"{NYC_URL}?{qs}", timeout=30)
            resp.raise_for_status()
        except requests.HTTPError as e:
            if resp.status_code == 500 and offset > 0:
                # Socrata deep-pagination limit hit; use what we have
                break
            raise
        batch = resp.json()
        if not batch:
            break
        for r in batch:
            loc = r.get('location', {})
            coords = loc.get('coordinates', [None, None]) if isinstance(loc, dict) else [None, None]
            if None in coords:
                continue
            latin, common = _parse_genusspecies(r.get('genusspecies'))
            records.append({
                'park': park_name,
                'lat': float(coords[1]),
                'lon': float(coords[0]),
                'species_latin': latin,
                'species_common': common,
                'dbh_cm': r.get('dbh'),
                'health': r.get('tpcondition'),
                'provenance': 'NYC Open Data - Forestry Tree Points',
                'date_observed': (r.get('createddate') or '')[:10],
            })
        if len(batch) < limit:
            break
        offset += limit
    return records


def fetch_osm(park_name, bounds):
    query = f"""
[out:json][timeout:90];
node[natural=tree]
  ({bounds['lat_min']},{bounds['lon_min']},{bounds['lat_max']},{bounds['lon_max']});
out body;
"""
    headers = {'User-Agent': 'belowground-nyc-tree-research/1.0'}
    for attempt in range(3):
        resp = requests.post(OVERPASS_URL, data={'data': query}, headers=headers, timeout=120)
        if resp.status_code == 429:
            time.sleep(10 * (attempt + 1))
            continue
        resp.raise_for_status()
        break
    records = []
    for el in resp.json().get('elements', []):
        tags = el.get('tags', {})
        latin = tags.get('species') or tags.get('taxon') or tags.get('genus') or ''
        common = tags.get('species:en') or tags.get('name') or ''
        records.append({
            'park': park_name,
            'lat': el['lat'],
            'lon': el['lon'],
            'species_latin': latin,
            'species_common': common,
            'dbh_cm': tags.get('circumference'),
            'health': None,
            'provenance': 'OpenStreetMap',
            'date_observed': (el.get('timestamp') or '')[:10],
        })
    return records


def fetch_inat(park_name, place_id):
    records, page, per_page = [], 1, 200
    while True:
        resp = requests.get(INAT_URL, params={
            'place_id': place_id,
            'iconic_taxa': 'Plantae',
            'quality_grade': 'research',
            'per_page': per_page,
            'page': page,
        }, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        results = data.get('results', [])
        if not results:
            break
        for obs in results:
            loc = obs.get('location') or ''
            if not loc:
                continue
            try:
                lat, lon = map(float, loc.split(','))
            except ValueError:
                continue
            taxon = obs.get('taxon') or {}
            records.append({
                'park': park_name,
                'lat': lat,
                'lon': lon,
                'species_latin': taxon.get('name', ''),
                'species_common': taxon.get('preferred_common_name', ''),
                'dbh_cm': None,
                'health': None,
                'provenance': 'iNaturalist',
                'date_observed': obs.get('observed_on', ''),
            })
        total = data.get('total_results', 0)
        if page * per_page >= min(total, 10000):
            break
        page += 1
        time.sleep(0.5)
    return records


# ── map ───────────────────────────────────────────────────────────────────────

def build_map(df):
    center = [df['lat'].mean(), df['lon'].mean()]
    m = folium.Map(location=center, zoom_start=12, tiles='CartoDB positron')

    for park_name in df['park'].unique():
        cluster = MarkerCluster(name=park_name).add_to(m)
        for _, row in df[df['park'] == park_name].iterrows():
            color = PROVENANCE_COLORS.get(row['provenance'], '#999999')
            folium.CircleMarker(
                location=[row['lat'], row['lon']],
                radius=5,
                color=color,
                fill=True,
                fill_color=color,
                fill_opacity=0.75,
                popup=folium.Popup(
                    f"<b>{row['species_latin'] or 'unknown'}</b><br>"
                    f"{row['species_common']}<br>"
                    f"Park: {row['park']}<br>"
                    f"Health: {row['health'] or 'N/A'}<br>"
                    f"DBH: {row['dbh_cm'] or 'N/A'} cm<br>"
                    f"Source: {row['provenance']}<br>"
                    f"Date: {row['date_observed'] or 'N/A'}",
                    max_width=220,
                ),
                tooltip=row['species_latin'] or 'unknown',
            ).add_to(cluster)

    legend_html = (
        '<div style="position:fixed;bottom:30px;left:30px;z-index:1000;'
        'background:white;padding:10px 14px;border-radius:6px;font-size:12px;'
        'box-shadow:2px 2px 6px rgba(0,0,0,0.3)">'
        '<b>Data source</b><br>'
    )
    for label, color in PROVENANCE_COLORS.items():
        legend_html += (
            f'<span style="background:{color};display:inline-block;'
            f'width:12px;height:12px;margin-right:6px;border-radius:50%">'
            f'</span>{label}<br>'
        )
    legend_html += '</div>'
    m.get_root().html.add_child(folium.Element(legend_html))
    folium.LayerControl().add_to(m)
    return m


# ── main ──────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    all_records = []

    for park_name, bounds in PARKS.items():
        print(f'\n{park_name}')

        if bounds['nyc_opendata']:
            print('  fetching NYC Open Data...')
            rows = fetch_nyc_opendata(park_name, bounds)
            print(f'  {len(rows)} trees')
            all_records.extend(rows)

        print('  fetching OpenStreetMap...')
        rows = fetch_osm(park_name, bounds)
        print(f'  {len(rows)} trees')
        all_records.extend(rows)

        print('  fetching iNaturalist...')
        rows = fetch_inat(park_name, bounds['inat_place_id'])
        print(f'  {len(rows)} observations')
        all_records.extend(rows)

    df = pd.DataFrame(all_records)
    print(f'\nTotal records: {len(df)}')
    print(df.groupby(['park', 'provenance']).size().to_string())

    df.to_csv('trees_combined.csv', index=False)
    print('\nSaved to trees_combined.csv')

    m = build_map(df)
    m.save('trees_combined_map.html')
    print('Saved to trees_combined_map.html')
