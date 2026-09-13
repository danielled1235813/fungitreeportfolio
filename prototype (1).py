#!/usr/bin/env python3
"""
Belowground prototype - produces three mycorrhizal network maps.

Map 1  map1_suitability.html     Habitat suitability heatmap
Map 2  map2_fungal_network.html  Inferred underground fungal connectivity
Map 3  map3_tree_network.html    Tree-to-tree network via shared fungi

Guild assignments and tree mycorrhizal types are approximated from a
hard-coded genus lookup table (FungalTraits join is a planned next step).
Field data (soil chemistry + microBIOMETER) is not yet available; when it
is, it feeds into Maps 2 and 3 as edge weight modifiers via the
suitability surface from Map 1.
"""

import math
import colorsys
import warnings
from collections import Counter
import pandas as pd
import numpy as np
import folium
from folium.plugins import HeatMap, FastMarkerCluster
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")

# ── Configuration ─────────────────────────────────────────────────────────────

FUNGI_CSV = "fungi_observations_nyc_parks.csv"
TREES_CSV = "trees_combined.csv"

ECM_RADIUS_M   = 30    # max ECM hyphal network extent (meters)
AM_RADIUS_M    = 8     # max AM hyphal network extent (meters)
TREE_RADIUS_MULT = 1.5 # tree-fungi search radius = fungal radius × this

MAP2_THRESHOLD = 0.15  # min edge probability to draw on Map 2
MAP3_THRESHOLD = 0.10  # min edge probability to draw on Map 3

MAX_TREE_MARKERS = 2000   # kept for reference; Map 1 now shows all host trees
MAX_TREES_NETWORK = 100000  # effectively uncapped: consider every tree per park/guild

MAP_CENTER = [40.72, -73.96]

# Basemaps: Esri gray canvas (keyless). CartoDB now stamps "API KEY REQUIRED"
# on its free tiles when loaded from a browser, so we avoid it.
ESRI_DARK  = "https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Dark_Gray_Base/MapServer/tile/{z}/{y}/{x}"
ESRI_LIGHT = "https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}"
ESRI_ATTR  = "Tiles &copy; Esri, HERE, Garmin, &copy; OpenStreetMap contributors"


def add_basemap(m, url):
    """Add a keyless Esri basemap. max_native_zoom=16 upscales past its top zoom."""
    folium.TileLayer(
        tiles=url, attr=ESRI_ATTR, name="Basemap",
        max_native_zoom=16, max_zoom=20, control=False,
    ).add_to(m)

# ── Guild lookup (placeholder for FungalTraits join) ──────────────────────────
# Source: known ecology of genera commonly observed in northeastern US parks.
# Replace with FungalTraits primary_lifestyle join once data is available.

ECM_GENERA = {
    "amanita", "boletus", "suillus", "russula", "lactarius", "cantharellus",
    "tricholoma", "cortinarius", "inocybe", "laccaria", "pisolithus",
    "scleroderma", "cenococcum", "hebeloma", "xerocomus", "tylopilus",
    "paxillus", "gyroporus", "boletellus", "strobilomyces", "chalciporus",
    "hygrophorus", "tomentella", "thelephora", "rhizopogon", "gautieria",
    "elaphomyces", "tuber", "hydnum", "sarcodon", "bankera", "clavulina",
    "sebacina", "wilcoxina", "amphinema", "piloderma", "cenococcum",
    "leccinum", "chroogomphus", "gomphidius", "truncocolumella",
}

AM_GENERA = {
    "glomus", "rhizophagus", "funneliformis", "claroideoglomus", "gigaspora",
    "scutellospora", "diversispora", "acaulospora", "ambispora", "archaeospora",
    "paraglomus", "redeckera", "septoglomus",
}

# ── Tree mycorrhizal type lookup (placeholder for FungalRoot join) ────────────

ECM_TREE_GENERA = {
    "quercus", "fagus", "betula", "pinus", "picea", "abies", "larix",
    "pseudotsuga", "tsuga", "carpinus", "corylus", "castanea",
    "alnus", "salix", "populus",
}

AM_TREE_GENERA = {
    "acer", "fraxinus", "ulmus", "prunus", "robinia", "platanus",
    "gleditsia", "gymnocladus", "liriodendron", "cercis", "liquidambar",
    "nyssa", "cornus", "amelanchier", "crataegus", "malus", "pyrus",
    "sorbus", "celtis", "morus", "catalpa", "ailanthus", "tilia",
    "koelreuteria", "styphnolobium", "sophora", "magnolia", "sassafras",
    "juglans", "carya", "maclura", "ginkgo", "metasequoia", "taxodium",
    "zelkova",
}

# Normalize park names across the two CSVs
PARK_NORM = {
    "Central Park":    "Central",
    "Prospect Park":   "Prospect",
    "Clove Lakes Park":"Clove_Lakes",
    "Central":         "Central",
    "Prospect":        "Prospect",
    "Clove_Lakes":     "Clove_Lakes",
}

# ── Helpers ───────────────────────────────────────────────────────────────────

def haversine(lat1, lon1, lat2, lon2):
    R = 6_371_000
    p1, p2 = math.radians(lat1), math.radians(lat2)
    a = (math.sin((p2 - p1) / 2) ** 2
         + math.cos(p1) * math.cos(p2) * math.sin(math.radians(lon2 - lon1) / 2) ** 2)
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

def to_xy(lats, lons, ref_lat, ref_lon):
    """Approximate lat/lon arrays to x/y metres from a reference point."""
    cos_ref = math.cos(math.radians(ref_lat))
    x = (np.asarray(lons) - ref_lon) * 111_320 * cos_ref
    y = (np.asarray(lats) - ref_lat) * 111_320
    return np.column_stack([x, y])

def assign_guild(taxon_name):
    if not isinstance(taxon_name, str):
        return "other"
    g = taxon_name.strip().split()[0].lower()
    if g in ECM_GENERA:
        return "ectomycorrhizal"
    if g in AM_GENERA:
        return "arbuscular_mycorrhizal"
    return "other"

def assign_tree_myco(species_latin):
    if not isinstance(species_latin, str) or not species_latin.strip():
        return "unresolved"
    g = species_latin.strip().split()[0].lower()
    if g in ECM_TREE_GENERA:
        return "ECM"
    if g in AM_TREE_GENERA:
        return "AM"
    return "unresolved"

def quality_w(q):
    return {"research": 1.0, "needs_id": 0.7, "casual": 0.4}.get(q, 0.5)

def taxon_base_prob(a, b):
    """Base connection probability from taxonomic similarity."""
    if not isinstance(a, str) or not isinstance(b, str):
        return 0.10
    a, b = a.strip().lower(), b.strip().lower()
    if a == b:
        return 0.75
    if a.split()[0] == b.split()[0]:
        return 0.45
    return 0.12

# ── Data loading ──────────────────────────────────────────────────────────────

def load_and_enrich():
    print("Loading data...")
    fungi = pd.read_csv(FUNGI_CSV, low_memory=False)
    trees = pd.read_csv(TREES_CSV, low_memory=False)

    fungi["latitude"]  = pd.to_numeric(fungi["latitude"],  errors="coerce")
    fungi["longitude"] = pd.to_numeric(fungi["longitude"], errors="coerce")
    trees["lat"] = pd.to_numeric(trees["lat"], errors="coerce")
    trees["lon"] = pd.to_numeric(trees["lon"], errors="coerce")

    fungi = fungi.dropna(subset=["latitude", "longitude"])
    trees = trees.dropna(subset=["lat", "lon"])

    fungi["park_key"] = fungi["park"].map(PARK_NORM).fillna(fungi["park"])
    trees["park_key"] = trees["park"].map(PARK_NORM).fillna(trees["park"])

    fungi["guild"]     = fungi["taxon_name"].apply(assign_guild)
    trees["myco_type"] = trees["species_latin"].apply(assign_tree_myco)

    print(f"  Fungi: {len(fungi):,}  |  {fungi['guild'].value_counts().to_dict()}")
    print(f"  Trees: {len(trees):,}  |  {trees['myco_type'].value_counts().to_dict()}")
    return fungi, trees

def network_palette(n):
    """n visually distinct hex colors, spread around the wheel by golden ratio."""
    cols = []
    for i in range(max(n, 1)):
        h = (i * 0.6180339887) % 1.0
        r, g, b = colorsys.hsv_to_rgb(h, 0.62, 0.88)
        cols.append("#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255)))
    return cols


# ── Map 1 - Habitat Suitability ───────────────────────────────────────────────

def build_map1(fungi, trees):
    print("\nMap 1 - Habitat Suitability Heatmap")

    myco_f  = fungi[fungi["guild"].isin(["ectomycorrhizal", "arbuscular_mycorrhizal"])]
    ecm_t   = trees[trees["myco_type"] == "ECM"]
    am_t    = trees[trees["myco_type"] == "AM"]
    print(f"  Mycorrhizal fungi: {len(myco_f):,}  ECM trees: {len(ecm_t):,}  AM trees: {len(am_t):,}")

    # prefer_canvas renders the full tree layer fast even at 20k+ points.
    m = folium.Map(location=MAP_CENTER, zoom_start=12, tiles=None, prefer_canvas=True)
    add_basemap(m, ESRI_DARK)

    # Suitability proxy: density of mycorrhizal fungal observations.
    heat_pts = myco_f[["latitude", "longitude"]].values.tolist()
    if heat_pts:
        HeatMap(
            heat_pts,
            name="Suitability proxy (mycorrhizal observation density)",
            min_opacity=0.25, radius=22, blur=18,
            gradient={0.2: "#0d47a1", 0.4: "#00c853", 0.65: "#ffeb3b",
                      0.85: "#ff6d00", 1.0: "#b71c1c"},
        ).add_to(m)

    # ALL ECM host trees, no sampling.
    ecm_layer = folium.FeatureGroup(name=f"ECM host trees ({len(ecm_t):,})", show=True)
    for _, r in ecm_t.iterrows():
        folium.CircleMarker(
            [r["lat"], r["lon"]], radius=2,
            color="#4fc3f7", fill=True, fill_color="#4fc3f7",
            fill_opacity=0.6, weight=0,
            tooltip=r.get("species_latin") or "ECM tree",
        ).add_to(ecm_layer)
    ecm_layer.add_to(m)

    # ALL AM host trees, no sampling.
    am_layer = folium.FeatureGroup(name=f"AM host trees ({len(am_t):,})", show=True)
    for _, r in am_t.iterrows():
        folium.CircleMarker(
            [r["lat"], r["lon"]], radius=2,
            color="#81c784", fill=True, fill_color="#81c784",
            fill_opacity=0.6, weight=0,
            tooltip=r.get("species_latin") or "AM tree",
        ).add_to(am_layer)
    am_layer.add_to(m)

    folium.LayerControl(collapsed=False).add_to(m)

    m.get_root().html.add_child(folium.Element(f"""
    <div style="position:fixed;bottom:30px;left:30px;z-index:9999;
                background:rgba(0,0,0,0.82);padding:14px 18px;border-radius:8px;
                color:white;font-family:sans-serif;font-size:13px;line-height:1.6">
      <b>Map 1 - Habitat Suitability</b><br>
      <span style="color:#b71c1c">&#9632;</span> High suitability<br>
      <span style="color:#ff6d00">&#9632;</span> Moderate<br>
      <span style="color:#ffeb3b">&#9632;</span> Low<br>
      <span style="color:#0d47a1">&#9632;</span> Minimal signal<br><br>
      <span style="color:#4fc3f7">&#9679;</span> ECM host trees ({len(ecm_t):,})<br>
      <span style="color:#81c784">&#9679;</span> AM host trees ({len(am_t):,})<br>
      <span style="font-size:11px;color:#cfcfcf">All {len(ecm_t) + len(am_t):,} host trees shown (no sampling).</span><br><br>
      <i style="font-size:11px">Proxy: iNaturalist mycorrhizal<br>
      observation density. Soil chemistry<br>
      + microBIOMETER data will replace<br>
      this surface in the next iteration.</i>
    </div>"""))

    m.save("map1_suitability.html")
    print(f"  -> map1_suitability.html  ({len(ecm_t) + len(am_t):,} host trees)")

# ── Map 2 - Fungal Network ────────────────────────────────────────────────────

def fungal_edges_for_park(df, guild, radius_m):
    """
    Compute probabilistic edges between fruiting body observations of the
    same guild within one park.  Returns a list of edge dicts.
    """
    if len(df) < 2:
        return []

    ref_lat, ref_lon = df["latitude"].mean(), df["longitude"].mean()
    xy   = to_xy(df["latitude"].values, df["longitude"].values, ref_lat, ref_lon)
    kdtree = cKDTree(xy)
    pairs  = kdtree.query_pairs(r=radius_m)

    recs   = df.reset_index(drop=True)
    edges  = []
    for i, j in pairs:
        a, b = recs.iloc[i], recs.iloc[j]
        dist = haversine(a["latitude"], a["longitude"], b["latitude"], b["longitude"])
        base = taxon_base_prob(a.get("taxon_name"), b.get("taxon_name"))
        decay = math.exp(-dist / (radius_m / 3))
        prob  = base * decay * math.sqrt(quality_w(a.get("quality_grade", "casual"))
                                         * quality_w(b.get("quality_grade", "casual")))
        if prob >= MAP2_THRESHOLD:
            edges.append({
                "lat_a": a["latitude"],  "lon_a": a["longitude"],
                "lat_b": b["latitude"],  "lon_b": b["longitude"],
                "probability": round(prob, 3),
                "taxon_a": a.get("taxon_name", ""),
                "taxon_b": b.get("taxon_name", ""),
                "distance_m": round(dist, 1),
                "guild": guild,
            })
    return edges

def build_map2(fungi):
    print("\nMap 2 - Underground Fungal Network")

    ecm = fungi[fungi["guild"] == "ectomycorrhizal"]
    am  = fungi[fungi["guild"] == "arbuscular_mycorrhizal"]
    print(f"  ECM observations: {len(ecm):,}  AM observations: {len(am):,}")

    all_edges = []
    for park in fungi["park_key"].unique():
        pe = ecm[ecm["park_key"] == park]
        pa = am[am["park_key"] == park]
        ecm_edges = fungal_edges_for_park(pe, "ectomycorrhizal", ECM_RADIUS_M)
        am_edges  = fungal_edges_for_park(pa, "arbuscular_mycorrhizal", AM_RADIUS_M)
        all_edges.extend(ecm_edges)
        all_edges.extend(am_edges)
        print(f"  {park}: {len(ecm_edges):,} ECM edges, {len(am_edges):,} AM edges")

    # Group fungi into networks (connected components) so color means network.
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        parent[find(a)] = find(b)

    def key(lat, lon):
        return (round(lat, 5), round(lon, 5))

    for e in all_edges:
        union(key(e["lat_a"], e["lon_a"]), key(e["lat_b"], e["lon_b"]))

    comp_size = Counter(find(n) for n in parent)
    roots_sorted = [root for root, _ in comp_size.most_common()]
    palette = network_palette(len(roots_sorted))
    root_color = {root: palette[i] for i, root in enumerate(roots_sorted)}
    n_networks = len(roots_sorted)
    connected = set(parent.keys())

    m = folium.Map(location=MAP_CENTER, zoom_start=12, tiles=None)
    add_basemap(m, ESRI_DARK)

    # Inferred connections, colored by the network they belong to.
    edge_layer = folium.FeatureGroup(name="Inferred fungal connections", show=True)
    for e in all_edges:
        col = root_color[find(key(e["lat_a"], e["lon_a"]))]
        folium.PolyLine(
            [[e["lat_a"], e["lon_a"]], [e["lat_b"], e["lon_b"]]],
            color=col, weight=2, opacity=min(e["probability"] + 0.3, 0.9),
            tooltip=(f"{e['taxon_a'] or '?'} to {e['taxon_b'] or '?'}  "
                     f"P={e['probability']}  {e['distance_m']}m"),
        ).add_to(edge_layer)
    edge_layer.add_to(m)

    # Fungi that belong to a network, filled with their network color.
    net_layer = folium.FeatureGroup(name="Fungi in a network (colored)", show=True)
    # Fungi with no inferred link, gray and off by default so networks stand out.
    iso_layer = folium.FeatureGroup(name="Unconnected fungi (gray)", show=False)
    for _, r in pd.concat([ecm, am]).iterrows():
        k = key(r["latitude"], r["longitude"])
        label = f"{r.get('taxon_name','')}  ({r.get('quality_grade','')})"
        if k in connected:
            root = find(k)
            folium.CircleMarker(
                [r["latitude"], r["longitude"]], radius=5,
                color="#ffffff", weight=0.5, fill=True,
                fill_color=root_color[root], fill_opacity=0.95,
                tooltip=f"{label} | network of {comp_size[root]} fungi",
            ).add_to(net_layer)
        else:
            folium.CircleMarker(
                [r["latitude"], r["longitude"]], radius=3,
                color="#888888", weight=0, fill=True,
                fill_color="#888888", fill_opacity=0.45,
                tooltip=f"{label} | not in a network",
            ).add_to(iso_layer)
    net_layer.add_to(m)
    iso_layer.add_to(m)

    folium.LayerControl(collapsed=False).add_to(m)

    if connected:
        lats = [c[0] for c in connected]
        lons = [c[1] for c in connected]
        m.fit_bounds([[min(lats), min(lons)], [max(lats), max(lons)]])

    m.get_root().html.add_child(folium.Element(f"""
    <div style="position:fixed;bottom:30px;left:30px;z-index:9999;
                background:rgba(0,0,0,0.82);padding:14px 18px;border-radius:8px;
                color:white;font-family:sans-serif;font-size:13px;line-height:1.6">
      <b>Map 2 - Underground Fungal Network</b><br>
      <i style="font-size:11px">Are these fruiting bodies connected below ground?</i><br><br>
      <b>Each color = one fungal network.</b><br>
      <span style="font-size:11px;color:#cfcfcf">Fungi of the same color are inferred to<br>
      share an underground network. Gray fungi<br>
      are not linked to any other (toggle on at right).</span><br><br>
      <b>{n_networks}</b> distinct networks<br>
      <b>{len(connected):,}</b> connected fungi of {len(ecm) + len(am):,} total<br>
      <b>{len(all_edges):,}</b> inferred links (P &gt;= {MAP2_THRESHOLD})<br>
      <span style="font-size:11px;color:#cfcfcf">Line opacity scales with probability.</span><br><br>
      <i style="font-size:11px">Speculative. Same-individual vs same-network<br>
      cannot be told apart without genetic analysis.</i>
    </div>"""))

    m.save("map2_fungal_network.html")
    print(f"  -> map2_fungal_network.html  ({len(all_edges):,} edges, {n_networks} networks)")

# ── Map 3 - Tree-to-Tree Network ──────────────────────────────────────────────

def tree_edges_for_park(park_trees, park_fungi, guild, radius_m):
    """
    For trees and guild-compatible fungi in one park, compute tree-to-tree
    edge probabilities using the complement-product formula:

        P(T1-T2 connected) = 1 - product(1 - p_i)

    where p_i = exp(-d(T1,F)/r) × exp(-d(T2,F)/r) × quality(F)
    for each shared fungus F within radius of both trees.
    """
    if len(park_trees) < 2 or len(park_fungi) < 1:
        return []

    # Downsample trees for performance
    if len(park_trees) > MAX_TREES_NETWORK:
        park_trees = park_trees.sample(MAX_TREES_NETWORK, random_state=42)

    ref_lat = park_trees["lat"].mean()
    ref_lon = park_trees["lon"].mean()

    fungi_xy = to_xy(park_fungi["latitude"].values, park_fungi["longitude"].values,
                     ref_lat, ref_lon)
    trees_xy = to_xy(park_trees["lat"].values, park_trees["lon"].values,
                     ref_lat, ref_lon)

    fungi_kd = cKDTree(fungi_xy)

    # tree_idx → set of nearby fungi indices
    tree_to_fungi = {}
    for t_idx, t_xy in enumerate(trees_xy):
        nearby = fungi_kd.query_ball_point(t_xy, r=radius_m * TREE_RADIUS_MULT)
        if nearby:
            tree_to_fungi[t_idx] = nearby

    # fungus_idx → set of adjacent tree indices
    fungi_to_trees: dict[int, list[int]] = {}
    for t_idx, f_list in tree_to_fungi.items():
        for f_idx in f_list:
            fungi_to_trees.setdefault(f_idx, []).append(t_idx)

    # Find all tree pairs sharing >=1 fungus
    tree_pair_fungi: dict[tuple, list[int]] = {}
    for f_idx, t_list in fungi_to_trees.items():
        t_sorted = sorted(t_list)
        for i in range(len(t_sorted)):
            for j in range(i + 1, len(t_sorted)):
                tree_pair_fungi.setdefault((t_sorted[i], t_sorted[j]), []).append(f_idx)

    trees_reset = park_trees.reset_index(drop=True)
    fungi_reset = park_fungi.reset_index(drop=True)

    edges = []
    for (ta_idx, tb_idx), shared in tree_pair_fungi.items():
        ta = trees_reset.iloc[ta_idx]
        tb = trees_reset.iloc[tb_idx]
        tree_dist = haversine(ta["lat"], ta["lon"], tb["lat"], tb["lon"])
        if tree_dist > 200:
            continue

        # 1 - product(1 - p_i) over all shared fungi
        p_disconnect = 1.0
        for f_idx in shared:
            f = fungi_reset.iloc[f_idx]
            da = haversine(ta["lat"], ta["lon"], f["latitude"], f["longitude"])
            db = haversine(tb["lat"], tb["lon"], f["latitude"], f["longitude"])
            pi = (math.exp(-da / radius_m)
                  * math.exp(-db / radius_m)
                  * quality_w(f.get("quality_grade", "casual")))
            p_disconnect *= (1.0 - pi)
        prob = 1.0 - p_disconnect

        if prob >= MAP3_THRESHOLD:
            edges.append({
                "lat_a": ta["lat"],       "lon_a": ta["lon"],
                "lat_b": tb["lat"],       "lon_b": tb["lon"],
                "probability":    round(prob, 3),
                "species_a":      ta.get("species_latin") or "",
                "species_b":      tb.get("species_latin") or "",
                "shared_fungi":   len(shared),
                "tree_dist_m":    round(tree_dist, 1),
                "guild":          guild,
            })
    return edges

def build_map3(fungi, trees):
    print("\nMap 3 - Tree-to-Tree Network")

    ecm_f = fungi[fungi["guild"] == "ectomycorrhizal"]
    am_f  = fungi[fungi["guild"] == "arbuscular_mycorrhizal"]
    ecm_t = trees[trees["myco_type"] == "ECM"]
    am_t  = trees[trees["myco_type"] == "AM"]

    all_edges = []
    for park in trees["park_key"].unique():
        print(f"  {park}...")
        for guild, pt, pf, r in [
            ("ECM", ecm_t[ecm_t["park_key"] == park], ecm_f[ecm_f["park_key"] == park], ECM_RADIUS_M),
            ("AM",  am_t[am_t["park_key"] == park],   am_f[am_f["park_key"] == park],   AM_RADIUS_M),
        ]:
            edges = tree_edges_for_park(pt, pf, guild, r)
            all_edges.extend(edges)
            print(f"    {guild}: {len(pt):,} trees, {len(pf):,} fungi -> {len(edges):,} tree edges")

    # Group connected trees into networks (union-find); one color per network.
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        parent[find(a)] = find(b)

    def key(lat, lon):
        return (round(lat, 5), round(lon, 5))

    for e in all_edges:
        union(key(e["lat_a"], e["lon_a"]), key(e["lat_b"], e["lon_b"]))

    comp_size = Counter(find(n) for n in parent)
    roots_sorted = [root for root, _ in comp_size.most_common()]
    palette = network_palette(len(roots_sorted))
    root_color = {root: palette[i] for i, root in enumerate(roots_sorted)}
    n_networks = len(roots_sorted)
    connected = set(parent.keys())

    # Light basemap keeps Map 3 distinct; canvas keeps the context layer fast.
    m = folium.Map(location=MAP_CENTER, zoom_start=12, tiles=None, prefer_canvas=True)
    add_basemap(m, ESRI_LIGHT)

    # CONTEXT: every host tree as a faint gray dot behind the network, so the
    # colored networks sit inside the full tree population. FastMarkerCluster
    # keeps the file small; it splits into individual dots once you zoom in.
    context = [[r["lat"], r["lon"]] for _, r in ecm_t.iterrows()]
    context += [[r["lat"], r["lon"]] for _, r in am_t.iterrows()]
    ctx_callback = (
        "function(row){return L.circleMarker([row[0],row[1]],"
        "{radius:2,color:'#9aa0a6',weight:0,fillColor:'#9aa0a6',fillOpacity:0.55});}"
    )
    ctx_cluster = (
        "function(cluster){var n=cluster.getChildCount();"
        "var s=n<100?26:(n<1000?34:44);"
        "return L.divIcon({html:'<div style=\"background:rgba(130,130,130,0.5);"
        "width:'+s+'px;height:'+s+'px;border-radius:50%;border:1px solid #eee;"
        "display:flex;align-items:center;justify-content:center;color:white;"
        "font-family:sans-serif;font-size:11px\">'+n+'</div>',"
        "className:'',iconSize:L.point(s,s)});}"
    )
    FastMarkerCluster(
        context, callback=ctx_callback, icon_create_function=ctx_cluster,
        name=f"All host trees (context, {len(context):,})", show=True,
        disableClusteringAtZoom=14, maxClusterRadius=50,
    ).add_to(m)

    # Links, colored by network, faint so the tree nodes read on top.
    edge_layer = folium.FeatureGroup(name="Fungal links between trees", show=True)
    for e in all_edges:
        folium.PolyLine(
            [[e["lat_a"], e["lon_a"]], [e["lat_b"], e["lon_b"]]],
            color=root_color[find(key(e["lat_a"], e["lon_a"]))],
            weight=1, opacity=0.3,
            tooltip=(f"{e['species_a'] or 'tree'} to {e['species_b'] or 'tree'}  "
                     f"P={e['probability']}  {e['shared_fungi']} shared fungi  "
                     f"{e['tree_dist_m']}m"),
        ).add_to(edge_layer)
    edge_layer.add_to(m)

    # Connected trees, filled with their network color.
    node_layer = folium.FeatureGroup(name="Connected trees (colored by network)", show=True)
    node_pts = []
    for _, r in pd.concat([ecm_t, am_t]).iterrows():
        k = key(r["lat"], r["lon"])
        if k in connected:
            root = find(k)
            node_pts.append((r["lat"], r["lon"], r["park_key"]))
            folium.CircleMarker(
                [r["lat"], r["lon"]], radius=5, color="#333333", weight=0.5,
                fill=True, fill_color=root_color[root], fill_opacity=0.9,
                tooltip=f"{r.get('species_latin') or 'tree'} | network of {comp_size[root]:,} trees",
            ).add_to(node_layer)
    node_layer.add_to(m)

    folium.LayerControl(collapsed=False).add_to(m)

    # Auto-zoom to the dense core (Prospect holds ~95% of the connected trees).
    core = [(lat, lon) for lat, lon, pk in node_pts if pk == "Prospect"]
    if not core:
        core = [(lat, lon) for lat, lon, pk in node_pts]
    if core:
        lats = [c[0] for c in core]
        lons = [c[1] for c in core]
        m.fit_bounds([[min(lats), min(lons)], [max(lats), max(lons)]])

    m.get_root().html.add_child(folium.Element(f"""
    <div style="position:fixed;bottom:30px;left:30px;z-index:9999;
                background:rgba(255,255,255,0.95);padding:14px 18px;border-radius:8px;
                color:#222;font-family:sans-serif;font-size:13px;line-height:1.6;
                border:1px solid #ccc;box-shadow:0 1px 6px rgba(0,0,0,0.2)">
      <b>Map 3 - Tree-to-Tree Network</b><br>
      <i style="font-size:11px">Only trees linked to another tree through a shared fungus</i><br><br>
      <b>Each color = one underground network.</b><br>
      <span style="font-size:11px;color:#555">Same color means those trees connect<br>
      to each other; different colors are separate networks.</span><br><br>
      <b>{n_networks}</b> distinct networks<br>
      <b>{len(node_pts):,}</b> connected trees<br>
      <b>{len(all_edges):,}</b> fungal links (P &gt;= {MAP3_THRESHOLD})<br>
      <span style="font-size:11px;color:#777">Gray dots = all {len(context):,} host trees<br>
      (context). Toggle at right. Colors sit on top.</span><br><br>
      <i style="font-size:11px">All links are ECM; AM fungi are rarely<br>
      seen above ground. Speculative until field data.</i>
    </div>"""))

    m.save("map3_tree_network.html")
    print(f"  -> map3_tree_network.html  ({len(all_edges):,} edges, "
          f"{len(node_pts):,} connected trees, {n_networks} networks, "
          f"{len(context):,} context trees)")

# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    fungi, trees = load_and_enrich()
    build_map1(fungi, trees)
    build_map2(fungi)
    build_map3(fungi, trees)
    print("\nDone.")
    print("  map1_suitability.html   - habitat suitability heatmap")
    print("  map2_fungal_network.html - inferred underground fungal connectivity")
    print("  map3_tree_network.html   - tree-to-tree network via shared fungi")
