#!/usr/bin/env python3
"""
Belowground prototype — produces three mycorrhizal network maps.

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
import warnings
import pandas as pd
import numpy as np
import folium
from folium.plugins import HeatMap
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

MAX_TREE_MARKERS = 2000   # max tree dots per guild on Map 1 (performance)
MAX_TREES_NETWORK = 4000  # max trees considered per park/guild in Map 3

MAP_CENTER = [40.72, -73.96]

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

# ── Map 1 — Habitat Suitability ───────────────────────────────────────────────

def build_map1(fungi, trees):
    print("\nMap 1 — Habitat Suitability Heatmap")

    myco_f  = fungi[fungi["guild"].isin(["ectomycorrhizal", "arbuscular_mycorrhizal"])]
    ecm_t   = trees[trees["myco_type"] == "ECM"]
    am_t    = trees[trees["myco_type"] == "AM"]
    print(f"  Mycorrhizal fungi: {len(myco_f):,}  ECM trees: {len(ecm_t):,}  AM trees: {len(am_t):,}")

    m = folium.Map(location=MAP_CENTER, zoom_start=12, tiles="CartoDB dark_matter")

    # Suitability proxy: density of mycorrhizal fungal observations
    heat_pts = myco_f[["latitude", "longitude"]].values.tolist()
    if heat_pts:
        HeatMap(
            heat_pts,
            name="Suitability proxy (mycorrhizal observation density)",
            min_opacity=0.25,
            radius=22,
            blur=18,
            gradient={0.2: "#0d47a1", 0.4: "#00c853", 0.65: "#ffeb3b", 0.85: "#ff6d00", 1.0: "#b71c1c"},
        ).add_to(m)

    # ECM host trees — sample for performance
    ecm_sample = ecm_t.sample(min(MAX_TREE_MARKERS, len(ecm_t)), random_state=42)
    ecm_layer = folium.FeatureGroup(name=f"ECM host trees (sample n={len(ecm_sample):,})", show=True)
    for _, r in ecm_sample.iterrows():
        folium.CircleMarker(
            [r["lat"], r["lon"]], radius=3,
            color="#4fc3f7", fill=True, fill_opacity=0.65, weight=0,
            tooltip=r.get("species_latin") or "ECM tree",
        ).add_to(ecm_layer)
    ecm_layer.add_to(m)

    # AM host trees — sample for performance
    am_sample = am_t.sample(min(MAX_TREE_MARKERS, len(am_t)), random_state=42)
    am_layer = folium.FeatureGroup(name=f"AM host trees (sample n={len(am_sample):,})", show=False)
    for _, r in am_sample.iterrows():
        folium.CircleMarker(
            [r["lat"], r["lon"]], radius=3,
            color="#81c784", fill=True, fill_opacity=0.65, weight=0,
            tooltip=r.get("species_latin") or "AM tree",
        ).add_to(am_layer)
    am_layer.add_to(m)

    folium.LayerControl().add_to(m)

    m.get_root().html.add_child(folium.Element("""
    <div style="position:fixed;bottom:30px;left:30px;z-index:9999;
                background:rgba(0,0,0,0.82);padding:14px 18px;border-radius:8px;
                color:white;font-family:sans-serif;font-size:13px;line-height:1.6">
      <b>Map 1 — Habitat Suitability</b><br>
      <span style="color:#b71c1c">■</span> High suitability<br>
      <span style="color:#ff6d00">■</span> Moderate<br>
      <span style="color:#ffeb3b">■</span> Low<br>
      <span style="color:#0d47a1">■</span> Minimal signal<br><br>
      <span style="color:#4fc3f7">●</span> ECM host trees<br>
      <span style="color:#81c784">●</span> AM host trees<br><br>
      <i style="font-size:11px">Proxy: iNaturalist mycorrhizal<br>
      observation density. Soil chemistry<br>
      + microBIOMETER data will replace<br>
      this surface in the next iteration.</i>
    </div>"""))

    m.save("map1_suitability.html")
    print("  → map1_suitability.html")

# ── Map 2 — Fungal Network ────────────────────────────────────────────────────

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
    print("\nMap 2 — Underground Fungal Network")

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

    m = folium.Map(location=MAP_CENTER, zoom_start=12, tiles="CartoDB dark_matter")

    # Edge layers
    ecm_edge_layer = folium.FeatureGroup(name="ECM connections", show=True)
    am_edge_layer  = folium.FeatureGroup(name="AM connections",  show=True)
    for e in all_edges:
        color = "#1565c0" if e["guild"] == "ectomycorrhizal" else "#2e7d32"
        layer = ecm_edge_layer if e["guild"] == "ectomycorrhizal" else am_edge_layer
        folium.PolyLine(
            [[e["lat_a"], e["lon_a"]], [e["lat_b"], e["lon_b"]]],
            color=color, weight=1.5,
            opacity=min(e["probability"], 0.9),
            tooltip=(f"{e['taxon_a'] or '?'} ↔ {e['taxon_b'] or '?'}  "
                     f"P={e['probability']}  {e['distance_m']}m"),
        ).add_to(layer)
    ecm_edge_layer.add_to(m)
    am_edge_layer.add_to(m)

    # Observation point layers
    ecm_pt_layer = folium.FeatureGroup(name="ECM fruiting bodies", show=True)
    am_pt_layer  = folium.FeatureGroup(name="AM fruiting bodies",  show=True)
    for _, r in ecm.iterrows():
        folium.CircleMarker(
            [r["latitude"], r["longitude"]], radius=4,
            color="#42a5f5", fill=True, fill_opacity=0.8, weight=0.5,
            tooltip=f"{r.get('taxon_name','')}  ({r.get('quality_grade','')})",
        ).add_to(ecm_pt_layer)
    for _, r in am.iterrows():
        folium.CircleMarker(
            [r["latitude"], r["longitude"]], radius=4,
            color="#66bb6a", fill=True, fill_opacity=0.8, weight=0.5,
            tooltip=f"{r.get('taxon_name','')}  ({r.get('quality_grade','')})",
        ).add_to(am_pt_layer)
    ecm_pt_layer.add_to(m)
    am_pt_layer.add_to(m)

    folium.LayerControl().add_to(m)

    m.get_root().html.add_child(folium.Element(f"""
    <div style="position:fixed;bottom:30px;left:30px;z-index:9999;
                background:rgba(0,0,0,0.82);padding:14px 18px;border-radius:8px;
                color:white;font-family:sans-serif;font-size:13px;line-height:1.6">
      <b>Map 2 — Underground Fungal Network</b><br>
      <i style="font-size:11px">Are these fruiting bodies connected below ground?</i><br><br>
      <span style="color:#42a5f5">●</span> ECM fruiting body<br>
      <span style="color:#66bb6a">●</span> AM fruiting body<br><br>
      <span style="color:#1565c0">—</span> ECM inferred connection<br>
      <span style="color:#2e7d32">—</span> AM inferred connection<br><br>
      Total edges: {len(all_edges):,}  (P ≥ {MAP2_THRESHOLD})<br>
      Line opacity = probability<br><br>
      <i style="font-size:11px">Speculative. Same-individual vs<br>
      same-network indistinguishable<br>
      without genetic analysis.</i>
    </div>"""))

    m.save("map2_fungal_network.html")
    print(f"  → map2_fungal_network.html  ({len(all_edges):,} edges total)")

# ── Map 3 — Tree-to-Tree Network ──────────────────────────────────────────────

def tree_edges_for_park(park_trees, park_fungi, guild, radius_m):
    """
    For trees and guild-compatible fungi in one park, compute tree-to-tree
    edge probabilities using the complement-product formula:

        P(T1–T2 connected) = 1 − ∏(1 − p_i)

    where p_i = exp(−d(T1,F)/r) × exp(−d(T2,F)/r) × quality(F)
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

    # Find all tree pairs sharing ≥1 fungus
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

        # 1 − ∏(1 − p_i) over all shared fungi
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
    print("\nMap 3 — Tree-to-Tree Network")

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
            print(f"    {guild}: {len(pt):,} trees, {len(pf):,} fungi → {len(edges):,} tree edges")

    # Collect coordinates of connected trees for node rendering
    connected = set()
    for e in all_edges:
        connected.add((round(e["lat_a"], 5), round(e["lon_a"], 5), e["guild"]))
        connected.add((round(e["lat_b"], 5), round(e["lon_b"], 5), e["guild"]))

    m = folium.Map(location=MAP_CENTER, zoom_start=12, tiles="CartoDB positron")

    # Edge layers
    ecm_edge_layer = folium.FeatureGroup(name="ECM tree connections", show=True)
    am_edge_layer  = folium.FeatureGroup(name="AM tree connections",  show=True)
    for e in all_edges:
        color = "#0d47a1" if e["guild"] == "ECM" else "#1b5e20"
        layer = ecm_edge_layer if e["guild"] == "ECM" else am_edge_layer
        folium.PolyLine(
            [[e["lat_a"], e["lon_a"]], [e["lat_b"], e["lon_b"]]],
            color=color, weight=2,
            opacity=min(e["probability"] * 1.4, 0.85),
            tooltip=(f"{e['species_a'] or 'tree'} ↔ {e['species_b'] or 'tree'}  "
                     f"P={e['probability']}  "
                     f"{e['shared_fungi']} shared fungi  "
                     f"{e['tree_dist_m']}m"),
        ).add_to(layer)
    ecm_edge_layer.add_to(m)
    am_edge_layer.add_to(m)

    # Node layers — only draw trees that appear in at least one edge
    ecm_node_layer = folium.FeatureGroup(name="ECM trees (connected)", show=True)
    am_node_layer  = folium.FeatureGroup(name="AM trees (connected)",  show=True)

    for _, r in ecm_t.iterrows():
        if (round(r["lat"], 5), round(r["lon"], 5), "ECM") in connected:
            folium.CircleMarker(
                [r["lat"], r["lon"]], radius=5,
                color="#1565c0", fill=True, fill_color="#42a5f5",
                fill_opacity=0.85, weight=1,
                tooltip=r.get("species_latin") or "ECM tree",
            ).add_to(ecm_node_layer)

    for _, r in am_t.iterrows():
        if (round(r["lat"], 5), round(r["lon"], 5), "AM") in connected:
            folium.CircleMarker(
                [r["lat"], r["lon"]], radius=5,
                color="#2e7d32", fill=True, fill_color="#66bb6a",
                fill_opacity=0.85, weight=1,
                tooltip=r.get("species_latin") or "AM tree",
            ).add_to(am_node_layer)

    ecm_node_layer.add_to(m)
    am_node_layer.add_to(m)
    folium.LayerControl().add_to(m)

    m.get_root().html.add_child(folium.Element(f"""
    <div style="position:fixed;bottom:30px;left:30px;z-index:9999;
                background:rgba(255,255,255,0.93);padding:14px 18px;border-radius:8px;
                color:#222;font-family:sans-serif;font-size:13px;line-height:1.6;
                border:1px solid #ccc">
      <b>Map 3 — Tree-to-Tree Network</b><br>
      <i style="font-size:11px">Trees connected via shared fungal networks</i><br><br>
      <span style="color:#1565c0">●</span> ECM tree (connected node)<br>
      <span style="color:#2e7d32">●</span> AM tree (connected node)<br><br>
      <span style="color:#0d47a1">—</span> ECM connection<br>
      <span style="color:#1b5e20">—</span> AM connection<br><br>
      Tree-tree edges: {len(all_edges):,}  (P ≥ {MAP3_THRESHOLD})<br>
      Connected nodes: {len(connected):,}<br><br>
      <i style="font-size:11px">Speculative. Hover edges for species,<br>
      probability, shared fungi count.<br>
      Field data will adjust edge weights<br>
      once integrated.</i>
    </div>"""))

    m.save("map3_tree_network.html")
    print(f"  → map3_tree_network.html  ({len(all_edges):,} edges, {len(connected):,} nodes)")

# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    fungi, trees = load_and_enrich()
    build_map1(fungi, trees)
    build_map2(fungi)
    build_map3(fungi, trees)
    print("\nDone.")
    print("  map1_suitability.html   — habitat suitability heatmap")
    print("  map2_fungal_network.html — inferred underground fungal connectivity")
    print("  map3_tree_network.html   — tree-to-tree network via shared fungi")
