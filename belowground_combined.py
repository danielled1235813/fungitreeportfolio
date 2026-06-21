#!/usr/bin/env python3
"""
belowground_combined.py
-----------------------
Single-file combined Belowground map.  Run:  python3 belowground_combined.py
Output: belowground_combined.html — open directly in a browser.

Data sources
  Trees : trees_combined.csv  (OSM + NYC Open Data + iNaturalist plants)
  Fungi : fungi_observations_nyc_parks.csv  (iNaturalist, kingdom Fungi)

Layers (toggleable via LayerControl)
  <park> — Trees   diamond markers, species-colored (green/blue hues)
  <park> — Fungi   circle markers, species-colored (orange/red hues)
  ECM fungal connections   static polylines (radius 30 m, P ≥ 0.15)
  AM  fungal connections   static polylines (radius 8 m,  P ≥ 0.15)
  [conn_0 .. conn_10 are slider-only, intentionally absent from LayerControl]

Slider (bottom-right)
  Controls "connection distance" = the fungal search radius r used when
  computing tree–tree edges via shared fungi.  At r=0 no edges exist.
  All 11 edge sets (r=0..10 m) are precomputed in Python and stored as
  separate Leaflet FeatureGroups; the slider only toggles visibility —
  no graph logic runs in the browser.

Distance metric
  haversine() computes all true distances (great-circle metres).
  cKDTree uses to_xy() (equirectangular metres) as a fast candidate filter
  only — an accepted approximation for small radii in a compact area.
  No change from existing code: haversine was already used throughout.

Species grouping decision
  Top 10 species/taxa by count → distinct colors; rest → "Other" (gray).
  Reason for top-10+Other rather than genus collapse:
    - Many fungi are already at genus level in iNaturalist data; collapsing
      further would lose more information than it saves.
    - Top tree species include ecologically distinct Quercus spp. worth
      keeping separate.
  Trees palette: green/blue/teal.  Fungi palette: orange/red/amber.
  Tree palette built only from rows where myco_type != "unresolved" to
  exclude non-tree iNaturalist plants (Artemisia vulgaris, Ficaria verna,
  etc.) that appear in trees_combined.csv.

Positional accuracy
  iNaturalist positional_accuracy (GPS uncertainty radius in metres) is
  fetched in tree_map.py.  If absent from the CSV, popups say "not in
  dataset."  Re-run tree_map.py to populate it.

Note on NetworkX
  The project description mentions a "NetworkX graph model" but none of
  these scripts import NetworkX.  Edges are computed with scipy cKDTree
  for spatial indexing and plain Python dicts for graph traversal.
"""

import math
import warnings
import pandas as pd
import numpy as np
import folium
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
FUNGI_CSV   = "fungi_observations_nyc_parks.csv"
TREES_CSV   = "trees_combined.csv"
OUTPUT_HTML = "belowground_combined.html"

# ── Network radii ─────────────────────────────────────────────────────────────
ECM_RADIUS_M     = 30    # ECM hyphal extent — used for static fungal edge layers
AM_RADIUS_M      = 8     # AM hyphal extent
TREE_RADIUS_MULT = 1.5   # tree-to-fungi search = slider_radius × this

# ── Edge probability thresholds ───────────────────────────────────────────────
MAP2_THRESHOLD = 0.15    # min P for fungal–fungal edge (static layers)
MAP3_THRESHOLD = 0.10    # min P for tree–tree edge (slider layers)

# ── Performance caps ──────────────────────────────────────────────────────────
MAX_TREE_MARKERS    = 2000   # per park; sampled randomly if exceeded
MAX_FUNGI_MARKERS   = 1500   # per park
MAX_EDGES_PER_LAYER = 3000   # keep highest-P tree–tree edges per threshold
MAX_TREES_NETWORK   = 4000   # max trees in edge computation per park/guild

# ── Slider ────────────────────────────────────────────────────────────────────
SLIDER_MIN     = 0
SLIDER_MAX     = 10
SLIDER_DEFAULT = 5   # connection distance shown on load

# ── Display ───────────────────────────────────────────────────────────────────
MAP_CENTER = [40.72, -73.96]   # NYC-wide so all three parks are visible
MAP_ZOOM   = 11

# ── Park name normalisation ───────────────────────────────────────────────────
PARK_NORM = {
    "Central Park":     "Central",
    "Prospect Park":    "Prospect",
    "Clove Lakes Park": "Clove_Lakes",
    "Central":          "Central",
    "Prospect":         "Prospect",
    "Clove_Lakes":      "Clove_Lakes",
}
PARK_DISPLAY = {
    "Central":     "Central Park",
    "Prospect":    "Prospect Park",
    "Clove_Lakes": "Clove Lakes",
}

# ── Species color palettes ────────────────────────────────────────────────────
# Trees: cool green–blue–teal family (clearly distinct from warm fungi colors)
TREE_PALETTE = [
    "#1b5e20",  # forest green
    "#0d47a1",  # navy
    "#00695c",  # dark teal
    "#33691e",  # olive
    "#0277bd",  # cerulean
    "#004d40",  # deep mint
    "#1565c0",  # royal blue
    "#2e7d32",  # medium green
    "#37474f",  # slate
    "#006064",  # dark cyan
]
# Fungi: warm orange–red–amber family (never confusable with tree colors)
FUNGI_PALETTE = [
    "#bf360c",  # burnt sienna
    "#c62828",  # crimson
    "#e65100",  # pumpkin orange
    "#b71c1c",  # dark red
    "#f57f17",  # amber
    "#6d1a00",  # dark brown-red
    "#d84315",  # deep orange
    "#880e4f",  # raspberry
    "#f9a825",  # golden
    "#4e0012",  # burgundy
]
OTHER_COLOR = "#777777"
TOP_N = 10  # species slots with distinct color; remainder → "Other"

# ── Guild / myco-type lookup (genus-based; placeholder for FungalTraits join) ─
ECM_GENERA = {
    "amanita", "boletus", "suillus", "russula", "lactarius", "cantharellus",
    "tricholoma", "cortinarius", "inocybe", "laccaria", "pisolithus",
    "scleroderma", "cenococcum", "hebeloma", "xerocomus", "tylopilus",
    "paxillus", "gyroporus", "boletellus", "strobilomyces", "chalciporus",
    "hygrophorus", "tomentella", "thelephora", "rhizopogon", "gautieria",
    "elaphomyces", "tuber", "hydnum", "sarcodon", "bankera", "clavulina",
    "sebacina", "wilcoxina", "amphinema", "piloderma",
    "leccinum", "chroogomphus", "gomphidius", "truncocolumella",
}
AM_GENERA = {
    "glomus", "rhizophagus", "funneliformis", "claroideoglomus", "gigaspora",
    "scutellospora", "diversispora", "acaulospora", "ambispora", "archaeospora",
    "paraglomus", "redeckera", "septoglomus",
}
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
    "juglans", "carya", "maclura", "ginkgo", "metasequoia", "taxodium", "zelkova",
}


# ── Assignment helpers ────────────────────────────────────────────────────────
def assign_guild(name):
    if not isinstance(name, str):
        return "other"
    g = name.strip().split()[0].lower()
    if g in ECM_GENERA:
        return "ectomycorrhizal"
    if g in AM_GENERA:
        return "arbuscular_mycorrhizal"
    return "other"


def assign_tree_myco(species):
    if not isinstance(species, str) or not species.strip():
        return "unresolved"
    g = species.strip().split()[0].lower()
    if g in ECM_TREE_GENERA:
        return "ECM"
    if g in AM_TREE_GENERA:
        return "AM"
    return "unresolved"


def quality_w(q):
    """iNaturalist quality grade → edge weight multiplier."""
    return {"research": 1.0, "needs_id": 0.7, "casual": 0.4}.get(str(q), 0.5)


def taxon_base_prob(a, b):
    """Base connection probability from taxonomic similarity of two fungal observations."""
    if not isinstance(a, str) or not isinstance(b, str):
        return 0.10
    a, b = a.strip().lower(), b.strip().lower()
    if a == b:
        return 0.75   # same species
    if a.split()[0] == b.split()[0]:
        return 0.45   # same genus
    return 0.12


# ── Geometry ──────────────────────────────────────────────────────────────────
def haversine(lat1, lon1, lat2, lon2):
    """Great-circle distance in metres between two lat/lon points."""
    R = 6_371_000
    p1, p2 = math.radians(lat1), math.radians(lat2)
    a = (math.sin((p2 - p1) / 2) ** 2
         + math.cos(p1) * math.cos(p2) * math.sin(math.radians(lon2 - lon1) / 2) ** 2)
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def to_xy(lats, lons, ref_lat, ref_lon):
    """Equirectangular metres from a reference point.
    Used for cKDTree candidate lookup only — true distances use haversine()."""
    cos_ref = math.cos(math.radians(ref_lat))
    x = (np.asarray(lons) - ref_lon) * 111_320 * cos_ref
    y = (np.asarray(lats) - ref_lat) * 111_320
    return np.column_stack([x, y])


# ── Data loading ──────────────────────────────────────────────────────────────
def load_and_enrich():
    """Load CSVs, coerce coords, assign guild/myco_type, report excluded rows."""
    print("Loading CSVs…")
    fungi = pd.read_csv(FUNGI_CSV, low_memory=False)
    trees = pd.read_csv(TREES_CSV, low_memory=False)

    fungi["latitude"]  = pd.to_numeric(fungi["latitude"],  errors="coerce")
    fungi["longitude"] = pd.to_numeric(fungi["longitude"], errors="coerce")
    trees["lat"] = pd.to_numeric(trees["lat"], errors="coerce")
    trees["lon"] = pd.to_numeric(trees["lon"], errors="coerce")

    f_before, t_before = len(fungi), len(trees)
    fungi = fungi.dropna(subset=["latitude", "longitude"])
    trees = trees.dropna(subset=["lat", "lon"])
    print(f"  Fungi: {f_before} → {len(fungi)}  ({f_before - len(fungi)} dropped, no coords)")
    print(f"  Trees: {t_before} → {len(trees)}  ({t_before - len(trees)} dropped, no coords)")

    # positional_accuracy was added to tree_map.py; may not exist in older CSVs.
    if "positional_accuracy" not in fungi.columns:
        fungi["positional_accuracy"] = None
        print("  Note: positional_accuracy absent from fungi CSV — re-run tree_map.py to fetch it.")

    fungi["park_key"]  = fungi["park"].map(PARK_NORM).fillna(fungi["park"])
    trees["park_key"]  = trees["park"].map(PARK_NORM).fillna(trees["park"])
    fungi["guild"]     = fungi["taxon_name"].apply(assign_guild)
    trees["myco_type"] = trees["species_latin"].apply(assign_tree_myco)

    unresolved = (trees["myco_type"] == "unresolved").sum()
    print(f"  Trees with unresolved myco type: {unresolved:,} "
          f"(non-tree plants + empty genus — shown gray, excluded from edge computation)")
    return fungi, trees


# ── Species color map ─────────────────────────────────────────────────────────
def make_colormap(series, palette):
    """Top-N taxa by count → palette color.  Rest → OTHER_COLOR."""
    top = series.dropna().value_counts().head(TOP_N).index.tolist()
    return {sp: palette[i] for i, sp in enumerate(top)}


# ── Marker shapes ─────────────────────────────────────────────────────────────
def tree_icon(color):
    """Rotated-square (diamond) DivIcon — instantly distinct from circle fungi."""
    return folium.DivIcon(
        html=(f'<div style="width:10px;height:10px;background:{color};'
              f'border:1.5px solid rgba(0,0,0,0.45);transform:rotate(45deg)"></div>'),
        icon_size=(14, 14),
        icon_anchor=(7, 7),
    )


# ── Edge computation — fungal–fungal (static layers) ─────────────────────────
def fungal_edges_for_park(df, guild, radius_m):
    """
    Probabilistic edges between fruiting-body observations of the same guild.

    1. Build cKDTree in equirectangular metres (fast candidate filter).
    2. For each candidate pair within radius_m, compute haversine distance.
    3. P = base_prob × exp(−d / (r/3)) × √(quality_a × quality_b)
    4. Keep edges where P ≥ MAP2_THRESHOLD.
    """
    if len(df) < 2:
        return []
    ref_lat, ref_lon = df["latitude"].mean(), df["longitude"].mean()
    xy    = to_xy(df["latitude"].values, df["longitude"].values, ref_lat, ref_lon)
    pairs = cKDTree(xy).query_pairs(r=radius_m)
    recs  = df.reset_index(drop=True)
    edges = []
    for i, j in pairs:
        a, b  = recs.iloc[i], recs.iloc[j]
        dist  = haversine(a["latitude"], a["longitude"], b["latitude"], b["longitude"])
        base  = taxon_base_prob(a.get("taxon_name"), b.get("taxon_name"))
        decay = math.exp(-dist / (radius_m / 3))
        prob  = (base * decay
                 * math.sqrt(quality_w(a.get("quality_grade", "casual"))
                             * quality_w(b.get("quality_grade", "casual"))))
        if prob >= MAP2_THRESHOLD:
            edges.append({
                "lat_a": a["latitude"], "lon_a": a["longitude"],
                "lat_b": b["latitude"], "lon_b": b["longitude"],
                "probability": round(prob, 3),
                "taxon_a": a.get("taxon_name", ""),
                "taxon_b": b.get("taxon_name", ""),
                "distance_m": round(dist, 1),
                "guild": guild,
            })
    return edges


# ── Edge computation — tree–tree via shared fungi (slider layers) ─────────────
def tree_edges_for_park(park_trees, park_fungi, guild, radius_m):
    """
    Compute tree–tree edges mediated by shared fungal observations within one park.
    Called once per (park, guild, radius_m); the threshold loop calls it 11 times
    per park/guild combination (r = 0, 1, …, 10 m).

    radius_m = 0  →  return [] immediately (no fungal proximity possible at r=0).

    Algorithm:
      1. Build cKDTree of fungal positions in equirectangular metres.
      2. For each tree, collect fungi within radius_m × TREE_RADIUS_MULT metres.
      3. Find all tree pairs sharing ≥1 fungus.
      4. Complement-product formula over shared fungi F:
           P(T1–T2) = 1 − ∏ (1 − pᵢ)
           pᵢ = exp(−d(T1,F)/r) × exp(−d(T2,F)/r) × quality(F)
         All d() calls use haversine (true great-circle metres).
      5. Keep edges where P ≥ MAP3_THRESHOLD and tree–tree distance ≤ 200 m.
    """
    if radius_m == 0 or len(park_trees) < 2 or len(park_fungi) < 1:
        return []

    if len(park_trees) > MAX_TREES_NETWORK:
        park_trees = park_trees.sample(MAX_TREES_NETWORK, random_state=42)

    ref_lat = park_trees["lat"].mean()
    ref_lon = park_trees["lon"].mean()
    fungi_xy = to_xy(park_fungi["latitude"].values, park_fungi["longitude"].values,
                     ref_lat, ref_lon)
    trees_xy = to_xy(park_trees["lat"].values, park_trees["lon"].values,
                     ref_lat, ref_lon)
    fungi_kd = cKDTree(fungi_xy)

    # tree index → list of nearby fungus indices
    tree_to_fungi: dict[int, list] = {}
    for t_idx, t_xy in enumerate(trees_xy):
        nearby = fungi_kd.query_ball_point(t_xy, r=radius_m * TREE_RADIUS_MULT)
        if nearby:
            tree_to_fungi[t_idx] = nearby

    # fungus index → list of adjacent tree indices
    fungi_to_trees: dict[int, list] = {}
    for t_idx, f_list in tree_to_fungi.items():
        for f_idx in f_list:
            fungi_to_trees.setdefault(f_idx, []).append(t_idx)

    # tree pair → list of shared fungus indices
    tree_pair_fungi: dict[tuple, list] = {}
    for f_idx, t_list in fungi_to_trees.items():
        t_sorted = sorted(t_list)
        for i in range(len(t_sorted)):
            for j in range(i + 1, len(t_sorted)):
                tree_pair_fungi.setdefault((t_sorted[i], t_sorted[j]), []).append(f_idx)

    trees_r = park_trees.reset_index(drop=True)
    fungi_r = park_fungi.reset_index(drop=True)
    edges: list = []
    for (ta_idx, tb_idx), shared in tree_pair_fungi.items():
        ta, tb    = trees_r.iloc[ta_idx], trees_r.iloc[tb_idx]
        tree_dist = haversine(ta["lat"], ta["lon"], tb["lat"], tb["lon"])
        if tree_dist > 200:
            continue
        p_disconnect = 1.0
        for f_idx in shared:
            f  = fungi_r.iloc[f_idx]
            da = haversine(ta["lat"], ta["lon"], f["latitude"], f["longitude"])
            db = haversine(tb["lat"], tb["lon"], f["latitude"], f["longitude"])
            pi = (math.exp(-da / radius_m)
                  * math.exp(-db / radius_m)
                  * quality_w(f.get("quality_grade", "casual")))
            p_disconnect *= (1.0 - pi)
        prob = 1.0 - p_disconnect
        if prob >= MAP3_THRESHOLD:
            edges.append({
                "lat_a": ta["lat"],  "lon_a": ta["lon"],
                "lat_b": tb["lat"],  "lon_b": tb["lon"],
                "probability":  round(prob, 3),
                "species_a":    ta.get("species_latin") or "",
                "species_b":    tb.get("species_latin") or "",
                "shared_fungi": len(shared),
                "tree_dist_m":  round(tree_dist, 1),
                "guild":        guild,
            })
    return edges


# ── Threshold loop: precompute all 11 slider edge sets ────────────────────────
def compute_all_thresholds(fungi, trees):
    """
    Run tree_edges_for_park at each integer radius r ∈ {0, 1, …, 10} m,
    for every park × guild combination.  This is the threshold loop.

    At r=0 the function returns [] immediately (short-circuit).
    At r>0 cKDTree finds fungi within r×TREE_RADIUS_MULT m of each tree,
    then the complement-product formula computes tree–tree probabilities.

    Results per threshold are sorted by probability descending and capped
    at MAX_EDGES_PER_LAYER so the HTML stays manageable.

    Returns: dict {radius_int: [edge_dict, …]}
    """
    ecm_f = fungi[fungi["guild"] == "ectomycorrhizal"]
    am_f  = fungi[fungi["guild"] == "arbuscular_mycorrhizal"]
    ecm_t = trees[trees["myco_type"] == "ECM"]
    am_t  = trees[trees["myco_type"] == "AM"]

    all_sets: dict[int, list] = {}
    for T in range(SLIDER_MIN, SLIDER_MAX + 1):
        print(f"  r={T:2d} m … ", end="", flush=True)
        edges: list = []
        for park in trees["park_key"].unique():
            for guild, pt, pf in [
                ("ECM", ecm_t[ecm_t["park_key"] == park], ecm_f[ecm_f["park_key"] == park]),
                ("AM",  am_t[am_t["park_key"] == park],   am_f[am_f["park_key"] == park]),
            ]:
                edges.extend(tree_edges_for_park(pt, pf, guild, T))
        edges.sort(key=lambda e: e["probability"], reverse=True)
        all_sets[T] = edges[:MAX_EDGES_PER_LAYER]
        print(f"{len(all_sets[T]):,} edges")
    return all_sets


# ── Map builder ───────────────────────────────────────────────────────────────
def build_map(fungi, trees, threshold_edges):
    """
    Assemble the Folium map.  Layer order matters for LayerControl:
      1. Tree FeatureGroups (one per park)
      2. Fungi FeatureGroups (one per park)
      3. ECM fungal connections (static)
      4. AM fungal connections (static)
      << LayerControl inserted here — conn_N layers below are excluded >>
      5. conn_0 .. conn_10 (slider-only, not in LayerControl)

    Returns (map, tree_colors, fungi_colors, trees_sampled_out, fungi_sampled_out)
    """
    # Build palettes; restrict tree palette to confirmed-tree genera so
    # non-tree iNat plants (Artemisia, Ficaria…) don't consume color slots.
    known_trees  = trees[trees["myco_type"] != "unresolved"]
    tree_colors  = make_colormap(known_trees["species_latin"], TREE_PALETTE)
    fungi_colors = make_colormap(fungi["taxon_name"], FUNGI_PALETTE)

    has_acc = fungi["positional_accuracy"].notna().any()

    m = folium.Map(location=MAP_CENTER, zoom_start=MAP_ZOOM, tiles="CartoDB positron")

    # ── Tree layers ───────────────────────────────────────────────────────────
    trees_sampled_out = 0
    for park_key, park_display in PARK_DISPLAY.items():
        park_df = trees[trees["park_key"] == park_key]
        if len(park_df) > MAX_TREE_MARKERS:
            trees_sampled_out += len(park_df) - MAX_TREE_MARKERS
            park_df = park_df.sample(MAX_TREE_MARKERS, random_state=42)

        layer = folium.FeatureGroup(
            name=f"{park_display} — Trees ({len(park_df):,})", show=True)

        for _, r in park_df.iterrows():
            sp    = r.get("species_latin") or ""
            color = tree_colors.get(sp, OTHER_COLOR)
            folium.Marker(
                location=[r["lat"], r["lon"]],
                icon=tree_icon(color),
                tooltip=sp or "unknown species",
                popup=folium.Popup(
                    f"<b>{sp or 'unknown'}</b><br>"
                    f"{r.get('species_common') or ''}<br>"
                    f"Myco type: {r.get('myco_type') or 'unresolved'}<br>"
                    f"Health: {r.get('health') or 'N/A'}<br>"
                    f"DBH: {r.get('dbh_cm') or 'N/A'} cm<br>"
                    f"Source: {r.get('provenance') or 'N/A'}<br>"
                    f"Date: {r.get('date_observed') or 'N/A'}",
                    max_width=240,
                ),
            ).add_to(layer)
        layer.add_to(m)

    # ── Fungi layers ──────────────────────────────────────────────────────────
    fungi_sampled_out = 0
    for park_key, park_display in PARK_DISPLAY.items():
        park_df = fungi[fungi["park_key"] == park_key]
        if len(park_df) > MAX_FUNGI_MARKERS:
            fungi_sampled_out += len(park_df) - MAX_FUNGI_MARKERS
            park_df = park_df.sample(MAX_FUNGI_MARKERS, random_state=42)

        layer = folium.FeatureGroup(
            name=f"{park_display} — Fungi ({len(park_df):,})", show=True)

        for _, r in park_df.iterrows():
            sp    = r.get("taxon_name") or ""
            color = fungi_colors.get(sp, OTHER_COLOR)

            # positional_accuracy: GPS uncertainty radius in metres.
            # Present after re-running tree_map.py; None in older CSVs.
            acc = r.get("positional_accuracy")
            if pd.notna(acc) and acc is not None:
                acc_str = f"{int(acc)} m GPS uncertainty"
            else:
                acc_str = "not in dataset — re-run tree_map.py"

            folium.CircleMarker(
                location=[r["latitude"], r["longitude"]],
                radius=5,
                color=color, fill=True, fill_color=color,
                fill_opacity=0.8, weight=0.5,
                tooltip=sp or "unknown taxon",
                popup=folium.Popup(
                    f"<b>{sp or 'unknown'}</b><br>"
                    f"{r.get('common_name') or ''}<br>"
                    f"Guild: {r.get('guild') or 'other'}<br>"
                    f"Quality: {r.get('quality_grade') or 'N/A'}<br>"
                    f"Date: {r.get('observed_on') or 'N/A'}<br>"
                    f"GPS accuracy: {acc_str}<br>"
                    f"<a href='{r.get('url') or '#'}' target='_blank'>iNat link</a>",
                    max_width=250,
                ),
            ).add_to(layer)
        layer.add_to(m)

    # ── Static fungal–fungal edge layers ──────────────────────────────────────
    print("  Computing static fungal–fungal edges…")
    ecm_f = fungi[fungi["guild"] == "ectomycorrhizal"]
    am_f  = fungi[fungi["guild"] == "arbuscular_mycorrhizal"]
    ecm_fe, am_fe = [], []
    for park_key in fungi["park_key"].unique():
        ecm_fe.extend(fungal_edges_for_park(
            ecm_f[ecm_f["park_key"] == park_key], "ectomycorrhizal", ECM_RADIUS_M))
        am_fe.extend(fungal_edges_for_park(
            am_f[am_f["park_key"] == park_key], "arbuscular_mycorrhizal", AM_RADIUS_M))
    print(f"    ECM fungal edges: {len(ecm_fe):,}  AM fungal edges: {len(am_fe):,}")

    ecm_fl = folium.FeatureGroup(
        name=f"ECM fungal connections ({len(ecm_fe):,})", show=False)
    am_fl  = folium.FeatureGroup(
        name=f"AM fungal connections ({len(am_fe):,})",  show=False)
    for e in ecm_fe:
        folium.PolyLine(
            [[e["lat_a"], e["lon_a"]], [e["lat_b"], e["lon_b"]]],
            color="#1565c0", weight=1, opacity=min(e["probability"], 0.85),
            tooltip=f"{e['taxon_a']} ↔ {e['taxon_b']}  P={e['probability']}  {e['distance_m']}m",
        ).add_to(ecm_fl)
    for e in am_fe:
        folium.PolyLine(
            [[e["lat_a"], e["lon_a"]], [e["lat_b"], e["lon_b"]]],
            color="#2e7d32", weight=1, opacity=min(e["probability"], 0.85),
            tooltip=f"{e['taxon_a']} ↔ {e['taxon_b']}  P={e['probability']}  {e['distance_m']}m",
        ).add_to(am_fl)
    ecm_fl.add_to(m)
    am_fl.add_to(m)

    # LayerControl added HERE — only layers registered before this point
    # appear in the control panel.  The conn_N layers added below are
    # managed exclusively by the slider and stay out of the panel.
    folium.LayerControl(collapsed=False).add_to(m)

    # ── Tree–tree connection layers (one per slider threshold) ────────────────
    # BUG FIX vs. original: the old code used eachLayer() to find conn_N layers.
    # Folium does NOT call map.addLayer() for show=False FeatureGroups, so
    # eachLayer() only found the one default-visible layer.  Fix: collect
    # Folium's generated JS variable name for every layer and inject them
    # directly into the slider script as a {threshold: layerObject} mapping.
    conn_layer_vars: dict[int, str] = {}

    for T, edges in sorted(threshold_edges.items()):
        layer = folium.FeatureGroup(
            name=f"conn_{T}",
            show=(T == SLIDER_DEFAULT),  # only the default radius visible on load
        )
        for e in edges:
            color = "#0d47a1" if e["guild"] == "ECM" else "#1b5e20"
            folium.PolyLine(
                [[e["lat_a"], e["lon_a"]], [e["lat_b"], e["lon_b"]]],
                color=color, weight=2,
                opacity=min(e["probability"] * 1.4, 0.85),
                tooltip=(
                    f"{e['species_a'] or 'tree'} ↔ {e['species_b'] or 'tree'}  "
                    f"P={e['probability']}  "
                    f"{e['shared_fungi']} shared fungi  "
                    f"{e['tree_dist_m']} m apart"
                ),
            ).add_to(layer)
        layer.add_to(m)
        # get_name() returns the JS variable name Folium generates, e.g.
        # "feature_group_a1b2c3".  We embed this directly in the script below.
        conn_layer_vars[T] = layer.get_name()

    # ── Slider HTML + JavaScript ──────────────────────────────────────────────
    # conn_layers_js: JS object literal mapping threshold int → FeatureGroup var.
    # Because Folium declares all FeatureGroup variables before html.add_child
    # elements are parsed, these references are valid when the IIFE runs.
    conn_layers_js = "{" + ", ".join(
        f"{T}: {varname}" for T, varname in sorted(conn_layer_vars.items())) + "}"

    # edge_counts_js: precomputed counts for the "N edges" label.
    edge_counts_js = "{" + ", ".join(
        f"{T}: {len(e)}" for T, e in sorted(threshold_edges.items())) + "}"

    acc_note = ("GPS accuracy shown in fungal popups."
                if has_acc else
                "Re-run tree_map.py to add GPS accuracy to fungal popups.")

    slider_html = f"""
<div id="slider_widget" style="
    position:fixed; bottom:30px; right:30px; z-index:9999;
    background:rgba(255,255,255,0.96); padding:13px 17px;
    border-radius:8px; border:1px solid #ccc;
    font-family:sans-serif; font-size:13px; min-width:240px;
    box-shadow:2px 2px 8px rgba(0,0,0,0.18)">
  <b>Tree–tree connections</b><br>
  <span style="color:#666;font-size:11px">via shared fungal network (speculative)</span><br><br>
  Connection distance:&nbsp;<b><span id="dist_label">{SLIDER_DEFAULT}</span> m</b><br>
  <input type="range" id="dist_slider"
         min="{SLIDER_MIN}" max="{SLIDER_MAX}" step="1" value="{SLIDER_DEFAULT}"
         style="width:100%;margin:6px 0 3px 0">
  <div id="edge_count_label" style="color:#888;font-size:11px;min-height:14px"></div>
  <div style="color:#aaa;font-size:10px;margin-top:5px">
    0 m = no edges &nbsp;·&nbsp; line opacity = P<br>
    {acc_note}
  </div>
</div>

<script>
(function() {{
  // ── Slider wiring ──────────────────────────────────────────────────────────
  // Step 1: find the Leaflet map instance.
  // Folium names the map variable "map_<hash>"; scan window for an L.Map.
  var leafletMap = null;
  for (var k in window) {{
    try {{
      if (k.startsWith('map_') && window[k] instanceof L.Map) {{
        leafletMap = window[k]; break;
      }}
    }} catch(e) {{}}
  }}
  if (!leafletMap) {{ console.warn('Belowground: Leaflet map not found'); return; }}

  // Step 2: build the threshold → FeatureGroup mapping.
  // These are direct references to the Folium-generated JS variables
  // (e.g. "feature_group_abc123").  We inject them from Python rather than
  // using eachLayer(), which only iterates layers currently on the map and
  // would miss show=False groups.
  var connLayers = {conn_layers_js};

  // Step 3: precomputed edge counts for the status label.
  var edgeCounts = {edge_counts_js};

  // Step 4: show layer t, hide all others.
  function updateConnections(t) {{
    for (var key in connLayers) {{
      var layer = connLayers[key];
      if (!layer) continue;
      if (parseInt(key, 10) === t) {{
        if (!leafletMap.hasLayer(layer)) leafletMap.addLayer(layer);
      }} else {{
        if (leafletMap.hasLayer(layer))  leafletMap.removeLayer(layer);
      }}
    }}
    document.getElementById('dist_label').textContent = t;
    var n = edgeCounts[t] || 0;
    document.getElementById('edge_count_label').textContent =
        n > 0
          ? n.toLocaleString() + ' edges (P ≥ {MAP3_THRESHOLD})'
          : 'No edges at this radius';
  }}

  // Step 5: wire the range input.
  document.getElementById('dist_slider').addEventListener('input', function() {{
    updateConnections(parseInt(this.value, 10));
  }});

  // Step 6: initialize labels and enforce show/hide state on page load.
  updateConnections({SLIDER_DEFAULT});
}})();
</script>
"""

    # ── Fixed legend ──────────────────────────────────────────────────────────
    def swatch_square(color):
        return (f'<span style="display:inline-block;width:9px;height:9px;'
                f'background:{color};transform:rotate(45deg);'
                f'margin:0 8px 0 1px;vertical-align:middle"></span>')

    def swatch_circle(color):
        return (f'<span style="display:inline-block;width:10px;height:10px;'
                f'border-radius:50%;background:{color};'
                f'margin-right:6px;vertical-align:middle"></span>')

    tree_rows = "".join(
        f'<div style="margin:2px 0">{swatch_square(c)}{sp}</div>'
        for sp, c in tree_colors.items())
    tree_rows += (f'<div style="margin:2px 0">{swatch_square(OTHER_COLOR)}'
                  f'Other / unresolved</div>')

    fungi_rows = "".join(
        f'<div style="margin:2px 0">{swatch_circle(c)}{sp}</div>'
        for sp, c in fungi_colors.items())
    fungi_rows += (f'<div style="margin:2px 0">{swatch_circle(OTHER_COLOR)}'
                   f'Other taxa</div>')

    sample_note = ""
    if trees_sampled_out:
        sample_note += (f'<div style="color:#b00;font-size:10px;margin-top:4px">'
                        f'{trees_sampled_out:,} trees sampled out (perf cap)</div>')
    if fungi_sampled_out:
        sample_note += (f'<div style="color:#b00;font-size:10px">'
                        f'{fungi_sampled_out:,} fungi sampled out</div>')

    legend_html = f"""
<div id="legend" style="
    position:fixed; top:70px; right:30px; z-index:9998;
    background:rgba(255,255,255,0.96); padding:12px 15px;
    border-radius:8px; border:1px solid #ccc;
    font-family:sans-serif; font-size:12px;
    max-height:calc(100vh - 200px); overflow-y:auto; min-width:200px;
    box-shadow:2px 2px 8px rgba(0,0,0,0.18)">
  <b>Belowground</b><br>
  <span style="font-size:10px;color:#666">Shape = kingdom &nbsp;&middot;&nbsp; Color = species</span>
  <hr style="margin:6px 0;border-color:#e0e0e0">
  <b>&#9670; Trees</b><br>
  {tree_rows}
  <hr style="margin:6px 0;border-color:#e0e0e0">
  <b>&#9679; Fungi</b><br>
  {fungi_rows}
  <hr style="margin:6px 0;border-color:#e0e0e0">
  <b>Connections</b><br>
  <div style="margin:2px 0"><span style="color:#0d47a1;font-weight:bold">&#8212;</span>&nbsp;ECM tree&#8211;tree (slider)</div>
  <div style="margin:2px 0"><span style="color:#1b5e20;font-weight:bold">&#8212;</span>&nbsp;AM tree&#8211;tree (slider)</div>
  <div style="margin:2px 0"><span style="color:#1565c0">&#8211; &#8211;</span>&nbsp;ECM fungal (toggle)</div>
  <div style="margin:2px 0"><span style="color:#2e7d32">&#8211; &#8211;</span>&nbsp;AM fungal (toggle)</div>
  {sample_note}
</div>
"""

    m.get_root().html.add_child(folium.Element(slider_html))
    m.get_root().html.add_child(folium.Element(legend_html))
    return m, tree_colors, fungi_colors, trees_sampled_out, fungi_sampled_out


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    fungi, trees = load_and_enrich()

    print("\nPrecomputing tree–tree edges for slider thresholds 0–10 m…")
    threshold_edges = compute_all_thresholds(fungi, trees)

    print("\nBuilding map…")
    m, tree_colors, fungi_colors, trees_out, fungi_out = build_map(
        fungi, trees, threshold_edges)

    m.save(OUTPUT_HTML)
    print(f"\nSaved → {OUTPUT_HTML}")

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\nSpecies grouping (top {TOP_N} + Other):")
    print("  Trees (green/blue palette):")
    for sp, c in tree_colors.items():
        n = (trees["species_latin"] == sp).sum()
        print(f"    {c}  {sp}  ({n:,})")
    print("  Fungi (orange/red palette):")
    for sp, c in fungi_colors.items():
        n = (fungi["taxon_name"] == sp).sum()
        print(f"    {c}  {sp}  ({n:,})")

    if trees_out:
        print(f"\n  Warning: {trees_out:,} tree markers sampled out "
              f"(cap = {MAX_TREE_MARKERS}/park)")
    if fungi_out:
        print(f"  Warning: {fungi_out:,} fungi markers sampled out "
              f"(cap = {MAX_FUNGI_MARKERS}/park)")

    print("\nWorkarounds / data issues:")
    print("  • trees_combined.csv includes non-tree iNaturalist plants "
          "(Artemisia vulgaris, Ficaria verna, etc.).")
    print("    Excluded from the tree species palette (assigned gray 'Other').")
    print("  • positional_accuracy not yet in fungi CSV — re-run tree_map.py to add it.")
    print("  • Central Park trees lack NYC Open Data source; many OSM trees have no species.")
    print(f"  • Slider radius 0–10 m is well below ECM_RADIUS_M={ECM_RADIUS_M} m used for")
    print("    static fungal layers.  Expect sparse tree–tree connections at low values.")
    print("  • No NetworkX used: edges computed with scipy cKDTree + plain Python dicts.")
    m.save('index.html')
