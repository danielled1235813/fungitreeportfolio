# Belowground: Mycorrhizal Network Mapping — Product Requirements Document

**Project:** Belowground  
**Date:** 2026-05-02  
**Status:** Phase 1 in progress

---

## Overview

Belowground maps the hidden fungal networks connecting trees in three NYC parks — Central Park, Prospect Park, and Clove Lakes Park. The system fuses publicly available tree and fungal observation data to build a probabilistic graph where nodes are trees and edges represent inferred mycorrhizal connections via shared fungal networks.

The outputs are interactive maps and a graph structure that can answer questions like: *which trees in Central Park are likely exchanging carbon and nutrients through the same fungal intermediary?*

---

## Modeling Approach

Per Dr. Anna Paltseva (Purdue, soil scientist) and consistent with the limitations of observation-based and biomass-proxy data, the three maps are built at increasing levels of inference and speculation. Map 1 is defensible; Maps 2 and 3 are explicitly speculative. Field measurements (soil chemistry + microBIOMETER biomass) serve as edge weight modifiers that make Maps 2 and 3 sensitive to measured soil conditions, and secondarily as a validation signal against model predictions.

**Limitations explicitly acknowledged:** this project does not measure mycorrhizal colonization directly, host identity at fine taxonomic resolution, or physical barriers below the site scale. Without SNP or microsatellite analysis, same-individual and same-network connections cannot be distinguished. iNaturalist observation density is spatially biased toward foot traffic and may not proxy true mycelial presence.

---

## Deliverables

Three maps, each answering a distinct question at increasing levels of inference:

**Map 1 — Fungal Network Habitat Suitability (defensible)**
A heatmap showing where mycorrhizal networks are more or less likely to persist across each park. Driven by soil chemistry (pH, phosphorus), measured fungal biomass (microBIOMETER F:B ratio), and guild-compatible host tree presence. Output is relative: "more likely here, less likely there." Does not claim to show physical connections.

**Map 2 — Underground Fungal Network Inference (speculative)**
Addresses the question: *do these fruiting bodies share an underground mycelial network — are they fruits of the same organism?* A mushroom is the visible tip of a potentially vast underground individual; two fruiting bodies in proximity may be the same genetic individual, different individuals connected by hyphal anastomosis, or entirely separate. This map infers which surface observations are likely connected below ground, using species identity, spatial proximity, guild compatibility, and field biomass as edge weight inputs. Edges are candidates, not confirmations — without SNP or microsatellite analysis, same-individual vs. same-network cannot be distinguished.

**Map 3 — Tree-to-Tree Connectivity Network (speculative)**
Collapses Map 2's fungal edges into a tree-level network. Nodes are trees; edges are inferred mycorrhizal connections via shared fungal networks, weighted by connection probability. Answers: *which trees in this park are likely exchanging carbon and nutrients through the same fungal intermediary?*

---

## Goals

1. Produce a spatially accurate co-location map of all tree and fungi observations across all three parks.
2. Infer probabilistic underground connections between fungal fruiting body observations (Map 2).
3. Collapse fungal edges into tree-to-tree connection probabilities (Map 3).
4. Integrate field measurements (soil chemistry + fungal biomass) as edge weight modifiers in Maps 2 and 3.
5. Surface data gaps (missing guild assignments, unresolved genera, sparse parks) and close them with curated external databases.

---

## Non-Goals

- Real-time data — this is a research/analysis pipeline, not a live system.
- Soil sampling or molecular (eDNA) data integration in the current scope.
- Parks outside the three NYC study sites.
- Modeling nutrient flux quantities; this is network topology, not physiology.

---

## Current State

### Data collected

| Dataset | Records | Source | Coverage |
|---|---|---|---|
| Fungi observations | 24,291 | iNaturalist API | All three parks, 2008–2026 |
| Tree locations | 47,293 | NYC Open Data + OSM + iNaturalist | All three parks |

### Maps produced
- `trees_combined_map.html` — 47k tree points colored by data provenance
- `trees_map.html` — street tree baseline (superseded)

### Scripts
- `script1_load_data.py` — fungi pull from iNaturalist
- `script3_combined_trees.py` — tree pull from three sources, unified schema, Folium map

### Known data quality issues
- Prospect Park and Clove Lakes tree records use bounding boxes, not actual park polygons — includes street trees at park edges.
- iNaturalist plant records include non-tree taxa (shrubs, herbs).
- Central Park has ~19,900 trees; current data covers ~13,900 (OSM + iNat only — no bulk public dataset exists for CPC inventory).
- Many fungi observations are at genus or family rank, not species — reduces host-specificity matching precision.
- No deduplication across sources — same tree may appear in NYC Open Data and iNaturalist as separate rows.

---

## Missing Data — Critical Gaps

### Gap 1: Fungal guild assignments

**What's missing:** The fungi CSV has 1,000+ unique taxa but no column indicating ecological role. Without guild, every fungus looks the same — but only mycorrhizal fungi (ectomycorrhizal or arbuscular mycorrhizal) form tree networks. Saprotrophic, parasitic, and lichen-forming fungi do not.

**What's needed:** A `guild` column on fungi records, one of:
- `ectomycorrhizal` (ECM) — forms networks with oaks, beeches, birches, pines
- `arbuscular_mycorrhizal` (AM) — forms networks with most urban deciduous trees
- `saprotrophic` — decomposes dead matter, no tree network
- `parasitic` — harms host, no mutualistic network
- `lichen` — lichenized, not relevant for tree networks
- `unresolved` — insufficient data

**Data source:** [FungalTraits](https://link.springer.com/article/10.1007/s13225-020-00466-2) (Põlme et al. 2020, *Fungal Diversity*) — assigns primary and secondary guilds to 10,210 fungal genera. Covers the genera present in this dataset. Available as a downloadable spreadsheet via the paper's supplementary materials.

**Implementation:** Join `fungi_observations_nyc_parks.csv` on `taxon_name` genus token → FungalTraits genus → `primary_lifestyle` field. Flag unmatched taxa as `unresolved`.

---

### Gap 2: Fungal genus → tree genus host pairings

**What's missing:** Even among mycorrhizal fungi, host specificity varies. *Suillus* associates almost exclusively with pines; *Amanita* is broadly ECM across many hardwoods; AM fungi are generalists. Without a pairing table, edge probabilities cannot reflect host specificity — a *Suillus* near an oak should get a near-zero edge weight.

**What's needed:** A lookup table mapping fungal genera to compatible tree genera/families, with a specificity score.

**Data sources:**
- **FungalTraits** (same paper) — includes `host_genera` annotations for many taxa.
- **FungalRoot** (Soudzilovskaia et al. 2020, *New Phytologist*) — global plant-mycorrhizal associations database; maps plant species → mycorrhizal type. Inverse lookup (which fungi associate with which plant genera) is derivable. Available at [fungalroot.org](https://fungalroot.org).
- **FUNGuild** (embedded in FungalTraits) — nutrient uptake guild assignments with host notes.
- Literature-derived table: well-characterized ECM host ranges for key genera present in NYC parks (oaks *Quercus*, beeches *Fagus*, birches *Betula*, maples *Acer*, ashes *Fraxinus*).

**Implementation:** Build a `fungal_host_affinities.csv` reference table with columns `fungal_genus`, `tree_genus`, `association_type` (ECM/AM), `specificity` (broad/moderate/strict), `source`. Use this table to weight edges in Phase 2.

---

### Gap 3: Tree mycorrhizal type

**What's missing:** The trees CSV has species names but no column for which mycorrhizal type each tree species supports. An oak and a maple in the same spot may share no fungal network even if the same fungi were observed nearby, because oaks are ECM hosts and most maples are AM hosts.

**What's needed:** A `mycorrhizal_type` column on tree records (ECM / AM / both / non-mycorrhizal).

**Data source:** FungalRoot maps plant species directly to mycorrhizal association. For unresolved species, family-level defaults apply (e.g., Fagaceae → ECM, most Sapindaceae → AM).

**Implementation:** Join `trees_combined.csv` on `species_latin` genus token → FungalRoot → `mycorrhizal_type`. Add a `confidence` flag (species-level match vs. genus-level default vs. family-level default).

---

### Gap 4: Unresolved tree taxa (iNaturalist and OSM records)

**What's missing:** Many iNaturalist plant observations and OSM nodes lack species-level identification (`species_latin` is null or genus-only). These records have coordinates but no guild-assignable species — they're present on the map but invisible to the network model.

**What's needed:** Imputation strategy for unresolved taxa.

**Options (in order of preference):**
1. **Nearest-neighbor imputation** — assign species from the nearest identified tree of the same genus within the same park.
2. **Park-level prior** — assign the most common species for that genus in that park's identified records.
3. **Mark as `unresolved`** and exclude from Phase 2 edge inference, but retain as map points.

---

## Three-Phase Implementation Plan

---

### Phase 1 — Spatial Map (partially complete)

**Goal:** All trees and fungi co-located on a single interactive map per park.

**Remaining work:**

| Task | Priority | Notes |
|---|---|---|
| Overlay fungi observations onto `trees_combined_map.html` | High | Adds a togglable fungi layer |
| Clip Prospect Park + Clove Lakes tree records to actual park polygon | High | Remove street trees outside boundaries; use NYC Parks Properties dataset |
| Filter iNaturalist plant records to tree-only taxa | Medium | Join against known tree genera list or GBIF tree checklist |
| Deduplicate trees across sources by spatial proximity (< 5m) | Medium | Flag likely duplicates; keep highest-quality record |
| Map color legend: add guild layer once Gap 1 is resolved | Low | Hold until FungalTraits join is complete |

**Definition of done:** Single Folium map per park with toggleable layers for trees (colored by species/family) and fungi (colored by guild), polygon-clipped to actual park boundaries.

---

### Phase 2 — Probabilistic Fungal Edges (produces Map 2)

**Biological question:** Do these two fruiting bodies share an underground mycelial network — are they fruits of the same organism?

A mushroom is the visible tip of a potentially vast underground individual. Two surface observations in proximity may be: (a) the same genetic individual expressing multiple fruiting bodies, (b) different individuals whose hyphae have fused via anastomosis into a shared network, or (c) entirely unrelated. Without SNP or microsatellite analysis, (a) and (b) cannot be distinguished. This phase infers whether two observations are *likely connected below ground*, not whether they are genetically identical.

**Edge model inputs:**
1. **Guild match** (Gap 1): Only ECM-ECM or AM-AM pairs are considered. Cross-guild pairs get probability 0.
2. **Taxon identity**: Same species → higher base probability. Same genus → lower. Same family → lowest non-zero prior.
3. **Spatial distance**: Mycelial networks have realistic range limits. Proposed decay function: probability falls off with distance, parameterized by guild (ECM networks can extend 10s of meters; AM networks are typically shorter-range).
4. **Shared host compatibility** (Gap 2): Both observations must be near at least one tree genus that supports their guild. An ECM fungus near only AM-host trees gets a penalized edge weight.
5. **Field data modifier**: Map 1 suitability score at each observation's location scales its edge weights — observations in low-suitability zones (high P, wrong pH, compacted soil, low biomass) receive a downward adjustment on all incident edges.
6. **Observation quality**: `research` grade observations weighted higher than `needs_id` or `casual`.

**Output:** Edge list `fungal_edges.csv` — columns `obs_id_a`, `obs_id_b`, `edge_probability`, `shared_taxon`, `distance_m`, `guild`, `suitability_modifier`, `model_version`.

**Approach:** Start with a simple spatial + guild filter (distance threshold + guild match = binary edge), then iterate toward a probabilistic model once field data is integrated.

---

### Phase 3 — Tree-to-Tree Network (produces Map 3)

**Goal:** Collapse the fungal edge graph into a tree-level network where edge weight between two trees reflects the probability that one or more fungal networks physically connect them.

**Edge model:**

For each pair of trees `(T1, T2)`:
1. Find all fungal observations `F1` within a spatial radius of `T1` that are guild-compatible with `T1`'s mycorrhizal type.
2. Find all fungal observations `F2` within a spatial radius of `T2` that are guild-compatible with `T2`'s mycorrhizal type.
3. For each pair `(F1, F2)` that has a non-trivial Phase 2 edge probability, compute contribution to the `T1–T2` edge.
4. Aggregate across all `(F1, F2)` pairs: `P(T1–T2 connected) = 1 − ∏(1 − P(F1–F2))` (complement of all paths being absent).

**Output:** Tree edge list `tree_edges.csv` — columns `tree_id_a`, `tree_id_b`, `connection_probability`, `supporting_fungal_pairs`, `guild`, `park`.

**Visualization:** Force-directed graph overlaid on park map; edge opacity = connection probability. Filter controls: guild type, probability threshold, species filter.

---

### Phase 4 — Field Data Integration and Validation

**Goal:** Collect soil chemistry (pH, phosphorus) and fungal biomass (microBIOMETER F:B ratio and total biomass) at sample sites across all three parks, and use these measurements in two ways:

**Primary use — edge weight modifier (model input):** Field measurements feed directly into Maps 2 and 3 via the Map 1 suitability surface. Sites with low soil phosphorus, appropriate pH, and high fungal biomass get an upward modifier on incident edge probabilities; sites with high phosphorus, mismatched pH, or low biomass get a downward modifier. This makes Maps 2 and 3 sensitive to measured soil conditions rather than treating all locations as equally likely to support active networks.

**Secondary use — validation:** Site-level predicted connectivity (sum of incident Map 2 edge probabilities at each sample site) is compared against measured fungal biomass via Spearman rank correlation. This tests whether the model's spatial predictions reflect real biology rather than iNaturalist observer foot-traffic patterns.

**Sampling design:** Sites stratified by Map 1 predictions — high-suitability, low-suitability, and zero-suitability zones — across all three parks. Currently 30 sites collected (10 per park, random placement). Targeted re-sampling near identified tree nodes recommended to enable direct linkage between soil measurements and tree records in `trees_combined.csv`.

**Failure modes:** Failure to correlate is a reportable finding, not a project failure. It identifies specific model assumptions that break down — e.g., observation-density bias masking low-biomass sites, or high-biomass zones with no iNaturalist coverage due to access constraints.

---

## External Data Dependencies Summary

| Database | What it provides | URL | License |
|---|---|---|---|
| FungalTraits | Guild (primary_lifestyle) per fungal genus; host annotations | [Põlme et al. 2020](https://link.springer.com/article/10.1007/s13225-020-00466-2) | Open supplementary data |
| FungalRoot | Plant species → mycorrhizal type; inverse gives tree–fungi affinities | [fungalroot.org](https://fungalroot.org) | Open |
| FUNGuild | Alternative guild lookup; embedded in FungalTraits | Included in FungalTraits supplementary | Open |
| NYC Parks Properties | Park polygons for spatial clipping | NYC Open Data | Public domain |

---

## Data Schema — Target State

### `fungi_observations_nyc_parks.csv` (add columns)

| New Column | Type | Source | Notes |
|---|---|---|---|
| `guild` | string | FungalTraits join | ECM / AM / saprotrophic / parasitic / lichen / unresolved |
| `guild_confidence` | string | join metadata | species / genus / family / unresolved |
| `host_genera` | string (list) | FungalTraits / literature | Compatible tree genera, semicolon-delimited |

### `trees_combined.csv` (add columns)

| New Column | Type | Source | Notes |
|---|---|---|---|
| `mycorrhizal_type` | string | FungalRoot join | ECM / AM / both / non-mycorrhizal / unresolved |
| `myco_confidence` | string | join metadata | species / genus / family / unresolved |
| `in_park_polygon` | bool | spatial join | True if inside actual park polygon boundary |
| `is_tree` | bool | taxa filter | True if species_latin matches known tree genera |

---

## Prioritized Next Steps

1. **Download FungalTraits supplementary table** and join to fungi CSV on genus — closes Gap 1, unlocks Phase 2.
2. **Download FungalRoot** and join to trees CSV on species/genus — closes Gap 3.
3. **Clip tree records to park polygons** using NYC Parks Properties dataset — cleans Phase 1 map.
4. **Filter iNaturalist plant records to tree taxa** — reduces noise in tree layer.
5. **Build Phase 2 edge model** (spatial + guild binary filter first, then probabilistic).
6. **Build Phase 3 tree network** once fungal edges are validated.

---

## Open Questions

- **Spatial radius for fungal network extent:** What distance threshold is ecologically defensible for ECM vs. AM networks in urban parks? Literature suggests ECM hyphal networks can extend 10–30m from a tree; AM networks typically 1–5m. These parameters will need tuning.
- **Observation density bias:** iNaturalist observations are citizen-science and spatially biased toward trails and visible specimens. Does observation density proxy for mycelial presence, or does it just reflect human foot traffic?
- **Multi-year observations:** The same fungus at the same GPS point observed in 2015 and 2023 — one fungal individual or two? Temporal clustering needed.
- **Non-mycorrhizal connections:** Some saprotrophic fungi decompose dead wood connecting living trees (e.g., *Trametes versicolor* on fallen logs). Should these be modeled as indirect tree connections?
- **Central Park data ceiling:** With no bulk tree inventory available, Central Park's network will always be sparser than Prospect Park. Is this a limitation to flag or a reason to prioritize Prospect Park for the full analysis?
