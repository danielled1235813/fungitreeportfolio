# NYC Park Tree Data Pull — Session Notes

## Objective

Collect individual tree locations and species data *inside* park boundaries for three NYC parks:
- Central Park (Manhattan)
- Prospect Park (Brooklyn)
- Clove Lakes Park (Staten Island)

Goal: produce a unified dataset with GPS coordinates, species names, health data, and provenance for each tree, for use in spatial analysis and mapping alongside the fungi observation data.

---

## What Was Done

### 1. Initial Attempt — NYC Street Tree Census (2015)

**File:** `script2_map_trees.py`

The first approach used the NYC 2015 Street Tree Census dataset (`uvpi-gqnh`) from NYC Open Data, which has ~683,788 trees citywide. Trees were fetched using bounding box filters per park and plotted on a Folium map colored by species.

**Problem:** The Street Tree Census only covers trees planted in sidewalk pits along city streets — not trees growing inside park boundaries. The resulting map showed rings of trees around each park perimeter rather than forest interiors.

**Output:** `trees_map.html` — street tree map (superseded by combined map)

---

### 2. Source Research — Interior Park Tree Datasets

A research pass identified the following sources with tree data inside park boundaries:

| Source | Parks Covered | GPS + Species | Access |
|---|---|---|---|
| NYC Open Data Forestry Tree Points (`hn5i-inap`) | Prospect, Clove Lakes | Yes | Socrata REST API |
| OpenStreetMap (Overpass API) | All three (partial) | Partial | Overpass API |
| iNaturalist | All three | Yes (research-grade) | REST API |

**Key finding:** Central Park's ~19,900 trees are managed by the Central Park Conservancy in a proprietary database that is not publicly available for bulk download. OSM and iNaturalist are the best available public sources for Central Park interior trees.

---

### 3. Combined Dataset — Three Sources

**File:** `script3_combined_trees.py`

#### Source 1: NYC Open Data Forestry Tree Points (`hn5i-inap`)

- **URL:** `https://data.cityofnewyork.us/resource/hn5i-inap.json`
- **Coverage:** Prospect Park and Clove Lakes only
- **Query method:** Socrata `within_box(location, NW_lat, NW_lon, SE_lat, SE_lon)` spatial filter with bounding box per park
- **Note:** Socrata URL-encodes `$` in query params — query string must be built manually to avoid 400 errors
- **Filters applied:** `tpstructure != 'Retired' AND tpstructure != 'Stump'` to exclude removed trees
- **Date field used:** `createddate` (inventory entry date; no field-inspection date available)
- **Species field:** `genusspecies` — format `"Acer nigrum - black maple"`, parsed into `species_latin` and `species_common`

#### Source 2: OpenStreetMap (Overpass API)

- **URL:** `https://overpass-api.de/api/interpreter`
- **Coverage:** All three parks (volunteer-mapped; coverage varies)
- **Query:** `node[natural=tree]` within bounding box per park
- **Note:** Requires `User-Agent` header to avoid 406 errors
- **Date field:** OSM element `timestamp` (last edit date, not observation date)
- **Species tags used (when present):** `species`, `taxon`, `genus` → `species_latin`; `species:en`, `name` → `species_common`

#### Source 3: iNaturalist

- **URL:** `https://api.inaturalist.org/v1/observations`
- **Coverage:** All three parks
- **Filters:** `iconic_taxa=Plantae`, `quality_grade=research`, paginated 200 per page
- **Place IDs used:**

| Park | iNaturalist Place ID |
|---|---|
| Central Park | 49955 |
| Prospect Park | 55174 |
| Clove Lakes Park | 125420 |

- **Cap:** 10,000 observations per park (API hard limit); Central Park and Prospect Park both hit this cap
- **Date field:** `observed_on` (actual observation date)
- **Note:** Includes all research-grade plant observations — not filtered to trees only; may include shrubs, herbaceous plants, etc.

---

### 4. Unified Schema

Every record across all three sources is normalized to the following columns:

| Column | Description |
|---|---|
| `park` | Park name (`Central`, `Prospect`, `Clove_Lakes`) |
| `lat` / `lon` | GPS coordinates |
| `species_latin` | Scientific name |
| `species_common` | Common name |
| `dbh_cm` | Trunk diameter at breast height (cm); available for NYC Open Data records only |
| `health` | Condition rating; available for NYC Open Data records only |
| `provenance` | Source dataset (see values below) |
| `date_observed` | Date of observation/inventory entry (ISO `YYYY-MM-DD`) |

**Provenance values:**
- `NYC Open Data - Forestry Tree Points`
- `OpenStreetMap`
- `iNaturalist`

---

### 5. Output Files

| File | Description |
|---|---|
| `trees_combined.csv` | Unified dataset, all sources, all parks |
| `trees_combined_map.html` | Folium map — points colored by provenance, clustered by park |

**Map colors:**
- Blue (`#2196F3`) — NYC Open Data
- Green (`#4CAF50`) — OpenStreetMap
- Orange (`#FF9800`) — iNaturalist

---

## Results Summary

| Park | Source | Records |
|---|---|---|
| Central | OpenStreetMap | 3,898 |
| Central | iNaturalist | 10,000 (capped) |
| Prospect | NYC Open Data | 18,184 |
| Prospect | iNaturalist | 10,000 (capped) |
| Clove Lakes | NYC Open Data | 3,129 |
| Clove Lakes | OpenStreetMap | 4 |
| Clove Lakes | iNaturalist | 1,465 |
| **Total** | | **47,293** |

---

## Known Limitations

- **Central Park:** No NYC Open Data coverage. OSM has ~3,900 of ~19,900 actual trees; species tagging is inconsistent. iNaturalist is capped at 10,000 and skews toward visible/flowering specimens.
- **Prospect Park bounding box:** The `within_box` spatial filter uses a rectangular bounding box, not the actual park polygon. Some records near the park boundary may be street trees in the surrounding neighborhood rather than interior trees.
- **iNaturalist Plantae observations:** Not filtered to trees specifically — includes shrubs and other plants. Filter `species_latin` against known tree genera to narrow if needed.
- **OSM dates:** The `timestamp` field is the date of the last OSM edit, not the date the tree was observed or planted.
- **No deduplication across sources:** The same tree may appear in both NYC Open Data and iNaturalist as separate rows with different provenance. This is intentional — provenance is tracked per row.

---

## Pending Work

- [ ] Clip Prospect Park and Clove Lakes records to actual park polygon boundaries (NYC Parks Properties dataset) to remove surrounding street trees
- [ ] Filter iNaturalist records to tree-only taxa
- [ ] Overlay tree locations with fungi observation data (`fungi_observations_nyc_parks.csv`)
