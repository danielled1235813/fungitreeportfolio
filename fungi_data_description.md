# iNaturalist Fungi Data Pull — Session Notes

## Objective

Pull real fungi observation data from iNaturalist for three NYC parks:
- Central Park (Manhattan)
- Prospect Park (Brooklyn)
- Clove Lakes Park (Staten Island)

---

## What Was Done

### 1. Replaced Stub Script with iNaturalist API Script

The original `script1_load_data.py` was a stub that used `numpy` to generate fake random data and was incomplete (syntax error mid-function). It was fully replaced with a working script that calls the live iNaturalist public API.

**File:** `script1_load_data.py`

### 2. iNaturalist API Query Design

- **API base URL:** `https://api.inaturalist.org/v1/observations`
- **Taxon:** Kingdom Fungi, taxon ID `47170`
- **Location method:** Geographic bounding boxes (`swlat`, `swlng`, `nelat`, `nelng`) defined per park
- **No authentication required** — iNaturalist's public API allows unauthenticated read access

**Park bounding boxes used:**

| Park | SW Corner | NE Corner |
|---|---|---|
| Central Park | 40.7644, -73.9816 | 40.8005, -73.9493 |
| Prospect Park | 40.6544, -73.9779 | 40.6804, -73.9573 |
| Clove Lakes Park | 40.6228, -74.1208 | 40.6352, -74.1074 |

### 3. Handling the 10,000-Result API Limit

The iNaturalist API enforces a hard limit of 50 pages × 200 results = **10,000 results per query**. Central Park alone had ~13,948 fungi observations, which caused a `403 Forbidden` error when the script attempted to fetch page 51.

**Fix:** Queries are chunked by **year** (2008–2026). Each year window is fetched independently. If any single year exceeds 10,000 results, it automatically falls back to **monthly chunks** for that year.

Duplicate observations are de-duplicated using a `seen_ids` set keyed on `observation_id`.

### 4. Data Fields Collected per Observation

| Field | Description |
|---|---|
| `observation_id` | Unique iNaturalist observation ID |
| `park` | Park name |
| `observed_on` | Date of observation |
| `time_observed` | Timestamp (if provided) |
| `latitude` / `longitude` | Coordinates |
| `place_guess` | Observer's free-text location description |
| `quality_grade` | `casual`, `needs_id`, or `research grade` |
| `num_id_agreements` | Number of agreeing identifications |
| `num_id_disagreements` | Number of disagreeing identifications |
| `taxon_id` | iNaturalist taxon ID |
| `taxon_name` | Scientific name |
| `common_name` | Common name (lowercased) |
| `taxon_rank` | Taxonomic rank (species, genus, family, etc.) |
| `iconic_taxon` | Broad category (always Fungi here) |
| `description` | Observer's notes |
| `observer_login` | iNaturalist username |
| `url` | Link to observation on iNaturalist |
| `photo_url` | URL of first photo (if any) |

### 5. Output File

**Path:** `/Users/danielledubov/belowground/fungi_observations_nyc_parks.csv`

All records are sorted by park then by observed date descending.

---

## Results Summary

| Park | Total Observations | Unique Taxa | Date Range |
|---|---|---|---|
| Central Park | 13,933 | 1,032 | 2009–2026 |
| Prospect Park | 10,250 | 1,094 | 2008–2026 |
| Clove Lakes Park | 108 | 80 | 2018–2026 |
| **Total** | **24,291** | — | — |

### Top Fungal Taxa — Central Park
| Taxon | Common Name | Observations |
|---|---|---|
| Candelaria concolor | candleflame lichen | 1,005 |
| Fungi | fungi including lichens | 916 |
| Erysiphe euonymicola | spindletree powdery mildew | 686 |
| Lecanoromycetes | common lichens | 481 |
| Physcia millegrana | rosette lichen | 378 |

### Top Fungal Taxa — Prospect Park
| Taxon | Common Name | Observations |
|---|---|---|
| Fungi | fungi including lichens | 606 |
| Agaricomycetes | mushrooms, bracket fungi, puffballs | 260 |
| Polyporales | shelf fungi | 240 |
| Trametes versicolor | turkey-tail | 197 |
| Auricularia | wood ear fungi | 168 |

### Top Fungal Taxa — Clove Lakes Park
| Taxon | Common Name | Observations |
|---|---|---|
| Fungi | fungi including lichens | 9 |
| Agaricomycetes | mushrooms, bracket fungi, puffballs | 7 |
| Amanita | amanita mushrooms | 3 |
| Auricularia | wood ear fungi | 3 |
| Pleurotus ostreatus | oyster mushroom | 3 |

---

## Pending Work

- [ ] Map fungi observations onto `trees_map.html`
