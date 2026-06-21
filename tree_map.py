import calendar
import requests
import pandas as pd
import time
from pathlib import Path

INAT_API = "https://api.inaturalist.org/v1"
FUNGI_TAXON_ID = 47170  # Kingdom Fungi

PARKS = {
    "Central Park": {
        "swlat": 40.7644, "swlng": -73.9816,
        "nelat": 40.8005, "nelng": -73.9493,
    },
    "Prospect Park": {
        "swlat": 40.6544, "swlng": -73.9779,
        "nelat": 40.6804, "nelng": -73.9573,
    },
    "Clove Lakes Park": {
        "swlat": 40.6228, "swlng": -74.1208,
        "nelat": 40.6352, "nelng": -74.1074,
    },
}


INAT_MAX_PAGE = 50  # API hard limit: 50 pages × 200 = 10,000 results per query


def _fetch_window(park_name: str, bounds: dict, d1: str, d2: str) -> tuple[list[dict], int]:
    """Fetch one date-windowed page of results; returns (records, total_in_window)."""
    records = []
    page = 1
    total = None

    while True:
        params = {
            "taxon_id": FUNGI_TAXON_ID,
            "swlat": bounds["swlat"],
            "swlng": bounds["swlng"],
            "nelat": bounds["nelat"],
            "nelng": bounds["nelng"],
            "d1": d1,
            "d2": d2,
            "per_page": 200,
            "page": page,
            "order": "desc",
            "order_by": "observed_on",
        }

        resp = requests.get(f"{INAT_API}/observations", params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()

        results = data.get("results", [])
        if total is None:
            total = data.get("total_results", 0)

        for obs in results:
            taxon = obs.get("taxon") or {}
            location = obs.get("location")
            lat, lng = (None, None)
            if location:
                parts = location.split(",")
                if len(parts) == 2:
                    lat, lng = float(parts[0]), float(parts[1])

            records.append({
                "observation_id":        obs.get("id"),
                "park":                  park_name,
                "observed_on":           obs.get("observed_on"),
                "time_observed":         obs.get("time_observed_at"),
                "latitude":              lat,
                "longitude":             lng,
                "place_guess":           obs.get("place_guess"),
                "quality_grade":         obs.get("quality_grade"),
                "num_id_agreements":     obs.get("num_identification_agreements"),
                "num_id_disagreements":  obs.get("num_identification_disagreements"),
                "taxon_id":              taxon.get("id"),
                "taxon_name":            taxon.get("name"),
                "common_name":           (taxon.get("preferred_common_name") or "").lower() or None,
                "taxon_rank":            taxon.get("rank"),
                "iconic_taxon":          obs.get("iconic_taxon_name"),
                "description":           obs.get("description"),
                "observer_login":        (obs.get("user") or {}).get("login"),
                "url":                   obs.get("uri"),
                "photo_url":             _first_photo(obs),
                "positional_accuracy":   obs.get("positional_accuracy"),
            })

        fetched = (page - 1) * 200 + len(results)
        if page >= INAT_MAX_PAGE or fetched >= total or not results:
            break

        page += 1
        time.sleep(0.5)

    return records, total


def fetch_observations(park_name: str, bounds: dict) -> list[dict]:
    """Fetch all fungi observations within a park's bounding box.

    Chunks by year to stay under the 10,000-result API limit per query.
    Falls back to monthly chunks for years that still exceed the limit.
    """
    all_records: list[dict] = []
    seen_ids: set[int] = set()
    current_year = 2026
    start_year = 2008  # iNaturalist launched publicly in 2008

    print(f"\n  {park_name}: fetching fungi observations (year-by-year)...")

    for year in range(start_year, current_year + 1):
        d1, d2 = f"{year}-01-01", f"{year}-12-31"
        records, total = _fetch_window(park_name, bounds, d1, d2)

        if total > INAT_MAX_PAGE * 200:
            # year exceeds limit — fall back to monthly chunks
            print(f"    {year}: {total} total, chunking by month...")
            records = []
            for month in range(1, 13):
                last_day = calendar.monthrange(year, month)[1]
                md1 = f"{year}-{month:02d}-01"
                md2 = f"{year}-{month:02d}-{last_day:02d}"
                mo_records, mo_total = _fetch_window(park_name, bounds, md1, md2)
                new = [r for r in mo_records if r["observation_id"] not in seen_ids]
                seen_ids.update(r["observation_id"] for r in new)
                records.extend(new)
                if mo_total > 0:
                    print(f"      {year}-{month:02d}: {mo_total} obs, fetched {len(mo_records)}")
                time.sleep(0.5)
        else:
            new = [r for r in records if r["observation_id"] not in seen_ids]
            seen_ids.update(r["observation_id"] for r in new)
            records = new
            if total > 0:
                print(f"    {year}: {total} obs, fetched {len(records)}")

        all_records.extend(records)
        time.sleep(0.5)

    return all_records


def _first_photo(obs: dict) -> str | None:
    photos = obs.get("photos") or []
    if photos:
        return photos[0].get("url")
    return None


def main():
    all_records = []

    for park_name, bounds in PARKS.items():
        records = fetch_observations(park_name, bounds)
        all_records.extend(records)
        print(f"    -> {len(records)} fungi observations collected for {park_name}")

    df = pd.DataFrame(all_records)
    df["observed_on"] = pd.to_datetime(df["observed_on"], errors="coerce")
    df = df.sort_values(["park", "observed_on"], ascending=[True, False])

    out_path = Path("fungi_observations_nyc_parks.csv")
    df.to_csv(out_path, index=False)
    print(f"\nSaved {len(df)} total observations to {out_path}")

    # quick summary
    summary = (
        df.groupby("park")
        .agg(
            total_obs=("observation_id", "count"),
            research_grade=("quality_grade", lambda x: (x == "research grade").sum()),
            unique_taxa=("taxon_name", "nunique"),
            earliest=("observed_on", "min"),
            latest=("observed_on", "max"),
        )
        .reset_index()
    )
    print("\nSummary by park:")
    print(summary.to_string(index=False))

    top_taxa = (
        df.groupby(["park", "taxon_name", "common_name"])
        .size()
        .reset_index(name="count")
        .sort_values(["park", "count"], ascending=[True, False])
        .groupby("park")
        .head(10)
    )
    print("\nTop 10 fungal taxa per park:")
    print(top_taxa.to_string(index=False))

    return df


if __name__ == "__main__":
    main()
