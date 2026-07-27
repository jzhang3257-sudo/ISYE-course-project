"""
preprocess.py — CitiBike data preprocessing for ISyE 524 course project.
Reads raw CitiBike trip CSV, aggregates net flow per station,
computes Haversine distance matrix, outputs clean CSV files for Julia.
"""
import csv
import math
from pathlib import Path
from collections import defaultdict

DATA_DIR = Path(__file__).parent.parent / "data"
RAW_CSV = DATA_DIR / "JC-202605-citibike-tripdata.csv"
OUT_STATIONS = DATA_DIR / "stations.csv"
OUT_DISTANCES = DATA_DIR / "distances.csv"


def haversine(lat1, lon1, lat2, lon2):
    """Great-circle distance in km between two lat/lon points."""
    R = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) ** 2 +
         math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) *
         math.sin(dlon / 2) ** 2)
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def main():
    # --- Step 1: Aggregate arrivals and departures per station ---
    arrivals = defaultdict(int)
    departures = defaultdict(int)
    station_info = {}  # station_id -> (name, lat, lon)

    raw_rows = 0
    valid_start_ids = 0
    valid_end_ids = 0
    valid_both_ids = 0

    with open(RAW_CSV, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            raw_rows += 1
            start_id = row["start_station_id"].strip()
            end_id = row["end_station_id"].strip()

            valid_start_ids += bool(start_id)
            valid_end_ids += bool(end_id)
            valid_both_ids += bool(start_id and end_id)

            if start_id:
                departures[start_id] += 1
                if start_id not in station_info:
                    try:
                        station_info[start_id] = (
                            row["start_station_name"].strip(),
                            float(row["start_lat"]),
                            float(row["start_lng"]),
                        )
                    except (ValueError, KeyError):
                        pass

            if end_id:
                arrivals[end_id] += 1
                if end_id not in station_info:
                    try:
                        station_info[end_id] = (
                            row["end_station_name"].strip(),
                            float(row["end_lat"]),
                            float(row["end_lng"]),
                        )
                    except (ValueError, KeyError):
                        pass

    # --- Step 2: Compute net flow ---
    all_stations = sorted(set(list(arrivals.keys()) + list(departures.keys())))
    station_data = []
    missing_coordinate_ids = []
    for sid in all_stations:
        # Do not assign missing stations to (0, 0), since that would create
        # unrealistic Haversine distances. Skip them and report the count.
        if sid not in station_info:
            missing_coordinate_ids.append(sid)
            continue

        net = arrivals.get(sid, 0) - departures.get(sid, 0)
        name, lat, lon = station_info[sid]
        station_data.append((sid, name, lat, lon, net))

    print(f"Raw rows: {raw_rows}")
    print(f"Trips with valid start station ID: {valid_start_ids}")
    print(f"Trips with valid end station ID: {valid_end_ids}")
    print(f"Trips with both station IDs: {valid_both_ids}")
    print(f"Stations skipped for missing coordinates: {len(missing_coordinate_ids)}")
    print(f"Total stations written: {len(station_data)}")
    surplus = sum(1 for s in station_data if s[4] > 0)
    deficit = sum(1 for s in station_data if s[4] < 0)
    balanced = sum(1 for s in station_data if s[4] == 0)
    print(f"Surplus: {surplus}, Deficit: {deficit}, Balanced: {balanced}")

    # --- Step 3: Write stations.csv ---
    with open(OUT_STATIONS, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["station_id", "name", "lat", "lon", "net_flow"])
        for sid, name, lat, lon, net in station_data:
            writer.writerow([sid, name, lat, lon, net])
    print(f"Wrote {OUT_STATIONS}")

    # --- Step 4: Compute and write distance matrix ---
    n = len(station_data)
    ids = [s[0] for s in station_data]
    lats = [s[2] for s in station_data]
    lons = [s[3] for s in station_data]

    with open(OUT_DISTANCES, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([""] + ids)  # header row
        for i in range(n):
            row = [ids[i]]
            for j in range(n):
                d = haversine(lats[i], lons[i], lats[j], lons[j])
                row.append(round(d, 4))
            writer.writerow(row)
    print(f"Wrote {OUT_DISTANCES} ({n}x{n} matrix)")


if __name__ == "__main__":
    main()
