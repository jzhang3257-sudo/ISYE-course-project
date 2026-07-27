"""
generate_subsets.py — Create small/medium problem instances for incremental testing.
Selects the top N stations by absolute net flow (most imbalanced stations),
then filters both stations.csv and distances.csv accordingly.

Usage: python generate_subsets.py
Output: data/subset_10/, data/subset_25/, data/subset_50/
"""
import csv
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"
STATIONS_FILE = DATA_DIR / "stations.csv"
DISTANCES_FILE = DATA_DIR / "distances.csv"

SIZES = [10, 25, 50]


def load_stations(path):
    """Load stations.csv, return list of dicts."""
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = []
        for row in reader:
            row["net_flow"] = int(row["net_flow"])
            rows.append(row)
    return rows


def load_distances(path):
    """Load full distance matrix, return (station_ids, matrix)."""
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        ids = header[1:]  # first col is empty string
        matrix = {}
        for row in reader:
            src = row[0]
            matrix[src] = {}
            for j, dst in enumerate(ids):
                matrix[src][dst] = float(row[j + 1])
    return ids, matrix


def select_stations(stations, n):
    """Select top N stations, balanced between surplus and deficit.

    These are selected stress-test instances, not random or geographically
    representative samples. Results across sizes should be interpreted with
    that limitation in mind.
    Takes the top n_surplus surplus stations (by net_flow descending) and
    top n_deficit deficit stations (by |net_flow| descending).
    """
    surplus_stations = [s for s in stations if s["net_flow"] > 0]
    deficit_stations = [s for s in stations if s["net_flow"] < 0]

    surplus_stations.sort(key=lambda s: s["net_flow"], reverse=True)
    deficit_stations.sort(key=lambda s: s["net_flow"])  # most negative first

    n_surplus = n // 2
    n_deficit = n - n_surplus

    # Don't take more than available
    n_surplus = min(n_surplus, len(surplus_stations))
    n_deficit = min(n_deficit, len(deficit_stations))

    selected = surplus_stations[:n_surplus] + deficit_stations[:n_deficit]
    return selected


def write_subset(subset_dir, selected, full_ids, full_dist_matrix):
    """Write filtered stations.csv and distances.csv for the subset."""
    subset_dir.mkdir(parents=True, exist_ok=True)
    selected_ids = {s["station_id"] for s in selected}

    # Write stations.csv
    with open(subset_dir / "stations.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["station_id", "name", "lat", "lon", "net_flow"])
        for s in selected:
            writer.writerow([s["station_id"], s["name"], s["lat"], s["lon"], s["net_flow"]])

    # Write distances.csv (subset × subset)
    with open(subset_dir / "distances.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        selected_id_list = [s["station_id"] for s in selected]
        writer.writerow([""] + selected_id_list)
        for s_src in selected:
            src_id = s_src["station_id"]
            row = [src_id]
            for s_dst in selected:
                dst_id = s_dst["station_id"]
                row.append(round(full_dist_matrix[src_id][dst_id], 4))
            writer.writerow(row)

    # Summary stats
    surplus = sum(1 for s in selected if s["net_flow"] > 0)
    deficit = sum(1 for s in selected if s["net_flow"] < 0)
    total_supply = sum(s["net_flow"] for s in selected if s["net_flow"] > 0)
    total_demand = sum(-s["net_flow"] for s in selected if s["net_flow"] < 0)
    print(f"  {subset_dir.name}: {len(selected)} stations "
          f"(surplus={surplus}, deficit={deficit}), "
          f"supply={total_supply}, demand={total_demand}")


def main():
    print("Loading full data...")
    stations = load_stations(STATIONS_FILE)
    full_ids, full_dist = load_distances(DISTANCES_FILE)
    print(f"Full dataset: {len(stations)} stations, {len(full_ids)}×{len(full_ids)} distance matrix")

    for n in SIZES:
        subset = select_stations(stations, n)
        subset_dir = DATA_DIR / f"subset_{n}"
        write_subset(subset_dir, subset, full_ids, full_dist)
        print(f"  Wrote {subset_dir}/")

    print("\nDone. Subsets ready for incremental testing.")


if __name__ == "__main__":
    main()