from pathlib import Path
import os

repo = Path(__file__).resolve().parents[1]
project = Path(os.environ.get("CGN_PROJECT_DIR", repo / "workspace")).expanduser().resolve()

required_first_run = [
    "GSE294965_processed_data.h5ad",
]

print("Working directory:", project)
print("\nRequired source input:")
ok = True
for name in required_first_run:
    p = project / name
    state = "OK" if p.exists() else "MISSING"
    print(f"  {state:7s} {p}")
    ok = ok and p.exists()

print("\nGenerated/intermediate files (may be absent before running notebooks):")
for name in [
    "roi_782_PC1_primary.h5ad",
    "roi782_cell_metadata.pkl",
    "figure4_driver_neighbor_data_smoothed.csv",
    "figure5_MAC_FIB_roi_mean_counts.csv",
    "figure5_frozen_module_scores.csv",
]:
    p = project / name
    print(f"  {'PRESENT' if p.exists() else 'not yet'}  {name}")

if not ok:
    raise SystemExit(
        "\nPlace the processed GSE294965 h5ad in the working directory "
        "or set CGN_PROJECT_DIR to the directory containing it."
    )

print("\nCore source input found.")
