"""Run seasonal calculations and save project results."""

import csv
from pathlib import Path

from src.task1 import run_task1
from src.condenser import compare_condensers
from src.decanter import estimate_condensate

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def save_csv(rows, filename):
    if not rows:
        return

    filename.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(
        key for row in rows for key in row
    ))

    with filename.open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.DictWriter(
            f, fieldnames=fields, extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    seasonal = run_task1()
    save_csv(seasonal, RESULTS / "tables/task1_seasonal.csv")

    condenser = compare_condensers()
    save_csv(condenser, RESULTS / "tables/task2_condenser.csv")

    condensate = []

    for row in seasonal:
        T = row["T_out_C"]
        RH = row["RH_out_percent"] / 100.0

        result = estimate_condensate(T, RH)
        result["condition"] = row["condition"]
        condensate.append(result)

    save_csv(
        condensate,
        RESULTS / "tables/task3_condensate.csv",
    )

    print("Project calculations finished.")
    print("Results directory:", RESULTS)
    print(
        "WARNING: preliminary model; verify the correlations "
        "and design constraints before submitting."
    )


if __name__ == "__main__":
    main()
