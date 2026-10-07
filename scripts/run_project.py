"""Single reproducible driver for ME3160 Project 2."""

from pathlib import Path
import sys
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.task1 import run_seasonal
from src.condenser import compare_condenser_types, size_condenser
from src.task1 import CONDITIONS


def main():
    results30 = run_seasonal(0.30)
    results40 = run_seasonal(0.40)

    out_table = ROOT / "results" / "tables"
    out_fig = ROOT / "results" / "figures"
    out_table.mkdir(parents=True, exist_ok=True)
    out_fig.mkdir(parents=True, exist_ok=True)

    df30 = pd.DataFrame(results30)
    df40 = pd.DataFrame(results40)

    df30.to_csv(out_table / "seasonal_30pct_fresh_air.csv", index=False)
    df40.to_csv(out_table / "seasonal_40pct_fresh_air.csv", index=False)

    # Sensitivity table.
    key = "coil_load_zone_kW"
    sens = df40[["condition", key, "Q_evap_kW", "Q_gen_kW",
                 "m_ref_kg_s", "m_fg_kg_s", "condensate_kg_h"]].copy()
    base = df30.set_index("condition")
    for col in sens.columns[1:]:
        sens[col + "_pct_change"] = [
            100.0 * (row[col] - base.loc[row["condition"], col]) /
            base.loc[row["condition"], col]
            for _, row in sens.iterrows()
        ]
    sens.to_csv(out_table / "fresh_air_sensitivity.csv", index=False)

    # Condenser comparison based on the calculated condenser duties.
    q_map = {r["condition"]: r["Q_evap_kW"] + r["Q_gen_kW"]
             for r in results30}
    cond_rows = compare_condenser_types(CONDITIONS, q_map)
    pd.DataFrame(cond_rows).to_csv(
        out_table / "condenser_comparison.csv", index=False
    )

    # Plot cooling load.
    plt.figure()
    plt.bar(df30["condition"], df30["Q_evap_kW"])
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Evaporator duty [kW]")
    plt.title("Seasonal evaporator duty — 30% fresh air")
    plt.tight_layout()
    plt.savefig(out_fig / "seasonal_evaporator_duty.png", dpi=200)
    plt.close()

    # Plot COP.
    plt.figure()
    plt.plot(df30["condition"], df30["COP"], marker="o")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("COP")
    plt.title("Seasonal absorption-chiller COP")
    plt.tight_layout()
    plt.savefig(out_fig / "seasonal_cop.png", dpi=200)
    plt.close()

    print("\nSeasonal results:")
    print(df30[[
        "condition", "Tmix_C", "T_ADP_C", "RH_room",
        "Q_evap_kW", "Q_gen_kW", "COP", "m_ref_kg_s",
        "m_fg_kg_s", "condensate_kg_h"
    ]].to_string(index=False))

    print("\nResults written to:")
    print(out_table)
    print(out_fig)


if __name__ == "__main__":
    main()
