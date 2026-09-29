"""Reproduce Article 1 numerical analyses and regression figures."""
import argparse
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd

from scientific_methods import (FEATURES, SELECTED_FEATURES, read_csv,
    build_specimens, specimen_summaries, pearson_tables, exhaustive_models,
    fit_selected, fit_univariate, top_frequency)

ROOT = Path(__file__).resolve().parents[1]


def hardness_summaries(data_dir, specimens):
    raw = read_csv(data_dir / "hardness.csv")
    grouped = raw.groupby(["condition", "alloy"])["HRC"]
    by_condition = grouped.agg(n="count", mean_HRC="mean", std_HRC="std").reset_index()
    by_specimen = (raw[raw["condition"] == "tempered"].groupby("specimen")["HRC"]
        .agg(n="count", mean_HRC="mean", std_HRC="std").reset_index()
        .rename(columns={"specimen": "Amostra"}))
    by_specimen = by_specimen.merge(specimens[["Amostra", "k_plot"]], on="Amostra", validate="one_to_one")
    return by_condition, by_specimen


def analyze(data_dir, output_dir, make_figures=True):
    fields = read_csv(data_dir / "fields.csv")
    medians = read_csv(data_dir / "field_particle_medians.csv")
    specimens = build_specimens(fields, medians)
    if not np.isfinite(specimens.select_dtypes("number").to_numpy()).all():
        raise ValueError("Specimen inputs contain non-finite values")
    summaries = specimen_summaries(fields, data_dir)
    pearson = pearson_tables(fields, summaries)
    ranking = exhaustive_models(specimens)
    selected = fit_selected(specimens)
    selected_name = " | ".join(SELECTED_FEATURES)
    selected_position = int(ranking.index[ranking["features"] == selected_name][0]) + 1
    regression_rows = []
    for feature in FEATURES + ["FVC_mean"]:
        linear, quadratic = fit_univariate(specimens, feature), fit_univariate(specimens, feature, True)
        regression_rows.append({"Feature": feature, "N": len(specimens),
            "R2_linear": linear.rsquared, "R2_quadratic": quadratic.rsquared,
            "linear_intercept": linear.params.iloc[0], "linear_slope": linear.params.iloc[1],
            "quadratic_intercept": quadratic.params.iloc[0],
            "quadratic_linear_term": quadratic.params.iloc[1],
            "quadratic_squared_term": quadratic.params.iloc[2]})
    regressions = pd.DataFrame(regression_rows)
    subset = specimens[~specimens["Amostra"].isin(["A2", "D2"])]
    morphology = specimens.assign(zum_gahr=specimens["Dm_median"]**1.5 * specimens["FVC_mean"] / specimens["lambda_median_median"])
    frequency = top_frequency(ranking)
    results = {
        "fields": len(fields), "specimens": len(specimens),
        "sections": int(fields["Número de carbonetos"].sum()),
        "models": len(ranking), "scale_um_per_px": sorted(fields["um_per_px"].unique().tolist()),
        "selected_model": {"features": SELECTED_FEATURES, "rank": selected_position,
            "R2": selected.rsquared, "Adj_R2": selected.rsquared_adj,
            "coefficients": selected.params.to_dict()},
        "highest_R2_model": {"features": ranking.iloc[0]["features"],
            "R2": float(ranking.iloc[0]["R2"]), "Adj_R2": float(ranking.iloc[0]["Adj_R2"])},
        "top50_frequency": frequency,
        "cvf_excluding_A2_D2_R2": [fit_univariate(subset, "FVC_percent").rsquared,
            fit_univariate(subset, "FVC_percent", True).rsquared],
        "zum_gahr_R2": fit_univariate(morphology, "zum_gahr").rsquared,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    hardness = hardness_summaries(data_dir, specimens) if (data_dir / "hardness.csv").exists() else None
    if hardness is not None:
        hardness[0].to_csv(output_dir / "hardness_by_condition.csv", index=False)
        hardness[1].to_csv(output_dir / "hardness_by_specimen.csv", index=False)
    specimens.to_csv(output_dir / "specimens.csv", index=False)
    for statistic, summary in summaries.items():
        summary.to_csv(output_dir / f"specimen_{statistic}.csv")
    pearson.to_csv(output_dir / "pearson_correlations.csv", index=False)
    ranking.to_csv(output_dir / "model_ranking.csv", index=False)
    regressions.to_csv(output_dir / "individual_regressions.csv", index=False)
    pd.DataFrame(list(frequency.items()), columns=["feature", "top50_count"]).to_csv(output_dir / "top50_frequency.csv", index=False)
    pd.DataFrame({"Amostra": specimens["Amostra"], "measured": specimens["k_plot"],
        "fitted": selected.fittedvalues, "residual": selected.resid}).to_csv(output_dir / "selected_model_predictions.csv", index=False)
    (output_dir / "results.json").write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    manifest = {"python": platform.python_version(),
        "packages": {name: importlib.metadata.version(name) for name in ["numpy", "pandas", "matplotlib", "scipy", "statsmodels"]},
        "input_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(data_dir.glob("*.csv"))}}
    (output_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    if make_figures:
        from plots import make_plots, hardness_plots
        make_plots(specimens, selected, output_dir / "figures")
        if hardness is not None:
            hardness_plots(*hardness, output_dir / "figures")
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs")
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()
    results = analyze(args.data_dir, args.output_dir, not args.no_figures)
    print(json.dumps({"fields": results["fields"], "sections": results["sections"],
        "models": results["models"], "selected_model": results["selected_model"]}, indent=2))


if __name__ == "__main__":
    main()
