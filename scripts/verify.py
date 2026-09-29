"""Compare recomputed outputs with the archived current scientific results."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analyze import ROOT, analyze
from scientific_methods import SELECTED_FEATURES, read_csv


def verify(data_dir, output_dir, reference_dir):
    results = analyze(data_dir, output_dir, make_figures=False)
    expected = json.loads((reference_dir / "results.json").read_text(encoding="utf-8"))
    checks = {}
    for name in ("specimens", "specimen_mean", "specimen_median", "model_ranking"):
        actual = read_csv(output_dir / f"{name}.csv")
        reference = read_csv(reference_dir / f"{name}.csv")
        if name == "model_ranking":
            pd.testing.assert_series_equal(actual["features"], reference["features"])
            checks["ranking_order_all_381"] = True
        else:
            pd.testing.assert_series_equal(actual["Amostra"], reference["Amostra"])
        numeric = [c for c in actual.select_dtypes("number").columns if c in reference]
        a, b = actual[numeric].to_numpy(), reference[numeric].to_numpy()
        np.testing.assert_allclose(a, b, rtol=1e-9, atol=1e-9, equal_nan=True)
        checks[f"{name}_max_absolute_error"] = float(np.nanmax(np.abs(a-b)))

    pearson = read_csv(output_dir / "pearson_correlations.csv")
    for row in pearson.itertuples():
        np.testing.assert_allclose(np.round([row.mean, row.median, row.std], 3),
            expected["pearson"][row.metric], rtol=0, atol=1e-12)
    checks["pearson_cells_at_published_precision"] = len(pearson) * 3
    regressions = read_csv(output_dir / "individual_regressions.csv")
    for row in regressions.itertuples():
        np.testing.assert_allclose(np.round([row.R2_linear, row.R2_quadratic], 3),
            expected["R2_lin_quad"][row.Feature], rtol=0, atol=1e-12)
    checks["individual_R2_values_at_published_precision"] = len(regressions) * 2
    assert results["top50_frequency"] == expected["top50_freq"]
    selected = results["selected_model"]
    assert selected["rank"] == 4 and selected["features"] == SELECTED_FEATURES
    assert results["fields"] == 1194 and results["sections"] == 489253
    assert results["specimens"] == 15 and results["models"] == 381
    np.testing.assert_allclose([round(selected["R2"], 4), round(selected["Adj_R2"], 4)],
        [expected["modelo_artigo1_no_dado"]["R2"], expected["modelo_artigo1_no_dado"]["Adj_R2"]], rtol=0, atol=1e-12)
    for feature, value in selected["coefficients"].items():
        assert round(value, 5) == expected["modelo_artigo1_no_dado"]["params"][feature]
    np.testing.assert_allclose(np.round(results["cvf_excluding_A2_D2_R2"], 3), expected["fig6_R2"], atol=1e-12, rtol=0)
    assert round(results["zum_gahr_R2"], 4) == expected["zumgahr_R2"]
    fields = read_csv(data_dir / "fields.csv")
    counts = read_csv(data_dir / "particle_counts.csv").set_index("Amostra")["sections"]
    np.testing.assert_array_equal(fields.groupby("Amostra")["Número de carbonetos"].sum().sort_index(), counts.sort_index())
    assert (fields["Amostra"] == "D2").sum() == 97 and int(counts.loc["D2"]) == 43781
    checks.update({"population_counts": True, "selected_model_rank": 4,
        "top50_frequency": True, "selected_model_coefficients": True,
        "cvf_subset_and_morphological_regressions": True, "status": "passed"})
    hardness = read_csv(data_dir / "hardness.csv")
    condition = read_csv(output_dir / "hardness_by_condition.csv")
    specimen = read_csv(output_dir / "hardness_by_specimen.csv")
    assert len(hardness) == 300 and len(condition) == 15 and len(specimen) == 15
    assert (specimen["n"] == 10).all()
    for row in condition.itertuples():
        values = hardness[(hardness["condition"] == row.condition) & (hardness["alloy"] == row.alloy)]["HRC"].to_numpy()
        assert len(values) == row.n
        np.testing.assert_allclose([row.mean_HRC, row.std_HRC],
            [np.mean(values), np.std(values, ddof=1)], rtol=1e-12, atol=1e-12)
    checks["hardness_measurements_and_summaries"] = 300
    (output_dir / "verification.json").write_text(json.dumps(checks, indent=2) + "\n", encoding="utf-8")
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs")
    parser.add_argument("--reference-dir", type=Path, default=ROOT / "reference")
    args = parser.parse_args()
    print(json.dumps(verify(args.data_dir, args.output_dir, args.reference_dir), indent=2))


if __name__ == "__main__":
    main()
