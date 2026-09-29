"""Specimen aggregation, Pearson correlations and exhaustive OLS analysis."""
from collections import Counter
from itertools import combinations

import numpy as np
import pandas as pd
import statsmodels.api as sm

FIELD_METRICS = [
    "KNN", "λmédio", "λmediano", "λcontagem", "λsoma", "λvariância",
    "Distância Delaunay", "Área de Voronoi (média)",
    "CV Área de Voronoi (imagem)", "Número de carbonetos",
    "Fração volumétrica de carbonetos",
]
PARTICLE_METRICS = [
    "Área", "Perímetro", "Largura (Width)", "Altura (Height)",
    "Eixo maior (Major)", "Eixo menor (Minor)", "Ângulo", "Circularidade",
    "Diâmetro de Feret", "Ângulo de Feret", "Feret mínimo",
    "Razão de aspecto (AR)", "Arredondamento", "Solidez",
    "Superfície específica", "Razão de Feret", "Diâmetro médio",
]
PEARSON_GROUPS = {
    "size_and_fraction": ["Fração volumétrica de carbonetos", "Número de carbonetos",
        "Área", "Largura (Width)", "Altura (Height)", "Feret mínimo",
        "Eixo maior (Major)", "Eixo menor (Minor)", "Diâmetro médio",
        "Diâmetro de Feret", "Perímetro"],
    "shape": ["Circularidade", "Superfície específica", "Arredondamento",
        "Razão de Feret", "Razão de aspecto (AR)", "Solidez"],
    "spatial": ["λmédio", "λmediano", "λcontagem", "λvariância",
        "Área de Voronoi (média)", "CV Área de Voronoi (imagem)",
        "Distância Delaunay", "KNN"],
}
FEATURES = [
    "Circ_median", "Dm_median", "Black_Median_Variance_Ratio",
    "lambda_median_median", "lambda_var_median", "FVC_percent",
    "Delaunay_mean", "Voronoi_mean", "Voronoi_CV_mean",
]
SELECTED_FEATURES = [
    "lambda_var_median", "FVC_percent", "Delaunay_mean",
    "Voronoi_mean", "Voronoi_CV_mean",
]


def read_csv(path):
    return pd.read_csv(path, encoding="utf-8-sig")


def build_specimens(fields, particle_fields):
    """Use median of field medians for size/shape and mean of fields for geometry."""
    if fields["ID"].duplicated().any() or particle_fields["ID"].duplicated().any():
        raise ValueError("Field IDs must be unique")
    if set(fields["ID"]) != set(particle_fields["ID"]):
        raise ValueError("Field and particle-summary IDs do not match")
    merged = fields.merge(particle_fields, on=["ID", "Amostra"], validate="one_to_one")
    if len(merged) != len(fields):
        raise ValueError("Specimen labels do not match between field tables")
    if fields.groupby("Amostra")["k"].nunique().ne(1).any():
        raise ValueError("Wear coefficient must be constant within each specimen")
    data = merged.groupby("Amostra", as_index=False).agg(
        k=("k", "first"), Circ_median=("Circ_img_median", "median"),
        Dm_median=("Dm_img_median", "median"),
        lambda_median_median=("λmediano", "median"),
        lambda_var_median=("λvariância", "median"),
        FVC_mean=("Fração volumétrica de carbonetos", "mean"),
        Delaunay_mean=("Distância Delaunay", "mean"),
        Voronoi_mean=("Área de Voronoi (média)", "mean"),
        Voronoi_CV_mean=("CV Área de Voronoi (imagem)", "mean"),
    )
    data["FVC_percent"] = data["FVC_mean"] * 100
    data["Black_Median_Variance_Ratio"] = data["lambda_median_median"] / data["lambda_var_median"]
    # The supplied wear coefficients already use the display scale of 10^-5 mm^3/(N m).
    data["k_plot"] = data["k"]
    return data.sort_values("Amostra").reset_index(drop=True)


def specimen_summaries(fields, data_dir):
    summaries = {}
    for statistic in ("mean", "median", "std"):
        field_summary = getattr(fields.groupby("Amostra")[FIELD_METRICS], statistic)()
        particle_summary = read_csv(data_dir / f"particle_{statistic}.csv").set_index("Amostra")
        if set(field_summary.index) != set(particle_summary.index):
            raise ValueError(f"Specimen IDs differ in particle_{statistic}.csv")
        summaries[statistic] = pd.concat([field_summary, particle_summary[PARTICLE_METRICS]], axis=1).sort_index()
    return summaries


def pearson_tables(fields, summaries):
    k = fields.groupby("Amostra")["k"].first()
    rows = []
    for group, metrics in PEARSON_GROUPS.items():
        for metric in metrics:
            row = {"group": group, "metric": metric}
            for statistic, summary in summaries.items():
                paired = pd.concat([summary[metric], k.rename("wear")], axis=1).dropna()
                row[statistic] = float(np.corrcoef(paired[metric], paired["wear"])[0, 1])
            rows.append(row)
    return pd.DataFrame(rows)


def exhaustive_models(data, max_features=5):
    """Enumerate untransformed subsets; fit an intercept by NumPy least squares."""
    rows = []
    y = data["k_plot"].to_numpy(float)
    for size in range(1, max_features + 1):
        for features in combinations(FEATURES, size):
            x = data[list(features)].to_numpy(float)
            valid = np.isfinite(y) & np.isfinite(x).all(axis=1)
            xv, yv = x[valid], y[valid]
            n = len(yv)
            if n < 10 or np.any(np.nanstd(xv, axis=0) <= 1e-12):
                continue
            design = np.column_stack([np.ones(n), xv])
            beta, _, rank, _ = np.linalg.lstsq(design, yv, rcond=None)
            sse = float(np.sum((yv - design @ beta) ** 2))
            sst = float(np.sum((yv - yv.mean()) ** 2))
            if sst <= 0:
                continue
            r2 = 1 - sse / sst
            row = dict(SearchType=f"RAW_{size}VARS", n_feat=size, R2=r2,
                Adj_R2=1-(1-r2)*(n-1)/(n-size-1), n_obs=n, rank=int(rank),
                SSE=sse, features=" | ".join(features), b_const=beta[0])
            for i in range(1, 6):
                row[f"f{i}"] = features[i-1] if i <= size else ""
                row[f"b_f{i}"] = beta[i] if i <= size else np.nan
            rows.append(row)
    return pd.DataFrame(rows).sort_values(["R2", "Adj_R2", "n_feat"],
        ascending=[False, False, True]).reset_index(drop=True)


def fit_selected(data):
    return sm.OLS(data["k_plot"], sm.add_constant(data[SELECTED_FEATURES])).fit()


def fit_univariate(data, feature, quadratic=False):
    x = data[feature]
    design = pd.DataFrame({"x": x, "x2": x*x}) if quadratic else x
    return sm.OLS(data["k_plot"], sm.add_constant(design)).fit()


def top_frequency(ranking, n=50):
    counts = Counter()
    for features in ranking.head(n)["features"]:
        counts.update(features.split(" | "))
    return dict(counts.most_common())
