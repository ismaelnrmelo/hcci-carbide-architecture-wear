"""Generate labelled scientific plots from specimen-level inputs."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.text import Annotation
import numpy as np

from scientific_methods import FEATURES, SELECTED_FEATURES, fit_univariate

LABELS = {
    "Circ_median": "Circularity [dimensionless]",
    "Dm_median": r"$D_m$ [$\mu$m]",
    "Black_Median_Variance_Ratio": r"$\lambda/\sigma^2_\lambda$ [$\mu$m$^{-1}$]",
    "lambda_median_median": r"$\lambda$ [$\mu$m]",
    "lambda_var_median": r"$\sigma^2_\lambda$ [$\mu$m$^2$]",
    "FVC_percent": "Carbide volume fraction [%]",
    "Delaunay_mean": r"$\mu_{Del}$ [$\mu$m]",
    "Voronoi_mean": r"$\mu_{Vor}$ [$\mu$m$^2$]",
    "Voronoi_CV_mean": "Voronoi area CV [dimensionless]",
    "zum_gahr": r"$D_m^{3/2}\,\mathrm{CVF}/\lambda$ [$\mu$m$^{1/2}$]",
}
MARKERS = {"A": "s", "B": "^", "C": "o", "D": "D", "E": "P"}
WEAR_LABEL = r"$k$ [$10^{-5}$ mm$^3$/(N m)]"


def save_figure(fig, output_dir, name):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in fig.axes:
        occupied = []
        for annotation in [text for text in ax.texts if isinstance(text, Annotation)]:
            for dx, dy, alignment in [(4, 4, "left"), (-5, -13, "right"),
                    (5, -13, "left"), (-5, 5, "right"), (4, 15, "left"), (-5, -24, "right")]:
                annotation.set_position((dx, dy))
                annotation.set_ha(alignment)
                box = annotation.get_window_extent(renderer).expanded(1.1, 1.15)
                if not any(box.overlaps(other) for other in occupied):
                    break
            occupied.append(box)
    for suffix in ("png", "pdf"):
        fig.savefig(Path(output_dir) / f"{name}.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def labelled_points(ax, data, x, y):
    for specimen, xv, yv in zip(data["Amostra"], x, y):
        ax.scatter(xv, yv, marker=MARKERS.get(specimen[0], "o"),
            color="0.35", edgecolor="black", s=48, zorder=3)
        ax.annotate(specimen, (xv, yv), xytext=(4, 4),
            textcoords="offset points", fontsize=8)
    ax.margins(0.12)
    ax.grid(alpha=0.2)


def regression_panel(ax, data, feature, quadratic=True):
    labelled_points(ax, data, data[feature], data["k_plot"])
    xs = np.linspace(data[feature].min(), data[feature].max(), 250)
    linear = fit_univariate(data, feature)
    ax.plot(xs, linear.params.iloc[0] + linear.params.iloc[1] * xs,
        "k--", lw=1.3, label=rf"Linear, $R^2={linear.rsquared:.3f}$")
    if quadratic:
        model = fit_univariate(data, feature, quadratic=True)
        ax.plot(xs, model.params.iloc[0] + model.params.iloc[1] * xs + model.params.iloc[2] * xs**2,
            "k-", lw=1.3, label=rf"Quadratic, $R^2={model.rsquared:.3f}$")
    ax.set_xlabel(LABELS[feature])
    ax.set_ylabel(WEAR_LABEL)
    ax.legend(fontsize=8, loc="best", framealpha=0.95)


def make_plots(data, selected, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "dejavuserif", "font.size": 10})
    for feature in FEATURES:
        fig, ax = plt.subplots(figsize=(7.2, 4.8), layout="constrained")
        regression_panel(ax, data, feature)
        save_figure(fig, output_dir, f"regression_{feature}")

    figure8_features = ["Circ_median", "Dm_median", "lambda_median_median",
        "lambda_var_median", "Delaunay_mean", "Voronoi_mean"]
    fig, axes = plt.subplots(3, 2, figsize=(12, 13), layout="constrained")
    for letter, feature, ax in zip("abcdef", figure8_features, axes.flat):
        regression_panel(ax, data, feature)
        ax.set_title(f"({letter})", loc="left")
    save_figure(fig, output_dir, "figure8_regressions")

    subset = data[~data["Amostra"].isin(["A2", "D2"])]
    fig, ax = plt.subplots(figsize=(7.2, 4.8), layout="constrained")
    regression_panel(ax, subset, "FVC_percent")
    save_figure(fig, output_dir, "cvf_regression_excluding_A2_D2")

    morph = data.assign(zum_gahr=data["Dm_median"]**1.5 * data["FVC_mean"] / data["lambda_median_median"])
    fig, ax = plt.subplots(figsize=(7.2, 4.8), layout="constrained")
    regression_panel(ax, morph, "zum_gahr", quadratic=False)
    save_figure(fig, output_dir, "figure9_zum_gahr")

    fig, ax = plt.subplots(figsize=(7.2, 5.2), layout="constrained")
    y = data["k_plot"]
    predicted = selected.fittedvalues
    labelled_points(ax, data, y, predicted)
    limits = [min(y.min(), predicted.min()) * 0.9, max(y.max(), predicted.max()) * 1.06]
    ax.plot(limits, limits, "k--", lw=1)
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_xlabel("Measured " + WEAR_LABEL)
    ax.set_ylabel("Fitted " + WEAR_LABEL)
    ax.text(0.04, 0.96, rf"$R^2={selected.rsquared:.3f}$" + "\n" +
        rf"Adjusted $R^2={selected.rsquared_adj:.3f}$", transform=ax.transAxes,
        va="top", bbox={"facecolor": "white", "edgecolor": "0.5", "alpha": 0.95})
    save_figure(fig, output_dir, "figure10_selected_model")


def hardness_plots(by_condition, by_specimen, output_dir):
    fig, ax = plt.subplots(figsize=(7.2, 4.8), layout="constrained")
    positions = np.arange(5)
    for offset, condition, shade in [(-0.25, "as_cast", "0.3"),
            (0, "annealed", "0.65"), (0.25, "tempered", "0.9")]:
        rows = by_condition[by_condition["condition"] == condition].sort_values("alloy")
        ax.bar(positions + offset, rows["mean_HRC"], yerr=rows["std_HRC"],
            width=0.24, capsize=3, color=shade, edgecolor="black", label=condition.replace("_", " "))
    ax.set_xticks(positions, list("ABCDE"))
    ax.set_xlabel("Alloy")
    ax.set_ylabel("Hardness [HRC]")
    ax.legend(title="Condition")
    save_figure(fig, output_dir, "hardness_by_condition")

    fig, ax = plt.subplots(figsize=(7.2, 4.8), layout="constrained")
    ax.errorbar(by_specimen["mean_HRC"], by_specimen["k_plot"],
        xerr=by_specimen["std_HRC"], fmt="none", ecolor="0.6", elinewidth=1, capsize=3)
    labelled_points(ax, by_specimen, by_specimen["mean_HRC"], by_specimen["k_plot"])
    ax.set_xlabel("Tempered hardness [HRC]")
    ax.set_ylabel(WEAR_LABEL)
    save_figure(fig, output_dir, "hardness_vs_wear")
