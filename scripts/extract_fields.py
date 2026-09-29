"""Recompute accepted-field descriptors from external ImageJ/LinearDistance tables.

Centroid coordinates and upstream mask acceptance must already be validated.
The three export subdirectories are EXTRATO, RESUMO and LINEAR_DISTANCE.
"""
import argparse
from pathlib import Path
import re

import numpy as np
import pandas as pd

from geometry import voronoi_areas, knn_mean_distance, delaunay_neighbor_distance
from scientific_methods import read_csv


def extract(exports_dir, accepted, source_scale, target_scale, width=640, height=480):
    records = []
    accepted = accepted.set_index("ID")
    seen = set()
    ratio = target_scale / source_scale
    for path in sorted((exports_dir / "EXTRATO").glob("*_EXTRATO.csv")):
        base = path.stem.removesuffix("_EXTRATO")
        match = re.search(r"([A-Z]\d+[A-Z]\d+)$", base)
        if not match or match.group(1) not in accepted.index:
            continue
        code = match.group(1)
        if code in seen:
            raise ValueError(f"Duplicate field export for {code}")
        seen.add(code)
        particles = read_csv(path)
        points = particles[["X", "Y"]].to_numpy(float)
        if not np.isfinite(points).all():
            raise ValueError(f"Non-finite centroid in {path.name}")
        if (points < 0).any() or (points[:, 0] > width*source_scale).any() or (points[:, 1] > height*source_scale).any():
            raise ValueError(f"Centroids exceed the physical field bounds in {path.name}")
        summary = read_csv(exports_dir / "RESUMO" / f"{base}_RESUMO.csv").iloc[0]
        ld = read_csv(exports_dir / "LINEAR_DISTANCE" / f"{base}_LD.csv").iloc[0]
        if int(summary["Count"]) != len(particles) or int(summary["Count"]) != int(accepted.loc[code, "Número de carbonetos"]):
            raise ValueError(f"Particle count differs for {code}")
        areas = voronoi_areas(points, (0, width*source_scale, 0, height*source_scale))
        fraction = float(summary["FVC"])
        if fraction > 1:
            fraction /= 100
        record = {"ID": code, "Liga": accepted.loc[code, "Liga"],
            "Amostra": accepted.loc[code, "Amostra"], "Região": accepted.loc[code, "Região"],
            "Imagem": int(accepted.loc[code, "Imagem"]),
            "KNN": knn_mean_distance(points) * ratio,
            "λmédio": float(ld["Black X and Y Mean"]) * ratio,
            "λmediano": float(ld["Black X and Y Median"]) * ratio,
            "λcontagem": float(ld["Black X and Y Number"]),
            "λsoma": float(ld["Black X and Y Sum"]) * ratio,
            "λvariância": float(ld["Black X and Y Variance"]) * ratio**2,
            "Distância Delaunay": delaunay_neighbor_distance(points) * ratio,
            "Área de Voronoi (média)": np.nanmean(areas) * ratio**2,
            "CV Área de Voronoi (imagem)": np.nanstd(areas, ddof=0) / np.nanmean(areas),
            "Número de carbonetos": int(summary["Count"]),
            "Fração volumétrica de carbonetos": fraction,
            "um_per_px": target_scale, "px_per_um": 1 / target_scale,
            "k": float(accepted.loc[code, "k"])}
        records.append(record)
    missing = set(accepted.index) - seen
    if missing:
        raise ValueError(f"Missing exports for {len(missing)} accepted fields")
    return pd.DataFrame(records).sort_values(["Amostra", "Região", "Imagem"]).reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exports-dir", type=Path, required=True)
    parser.add_argument("--accepted-fields", type=Path, required=True)
    parser.add_argument("--specimen", help="Optionally process only one accepted specimen")
    parser.add_argument("--source-scale", type=float, required=True,
        help="Source micrometres per pixel; X and Y must already be in source micrometres")
    parser.add_argument("--target-scale", type=float, default=0.4453125)
    parser.add_argument("--width", type=int, default=640)
    parser.add_argument("--height", type=int, default=480)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if min(args.source_scale, args.target_scale, args.width, args.height) <= 0:
        parser.error("Scales and field dimensions must be positive")
    accepted = read_csv(args.accepted_fields)
    if args.specimen:
        accepted = accepted[accepted["Amostra"] == args.specimen]
    if accepted.empty or accepted["ID"].duplicated().any():
        parser.error("Accepted field IDs must be nonempty and unique")
    result = extract(args.exports_dir, accepted, args.source_scale, args.target_scale, args.width, args.height)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f"Wrote {len(result)} accepted fields")


if __name__ == "__main__":
    main()
