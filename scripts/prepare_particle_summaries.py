"""Prepare compact particle summaries from an export and optional ImageJ replacements.

The historical export and replacement CSVs must share one physical scale.
This script consumes measured particle tables; it does not segment images.
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from scientific_methods import PARTICLE_METRICS, read_csv

IMAGEJ_COLUMNS = {
    "Area": "Área", "Perim.": "Perímetro", "Width": "Largura (Width)",
    "Height": "Altura (Height)", "Major": "Eixo maior (Major)",
    "Minor": "Eixo menor (Minor)", "Angle": "Ângulo", "Circ.": "Circularidade",
    "Feret": "Diâmetro de Feret", "FeretAngle": "Ângulo de Feret",
    "MinFeret": "Feret mínimo", "AR": "Razão de aspecto (AR)",
    "Round": "Arredondamento", "Solidity": "Solidez",
}
LENGTHS = ["Perímetro", "Largura (Width)", "Altura (Height)",
    "Eixo maior (Major)", "Eixo menor (Minor)", "Diâmetro de Feret",
    "Feret mínimo", "Diâmetro médio"]


def prepare(source, replacement_dir, replace_specimen, source_scale, target_scale, fields):
    particles = pd.read_parquet(source)
    if replacement_dir is not None:
        replacements = []
        for path in sorted(replacement_dir.glob("*_EXTRATO.csv")):
            match = re.search(r"([A-Z]\d+[A-Z]\d+)_EXTRATO$", path.stem)
            if match is None:
                raise ValueError(f"Cannot parse field ID from {path.name}")
            field_id = match.group(1)
            specimen = re.match(r"[A-Z]\d+", field_id).group(0)
            if specimen != replace_specimen:
                raise ValueError("Replacement directory contains another specimen")
            raw = read_csv(path)
            table = raw.rename(columns=IMAGEJ_COLUMNS).copy()
            table["ID"], table["Amostra"] = field_id, specimen
            table["Diâmetro médio"] = 2 * np.sqrt(table["Área"] / np.pi)
            table["Superfície específica"] = table["Perímetro"] / table["Área"]
            table["Razão de Feret"] = table["Diâmetro de Feret"] / table["Feret mínimo"]
            replacements.append(table[["ID", "Amostra"] + PARTICLE_METRICS])
        if not replacements:
            raise ValueError("No *_EXTRATO.csv replacement files found")
        particles = pd.concat([particles[particles["Amostra"] != replace_specimen],
            *replacements], ignore_index=True)
    particles = particles[["ID", "Amostra"] + PARTICLE_METRICS].copy()
    ratio = target_scale / source_scale
    particles[LENGTHS] *= ratio
    particles["Área"] *= ratio * ratio
    particles["Superfície específica"] /= ratio
    expected = fields.set_index("ID")["Número de carbonetos"].sort_index()
    observed = particles.groupby("ID").size().sort_index()
    if not expected.index.equals(observed.index) or not np.array_equal(expected.to_numpy(), observed.to_numpy()):
        raise ValueError("Particle counts do not match the supplied field table")
    return particles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--particles", type=Path, required=True)
    parser.add_argument("--fields", type=Path, required=True)
    parser.add_argument("--replacement-dir", type=Path)
    parser.add_argument("--replace-specimen", default="D2")
    parser.add_argument("--source-scale", type=float, required=True, help="Micrometres per pixel used in the source tables")
    parser.add_argument("--target-scale", type=float, default=0.4453125)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if min(args.source_scale, args.target_scale) <= 0:
        parser.error("Scales must be positive")
    particles = prepare(args.particles, args.replacement_dir, args.replace_specimen,
        args.source_scale, args.target_scale, read_csv(args.fields))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for statistic in ("mean", "median", "std"):
        summary = getattr(particles.groupby("Amostra")[PARTICLE_METRICS], statistic)()
        summary.to_csv(args.output_dir / f"particle_{statistic}.csv", encoding="utf-8", index=True)
    field_medians = particles.groupby(["ID", "Amostra"], as_index=False).agg(
        Circ_img_median=("Circularidade", "median"), Dm_img_median=("Diâmetro médio", "median"))
    field_medians.to_csv(args.output_dir / "field_particle_medians.csv", index=False)
    counts = particles.groupby("Amostra").size().rename("sections").reset_index()
    counts.to_csv(args.output_dir / "particle_counts.csv", index=False)
    print(json.dumps({"fields": len(field_medians), "sections": len(particles),
        "source_scale_um_per_px": args.source_scale, "target_scale_um_per_px": args.target_scale}))


if __name__ == "__main__":
    main()
