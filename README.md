# Carbide architecture and wear in high chromium cast irons

Portable numerical analysis for the current Article 1 dataset, revised on 25 September 2026. The repository contains scientific data, specimen aggregation, Pearson tables, exhaustive ordinary least squares analysis and regression plots.

The supplied population comprises 1,194 accepted optical fields from 15 specimens and 489,253 carbide sections. D2 contributes 97 fields and 43,781 sections. The calibrated field scale is 0.4453125 µm per pixel.

## Run

Python 3.10 or later is required. From the repository directory, create a virtual environment and install the analysis dependencies.

```sh
python -m venv .venv
# Linux or macOS
source .venv/bin/activate
# Windows PowerShell
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python scripts/analyze.py
python scripts/verify.py
```

Use the activation command for your operating system. The scripts also work from another working directory, because default data and output locations are resolved relative to the scripts.

```sh
python scripts/analyze.py --data-dir data --output-dir outputs
python scripts/analyze.py --no-figures
python scripts/verify.py --output-dir outputs
```

`analyze.py` writes numerical outputs and PNG/PDF plots to `outputs/`. `verify.py` independently recomputes the numerical outputs and compares them against the archived current results in `reference/`. A failed comparison raises an error and exits unsuccessfully. Verification is specific to this dataset; custom inputs should use `analyze.py`.

## Included inputs

| File | Contents |
|---|---|
| `data/fields.csv` | One row per accepted field, including spatial descriptors, carbide fraction, section count, calibration and the specimen wear coefficient |
| `data/field_particle_medians.csv` | Field medians of particle circularity and equivalent-circle diameter |
| `data/particle_mean.csv` | Specimen means across all measured particles in that specimen |
| `data/particle_median.csv` | Specimen medians across all measured particles in that specimen |
| `data/particle_std.csv` | Specimen standard deviations across particles, with sample correction `ddof=1` |
| `data/particle_counts.csv` | Section counts per specimen |
| `data/hardness.csv` | The 300 preserved HRC measurements across as-cast, annealed and tempered conditions |
| `reference/specimens.csv` | Archived regression input from the current analysis |
| `reference/specimen_mean.csv`, `reference/specimen_median.csv` | Archived current specimen summaries |
| `reference/model_ranking.csv` | Archived complete ranking of the 381 candidate models |
| `reference/results.json` | Archived current Pearson values, regression coefficients and rounded published results |

CSV files use UTF-8, comma delimiters and decimal points. Original scientific column names are retained to preserve correspondence with the source tables. `Amostra` identifies the specimen, `ID` the field, and `Liga` the alloy group. See [DATA_DICTIONARY.md](DATA_DICTIONARY.md) for units and aggregation rules. Input and source-export checksums are in [data/provenance.json](data/provenance.json).

## Computation

The code implements the calculations used by the original scientific programs `tratamento_dados_TransUNET.py`, `graficos_ingles.py` and `testes_combinacoes_atributos.py`, with portable input/output handling.

For the Pearson tables, field descriptors are aggregated by specimen using their mean, median and sample standard deviation. Particle descriptors use the supplied pooled particle summaries. Each resulting specimen descriptor is correlated with the specimen wear coefficient over 15 specimens. The output contains all 75 Pearson entries from the size/fraction, shape and spatial tables.

The regression dataset uses the median across field-level particle medians for circularity and equivalent diameter. It uses the median across fields for λ and its variance, and the mean across fields for carbide fraction, Delaunay distance, Voronoi area and Voronoi area CV. Thus the regression diameter is not the pooled particle median used in the Pearson tables.

The exhaustive search evaluates all subsets of one through five predictors from nine candidate columns, giving 381 models. Each model has an intercept and is fitted by ordinary least squares. Models are ordered by descending R², descending adjusted R² and ascending predictor count. The reported fits are in-sample fits; this search does not provide independent predictive validation.

The article's selected model remains fixed to λ variance, carbide volume fraction, mean Delaunay distance, mean Voronoi area and mean Voronoi area CV. Its R² is 0.7670221369811973, adjusted R² is 0.6375899908596402, and it ranks fourth. The highest-R² model is retained in the ranking without replacing the selected article model.

The carbide-fraction regression excludes A2 and D2, matching the source analysis. Other individual regressions use all 15 specimens. The morphological parameter is `Dm_median**1.5 * FVC_mean / lambda_median_median`, with carbide fraction expressed as a fraction, not a percentage.

## Outputs

Numerical outputs include regenerated specimen inputs and summaries, the three Pearson groups, the complete model ranking, frequencies among the top 50 models, individual linear/quadratic coefficients and R² values, selected-model predictions and residuals, and `results.json`. The run manifest records package versions and input hashes.

Figures include individual regressions for all nine predictors, the six-panel Figure 8, the carbide-fraction subset regression, the morphological-parameter Figure 9 and the selected-model Figure 10. Figure 10 uses the selected five-predictor model. Plot layout is regenerated and is not claimed to be pixel-identical to the manuscript artwork.

Hardness outputs report means and sample standard deviations by alloy/condition and by tempered specimen. Two further plots show hardness by condition and tempered hardness against wear. The source is the preserved measurement arrays in `testes_estatisticos_durezas.py`. Ten consecutive tempered measurements identify each specimen; the as-cast and annealed arrays are retained at alloy level. Error bars are sample standard deviations, not standard errors. This package does not add hardness hypothesis tests.

## Particle-summary reconstruction

The compact particle summaries were reconstructed from the preserved historical particle export by replacing all D2 rows with the adopted D2 ImageJ particle measurements and rescaling physical lengths from 0.5780346820809248 to 0.4453125 µm per pixel. Areas scale by the squared ratio; perimeter/area scales inversely. Dimensionless shape descriptors are unchanged.

The optional preparation program applies the same transformation when those external scientific inputs are available. It checks every field's particle count against the supplied current field table. The example paths below are placeholders for external scientific inputs.

```sh
python -m pip install -r requirements-optional.txt
python scripts/prepare_particle_summaries.py \
  --particles local_inputs/historical_particles.parquet \
  --replacement-dir local_inputs/D2/EXTRATO \
  --fields data/fields.csv \
  --source-scale 0.5780346820809248 \
  --target-scale 0.4453125 \
  --output-dir outputs/rebuilt_particle_summaries
```

In PowerShell, enter that command on one line or use PowerShell line continuation. The replacement directory must contain only D2 files named `*_EXTRATO.csv`, with field IDs such as `0646-D2C1_EXTRATO.csv`. To summarize an already corrected particle export, omit the replacement directory and set the source and target scales to the same value. This program consumes measured particle tables and does not run image segmentation.

## Reproduction boundary

The included compact data are sufficient to recompute the article's Pearson tables, specimen regression inputs, 381-model ranking, top-50 frequencies and the supplied regression plots. Field geometry is supplied as measured input. Raw micrographs, segmentation masks, full particle tables, segmentation training, ImageJ extraction, LinearDistance extraction and image-overlay generation are outside this package. Optional reconstruction requires external exports; the default analysis does not.

`scripts/geometry.py` also provides the original clipped Voronoi-area and Delaunay-edge algorithms and the mean three-nearest-neighbour distance. `scripts/extract_fields.py` can regenerate accepted-field descriptors from externally supplied `EXTRATO`, `RESUMO` and `LINEAR_DISTANCE` directories. It requires previously validated centroid coordinates in physical units, verifies within-field coordinates and counts, and uses the accepted-field manifest to exclude other fields. It does not repeat centroid-to-mask correspondence or segmentation quality checks.

```sh
python scripts/extract_fields.py --exports-dir local_inputs/D2 --accepted-fields data/fields.csv --specimen D2 --source-scale 0.5780346820809248 --target-scale 0.4453125 --output outputs/rebuilt_D2_fields.csv
```

The geometry implementation computes descriptors in the source scale, then rescales lengths and areas, matching the existing analysis. `X` and `Y` must already use the source physical calibration. Mixing pixel coordinates with physical coordinates is unsupported. The provided D2 export schema has one summary and one LinearDistance record per field.

The archived specimen mean/median tables also contain pooled particle Voronoi area. That column is not used by the article's Pearson or regression analysis and is not recomputed from the compact particle summaries. The field-level Voronoi descriptors used by those analyses are fully included.

Verification checks numerical agreement with the archived current results. It does not establish image segmentation accuracy, calibration traceability, causation or out-of-sample model performance. No raw measurement, scientific claim or model choice is changed by the portability work.

## Verified environment

Numerical verification passed with Python 3.10.11, NumPy 2.2.6, pandas 2.3.3, SciPy 1.15.3, Matplotlib 3.10.9 and statsmodels 0.15.0. The optional reconstruction used PyArrow 25.0.1. All 381 feature combinations retain their archived order, all 75 Pearson entries and 20 individual-regression R² values agree at the archived precision, and the selected-model coefficients and top-50 frequencies agree. Floating-point tolerances are specified in `scripts/verify.py`.

The optional field-extraction program was also executed against all 97 adopted D2 exports. All 11 field descriptors matched the current D2 CSV values exactly after CSV parsing. This comparison covers the adopted D2 exports; it does not claim a new segmentation validation. The numerical verification command also passed when run from outside the repository directory.
