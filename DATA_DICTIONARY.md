# Data dictionary

The physical scale is 0.4453125 µm per pixel. The analysis uses 640 × 480 pixel optical fields, corresponding to 285 × 213.75 µm. One specimen has one wear coefficient and multiple accepted fields. Field counts differ across specimens.

## Field columns

| Column | Meaning and unit |
|---|---|
| `ID`, `Liga`, `Amostra`, `Região`, `Imagem` | Field ID, alloy group, specimen ID, region code and field index |
| `KNN` | Mean distance to the three nearest particle centroids, µm |
| `λmédio`, `λmediano` | Mean and median matrix free path from LinearDistance, µm |
| `λcontagem` | Number of measured free-path segments |
| `λsoma` | Sum of measured free paths, µm |
| `λvariância` | Variance of measured free paths, µm² |
| `Distância Delaunay` | Mean centroid-neighbour distance on the Delaunay graph, µm |
| `Área de Voronoi (média)` | Mean area of centroid Voronoi cells clipped to the field, µm² |
| `CV Área de Voronoi (imagem)` | Field Voronoi-area population standard deviation divided by mean, dimensionless |
| `Número de carbonetos` | Count of accepted carbide sections in the field |
| `Fração volumétrica de carbonetos` | Measured planar area fraction used as the carbide volume-fraction estimator, fraction from 0 to 1 |
| `um_per_px`, `px_per_um` | Current analysis calibration and reciprocal |
| `centroid_cols`, `centroid_raw_units`, `centroid_swap_xy` | Source centroid-column and coordinate-conversion metadata |
| `centroid_mask_score`, `centroid_inside_frac` | Stored centroid/mask correspondence and within-field fractions from extraction |
| `um_per_px_tiff` | Original TIFF calibration metadata when present; may be empty |
| `k` | Specimen wear coefficient expressed numerically in units of 10⁻⁵ mm³/(N m); these values are already scaled for plotting |

## Particle columns

`Área` is section area in µm². `Perímetro`, `Largura (Width)`, `Altura (Height)`, `Eixo maior (Major)`, `Eixo menor (Minor)`, `Diâmetro de Feret` and `Feret mínimo` are lengths in µm. `Ângulo` and `Ângulo de Feret` are angles in degrees.

`Circularidade`, `Razão de aspecto (AR)`, `Arredondamento` and `Solidez` retain their ImageJ descriptors. `Superfície específica` is perimeter divided by section area in µm⁻¹. `Razão de Feret` is maximum divided by minimum Feret diameter. `Diâmetro médio` is the equivalent-circle diameter `2 * sqrt(area / pi)` in µm. Pooled particle means, medians and standard deviations are computed separately for each specimen. The standard deviation uses `ddof=1`.

## Regression columns

| Column | Aggregation and unit |
|---|---|
| `Circ_median` | Median across fields of each field's particle-median circularity |
| `Dm_median` | Median across fields of each field's particle-median equivalent diameter, µm |
| `lambda_median_median` | Median across fields of field median free path, µm |
| `lambda_var_median` | Median across fields of field free-path variance, µm² |
| `FVC_mean`, `FVC_percent` | Mean field carbide fraction; the second column is that value multiplied by 100 |
| `Delaunay_mean` | Mean field Delaunay distance, µm |
| `Voronoi_mean` | Mean field mean Voronoi area, µm² |
| `Voronoi_CV_mean` | Mean field Voronoi-area coefficient of variation |
| `Black_Median_Variance_Ratio` | `lambda_median_median / lambda_var_median`, µm⁻¹ |
| `k_plot` | Equal to the supplied already-scaled `k` |

The matrix `rank` column in `model_ranking.csv` is the linear-algebra rank of the fitted design matrix. The ranking position is the CSV row order, starting at one. The selected model occupies position four.

## Hardness columns

`hardness.csv` contains condition, alloy group, specimen where identified, sequence number within the alloy/condition array, and `HRC`. The conditions are `as_cast`, `annealed` and `tempered`, corresponding to the preserved source arrays `BRUTA`, `RECOZ` and `TEMP`. As-cast and annealed conditions each have 15 measurements per alloy; tempered conditions have 10 per specimen, or 30 per alloy. Specimen labels for as-cast and annealed records are left empty because only alloy-level identification is used here.
