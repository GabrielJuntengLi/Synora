# **Synora: Spatial Tissue Architecture Analysis for Omics Data**

Synora — named after the Greek term for *boundary* — is an R package that resolves tissue architecture in spatial omics data. Rather than simply classifying cell types, Synora quantifies biologically meaningful interfaces, enabling downstream analyses that depend on knowing **where** a cell is relative to a tissue boundary, not just **what type** it is.


By introducing a novel vector-based **Orientedness** metric, Synora distinguishes true structured boundaries from random cellular infiltration — a distinction that conventional heterogeneity-only methods cannot make.

## Installation

```r
# Install from GitHub (requires remotes)
# Stable version
remotes::install_github("lzxlab/Synora")

# Development version (experimental features)
remotes::install_github("GabrielJuntengLi/Synora")
```

---

## What Synora Enables

### 🔬 Boundary-Aware Gene & Protein Expression
By assigning every cell a precise distance to the nearest tissue boundary, Synora enables:
- Identification of spatially graded gene/protein expression across tissue compartments
- Discovery of boundary-enriched or boundary-depleted molecular signatures
- Pseudospace analysis: treat distance-to-boundary as a continuous axis (analogous to pseudotime)

### 🧬 Tissue Interface Characterization
- Quantify how sharp or diffuse a boundary is across samples or conditions
- Compare boundary morphology between patient groups, treatment conditions, or tissue regions
- Extract shape metrics (e.g., Boundary-to-Nest Ratio) as quantitative features for downstream modeling

### 🗺️ Spatial Stratification of Cells
- Move beyond bulk or cell-type-level comparisons — stratify cells by their **spatial niche** (Boundary / Nest / Outside)
- Enables spatially-resolved differential expression, deconvolution, or cell-cell interaction analyses
- Compatible with any downstream R-based spatial analysis framework

### 🧫 Broad Tissue Biology Applications
Synora is not limited to tumor biology. Any tissue with a meaningful spatial interface benefits:
- **Cancer**: tumor–stroma interface, invasion front characterization
- **Development**: organ boundary formation, tissue compartmentalization
- **Immunology**: inflammatory lesion margins, lymphoid tissue organization
- **Organ biology**: zonation patterns (e.g., liver zones, cortex-medulla boundaries)


## Key Features
- **Precise Boundary Detection**: Introduces "Orientedness" metric to distinguish structured boundaries from random infiltration
- **Minimal Input Requirements**: Only needs cell coordinates and binary annotations
- **Robust Performance**: Maintains accuracy with 50% missing cells and 25% infiltration
- **Platform Agnostic**: Works with any spatial omics platform (Xenium, CosMx, Visium HD, CODEX, MIBI-TOF, etc.)
- **Comprehensive Analysis**: Three modular functions for boundary detection, distance calculation, and shape metrics


## How Synora Works

![Rationale of Synora](https://github.com/GabrielJuntengLi/Synora/blob/main/Synora_Figure.png)

### The Challenge
Traditional methods using cellular heterogeneity (Mixedness) alone cannot distinguish between:
- **True boundaries**: spatially segregated cell types (e.g., tumor–stroma interface)
- **Infiltration**: randomly mixed cell types (e.g., immune cells within tumor nests)
  
### Synora's Solution
- **Mixedness**: quantifies local cellular diversity (0 = homogeneous, 1 = maximum diversity)
- **Orientedness**: quantifies directional spatial bias using vector calculus
- **BoundaryScore**: Mixedness × Orientedness identifies true boundaries

---

## Dependencies

Requires R ≥ 4.1.0. Run the following to check and install any missing dependencies:

```r
pkgs <- c(
  "tidyverse",       # dplyr, tidyr, purrr, tibble, forcats, magrittr
  "dbscan",
  "concaveman",
  "sf",
  "spatstat.geom",
  "spatstat.explore",
  "deldir",
  "tidygraph"
)

missing <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing)
} else {
  message("All dependencies are already installed.")
}
```

---

## Quick Start

```r
library(Synora)
library(dplyr)
library(ggplot2)
library(patchwork)

data("DummyData")
# DummyData — named list of 5 data frames (one per Perlin noise frequency).
# Columns: Cell_ID, X, Y, CT (0/1 annotation), Distance2TrueContour
```

### Step 1 — Detect boundaries

```r
BoundaryResultList <- GetBoundary(
  INPUT                = DummyData,
  CELL_ID_COLUMN       = "Cell_ID",
  X_POSITION           = "X",
  Y_POSITION           = "Y",
  ANNO_COLUMN          = "CT",
  RADIUS               = 20,
  NEST_SPECIFICITY     = 0.25,
  BOUNDARY_SPECIFICITY = 0.05
)
# Each cell is labelled: Nest | Boundary_Inner | Boundary_Outer | Outside | Noise
```

### Step 2 — Measure distance to boundary

```r
DistanceResultList <- GetDist2Boundary(
  INPUT            = BoundaryResultList,
  CELL_ID_COLUMN   = "Cell_ID",
  X_POSITION       = "X",
  Y_POSITION       = "Y",
  ANNO_COLUMN      = "SynoraAnnotation",
  ANNO_OF_BOUNDARY = "Boundary"
)
# Adds Distance2Boundary (signed): negative = Outside, 0 = Boundary, positive = Nest
```

### Step 3 — Extract shape metrics

#### 3a. Sample-level global metrics

```r
GlobalMetricsList <- GetShapeMetrics(
  INPUT          = BoundaryResultList,
  CELL_ID_COLUMN = "Cell_ID",
  X_POSITION     = "X",
  Y_POSITION     = "Y",
  ANNO_COLUMN    = "SynoraAnnotation",
  SHAPE_METRICS  = c("Boundary2NestRatio", "NestFragmentation", "NestDispersion")
)

# Tidy into one data frame
GlobalMetricsDF <- purrr::imap_dfr(GlobalMetricsList, ~ tibble::tibble(
  Sample             = .y,
  Boundary2NestRatio = .x$Global$Boundary2NestRatio,
  NestFragmentation  = .x$Global$NestFragmentation,
  NestDispersion     = .x$Global$NestDispersion
))
```

#### 3b. Per-nest metrics with closed-nest detection

```r
ShapeResultList <- GetShapeMetrics(
  INPUT          = BoundaryResultList,
  CELL_ID_COLUMN = "Cell_ID",
  X_POSITION     = "X",
  Y_POSITION     = "Y",
  ANNO_COLUMN    = "SynoraAnnotation",
  SHAPE_METRICS  = c(
    "Boundary2NestRatio",                                    # global
    "NestSolidity", "NestCompactness", "NestElongation",     # all nests
    "BoundaryIntegrity", "NestBoundaryBalance",              # closed nests only
    "NestEnclosure",     "NestBoundaryDensity"
  ),
  SEPARATE_NESTS   = TRUE,
  DIST_THRESHOLD   = 60,
  MIN_NEST_SIZE    = 10,
  CONCAVITY        = 1,
  MARGIN_DIST      = 50,    # units within this distance of tissue edge = "touching margin"
  MARGIN_TOLERANCE = 0.05,  # ≤ 5 % of hull vertices touching → nest classified as closed
  BOUNDARY_RADIUS  = 50     # search radius for boundary-cell association
)

# Each element contains:
#   $Global   — one-row data frame of sample-level metrics
#   $PerNest  — one row per nest: NestSize, IsClosedNest, ClosedFraction, all metrics
#   $HullMap  — hull polygon vertices with per-nest metrics joined (ready for ggplot2)
#   $CellMap  — cell-to-nest ID mapping
```

### Step 4 — Visualize

#### 4a. Boundary annotation + distance maps

```r
PlotList <- purrr::imap(DummyData, function(raw_df, nm) {

  p1 <- ggplot2::ggplot(raw_df, ggplot2::aes(X, Y, color = as.factor(CT))) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_color_manual(
      name   = "Cell Type",
      values = c(`0` = "#B9C7E2", `1` = "#DE7424"),
      labels = c("Non-tumor", "Tumor")
    ) +
    ggplot2::theme_void() + ggplot2::coord_equal() +
    ggplot2::labs(title = nm)

  p2 <- ggplot2::ggplot(
    DistanceResultList[[nm]],
    ggplot2::aes(X, Y, color = factor(SynoraAnnotation))
  ) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_color_manual(
      values = c(Boundary_Inner = "#5B6530", Boundary_Outer = "#F1C100",
                 Nest = "#ECAB99", Outside = "#BFB5D0"),
      name   = "Synora Annotation"
    ) +
    ggplot2::theme_void() + ggplot2::coord_equal()

  p3 <- DistanceResultList[[nm]] %>%
    dplyr::mutate(
      D = scales::rescale_mid(Distance2Boundary, to = c(-1, 1), mid = 0)
    ) %>%
    ggplot2::ggplot(ggplot2::aes(X, Y, color = D)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_color_gradient2(
      low = "#046C9A", mid = "#FFFFFF", high = "#CB2314",
      midpoint = 0, limits = c(-1, 1), name = "Distance\n(scaled)"
    ) +
    ggplot2::theme_void() + ggplot2::coord_equal()

  patchwork::wrap_plots(p1, p2, p3, ncol = 1)
})

patchwork::wrap_plots(PlotList, nrow = 1, guides = "collect")
```

#### 4b. Shape metrics — `PlotShapeMetrics`

```r
# Single sample: all three panels (spatial + per-nest scatter + summary bar)
PlotShapeMetrics(
  SPATIAL_INPUT  = BoundaryResultList[[1]],
  METRICS_RESULT = ShapeResultList[[1]],
  SAMPLE_NAME    = names(ShapeResultList)[[1]]
)

# All samples at once — returns a named list of patchwork plots
ShapePlotList <- PlotShapeMetrics(
  SPATIAL_INPUT    = BoundaryResultList,
  METRICS_RESULT   = ShapeResultList,
  HULL_TYPE        = "concave",
  HIGHLIGHT_CLOSED = TRUE
)
patchwork::wrap_plots(ShapePlotList, nrow = 1)

# Spatial panel only
PlotShapeMetrics(
  SPATIAL_INPUT  = BoundaryResultList[[3]],
  METRICS_RESULT = ShapeResultList[[3]],
  PLOT_TYPE      = "spatial",
  SAMPLE_NAME    = "Perlin_0.005"
)

# Cross-sample summary comparison
purrr::imap(ShapeResultList, ~ PlotShapeMetrics(
  SPATIAL_INPUT  = BoundaryResultList[[.y]],
  METRICS_RESULT = .x,
  PLOT_TYPE      = "summary",
  SAMPLE_NAME    = .y
)) %>%
  patchwork::wrap_plots(nrow = 1)
```

#### 4c. Hull polygon visualisation — `PlotNestHulls`

```r
# Closure colouring: closed nests in red, open nests in grey
PlotNestHulls(
  METRICS_RESULT    = ShapeResultList[[1]],
  FILL_MODE         = "closure",
  BACKGROUND_POINTS = BoundaryResultList[[1]],
  SAMPLE_NAME       = names(ShapeResultList)[[1]]
)

# Heatmap fill by a per-nest metric
PlotNestHulls(
  METRICS_RESULT = ShapeResultList[[1]],
  FILL_MODE      = "metric",
  FILL_METRIC    = "NestCompactness",
  FILL_PALETTE   = "viridis"
)

# All samples — returns a named list of ggplots
HullPlotList <- PlotNestHulls(
  METRICS_RESULT    = ShapeResultList,
  FILL_MODE         = "closure",
  BACKGROUND_POINTS = BoundaryResultList
)
patchwork::wrap_plots(HullPlotList, nrow = 1)
```

## Function Overview

| Function | Input | Purpose |
|---|---|---|
| `GetBoundary` | Coordinates + binary annotation | Detect and label boundary cells |
| `GetDist2Boundary` | `GetBoundary` output | Signed distance of every cell to the nearest boundary |
| `GetShapeMetrics` | `GetBoundary` output | Global and per-nest morphology metrics |
| `PlotShapeMetrics` | `GetBoundary` + `GetShapeMetrics` output | Multi-panel diagnostic plots |
| `PlotNestHulls` | `GetShapeMetrics` output | Hull polygon visualisation with flexible fill modes |

All five functions accept **three equivalent input modes**:

| Mode | Input | Returns |
|---|---|---|
| List | Named list of data frames | Named list of results |
| Data frame + `SAMPLE_COLUMN` | Single data frame with a sample ID column | Single reassembled data frame / result |
| Single data frame | One data frame, no sample column | Single result |

---

## Output Reference

| Function | Key output columns / slots |
|---|---|
| `GetBoundary` | `Nb_Count` · `Anno_Midpoint` · `Mixedness` · `Orientedness` · `BoundaryScore` · `SynoraAnnotation` |
| `GetDist2Boundary` | All `GetBoundary` columns + `Distance2Boundary` (signed) · `Boundary_kNN_IDs` |
| `GetShapeMetrics` | `$Global` · `$PerNest` · `$HullMap` · `$CellMap` |
| `PlotShapeMetrics` | `patchwork` object or named list thereof |
| `PlotNestHulls` | `ggplot` object or named list thereof |

---

## Metric Reference

### Global metrics (`SEPARATE_NESTS = FALSE`)

| Metric | Description | Range |
|---|---|---|
| `Boundary2NestRatio` | Boundary cells ÷ nest cells | ≥ 0 |
| `NestFragmentation` | Disconnected components ÷ total nest cells | 0–1 |
| `NestDispersion` | Clark–Evans index (mean NND ÷ expected NND) | > 0 |
| `NestFractalDimension` | Box-counting fractal dimension of all nest cells | 1–2 |

### Per-nest metrics — all nests (`SEPARATE_NESTS = TRUE`)

| Metric | Description | Range |
|---|---|---|
| `NestSolidity` | area(concave hull) ÷ area(convex hull) | 0–1 |
| `BoundaryRoughness` | perimeter(concave) ÷ perimeter(convex) | ≥ 1 |
| `NestCompactness` | Polsby–Popper: 4π · area ÷ perimeter² | 0–1 |
| `NestElongation` | PCA minor ÷ major axis ratio | 0–1 |
| `NestFractalDimension` | Box-counting fractal dimension per nest | 1–2 |

### Per-nest metrics — closed nests only

| Metric | Description | Range |
|---|---|---|
| `BoundaryIntegrity` | Fraction of nest cells with ≥ 1 boundary cell within `BOUNDARY_RADIUS` | 0–1 |
| `NestBoundaryBalance` | Normalised angular entropy of boundary cells around nest centroid | 0–1 |
| `NestEnclosure` | Mean distance of nest cells to tissue margin | ≥ 0 |
| `NestBoundaryDensity` | Boundary cells near nest ÷ concave hull perimeter | ≥ 0 |

---

## Citation

If you use Synora in your research, please cite:

> Li J-T, et al.
> **Synora: vector-based boundary detection for spatial omics.**
> *bioRxiv* 2026. doi: [10.64898/2026.02.26.708395](https://doi.org/10.64898/2026.02.26.708395)

---
