# **Synora: Spatial Tissue Architecture Analysis for Omics Data**

Synora — named after the Greek term for boundary — is an R package that resolves tissue architecture in spatial omics data. Rather than simply classifying cell types, Synora quantifies biologically meaningful interfaces, enabling downstream analyses that depend on knowing *where* a cell is relative to a tissue boundary, not just *what type* it is.

By accurately identifying tissue interfaces using a novel vector-based approach, Synora unlocks a new axis of biological inquiry: **spatial context**.

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
library(purrr)
library(ggplot2)
library(patchwork)

data("DummyData")
# DummyData is a named list of 5 data frames, one per Perlin noise frequency.
# Each data frame contains: Cell_ID, X, Y, CT (0/1), Distance2TrueContour
```

### Step 1 — Detect boundaries

```r
BoundaryResultList <- DummyData %>%
  GetBoundary(
    CELL_ID_COLUMN       = "Cell_ID",
    X_POSITION           = "X",
    Y_POSITION           = "Y",
    ANNO_COLUMN          = "CT",
    RADIUS               = 20,
    NEST_SPECIFICITY     = 0.25,
    BOUNDARY_SPECIFICITY = 0.05
  )
```

`GetBoundary` annotates every cell as **Nest**, **Boundary**,
**Outside**, or **Noise**.

---

### Step 2 — Measure distance to boundary

```r
DistanceResultList <- BoundaryResultList %>%
  GetDist2Boundary(
    CELL_ID_COLUMN   = "Cell_ID",
    X_POSITION       = "X",
    Y_POSITION       = "Y",
    ANNO_COLUMN      = "SynoraAnnotation",
    ANNO_OF_BOUNDARY = "Boundary"
  )
```

Each cell receives a **signed** `Distance2Boundary`: negative for Outside
cells, positive for Nest cells, and 0 for Boundary cells.

---

### Step 3 — Extract shape metrics

#### 3a. Global metrics (no nest separation)

```r
GlobalMetricsList <- BoundaryResultList %>%
  GetShapeMetrics(
    CELL_ID_COLUMN = "Cell_ID",
    X_POSITION     = "X",
    Y_POSITION     = "Y",
    ANNO_COLUMN    = "SynoraAnnotation",
    SHAPE_METRICS  = c("Boundary2NestRatio", "NestFragmentation",
                       "NestDispersion")
  )

# Tidy into a data frame
GlobalMetricsDF <- GlobalMetricsList %>%
  purrr::imap(~ tibble::tibble(
    Sample             = .y,
    Boundary2NestRatio = .x$Global$Boundary2NestRatio,
    NestFragmentation  = .x$Global$NestFragmentation,
    NestDispersion     = .x$Global$NestDispersion
  )) %>%
  dplyr::bind_rows()
```

#### 3b. Per-nest metrics with closed-nest detection

```r
ShapeResultList <- BoundaryResultList %>%
  GetShapeMetrics(
    CELL_ID_COLUMN   = "Cell_ID",
    X_POSITION       = "X",
    Y_POSITION       = "Y",
    ANNO_COLUMN      = "SynoraAnnotation",
    SHAPE_METRICS    = c(
      # Global
      "Boundary2NestRatio",
      # All-nest
      "NestSolidity", "NestCompactness", "NestElongation",
      # Closed-nest only
      "BoundaryIntegrity", "NestBoundaryBalance",
      "NestEnclosure",     "NestBoundaryDensity"
    ),
    SEPARATE_NESTS   = TRUE,
    DIST_THRESHOLD   = 60,
    MIN_NEST_SIZE    = 10,
    MARGIN_DIST      = 50,    # cells within 50 units of tissue edge = "touching margin"
    MARGIN_TOLERANCE = 0.05,  # ≤5% touching → nest is closed
    BOUNDARY_RADIUS  = 50     # search radius for boundary-cell association
  )

# ShapeResultList[[i]] has three slots:
#   $Global   — scalar global metrics
#   $PerNest  — one entry per nest (NestSize, NestIDs, IsClosedNest,
#               ClosedFraction, and all requested metrics)
#   $Summary  — mean/median/SD across AllNests and ClosedNests
```

---

### Step 4 — Visualize

#### 4a. Boundary annotation + distance maps (manual)

```r
PlotList <- vector("list", length(DummyData))

for (i in seq_along(DummyData)) {

  p1 <- DummyData[[i]] %>%
    ggplot2::ggplot(ggplot2::aes(X, Y, color = as.factor(CT))) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_color_manual(
      name   = "Cell Type",
      values = c(`0` = "#e9c46a", `1` = "#046C9A"),
      labels = c("Non-tumor cell", "Tumor cell")
    ) +
    ggplot2::theme_void() +
    ggplot2::labs(title = names(DummyData)[[i]]) +
    ggplot2::coord_equal()

  p2 <- DistanceResultList[[i]] %>%
    ggplot2::ggplot(ggplot2::aes(X, Y, color = factor(SynoraAnnotation))) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_color_manual(
      values = c(Boundary = "#4DAF4AFF", Nest = "#377EB8FF",
                 Outside  = "#984EA3FF"),
      name   = "Synora Annotation"
    ) +
    ggplot2::theme_void() +
    ggplot2::coord_equal()

  p3 <- DistanceResultList[[i]] %>%
    dplyr::mutate(
      Distance2Boundary = scales::rescale_mid(Distance2Boundary, c(-1, 1), mid = 0)
    ) %>%
    ggplot2::ggplot(ggplot2::aes(X, Y, color = Distance2Boundary)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_color_gradient2(
      low      = "#046C9A",
      mid      = "#FFFFFF",
      high     = "#CB2314",
      midpoint = 0,
      limits   = c(-1, 1),
      name     = "Distance to Boundary"
    ) +
    ggplot2::theme_void() +
    ggplot2::coord_equal()

  PlotList[[i]] <- patchwork::wrap_plots(ncol = 1, list(p1, p2, p3))
}

FinalPlot <- patchwork::wrap_plots(PlotList, nrow = 1,
                                   guides = "collect",
                                   axis_titles = "collect")
print(FinalPlot)
```

#### 4b. Shape metrics visualization (via `PlotShapeMetrics`)

```r
# Single sample — all three panels (spatial + per-nest + summary)
PlotShapeMetrics(
  SPATIAL_INPUT  = BoundaryResultList[[1]],
  METRICS_RESULT = ShapeResultList[[1]],
  SAMPLE_NAME    = names(ShapeResultList)[[1]]
)

# All samples — returns a named list of patchwork plots
ShapePlotList <- PlotShapeMetrics(
  SPATIAL_INPUT    = BoundaryResultList,
  METRICS_RESULT   = ShapeResultList,
  HULL_TYPE        = "concave",
  HIGHLIGHT_CLOSED = TRUE
)

# Print one sample
print(ShapePlotList[["Perlin_0.005"]])

# Spatial panel only
PlotShapeMetrics(
  SPATIAL_INPUT  = BoundaryResultList[[3]],
  METRICS_RESULT = ShapeResultList[[3]],
  PLOT_TYPE      = "spatial",
  SAMPLE_NAME    = "Perlin_0.005"
)

# Summary panel only — good for cross-sample comparison when combined
SummaryPlots <- purrr::imap(ShapeResultList, ~ PlotShapeMetrics(
  SPATIAL_INPUT  = BoundaryResultList[[.y]],
  METRICS_RESULT = .x,
  PLOT_TYPE      = "summary",
  SAMPLE_NAME    = .y
))
patchwork::wrap_plots(SummaryPlots, nrow = 1)
```

---

## Output structure at a glance

| Function | Returns |
|---|---|
| `GetBoundary` | Original data frame + `Nb_Count`, `Anno_Midpoint`, `Mixedness`, `Orientedness`, `BoundaryScore`, `SynoraAnnotation` |
| `GetDist2Boundary` | All columns from `GetBoundary` + `Distance2Boundary` (signed), `Boundary_kNN_IDs` |
| `GetShapeMetrics` | `list(Global, PerNest, Summary)` — see §3b above |
| `PlotShapeMetrics` | `patchwork` plot or named list thereof |

---

## Metric reference

### Global metrics (`SEPARATE_NESTS = FALSE`)

| Metric | Description | Range |
|---|---|---|
| `Boundary2NestRatio` | Boundary cells / nest cells | ≥ 0 |
| `NestFragmentation` | Disconnected components / total nest cells | 0–1 |
| `NestDispersion` | Clark-Evans index (mean NND / expected NND) | > 0 |

### Per-nest metrics — all nests (`SEPARATE_NESTS = TRUE`)

| Metric | Description | Range |
|---|---|---|
| `NestSolidity` | area(concave hull) / area(convex hull) | 0–1 |
| `BoundaryRoughness` | perimeter(concave) / perimeter(convex) | ≥ 1 |
| `NestCompactness` | Polsby-Popper: 4π·area / perimeter² | 0–1 |
| `NestElongation` | PCA minor / major axis ratio | 0–1 |
| `NestFractalDimension` | Box-counting fractal dimension | 1–2 |

### Per-nest metrics — closed nests only

| Metric | Description | Range |
|---|---|---|
| `BoundaryIntegrity` | Fraction of nest cells covered by ≥1 boundary cell within `BOUNDARY_RADIUS` | 0–1 |
| `NestBoundaryBalance` | Normalized angular entropy of boundary cells around nest centroid | 0–1 |
| `NestEnclosure` | Mean distance of nest cells to tissue margin | ≥ 0 |
| `NestBoundaryDensity` | Boundary cells near nest / concave hull perimeter | ≥ 0 |

---






