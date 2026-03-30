# **Synora: Spatial Tissue Architecture Analysis for Omics Data**

Synora is an R package that resolves tissue architecture in spatial omics data — enabling downstream analyses that depend on knowing *where* a cell is relative to a tissue boundary, not just *what type* it is.

By accurately identifying tissue interfaces using a novel vector-based approach, Synora unlocks a new axis of biological inquiry: **spatial context**.

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

---

## Key Features
- **Precise Boundary Detection**: Introduces "Orientedness" metric to distinguish structured boundaries from random infiltration
- **Minimal Input Requirements**: Only needs cell coordinates and binary annotations
- **Robust Performance**: Maintains accuracy with 50% missing cells and 25% infiltration
- **Platform Agnostic**: Works with any spatial omics platform (Visium HD, CODEX, MIBI-TOF, etc.)
- **Comprehensive Analysis**: Three modular functions for boundary detection, distance calculation, and shape metrics

---

## Installation

```r
# Stable version
devtools::install_github("lzxlab/Synora")

# Development version (experimental features)
devtools::install_github("GabrielJuntengLi/Synora")
```

---

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

## Quick Start

```r
library(Synora)
data("DummyData")

# 1. Detect boundaries
BoundaryResultList <- DummyData %>% 
    purrr::map(.progress = T, ~ Synora::GetBoundary(
        INPUT = .x,
        CELL_ID_COLUMN = 'Cell_ID',
        X_POSITION = 'X',
        Y_POSITION = 'Y',
        ANNO_COLUMN = 'CT',
        RADIUS = 20,
        NEST_SPECIFICITY = 0.25,
        BOUNDARY_SPECIFICITY = 0.05
    ))
    
# 2. Calculate distances
DistanceResultList <- BoundaryResultList %>% 
    purrr::map(.progress = T, ~ Synora::GetDist2Boundary(
        INPUT = .x,
        CELL_ID_COLUMN = 'Cell_ID',
        X_POSITION = 'X',
        Y_POSITION = 'Y',
        ANNO_COLUMN = 'SynoraAnnotation',
        ANNO_OF_BOUNDARY = 'Boundary'
    ))

# 3. Extract shape metrics
ShapeResultList <- BoundaryResultList %>% 
    purrr::map(.progress = T, ~ Synora::GetShapeMetrics(
        INPUT = .x,
        CELL_ID_COLUMN = 'Cell_ID',
        X_POSITION = 'X',
        Y_POSITION = 'Y',
        ANNO_COLUMN = 'SynoraAnnotation',
        SHAPE_METRICS = c('Boundary2NestRatio')
    )) %>% 
    tibble::enframe(name = 'PerlinFrequency', value = 'BNR') %>% 
    tidyr::unnest(BNR) %>% 
    tidyr::unnest(BNR)

# 4. Visualization
PlotList <- list()
for (i in 1:length(DummyData)) {
  p1 <- DummyData[[i]] %>%
    ggplot2::ggplot(ggplot2::aes(X, Y, color = as.factor(CT))) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_color_manual(name = 'Cell Type',
                                values = c(`0` = '#e9c46a', `1` = '#046C9A'),
                                labels = c('Non-tumor cell', 'Tumor cell')) +
    ggplot2::theme_void() +
    ggplot2::labs(title = names(DummyData)[[i]]) +
    # ggplot2::theme(axis.title = ggplot2::element_text()) +
    ggplot2::coord_equal()
  p2 <- DistanceResultList[[i]] %>%
    ggplot2::ggplot(ggplot2::aes(X, Y, color = factor(SynoraAnnotation))) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_color_manual(values = c(Boundary = '#4DAF4AFF', Nest = '#377EB8FF', Outside = '#984EA3FF'),
                                name = "Synora Annotation") +
    ggplot2::theme_void() +
    ggplot2::coord_equal()
  
  p3 <- DistanceResultList[[i]] %>%
    dplyr::mutate(Distance2Boundary = Distance2Boundary %>% scales::rescale_mid(c(-1, 1), mid = 0)) %>%
    ggplot2::ggplot(ggplot2::aes(X, Y, color = Distance2Boundary)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_color_gradient2(low = '#046C9A',
                                   mid = '#FFFFFF',
                                   high = '#CB2314',
                                   midpoint = 0,
                                   limits = c(-1, 1),
                                   name = "Distance to Boundary") +
    ggplot2::theme_void() +
    ggplot2::coord_equal()
  PlotList[[i]] <- patchwork::wrap_plots(ncol = 1, list(p1, p2, p3))
}
FinalPlot <- PlotList %>% 
  patchwork::wrap_plots(nrow = 1, guides = 'collect', axis_titles = 'collect')
print(FinalPlot)

```

---


## Dependencies
- R ≥ 4.0
- tidyverse
- dbscan
- concaveman
- sf
- spatstat

