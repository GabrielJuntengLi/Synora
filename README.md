# **Synora: A Vector-based Boundary Detection Tool for Spatial Omics**

Synora is an R package for accurate identification and characterization of tissue boundaries in spatial omics data. It uses a novel vector-based approach to distinguish true tissue interfaces from regions of cellular infiltration, requiring only cell coordinates and binary annotations.

## Key Features
- **Precise Boundary Detection**: Introduces "Orientedness" metric to distinguish structured boundaries from random infiltration
- **Minimal Input Requirements**: Only needs cell coordinates and tumor/non-tumor annotations
- **Robust Performance**: Maintains accuracy with 50% missing cells and 25% infiltration
- **Platform Agnostic**: Works with any spatial omics platform (Visium HD, CODEX, MIBI-TOF, etc.)
- **Comprehensive Analysis**: Three modular functions for boundary detection, distance calculation, and shape metrics



## Installation
```r
# Install from GitHub
devtools::install_github("GabrielJuntengLi/Synora")
```

## How Synora Works & Applications
![Rationale of Synora](https://github.com/GabrielJuntengLi/Synora/blob/main/Synora_Figure.jpg)
### The Challenge
Traditional methods using cellular heterogeneity (Mixedness) alone cannot distinguish between:
- True boundaries: spatially segregated cell types (e.g., tumor-stroma interface)
- Infiltration: randomly mixed cell types (e.g., immune cells within tumor nests)
  
### Synora's Solution
- Mixedness: quantifies local cellular diversity (0 = homogeneous, 1 = maximum diversity)
- Orientedness: quantifies directional spatial bias using vector calculus
- BoundaryScore: Mixedness × Orientedness identifies true boundaries

### Applications
Synora is broadly applicable to any tissue interface, enabling:
- Cancer Biology:
  - Tumor-stroma interface characterization
  - Boundary-enriched gene/protein identification
  - Morphological feature extraction
- Beyond Cancer:
  - Developmental boundaries
  - Organ zonation patterns
  - Inflammatory lesion interfaces
  - Spatial gradient analysis across any tissue compartment


 
## Quick Start
```r
library(Synora)
data("DummyData")

# 1. Detect boundaries
BoundaryResultList <- DummyData %>% 
    purrr::map(.progress = T, ~ GetBoundary(
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
    purrr::map(.progress = T, ~ GetDist2Boundary(
        INPUT = .x,
        CELL_ID_COLUMN = 'Cell_ID',
        X_POSITION = 'X',
        Y_POSITION = 'Y',
        ANNO_COLUMN = 'SynoraAnnotation',
        ANNO_OF_BOUNDARY = 'Boundary'
    ))

# 3. Extract shape metrics
ShapeResultList <- BoundaryResultList %>% 
    purrr::map(.progress = T, ~ GetShapeMetrics(
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

## Dependencies
- R ≥ 4.0
- tidyverse
- dbscan
- concaveman
- sf
- spatstat

