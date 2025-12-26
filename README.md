**Synora: A Vector-based Boundary Detection Tool for Spatial Omics**.

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
![Rationale of Synora](https://github.com/GabrielJuntengLi/Synora/blob/main/Synora_Figure.png)
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


 
# Quick Start
```r
library(Synora)

# 1. Detect boundaries
boundary_result <- GetBoundary(
  INPUT = your_data,
  X_POSITION = 'X',
  Y_POSITION = 'Y',
  ANNO_COLUMN = 'CellType',
  CELL_ID_COLUMN = 'Cell_ID',
  RADIUS = 20,
  NEST_SPECIFICITY = 0.25,
  BOUNDARY_SPECIFICITY = 0.05
)

# 2. Calculate distances
distance_result <- GetDist2Boundary(
  INPUT = boundary_result,
  CELL_ID_COLUMN = 'Cell_ID',
  X_POSITION = 'X',
  Y_POSITION = 'Y',
  ANNO_COLUMN = 'SynoraAnnotation',
  ANNO_OF_BOUNDARY = 'Boundary'
)

# 3. Extract shape metrics
shape_metrics <- GetShapeMetrics(
  INPUT = boundary_result,
  CELL_ID_COLUMN = 'Cell_ID',
  X_POSITION = 'X',
  Y_POSITION = 'Y',
  ANNO_COLUMN = 'SynoraAnnotation',
  SHAPE_METRICS = c('Boundary2NestRatio')
)

```

## Dependencies
- R ≥ 4.0
- tidyverse
- dbscan
- concaveman
- sf
- spatstat

