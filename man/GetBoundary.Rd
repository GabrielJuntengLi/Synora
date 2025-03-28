% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetBoundary.R
\name{GetBoundary}
\alias{GetBoundary}
\title{Get Boundary}
\usage{
GetBoundary(
  INPUT,
  X_POSITION,
  Y_POSITION,
  ANNO_COLUMN,
  CELL_ID_COLUMN,
  CELL_ID_PREFIX,
  ANNO_RANGE = c(0, 1),
  ANNO_MIDPOINT = 0.5,
  RADIUS = "auto",
  MEDIAN_NB = 15,
  NEST_MIN_SIZE = 5,
  NEST_SPECIFICITY = 0.25,
  BOUNDARY_SPECIFICITY = 0.05
)
}
\arguments{
\item{INPUT}{A data frame containing cell coordinates and annotations, where each row represents a cell in one single image.}

\item{X_POSITION}{The name of the column containing the x-coordinates of the cells.}

\item{Y_POSITION}{The name of the column containing the y-coordinates of the cells.}

\item{ANNO_COLUMN}{The name of the column containing the cell annotations, where 1 indicates a tumor/nest cell and 0 indicates a non-tumor/non-nest cell. Can be binary or continuous.}

\item{CELL_ID_COLUMN}{(optional) The name of the column containing unique cell IDs. If not provided, row indices will be used as cell IDs. Default is NULL.}

\item{CELL_ID_PREFIX}{(optional) The prefix used for creating cell IDs. Default is NULL.}

\item{ANNO_RANGE}{(optional) The range of the cell annotations. Default is c(0, 1).}

\item{ANNO_MIDPOINT}{(optional) The midpoint of the cell annotation range. Default is 0.5.}

\item{RADIUS}{(optional) The radius used for determination of neighbors. If set to 'auto', the radius will be determined automatically. Default is 'auto'.}

\item{MEDIAN_NB}{(optional) Used when RADIUS is set to 'auto' to determine the minimal radius that suffices the profiling of cell neighborhood. Default 15.}

\item{NEST_MIN_SIZE}{(optional) The minimum size of cell nests, recommended between 5 and 10. Default is 5.}

\item{NEST_SPECIFICITY}{(optional) The specificity of nest cell detection, recommended between 0.1 and 0.4. Default is 0.25. This is the minimal frequency of neighboring tumor cells for the definition of a nest cell. Higher nest specificity results in less nest cells detected.}

\item{BOUNDARY_SPECIFICITY}{(optional) The specificity of boundary cell detection, recommended between 0.01 and 0.1. Default is 0.05. Higher boundary specificity results in less boundary cells detected.}
}
\value{
A data frame with the following columns:
\describe{
\item{CELL_ID}{The unique identifier for each cell.}
\item{X_POSITION}{The x-coordinate of the cell.}
\item{Y_POSITION}{The y-coordinate of the cell.}
\item{ANNOTATION}{The original cell annotation (1 for tumor, 0 for non-tumor).}
\item{BOUNDARY}{A logical value indicating whether the cell is part of the tumor boundary (TRUE) or not (FALSE).}
}
}
\description{
This function takes a data frame of cell coordinates and annotations and returns a data frame with the tumor boundary annotation.
}
\examples{
data <- data.frame(
  Cell_ID = 1:100,
  X_position = rnorm(100, 50, 10),
  Y_position = rnorm(100, 50, 10),
  CT0 = sample(c(0, 1), 100, replace = TRUE)
)

boundary_df <- GetBoundary(
  INPUT = data,
  X_POSITION = 'X_position',
  Y_POSITION = 'Y_position',
  ANNO_COLUMN = 'CT0',
  CELL_ID_COLUMN = 'Cell_ID',
  RADIUS = 'auto',
  NEST_SPECIFICITY = 0.25,
  BOUNDARY_SPECIFICITY = 0.05
)

# Visualize the tumor boundary
plot(boundary_df$X_POSITION, boundary_df$Y_POSITION, col = ifelse(boundary_df$BOUNDARY, 'red', 'black'), pch = 16)
}
