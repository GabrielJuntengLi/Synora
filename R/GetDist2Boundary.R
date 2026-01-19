#' @title GetDist2Boundary
#' @description Calculate Distance to Boundary with Directionality
#' @param INPUT Input data frame containing cell coordinates and boundary annotations, where each row represents a cell in one single image.
#' @param X_POSITION Name of X coordinate column.
#' @param Y_POSITION Name of Y coordinate column.
#' @param ANNO_COLUMN Name of annotation column with boundary
#' @param CELL_ID_COLUMN (optional) Cell ID column. If not provided, row indices will be used as cell IDs. Default NULL.
#' @param CELL_ID_PREFIX (optional) Cell ID prefix used for creating cell IDs. Default NULL.
#' @param ANNO_COLUMN Name of annotation column containing boundary annotation. Default is "SynoraAnnotation".
#' @param ANNO_OF_BOUNDARY Annotation value representing boundaries. Default is "Boundary".
#' @param K Number of nearest boundary neighbors to consider. Default is 5.
#' @param DIRECTION Named list specifying direction mapping with:
#'   \itemize{
#'     \item \code{negative}: Annotation value for negative direction (default: 'Outside')
#'     \item \code{positive}: Annotation value for positive direction (default: 'Nest')
#'   }
#'   Set to \code{NULL} to ignore directionality
#' @importFrom dplyr transmute mutate nest_by arrange rename select
#' @importFrom tidyr unnest
#' @return A data frame with the following columns:
#' \describe{
#'   \item{CELL_ID}{Cell ID column.}
#'   \item{X_POSITION}{X coordinate column.}
#'   \item{Y_POSITION}{Y coordinate column.}
#'   \item{ANNO_COLUMN}{Cell annotation column containing boundary cells.}
#'   \item{Distance2Boundary}{Euclidean distance to boundary cell(s) (could be negative if DIRECTION applied).}
#' }
#' @export
#' @importFrom magrittr `%>%`
#' @examples
#' library(tidyverse)
#' library(patchwork)
#' library(Synora)
#'
#' # Generate Dummy Data
#' set.seed(123)
#' DummyData <- tidyr::expand_grid(X = seq(-1, 1, length.out = 65),
#'                                 Y = seq(-1, 1, length.out = 65)) %>%
#'   dplyr::mutate(R = 0.8,
#'                 Sinusoid = 0.4,
#'                 theta = 2 * atan(Y / (X + sqrt(X ^ 2 + Y ^ 2))),
#'                 theta = ifelse(is.na(theta), pi, theta),
#'                 CT = (X / R - Sinusoid * cos(theta) * sin(12 * theta)) ^ 2 +
#'                   (Y / R - Sinusoid * sin(theta) * sin(12 * theta)) ^ 2 < 1) %>%
#'   dplyr::mutate(CT = ifelse(CT & (as.logical(runif(n = nrow(.)) %/% 0.75)), !CT, CT) %>%
#'                   as.numeric()) %>%
#'   dplyr::mutate(X = 250 * X + rnorm(n = nrow(.)),
#'                 Y = 250 * Y + rnorm(n = nrow(.))) %>%
#'   dplyr::mutate(X = X %>% pmax(-250) %>% pmin(250),
#'                 Y = Y %>% pmax(-250) %>% pmin(250)) %>%
#'   dplyr::mutate(Cell_ID = paste0('Cell_', sprintf('%05.f', dplyr::row_number())))
#' set.seed(NULL)
#'
#' # Run GetBoundary
#' BoundaryResult <- Synora::GetBoundary(
#'   INPUT = DummyData,
#'   X_POSITION = 'X',
#'   Y_POSITION = 'Y',
#'   ANNO_COLUMN = 'CT',
#'   CELL_ID_COLUMN = 'Cell_ID',
#'   RADIUS = 20,
#'   NEST_SPECIFICITY = 0.25,
#'   BOUNDARY_SPECIFICITY = 0.05
#' )
#'
#' # Run GetDist2Boundary
#' LayerResult <- Synora::GetDist2Boundary(
#'   INPUT = BoundaryResult,
#'   CELL_ID_COLUMN = 'Cell_ID',
#'   X_POSITION = 'X',
#'   Y_POSITION = 'Y',
#'   ANNO_COLUMN = 'SynoraAnnotation',
#'   ANNO_OF_BOUNDARY = 'Boundary'
#' )
#'
#' # Visualization
#' p <- patchwork::wrap_plots(
#'   nrow = 2,
#'   DummyData %>%
#'     ggplot(aes(X, Y, color = as.factor(CT))) +
#'     geom_point(size = 1) +
#'     scale_color_manual(name = 'Cell Type',
#'                        values = c(`0` = '#e9c46a', `1` = '#046C9A'),
#'                        labels = c('Non-tumor cell', 'Tumor cell')) +
#'     labs(title = 'Cell Type') +
#'     theme_void() +
#'     coord_equal(),
#'   BoundaryResult %>%
#'     dplyr::mutate(BoundaryScore = pmin(BoundaryScore, 0.5)) %>%
#'     ggplot(aes(X, Y, color = BoundaryScore)) +
#'     geom_point(size = 1) +
#'     scale_color_gradient2(low = '#046C9A',
#'                           mid = '#FFFFFF',
#'                           high = '#CB2314',
#'                           midpoint = 0.25,
#'                           limits = c(0, 0.5),
#'                           name = "Boundary Score") +
#'     labs(title = "Boundary Score") +
#'     theme_void() +
#'     coord_equal(),
#'   BoundaryResult %>%
#'     ggplot(aes(X, Y, color = factor(SynoraAnnotation))) +
#'     geom_point(size = 1) +
#'     scale_color_manual(values = c(Boundary = '#4DAF4AFF', Nest = '#377EB8FF', Outside = '#984EA3FF'),
#'                        name = "Synora Annotation") +
#'     labs(title = "Synora Annotation") +
#'     theme_void() +
#'     coord_equal(),
#'   LayerResult %>%
#'     dplyr::mutate(Distance2Boundary = Distance2Boundary %>% scales::rescale_mid(c(-1, 1), mid = 0)) %>%
#'     ggplot(aes(X, Y, color = Distance2Boundary)) +
#'     geom_point(size = 1) +
#'     scale_color_gradient2(low = '#046C9A',
#'                           mid = '#FFFFFF',
#'                           high = '#CB2314',
#'                           midpoint = 0,
#'                           limits = c(-1, 1),
#'                           name = "Distance to Boundary") +
#'     labs(title = "Distance to Boundary") +
#'     theme_void() +
#'     coord_equal()) +
#'   patchwork::plot_layout(guides = "collect")
#' print(p)

GetDist2Boundary <- function(INPUT, X_POSITION, Y_POSITION,
                             CELL_ID_COLUMN, CELL_ID_PREFIX,
                             ANNO_COLUMN = 'SynoraAnnotation',
                             ANNO_OF_BOUNDARY = 'Boundary',
                             DIRECTION = list(negative = 'Outside', positive = 'Nest'),
                             K = 5) {
  if (missing(INPUT)) stop("INPUT data frame must be provided")
  if (!is.data.frame(INPUT)) stop("INPUT must be a data frame")
  if (missing(X_POSITION) || missing(Y_POSITION)) stop("Both X_POSITION and Y_POSITION column names must be provided")
  if (!X_POSITION %in% names(INPUT)) stop("X_POSITION: `\033[1;4;41m", X_POSITION, "\033[0m` not found in INPUT")
  if (!Y_POSITION %in% names(INPUT)) stop("Y_POSITION: `\033[1;4;41m", Y_POSITION, "\033[0m` not found in INPUT")
  if (missing(CELL_ID_COLUMN)) {
    message('Creating Cell_ID...')
    if (missing(CELL_ID_PREFIX)) {
      INPUT <- INPUT %>% dplyr::mutate(Cell_ID = dplyr::row_number())
    } else {
      INPUT <- INPUT %>% dplyr::mutate(Cell_ID = paste0(CELL_ID_PREFIX, '_', dplyr::row_number()))
    }
    CELL_ID_COLUMN <- 'Cell_ID'
  } else if (!CELL_ID_COLUMN %in% names(INPUT)) {
    stop("Specified CELL_ID_COLUMN not found in INPUT")
  }
  if (missing(ANNO_COLUMN)) stop("ANNO_COLUMN must be provided")
  if (!ANNO_COLUMN %in% names(INPUT)) stop("ANNO_COLUMN: `\033[1;4;41m", ANNO_COLUMN, "\033[0m` not found in INPUT")

  if (!missing(DIRECTION) && !is.null(DIRECTION) && !all(is.na(DIRECTION))) {
    if (!is.list(DIRECTION)) {
      stop("DIRECTION must be a named list")
    }
    if (!all(c("negative", "positive") %in% names(DIRECTION))) {
      stop("DIRECTION must contain both 'negative' and 'positive' elements")
    }
    if (identical(DIRECTION$negative, DIRECTION$positive)) {
      stop("DIRECTION$negative and DIRECTION$positive cannot be the same value")
    }
  }
  RESULT <- INPUT %>%
    dplyr::transmute(!!as.name(CELL_ID_COLUMN),
                     !!as.name(X_POSITION),
                     !!as.name(Y_POSITION),
                     !!as.name(ANNO_COLUMN),
                     Distance2Boundary = NA)

  if (sum(INPUT[[ANNO_COLUMN]] == ANNO_OF_BOUNDARY, na.rm = T) > 1) {
    RESULT1 <- INPUT %>%
      dplyr::transmute(Cell_ID = !!as.name(CELL_ID_COLUMN),
                       X = !!as.name(X_POSITION),
                       Y = !!as.name(Y_POSITION),
                       Anno = !!as.name(ANNO_COLUMN)) %>%
      dplyr::mutate(Query = Anno != ANNO_OF_BOUNDARY) %>%
      dplyr::nest_by(Query) %>%
      dplyr::mutate(INPUT_kNN = data %>% dplyr::transmute(X, Y) %>% list())
    
    
    kNN_RESULT <- dbscan::kNN(
      x = RESULT1$INPUT_kNN[[1]],
      query = RESULT1$INPUT_kNN[[2]],
      k = min(K, nrow(RESULT1$INPUT_kNN[[1]]) - 1))
    
    kNN_ID <- kNN_RESULT$id
    kNN_ID[] <- RESULT1$data[[1]]$Cell_ID[c(kNN_ID)]

    RESULT <- RESULT1 %>% 
      dplyr::mutate(Distance2Boundary = kNN_RESULT$dist %>% rowMeans() %>% list(),
                    kNN_ID = kNN_ID %>% {unname(as.list(as.data.frame(t(.))))} %>% list(),
                    kNN_dist = kNN_RESULT$dist %>% {unname(as.list(as.data.frame(t(.))))} %>% list()) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(data,
                       Distance2Boundary = ifelse(Query, Distance2Boundary, 0),
                       ) %>%
      tidyr::unnest(cols = c(data, Distance2Boundary)) %>%
      dplyr::arrange(Cell_ID) %>%
      dplyr::rename(!!as.name(CELL_ID_COLUMN) := Cell_ID,
                    !!as.name(X_POSITION) := X,
                    !!as.name(Y_POSITION) := Y,
                    !!as.name(ANNO_COLUMN) := Anno)
  }
  if (!is.null(DIRECTION) && !all(is.na(DIRECTION))) {
    RESULT <- RESULT %>%
      dplyr::mutate(
        Distance2Boundary = dplyr::case_when(
          .data[[ANNO_COLUMN]] == DIRECTION$negative ~ -Distance2Boundary,
          .data[[ANNO_COLUMN]] == DIRECTION$positive ~ Distance2Boundary,
          .data[[ANNO_COLUMN]] == ANNO_OF_BOUNDARY ~ Distance2Boundary,
          T ~ Distance2Boundary
        )
      )
  }
  return(RESULT)
}

