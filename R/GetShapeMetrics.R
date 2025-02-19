#' @title GetShapeMetrics
#'
#' @description This function takes a name and returns a greeting message.
#'
#' @param INPUT A data frame containing cell coordinates and boundary annotations, where each row represents a cell in one single image.
#' @param CELL_ID_COLUMN The name of the column containing unique cell IDs.
#' @param X_POSITION The name of the column containing the x-coordinates of the cells.
#' @param Y_POSITION The name of the column containing the y-coordinates of the cells.
#' @param ANNO_COLUMN The name of the column containing the boundary and nest annotations.
#' @param ANNO_OF_BOUNDARY The name of boundary annotation derived from GetBoundary. Default is "Boundary".
#' @param ANNO_OF_NEST The name of nest annotation derived from GetBoundary. Default is "Nest".
#' @param SHAPE_METRICS A vector of metrics to calculate. Default is c("Boundary2NestRatio", "Dist2BoundaryVariance", "Convexity").
#' @importFrom dplyr transmute mutate nest_by arrange rename select
#' @importFrom tidyr unnest
#' @importFrom concaveman concaveman
#' @importFrom sf st_polygon st_area
#' @return A dataframe that contain the microenvironment inside and outside the tumor divided by the boundary.
#' @export
#' @examples
#' GetShapeMetrics(CELL_ID_COLUMN = 'Cell_ID', X_POSITION = 'X_position', Y_POSITION = 'Y_position', ANNO_COLUMN = 'CT0_2', ANNO_OF_BOUNDARY = 'Boundary')

GetShapeMetrics <- function(INPUT, CELL_ID_COLUMN, X_POSITION, Y_POSITION, ANNO_COLUMN,
                            ANNO_OF_BOUNDARY = 'Boundary', ANNO_OF_NEST = 'Nest',
                            SHAPE_METRICS = c('Boundary2NestRatio', 'Nest2BoundaryRatio', 'NestSolidity', 'Dist2BoundaryVariance')) {

  INPUT <- INPUT %>%
    dplyr::filter(!!as.name(ANNO_COLUMN) %in% c(ANNO_OF_BOUNDARY, ANNO_OF_NEST)) %>%
    dplyr::transmute(Cell_ID = !!as.name(CELL_ID_COLUMN),
                     X = !!as.name(X_POSITION),
                     Y = !!as.name(Y_POSITION),
                     Anno = !!as.name(ANNO_COLUMN) %>%
                       forcats::fct_recode(Boundary = ANNO_OF_BOUNDARY, Nest = ANNO_OF_NEST) %>%
                       forcats::fct_drop()) %>%
    dplyr::arrange(Anno)
  BoundaryCount <- INPUT %>%
    dplyr::filter(Anno == 'Boundary') %>%
    nrow()
  NestCount <- INPUT %>%
    dplyr::filter(Anno == 'Nest') %>%
    nrow()
  RESULT <- SHAPE_METRICS %>%
    purrr::set_names() %>%
    purrr::map(function(X) {
      if (X == 'Boundary2NestRatio') {
        if (BoundaryCount > 0 | NestCount > 0) {
          sum(INPUT$Anno == 'Boundary', na.rm = T) / sum(INPUT$Anno == 'Nest', na.rm = T)
        } else {
          NA
        }
      } else if (X == 'Nest2BoundaryRatio') {
        if (BoundaryCount > 0 | NestCount > 0) {
          sum(INPUT$Anno == 'Nest', na.rm = T) / sum(INPUT$Anno == 'Boundary', na.rm = T)
        } else {
          NA
        }
      } else if (X == 'NestSolidity') {
        if (BoundaryCount > 0 | NestCount > 0) {
          c(1, 10000) %>%
            purrr::map(function(CONCAVITY) {
              INPUT %>%
                dplyr::filter(Anno == 'Nest') %>%
                dplyr::select(X, Y) %>%
                as.matrix() %>%
                concaveman::concaveman(concavity = CONCAVITY) %>%
                list() %>%
                sf::st_polygon() %>%
                sf::st_area()
            }) %>%
            {.[[1]] / .[[2]]}
        } else {
          NA
        }
      } else if (X == 'Dist2BoundaryVariance') {
        if (BoundaryCount > 0 & NestCount > 0) {
          stats::var(
            dbscan::kNN(x = INPUT %>%
                          dplyr::filter(Anno == 'Boundary') %>%
                          dplyr::select(X, Y) %>%
                          as.matrix(),
                        query = INPUT %>%
                          dplyr::filter(Anno == 'Nest') %>%
                          dplyr::select(X, Y) %>%
                          as.matrix(),
                        k = 1)$dist[,1]
          )
        } else {
          NA
        }
      }
    })
  return(RESULT)
}

