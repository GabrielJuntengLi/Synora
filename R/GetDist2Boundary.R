#' @title GetDist2Boundary
#'
#' @description This function takes a name and returns a greeting message.
#'
#' @param INPUT A data frame containing cell coordinates and boundary annotations, where each row represents a cell in one single image.
#' @param CELL_ID_COLUMN The name of the column containing unique cell IDs.
#' @param X_POSITION The name of the column containing the x-coordinates of the cells.
#' @param Y_POSITION The name of the column containing the y-coordinates of the cells.
#' @param ANNO_COLUMN The name of the column containing the boundary annotations.
#' @param ANNO_OF_BOUNDARY The name of boundary annotation derived from GetBoundary. Default is "Boundary".
#' @param K The number of neighboring boundary cells, used for calculating MeanDistance. Default is 5.
#' @importFrom dplyr transmute mutate nest_by arrange rename select
#' @importFrom tidyr unnest
#' @return A dataframe that contain the microenvironment inside and outside the tumor divided by the boundary.
#' @export
#' @examples
#' GetDist2Boundary(CELL_ID_COLUMN = 'Cell_ID', X_POSITION = 'X_position', Y_POSITION = 'Y_position', ANNO_COLUMN = 'CT0_2', ANNO_OF_BOUNDARY = 'Boundary')

GetDist2Boundary <- function(INPUT, CELL_ID_COLUMN, X_POSITION, Y_POSITION,
                             ANNO_COLUMN, ANNO_OF_BOUNDARY = 'Boundary', K = 5) {
  if (sum(INPUT[[ANNO_COLUMN]] == ANNO_OF_BOUNDARY, na.rm = T) > 1) {
    RESULT <- INPUT %>%
      dplyr::transmute(Cell_ID = !!as.name(CELL_ID_COLUMN),
                       X = !!as.name(X_POSITION),
                       Y = !!as.name(Y_POSITION),
                       Anno = !!as.name(ANNO_COLUMN)) %>%
      dplyr::mutate(Query = Anno != ANNO_OF_BOUNDARY) %>%
      dplyr::nest_by(Query) %>%
      dplyr::mutate(INPUT_kNN = data %>% dplyr::transmute(X, Y) %>% list()) %T>%
      {kNN_RESULT <<- dbscan::kNN(
        x = .$INPUT_kNN[[1]],
        query = .$INPUT_kNN[[2]],
        k = min(K, nrow(.$INPUT_kNN[[1]]) - 1))} %T>%
      {kNN_ID <<- kNN_RESULT$id; kNN_ID[] <<- .$data[[1]]$Cell_ID[c(kNN_ID)]} %>%
      dplyr::mutate(MinDistance = kNN_RESULT$dist[,1] %>% list(),
                    MeanDistance = kNN_RESULT$dist %>% rowMeans() %>% list(),
                    kNN_ID = kNN_ID %>% {unname(as.list(as.data.frame(t(.))))} %>% list(),
                    kNN_dist = kNN_RESULT$dist %>% {unname(as.list(as.data.frame(t(.))))} %>% list()) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(data,
                       MinDistance = ifelse(Query, MinDistance, 0),
                       MeanDistance = ifelse(Query, MeanDistance, 0),
                       kNN_ID = ifelse(Query, kNN_ID, NA),
                       kNN_dist = ifelse(Query, kNN_dist, NA)) %>%
      tidyr::unnest(cols = c(data, MinDistance, MeanDistance, kNN_ID, kNN_dist)) %>%
      dplyr::arrange(Cell_ID) %>%
      dplyr::rename(!!as.name(CELL_ID_COLUMN) := Cell_ID,
                    !!as.name(X_POSITION) := X,
                    !!as.name(Y_POSITION) := Y,
                    !!as.name(ANNO_COLUMN) := Anno)
  } else {
    RESULT <- INPUT %>%
      dplyr::select(dplyr::all_of(c(CELL_ID_COLUMN, X_POSITION, Y_POSITION, ANNO_COLUMN))) %>%
      dplyr::mutate(MinDistance = NA,
                    MeanDistance = NA,
                    kNN_ID = NULL,
                    kNN_dist = NULL)
  }
  return(RESULT)
}

