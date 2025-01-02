#' @title GetMO
#'
#' @description
#' This function takes a name and returns a greeting message.
#'
#' @param INPUT A data frame that must contain the image ID, cell ID, the coordinates of the cell in the image, and the cell type (tumor or non-tumor).
#' @param CELL_ID_COLUMN The cell ID.
#' @param X_POSITION The x-axis coordinate of cell in image.
#' @param Y_POSITION The y-axis coordinate of cell in image.
#' @param ANNO_COLUMN The cell type, binary classification (tumor -> 1 or non-tumor -> 0).
#' @param RADIUS Default 'auto'.
#' @param MEDIAN_NB Default 15. Recommended 15 to 30.
#' @param ORIENTEDNESS TRUE or FALSE. Default T.
#' @param VERBOSE Default 1.
#' @importFrom tidyverse
#' @importFrom data.table
#' @return A dataframe that contain tumor boundary and the microenvironment inside and outside the tumor divided by the boundary.
#' @examples
#' GetMO(CELL_ID_COLUMN = 'Cell_ID', X_POSITION = 'X_position', Y_POSITION = 'Y_position', ANNO_COLUMN = 'CT0', RADIUS = 'auto', MEDIAN_NB = 15, ORIENTEDNESS = F, VERBOSE = 1)

GetMO <- function(INPUT, CELL_ID_COLUMN, X_POSITION, Y_POSITION,
                  ANNO_COLUMN, RADIUS = 'auto', MEDIAN_NB = 15,
                  ORIENTEDNESS = T, VERBOSE = 1) {
  
  if(!is.numeric(INPUT[[ANNO_COLUMN]])) stop("ANNO_COLUMN must be a numeric")
  if(!(RADIUS == 'auto' | (is.numeric(RADIUS) & RADIUS > 0))) stop("RADIUS must be positive number")
  
  INPUT <- INPUT %>%
    dplyr::transmute(Cell_ID = !!as.name(CELL_ID_COLUMN),
                     X = !!as.name(X_POSITION),
                     Y = !!as.name(Y_POSITION),
                     Anno = !!as.name(ANNO_COLUMN))
  
  if (RADIUS == 'auto') {
    if (VERBOSE > 0) {
      message('Finding Radius Threshold...')
    }
    
    for (THR in c(seq(10, 100, by = 10),
                  seq(120, 200, by = 20),
                  seq(250, 500, by = 50))) {
      MEDIAN_NB_TEST <- INPUT %>%
        dplyr::select(X, Y) %>%
        dbscan::frNN(sort = F, eps = THR) %>%
        .$id %>%
        purrr::map_dbl(length) %>%
        summary() %>%
        .[['Median']]
      
      if (MEDIAN_NB_TEST >= MEDIAN_NB) {
        RADIUS <- THR
        break
      }
    }
    
    if (RADIUS == 'auto') {
      RADIUS <- 500
      message('Cell density is too low.')
    }
    
    if (VERBOSE > 0) {
      message('Choosing Radius of ', RADIUS, ' as threshold.')
    }
  }
  
  FRNN <- INPUT %>%
    dplyr::select(X, Y) %>%
    dbscan::frNN(sort = F, eps = RADIUS)
  
  if (ORIENTEDNESS) {
    RESULT <- INPUT %>%
      dplyr::mutate(ID = FRNN$id,
                    Dist = FRNN$dist) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(X_u = list((.$X[ID] - X) / Dist),
                    Y_u = list((.$Y[ID] - Y) / Dist),
                    Nb_Anno = list(.$Anno[ID])) %>%
      dplyr::mutate(Nb_Count = length(Nb_Anno),
                    Mean_Anno = mean(Nb_Anno),
                    X_u_SumBg = sum(X_u),
                    Y_u_SumBg = sum(Y_u),
                    X_u_Sum = sum(Nb_Anno * X_u),
                    Y_u_Sum = sum(Nb_Anno * Y_u)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Mixedness = 1 -(Mean_Anno) ^ 2) %>%
      dplyr::mutate(Orientedness = ((sqrt(X_u_Sum ^ 2 + Y_u_Sum ^ 2) - sqrt(X_u_SumBg ^ 2 + Y_u_SumBg ^ 2)) / Nb_Count)) %>%
      dplyr::mutate(Orientedness = ifelse(Orientedness < 0, 0, Orientedness)) %>%
      dplyr::mutate(Radius = RADIUS) %>%
      dplyr::select(Cell_ID, X, Y, Anno, Radius, Nb_Count, Mean_Anno, Mixedness, Orientedness) %>%
      dplyr::rename(!!as.name(CELL_ID_COLUMN) := Cell_ID,
                    !!as.name(X_POSITION) := X,
                    !!as.name(Y_POSITION) := Y,
                    !!as.name(ANNO_COLUMN) := Anno)
    
  } else {
    RESULT <- INPUT %>%
      dplyr::mutate(ID = FRNN$id,
                    Dist = FRNN$dist) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Nb_Anno = list(.$Anno[ID])) %>%
      dplyr::mutate(Nb_Count = length(Nb_Anno),
                    Mean_Anno = mean(Nb_Anno)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Mixedness = 1 -(Mean_Anno) ^ 2) %>%
      dplyr::mutate(Radius = RADIUS) %>%
      dplyr::select(Cell_ID, X, Y, Anno, Radius, Nb_Count, Mean_Anno, Mixedness) %>%
      dplyr::rename(!!as.name(CELL_ID_COLUMN) := Cell_ID,
                    !!as.name(X_POSITION) := X,
                    !!as.name(Y_POSITION) := Y,
                    !!as.name(ANNO_COLUMN) := Anno)
  }
  return(RESULT)
}

