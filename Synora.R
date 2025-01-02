#' Get Tumor Boundary
#'
#' This function takes a name and returns a greeting message.
#' 
#' @param INPUT A data frame that must contain the image ID, cell ID, the coordinates of the cell in the image, and the cell type (tumor or non-tumor).
#' @param CELL_ID_COLUMN The cell ID.
#' @param X_POSITION The x-axis coordinate of cell in image.
#' @param Y_POSITION The y-axis coordinate of cell in image.
#' @param ANNO_COLUMN The cell type, binary classification (tumor -> 1 or non-tumor -> 0).
#' @param ANNO_OF_INTEREST The cell type interested, binary classification ('Cancer' or non-tumor).
#' @param ANNO_MAX Default 1.
#' @param ANNO_MIN Default 0.
#' @param ANNO_MIDPOINT Default 0.5.
#' @param RADIUS Default 'auto'.
#' @param MEDIAN_NB Default 15. Recommended 15 to 30.
#' @param NEST_MIN_SIZE Default 5. Recommended 5 to 10.
#' @param NEST_SPECIFICITY Default 0.25. Recommended 0.2 to 0.4.
#' @param BOUNDARY_SPECIFICITY Default 0.25. Recommended 0.2 to 0.4.
#' @return A dataframe that contain tumor boundary and the microenvironment inside and outside the tumor divided by the boundary.
#' @examples
#' greet("Alice")

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


GetBoundary <- function(INPUT, CELL_ID_COLUMN, X_POSITION, Y_POSITION, 
                        ANNO_COLUMN, ANNO_OF_INTEREST, 
                        ANNO_MAX = 1, ANNO_MIN = 0, ANNO_MIDPOINT = 0.5,
                        RADIUS = 'auto',
                        MEDIAN_NB = 15,                 # 15-30
                        NEST_MIN_SIZE = 5,              # 5-10
                        NEST_SPECIFICITY = 0.25,        # 0.2-0.4
                        BOUNDARY_SPECIFICITY = 0.25     # 0.2-0.4
                        # BOUNDARY_SPECIFICITY = 0.02     # 0.01-0.03
) {
  INPUT <- INPUT %>% 
    dplyr::mutate(
      !!as.name(ANNO_COLUMN) := ifelse(!!as.name(ANNO_COLUMN) >= ANNO_MIDPOINT, 
                                       (!!as.name(ANNO_COLUMN) - ANNO_MIDPOINT) / (ANNO_MAX - ANNO_MIDPOINT), 
                                       (!!as.name(ANNO_COLUMN) - ANNO_MIDPOINT) / (ANNO_MIDPOINT - ANNO_MIN)))
  
  RESULT_1 <- INPUT %>% 
    GetMO(CELL_ID_COLUMN = CELL_ID_COLUMN, X_POSITION = X_POSITION, Y_POSITION = Y_POSITION, 
          ANNO_COLUMN = ANNO_COLUMN, RADIUS = RADIUS, MEDIAN_NB = MEDIAN_NB, 
          ORIENTEDNESS = F, VERBOSE = 1) %>% 
    dplyr::mutate(Nest = (0.5 * Mean_Anno + 0.5) >= NEST_SPECIFICITY)
  
  NOISE_COUNT <- c()
  i <- 1
  while (T) {
    RESULT_1 <- RESULT_1 %>% 
      dplyr::nest_by(Nest, Radius) %>% 
      dplyr::mutate(!!as.name(paste0('Noise_', i)) := data %>%
                      dplyr::select(all_of(c(X_POSITION, Y_POSITION))) %>% 
                      dbscan::dbscan(eps = Radius, minPts = NEST_MIN_SIZE) %>% 
                      .$cluster %>% 
                      `==`(0) %>% 
                      list()) %>% 
      tidyr::unnest(cols = c(data, !!as.name(paste0('Noise_', i))))  %>% 
      dplyr::ungroup() %>% 
      dplyr::bind_rows() %>% 
      dplyr::mutate(Nest = ifelse(!!as.name(paste0('Noise_', i)), !Nest, Nest))
    NOISE_COUNT <- c(NOISE_COUNT, sum(RESULT_1[[paste0('Noise_', i)]]))
    
    if (ifelse(length(NOISE_COUNT) == 1, NOISE_COUNT == 0, diff(tail(NOISE_COUNT, 2)) == 0)) {
      RESULT_1 <- RESULT_1 %>% 
        dplyr::mutate(!!as.name(paste0(ANNO_COLUMN, '_1')) := dplyr::case_when(
          !(!!as.name(paste0('Noise_', i))) & Nest ~ 1, 
          !(!!as.name(paste0('Noise_', i))) & !Nest ~ -1, 
          T ~ NA, 
        )# %>% 
        # forcats::fct_expand('Nest', 'Outside', 'Noise') %>% 
        # forcats::fct_relevel('Nest', 'Outside', 'Noise')
        )
      break
    } else {
      i <- i + 1
    }
  }
  RESULT <- RESULT_1 %>% 
    GetMO(CELL_ID_COLUMN = CELL_ID_COLUMN, X_POSITION = X_POSITION, Y_POSITION = Y_POSITION,
          ANNO_COLUMN = paste0(ANNO_COLUMN, '_1'), RADIUS = RADIUS, MEDIAN_NB = MEDIAN_NB,
          ORIENTEDNESS = T, VERBOSE = 0) %>%
    dplyr::mutate(Boundary_Score = Mixedness * Orientedness) %>%
    dplyr::mutate(!!as.name(paste0(ANNO_COLUMN, '_2')) := dplyr::case_when(
      Boundary_Score >= BOUNDARY_SPECIFICITY ~  'Boundary', 
      Boundary_Score < BOUNDARY_SPECIFICITY & !!as.name(paste0(ANNO_COLUMN, '_1')) == 1 ~  'Nest', 
      Boundary_Score < BOUNDARY_SPECIFICITY & !!as.name(paste0(ANNO_COLUMN, '_1')) == -1 ~  'Outside', 
      T ~  'Noise') %>%
        forcats::fct_expand('Boundary', 'Nest', 'Outside', 'Noise') %>%
        forcats::fct_relevel('Boundary', 'Nest', 'Outside', 'Noise')) %>% 
    dplyr::left_join(INPUT, ., by = c(CELL_ID_COLUMN, X_POSITION, Y_POSITION)) %>%
    dplyr::relocate(all_of(colnames(INPUT)), Radius, Nb_Count, paste0(ANNO_COLUMN, '_1'), paste0(ANNO_COLUMN, '_2'))
  
  return(RESULT)
}


GetGrid <- function(INPUT, CELL_ID_COLUMN, X_POSITION, Y_POSITION, 
                    ANNO_COLUMN, ANNO_OF_BOUNDARY, K = 5) {
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
      dplyr::select(dplyr::all_of(c(CELL_ID_COLUMN, X_POSITION, Y_POSITION, ANNO_COLUMN)))
  }
  return(RESULT)
}

