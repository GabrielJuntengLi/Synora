#' @title Get Boundary
#'
#' @description This function takes a name and returns a greeting message.
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
#' @importFrom dplyr mutate nest_by ungroup case_when left_join relocate
#' @importFrom dbscan dbscan
#' @importFrom tidyr unnest
#' @importFrom forcats fct_expand fct_relevel
#' @return A dataframe that contain the tumor boundary annotation.
#' @examples
#' GetBoundary(CELL_ID_COLUMN = 'Cell_ID', X_POSITION = 'X_position', Y_POSITION = 'Y_position', ANNO_COLUMN = 'CT0', ANNO_OF_INTEREST = 'Cancer', RADIUS = 'auto', NEST_SPECIFICITY = 0.4, BOUNDARY_SPECIFICITY = 0.01)

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

