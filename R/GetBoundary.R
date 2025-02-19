#' @title Get Boundary
#'
#' @description This function takes a data frame of cell coordinates and annotations and returns a data frame with the tumor boundary annotation.
#'
#' @param INPUT A data frame containing cell coordinates and annotations, where each row represents a cell in one single image.
#' @param X_POSITION The name of the column containing the x-coordinates of the cells.
#' @param Y_POSITION The name of the column containing the y-coordinates of the cells.
#' @param ANNO_COLUMN The name of the column containing the cell annotations, where 1 indicates a tumor/nest cell and 0 indicates a non-tumor/non-nest cell. Can be binary or continuous.
#' @param CELL_ID_COLUMN (optional) The name of the column containing unique cell IDs. If not provided, row indices will be used as cell IDs. Default is NULL.
#' @param CELL_ID_PREFIX (optional) The prefix used for creating cell IDs. Default is NULL.
#' @param ANNO_RANGE (optional) The range of the cell annotations. Default is c(0, 1).
#' @param ANNO_MIDPOINT (optional) The midpoint of the cell annotation range. Default is 0.5.
#' @param RADIUS (optional) The radius used for determination of neighbors. If set to 'auto', the radius will be determined automatically. Default is 'auto'.
#' @param MEDIAN_NB (optional) Used when RADIUS is set to 'auto' to determine the minimal radius that suffices the profiling of cell neighborhood. Default 15.
#' @param NEST_SPECIFICITY (optional) The specificity of nest cell detection, recommended between 0.1 and 0.4. Default is 0.25. This is the minimal frequency of neighboring tumor cells for the definition of a nest cell. Higher nest specificity results in less nest cells detected.
#' @param BOUNDARY_SPECIFICITY (optional) The specificity of boundary cell detection, recommended between 0.01 and 0.1. Default is 0.05. Higher boundary specificity results in less boundary cells detected.
#' @param NEST_MIN_SIZE (optional) The minimum size of cell nests, recommended between 5 and 10. Default is 5.
#' @importFrom dplyr mutate nest_by ungroup case_when left_join relocate
#' @importFrom dbscan dbscan
#' @importFrom tidyr unnest
#' @importFrom forcats fct_expand fct_relevel
#' @importFrom magrittr %>%
#' @return A data frame with the following columns:
#' \describe{
#'   \item{CELL_ID}{The unique identifier for each cell.}
#'   \item{X_POSITION}{The x-coordinate of the cell.}
#'   \item{Y_POSITION}{The y-coordinate of the cell.}
#'   \item{ANNOTATION}{The original cell annotation (1 for tumor, 0 for non-tumor).}
#'   \item{BOUNDARY}{A logical value indicating whether the cell is part of the tumor boundary (TRUE) or not (FALSE).}
#' }
#' @export
#' @examples
#' data <- data.frame(
#'   Cell_ID = 1:100,
#'   X_position = rnorm(100, 50, 10),
#'   Y_position = rnorm(100, 50, 10),
#'   CT0 = sample(c(0, 1), 100, replace = TRUE)
#' )
#'
#' boundary_df <- GetBoundary(
#'   INPUT = data,
#'   X_POSITION = 'X_position',
#'   Y_POSITION = 'Y_position',
#'   ANNO_COLUMN = 'CT0',
#'   CELL_ID_COLUMN = 'Cell_ID',
#'   RADIUS = 'auto',
#'   NEST_SPECIFICITY = 0.25,
#'   BOUNDARY_SPECIFICITY = 0.05
#' )
#'
#' # Visualize the tumor boundary
#' plot(boundary_df$X_POSITION, boundary_df$Y_POSITION, col = ifelse(boundary_df$BOUNDARY, 'red', 'black'), pch = 16)

GetBoundary <- function(INPUT, X_POSITION, Y_POSITION,
                        ANNO_COLUMN, CELL_ID_COLUMN, CELL_ID_PREFIX,
                        ANNO_RANGE = c(0, 1), ANNO_MIDPOINT = 0.5,
                        RADIUS = 'auto',
                        MEDIAN_NB = 15,                 # 15-30
                        NEST_MIN_SIZE = 5,              # 5-10
                        NEST_SPECIFICITY = 0.25,        # 0.2-0.4
                        BOUNDARY_SPECIFICITY = 0.05     # 0.2-0.4
) {
  INPUT <- INPUT %>%
    dplyr::mutate(
      !!as.name(ANNO_COLUMN) := ifelse(!!as.name(ANNO_COLUMN) >= ANNO_MIDPOINT,
                                       (!!as.name(ANNO_COLUMN) - ANNO_MIDPOINT) / (ANNO_RANGE[2] - ANNO_MIDPOINT),
                                       (!!as.name(ANNO_COLUMN) - ANNO_MIDPOINT) / (ANNO_MIDPOINT - ANNO_RANGE[1])))
  if (missing(CELL_ID_COLUMN) & missing(CELL_ID_PREFIX)) {
    INPUT <- INPUT %>% dplyr::mutate(Cell_ID = dplyr::row_number())
    CELL_ID_COLUMN <- 'Cell_ID'
  } else if (missing(CELL_ID_COLUMN) & !missing(CELL_ID_PREFIX)) {
    INPUT <- INPUT %>% dplyr::mutate(Cell_ID = paste0(CELL_ID_PREFIX, '_', dplyr::row_number()))
    CELL_ID_COLUMN <- 'Cell_ID'
  }

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

