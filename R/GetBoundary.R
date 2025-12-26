#' @title GetBoundary
#' @description Detect Boundary Cells
#' @param INPUT Input data frame containing cell coordinates and annotations, where each row represents a cell in one single image.
#' @param X_POSITION Name of X coordinate column.
#' @param Y_POSITION Name of Y coordinate column.
#' @param ANNO_COLUMN Name of annotation column, where 1 indicates cells which constitute nest (i.e. tumor cells if tumor nests are to be identified) and 0 indicates other cells. Can be binary or continuous.
#' @param CELL_ID_COLUMN (optional) Cell ID column. If not provided, row indices will be used as cell IDs. Default NULL.
#' @param CELL_ID_PREFIX (optional) Cell ID prefix used for creating cell IDs. Default NULL.
#' @param ANNO_RANGE Annotation range. Default [0, 1].
#' @param ANNO_MIDPOINT Annotation midpoint. Either numeric value or "auto" to automatically detect using EM-based cutoff determination. Default 0.5.
#' @param NEIGHBOR_METHOD Neighborhood detection method: "radius" (fixed radius), "knn" (k-nearest neighbors), or "hybrid" (kNN within radius). Default "radius".
#' @param RADIUS Neighborhood radius (required for "radius" and "hybrid" methods).
#' @param KNN_K Number of nearest neighbors (required for "knn" and "hybrid" methods).
#' @param DENOISE Whether to perform noise removal (TRUE) or keep all cells (FALSE). Default TRUE.
#' @param NEST_MIN_SIZE Minimum nest size. Default 5.
#' @param NEST_SPECIFICITY Nest specificity. Default 0.25, recommended between 0.1 and 0.4.
#' @param BOUNDARY_SPECIFICITY Boundary specificity. Default 0.05, recommended between 0.01 and 0.1.
#' @return A data frame with the following columns:
#' \describe{
#'   \item{CELL_ID}{Cell ID column.}
#'   \item{X_POSITION}{X coordinate column.}
#'   \item{Y_POSITION}{Y coordinate column.}
#'   \item{ANNO_COLUMN}{Original cell annotation column.}
#'   \item{Nb_Count}{Number of neighbors considered.}
#'   \item{BoundaryScore}{Boundary Score}
#'   \item{SynoraAnnotation}{Synora annotation column containing 'Boundary', 'Nest', 'Outside' and 'Noise'}
#' }
#' @export
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
#' # Visualization
#' p <- patchwork::wrap_plots(
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
#'     coord_equal()) +
#'   patchwork::plot_layout(guides = "collect")
#' print(p)

GetBoundary <- function(INPUT, X_POSITION, `Y_POSITION`,
                        ANNO_COLUMN, CELL_ID_COLUMN, CELL_ID_PREFIX,
                        ANNO_RANGE = c(0, 1), ANNO_MIDPOINT = 0.5,
                        NEIGHBOR_METHOD = c("radius", "knn", "hybrid"),
                        RADIUS, KNN_K,
                        DENOISE = TRUE,
                        NEST_MIN_SIZE = 5, NEST_SPECIFICITY = 0.25,
                        BOUNDARY_SPECIFICITY = 0.05,
                        VERBOSE = FALSE) {
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
  if (!is.numeric(INPUT[[ANNO_COLUMN]])) stop("ANNO_COLUMN must be a numeric. Convert to numeric first")
  if (any(is.na(INPUT[[ANNO_COLUMN]]) | is.nan(INPUT[[ANNO_COLUMN]]) | is.infinite(INPUT[[ANNO_COLUMN]]))) {
    stop("ANNO_COLUMN contains NA, NaN, or Inf values. Please clean the data first")
  }

  NEIGHBOR_METHOD <- match.arg(NEIGHBOR_METHOD)

  if (NEIGHBOR_METHOD %in% c("radius", "hybrid") && missing(RADIUS))
    stop("RADIUS required for this method")
  if (NEIGHBOR_METHOD %in% c("knn", "hybrid") && missing(KNN_K))
    stop("KNN_K required for this method")
  if (NEIGHBOR_METHOD == "radius") {
    if (!is.numeric(RADIUS) || RADIUS <= 0)
      stop("RADIUS must be positive numeric")
  }
  if (NEIGHBOR_METHOD %in% c("knn", "hybrid")) {
    if (!is.numeric(KNN_K) || KNN_K <= 0 || KNN_K != as.integer(KNN_K))
      stop("KNN_K must be positive integer")
  }

  neighbor_params <- list(
    method = NEIGHBOR_METHOD,
    radius = if (NEIGHBOR_METHOD %in% c("radius", "hybrid")) RADIUS else NULL,
    knn_k = if (NEIGHBOR_METHOD %in% c("knn", "hybrid")) KNN_K else NULL
  )


  if (!identical(ANNO_RANGE, "auto")) {
    if (!is.numeric(ANNO_RANGE) || length(ANNO_RANGE) != 2 || ANNO_RANGE[1] >= ANNO_RANGE[2]) {
      stop("ANNO_RANGE must be 'auto' or numeric vector of length 2 with min < max")
    }
  } else {
    if (!is.numeric(INPUT[[ANNO_COLUMN]])) stop("ANNO_COLUMN must be numeric for automatic range detection")
    ANNO_RANGE <- c(min(INPUT[[ANNO_COLUMN]], na.rm = TRUE), max(INPUT[[ANNO_COLUMN]], na.rm = TRUE))

    if (diff(ANNO_RANGE) == 0) {
      warning("ANNO_COLUMN has zero range, adjusting with epsilon")
      ANNO_RANGE <- ANNO_RANGE + c(-1e-6, 1e-6)
    }
  }

  if (!(is.numeric(ANNO_MIDPOINT) && length(ANNO_MIDPOINT) == 1) &&
      ANNO_MIDPOINT != "auto") {
    stop("ANNO_MIDPOINT must be either numeric or 'auto'")
  }

  if (ANNO_MIDPOINT == "auto") {
    captured <- capture.output({
      ANNO_MIDPOINT <- .FindCutoff(INPUT[[ANNO_COLUMN]])
    }, type = "output")

    auto_msg <- paste0(
      if (length(captured) > 0) paste0(paste0(captured, collapse = " "), " | "),
      "Automatically detected ANNO_MIDPOINT: ",
      round(ANNO_MIDPOINT, 3),
      ' [', round(ANNO_RANGE[1], 3), ',', round(ANNO_RANGE[2], 3), ']'
    )

    if (VERBOSE) message(auto_msg)

    if (!dplyr::between(ANNO_MIDPOINT, ANNO_RANGE[1], ANNO_RANGE[2])) {
      warning("Auto-detected midpoint outside ANNO_RANGE, using range midpoint")
      ANNO_MIDPOINT <- mean(ANNO_RANGE)
    }
  }

  INPUT <- INPUT %>%
    dplyr::mutate(
      !!as.name(paste0(ANNO_COLUMN, '_scaled')) := ifelse(
        !!as.name(ANNO_COLUMN) >= ANNO_MIDPOINT,
        (!!as.name(ANNO_COLUMN) - ANNO_MIDPOINT) / (ANNO_RANGE[2] - ANNO_MIDPOINT),
        (!!as.name(ANNO_COLUMN) - ANNO_MIDPOINT) / (ANNO_MIDPOINT - ANNO_RANGE[1])
      )
    )


  if (DENOISE) {
    RESULT_1 <-INPUT %>%
      .GetMO(CELL_ID_COLUMN = CELL_ID_COLUMN,
             X_POSITION = X_POSITION,
             Y_POSITION = Y_POSITION,
             ANNO_COLUMN = paste0(ANNO_COLUMN, '_scaled'),
             NEIGHBOR_METHOD = neighbor_params$method,
             RADIUS = neighbor_params$radius,
             KNN_K = neighbor_params$knn_k,
             ORIENTEDNESS = FALSE) %>%
      dplyr::mutate(Nest = (0.5 * Mean_Anno + 0.5) >= NEST_SPECIFICITY)

    NOISE_COUNT <- c()
    i <- 1
    while (T) {
      RESULT_1 <- RESULT_1 %>%
        dplyr::nest_by(Nest, NeighborhoodRadius) %>%
        dplyr::mutate(!!as.name(paste0('Noise_', i)) := data %>%
                        dplyr::select(all_of(c(X_POSITION, Y_POSITION))) %>%
                        dbscan::dbscan(eps = NeighborhoodRadius, minPts = NEST_MIN_SIZE) %>%
                        .$cluster %>%
                        {. == 0} %>%
                        list()) %>%
        tidyr::unnest(cols = c(data, !!as.name(paste0('Noise_', i))))  %>%
        dplyr::ungroup() %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(Nest = ifelse(!!as.name(paste0('Noise_', i)), !Nest, Nest))
      NOISE_COUNT <- c(NOISE_COUNT, sum(RESULT_1[[paste0('Noise_', i)]]))

      if (ifelse(length(NOISE_COUNT) == 1, NOISE_COUNT == 0, diff(tail(NOISE_COUNT, 2)) == 0)) {
        RESULT_1 <- RESULT_1 %>%
          dplyr::mutate(!!as.name(paste0(ANNO_COLUMN, '_denoised')) := dplyr::case_when(
            !(!!as.name(paste0('Noise_', i))) & Nest ~ 1,
            !(!!as.name(paste0('Noise_', i))) & !Nest ~ -1,
            T ~ NA,
          ))
        break
      } else {
        i <- i + 1
      }
    }
    RESULT <- RESULT_1 %>%
      .GetMO(CELL_ID_COLUMN = CELL_ID_COLUMN,
             X_POSITION = X_POSITION,
             Y_POSITION = Y_POSITION,
             ANNO_COLUMN = paste0(ANNO_COLUMN, '_denoised'),
             NEIGHBOR_METHOD = neighbor_params$method,
             RADIUS = neighbor_params$radius,
             KNN_K = neighbor_params$knn_k,
             ORIENTEDNESS = T) %>%
      dplyr::mutate(BoundaryScore = Mixedness * Orientedness) %>%
      dplyr::mutate(SynoraAnnotation := dplyr::case_when(
        BoundaryScore >= BOUNDARY_SPECIFICITY ~  'Boundary',
        BoundaryScore < BOUNDARY_SPECIFICITY & !!as.name(paste0(ANNO_COLUMN, '_denoised')) == 1 ~  'Nest',
        BoundaryScore < BOUNDARY_SPECIFICITY & !!as.name(paste0(ANNO_COLUMN, '_denoised')) == -1 ~  'Outside',
        T ~  'Noise') %>%
          forcats::fct_expand('Boundary', 'Nest', 'Outside', 'Noise') %>%
          forcats::fct_relevel('Boundary', 'Nest', 'Outside', 'Noise')) %>%
      dplyr::arrange(!!as.name(CELL_ID_COLUMN)) %>%
      dplyr::left_join(INPUT, ., by = c(CELL_ID_COLUMN, X_POSITION, Y_POSITION)) %>%
      dplyr::transmute(!!as.name(CELL_ID_COLUMN),
                       !!as.name(X_POSITION),
                       !!as.name(Y_POSITION),
                       !!as.name(ANNO_COLUMN),
                       Nb_Count,
                       Anno_Midpoint = ANNO_MIDPOINT,
                       Mixedness = ifelse(is.na(Mixedness), 0, Mixedness),
                       Orientedness = ifelse(is.na(Orientedness), 0, Orientedness),
                       BoundaryScore = ifelse(is.na(BoundaryScore), 0, BoundaryScore),
                       SynoraAnnotation = SynoraAnnotation)
  } else {
    RESULT <- INPUT %>%
      .GetMO(CELL_ID_COLUMN = CELL_ID_COLUMN,
             X_POSITION = X_POSITION,
             Y_POSITION = Y_POSITION,
             ANNO_COLUMN = paste0(ANNO_COLUMN, '_scaled'),
             NEIGHBOR_METHOD = neighbor_params$method,
             RADIUS = neighbor_params$radius,
             KNN_K = neighbor_params$knn_k,
             ORIENTEDNESS = T) %>%
      dplyr::mutate(BoundaryScore = Mixedness * Orientedness) %>%
      dplyr::mutate(SynoraAnnotation := dplyr::case_when(
        BoundaryScore >= BOUNDARY_SPECIFICITY ~  'Boundary',
        (0.5 * Mean_Anno + 0.5) >= NEST_SPECIFICITY ~ 'Nest',
        T ~  'Outside') %>%
          forcats::fct_expand('Boundary', 'Nest', 'Outside', 'Noise') %>%
          forcats::fct_relevel('Boundary', 'Nest', 'Outside', 'Noise')) %>%
      dplyr::arrange(!!as.name(CELL_ID_COLUMN)) %>%
      dplyr::left_join(INPUT, ., by = c(CELL_ID_COLUMN, X_POSITION, Y_POSITION)) %>%
      dplyr::transmute(!!as.name(CELL_ID_COLUMN),
                       !!as.name(X_POSITION),
                       !!as.name(Y_POSITION),
                       !!as.name(ANNO_COLUMN),
                       Nb_Count,
                       Anno_Midpoint = ANNO_MIDPOINT,
                       Mixedness = ifelse(is.na(Mixedness), 0, Mixedness),
                       Orientedness = ifelse(is.na(Orientedness), 0, Orientedness),
                       BoundaryScore = ifelse(is.na(BoundaryScore), 0, BoundaryScore),
                       SynoraAnnotation = SynoraAnnotation)
  }
  return(RESULT)
}

#' @inheritParams GetBoundary
#' @keywords internal
.GetMO <- function(INPUT, CELL_ID_COLUMN, X_POSITION, Y_POSITION,
                   ANNO_COLUMN, NEIGHBOR_METHOD = c("radius", "knn", "hybrid"),
                   RADIUS = NULL, KNN_K = NULL, ORIENTEDNESS = TRUE) {

  NEIGHBOR_METHOD <- match.arg(NEIGHBOR_METHOD)

  INPUT <- INPUT %>%
    dplyr::transmute(Cell_ID = !!as.name(CELL_ID_COLUMN),
                     X = !!as.name(X_POSITION),
                     Y = !!as.name(Y_POSITION),
                     Anno = !!as.name(ANNO_COLUMN))

  if (NEIGHBOR_METHOD == "radius") {
    NN <- INPUT %>%
      dplyr::select(X, Y) %>%
      dbscan::frNN(eps = RADIUS, sort = FALSE)
    NeighborhoodRadius <- RADIUS
  }
  else if (NEIGHBOR_METHOD == "knn") {
    NN <- INPUT %>%
      dplyr::select(X, Y) %>%
      dbscan::kNN(k = KNN_K, sort = FALSE)
    NeighborhoodRadius <- stats::quantile(NN$dist[, KNN_K], 0.75)
    NN <- list(id = split(NN$id, seq(nrow(NN$id))),
               dist = split(NN$dist, seq(nrow(NN$dist))))
  }
  else if (NEIGHBOR_METHOD == "hybrid") {
    NN <- INPUT %>%
      dplyr::select(X, Y) %>%
      dbscan::kNN(k = KNN_K, sort = FALSE)
    valid_mask <- NN$dist <= RADIUS
    NN <- list(
      id = lapply(1:nrow(NN$id), \(i) NN$id[i, valid_mask[i, ]]),
      dist = lapply(1:nrow(NN$dist), \(i) NN$dist[i, valid_mask[i, ]])
    )
    NeighborhoodRadius <- RADIUS
  }

  if (ORIENTEDNESS) {
    RESULT <- INPUT %>%
      dplyr::mutate(ID = NN$id,
                    Dist = NN$dist) %>%
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
      dplyr::mutate(NeighborhoodRadius = NeighborhoodRadius) %>%
      dplyr::select(Cell_ID, X, Y, Anno, Nb_Count, NeighborhoodRadius, Mean_Anno, Mixedness, Orientedness) %>%
      dplyr::rename(!!as.name(CELL_ID_COLUMN) := Cell_ID,
                    !!as.name(X_POSITION) := X,
                    !!as.name(Y_POSITION) := Y,
                    !!as.name(ANNO_COLUMN) := Anno)
  } else {
    RESULT <- INPUT %>%
      dplyr::mutate(ID = NN$id) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Nb_Anno = list(.$Anno[ID])) %>%
      dplyr::mutate(Nb_Count = length(Nb_Anno),
                    Mean_Anno = mean(Nb_Anno)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Mixedness = 1 - (Mean_Anno) ^ 2) %>%
      dplyr::mutate(NeighborhoodRadius = NeighborhoodRadius) %>%
      dplyr::select(Cell_ID, X, Y, Anno, Nb_Count, NeighborhoodRadius, Mean_Anno, Mixedness) %>%
      dplyr::rename(!!as.name(CELL_ID_COLUMN) := Cell_ID,
                    !!as.name(X_POSITION) := X,
                    !!as.name(Y_POSITION) := Y,
                    !!as.name(ANNO_COLUMN) := Anno)
  }
  return(RESULT)
}

.FindCutoff <- function(data, maxit = 1000, eps = 1e-6, seed = 123) {
  tryCatch({
    set.seed(seed)

    if (length(unique(data)) < 5 || sd(data) < 1e-3) {
      message("Insufficient variation for bimodal detection")
      return(median(data))
    }

    suppressWarnings({
      mix <- mixtools::normalmixEM(
        data, k = 2,
        maxit = maxit,
        epsilon = eps
      )
    })

    if (mix$mu[1] > mix$mu[2]) {
      mu <- mix$mu[2:1]
      sigma <- mix$sigma[2:1]
    } else {
      mu <- mix$mu
      sigma <- mix$sigma
    }

    search_min <- max(min(data), mu[1] - 2 * sigma[1])
    search_max <- min(max(data), mu[2] + 2 * sigma[2])

    if (search_min >= search_max) {
      message("Overlapping search range, using midpoint between modes")
      return(mean(mu))
    }

    f <- function(x) dnorm(x, mu[1], sigma[1]) - dnorm(x, mu[2], sigma[2])

    if (f(search_min) * f(search_max) >= 0) {
      message("No sign change detected, using weighted average")
      return((mu[1] * sigma[2] + mu[2] * sigma[1]) / (sigma[1] + sigma[2]))
    }

    uniroot(f, interval = c(search_min, search_max),
            extendInt = "yes", tol = .Machine$double.eps^0.5)$root
  }, error = function(e) {
    message("Fallback to median: ", e$message)
    median(data)
  })
}


