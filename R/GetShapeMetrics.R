#' @title GetShapeMetrics
#'
#' @description Analyzes spatial organization of cell nests using metrics like nest/boundary ratios, solidity, and boundary distance variation.
#'
#' @param INPUT A data frame containing cell coordinates and boundary annotations, where each row represents a cell in one single image.
#' @param CELL_ID_COLUMN The name of the column containing unique cell IDs.
#' @param X_POSITION The name of the column containing the x-coordinates of the cells.
#' @param Y_POSITION The name of the column containing the y-coordinates of the cells.
#' @param ANNO_COLUMN The name of the column containing the boundary and nest annotations.
#' @param ANNO_OF_BOUNDARY The name of boundary annotation derived from GetBoundary. Default is "Boundary".
#' @param ANNO_OF_NEST The name of nest annotation derived from GetBoundary. Default is "Nest".
#' @param SHAPE_METRICS A vector of metrics to calculate. Options:
#' \itemize{
#'   \item `Boundary2NestRatio`: Count ratio of boundary cells to nest cells (collective metric, only used when `SEPARATE_NESTS = FALSE`)
#'   \item `NestSolidity`: Area ratio of concave hull to convex hull of nest cells
#'   \item `Nest2BoundaryDistCoV`: Coefficient of variation of distances from nest cells to nearest boundary cells
#'   \item `NestFractalDimension`: Slope of log(N(s)) ~ log(1/s) regression with box-counting method (1 = simple shape, 2 = complex shape)
#'   \item `NestClarkEvansIndex`: Observed nearest neighbor distance divided by expected distance in Poisson distribution (<1 = clustered, 1 = random, >1 = dispersed)
#' }
#' Default is "Boundary2NestRatio".
#' @param SEPARATE_NESTS Analyze nests individually? When TRUE:
#' \itemize{
#'   \item Performs nest separation using Voronoi tessellation
#'   \item Returns metrics per nest with `NestSize` and `NestIDs`
#'   \item Excludes `Boundary2NestRatio`
#'   \item Requires `DIST_THRESHOLD` and `MIN_NEST_SIZE` parameters
#' }
#' Default is FALSE
#'
#' @param DIST_THRESHOLD Maximum connection distance for nest clustering (only when `SEPARATE_NESTS = TRUE`).
#' \itemize{
#'   \item Used only when `SEPARATE_NESTS = TRUE`
#'   \item Cells within this distance are considered connected
#'   \item Set to ~2× average cell spacing for biological tissues
#' }
#' @param MIN_NEST_SIZE Minimum number of cells required to qualify as a nest  (only when `SEPARATE_NESTS = TRUE`). Smaller clusters are ignored.
#'
#' @return A named list of calculated metrics. When `SEPARATE_NESTS = TRUE`, returns a list of lists containing metrics for each nest.
#' @export
#' @importFrom dplyr filter transmute mutate select arrange as_tibble
#' @importFrom purrr map set_names
#' @importFrom magrittr `%>%`
#' @importFrom concaveman concaveman
#' @importFrom sf st_polygon st_area
#' @importFrom dbscan kNN
#' @importFrom stats var sd
#' @examples
#' library(tidyverse)
#' library(patchwork)
#' library(Synora)
#'
#' # Generate Dummy Data
#' set.seed(123)
#' DummyData <- tidyr::expand_grid(X = seq(-1, 1, length.out = 65),
#'                                 Y = seq(-1, 1, length.out = 65)) %>%
#'   dplyr::nest_by(.key = 'Input') %>%
#'   tidyr::expand_grid(R = c(0.6, 0.7, 0.8, 0.9),
#'                      Sinusoid = c(0, 0.1, 0.2, 0.3, 0.4)) %>%
#'   dplyr::rowwise() %>%
#'   dplyr::mutate(Input = Input %>%
#'                   dplyr::mutate(
#'                     theta = 2 * atan(Y / (X + sqrt(X ^ 2 + Y ^ 2))),
#'                     theta = ifelse(is.na(theta), pi, theta),
#'                     CT = (X / R - Sinusoid * cos(theta) * sin(12 * theta)) ^ 2 +
#'                       (Y / R - Sinusoid * sin(theta) * sin(12 * theta)) ^ 2 <
#'                       (1 - Sinusoid) ^ 2) %>%
#'                   dplyr::mutate(CT = ifelse(CT & (as.logical(runif(n = nrow(.)) %/% 0.75)), !CT, CT) %>%
#'                                   as.numeric()) %>%
#'                   dplyr::mutate(X = 250 * X + rnorm(n = nrow(.)),
#'                                 Y = 250 * Y + rnorm(n = nrow(.))) %>%
#'                   dplyr::mutate(X = X %>% pmax(-250) %>% pmin(250),
#'                                 Y = Y %>% pmax(-250) %>% pmin(250)) %>%
#'                   dplyr::mutate(Cell_ID = paste0('Cell_', sprintf('%05.f', dplyr::row_number()))) %>%
#'                   list())
#' set.seed(NULL)
#'
#' # Run GetBoundary
#' Results <- DummyData %>%
#'   dplyr::mutate(BoundaryResult =
#'                   Synora::GetBoundary(
#'                     INPUT = Input,
#'                     X_POSITION = 'X',
#'                     Y_POSITION = 'Y',
#'                     ANNO_COLUMN = 'CT',
#'                     CELL_ID_COLUMN = 'Cell_ID',
#'                     RADIUS = 20,
#'                     NEST_SPECIFICITY = 0.25,
#'                     BOUNDARY_SPECIFICITY = 0.05
#'                   ) %>%
#'                   list()
#'   )
#' # Run GetShapeMetrics
#' Results <- Results %>%
#'   dplyr::mutate(ShapeMetrics =
#'                   Synora::GetShapeMetrics(
#'                     INPUT = BoundaryResult,
#'                     CELL_ID_COLUMN = 'Cell_ID',
#'                     X_POSITION = 'X',
#'                     Y_POSITION = 'Y',
#'                     ANNO_COLUMN = 'SynoraAnnotation',
#'                     ANNO_OF_BOUNDARY = 'Boundary',
#'                     ANNO_OF_NEST = 'Nest',
#'                     SHAPE_METRICS = c('Boundary2NestRatio', 'Nest2BoundaryDistCoV', 'NestSolidity', 'NestFractalDimension', 'NestClarkEvansIndex')
#'                   ) %>%
#'                   list()
#'   )
#'
#' # Visualization
#'
#' p <- patchwork::wrap_plots(
#'   nrow = 2,
#'   Results %>%
#'     dplyr::select(Input, R, Sinusoid) %>%
#'     tidyr::unnest(Input) %>%
#'     ggplot(aes(X, Y, color = as.factor(CT))) +
#'     facet_grid(R ~ Sinusoid) +
#'     geom_point(size = 1) +
#'     scale_color_manual(name = 'Cell Type',
#'                        values = c(`0` = '#e9c46a', `1` = '#046C9A'),
#'                        labels = c('Non-tumor cell', 'Tumor cell')) +
#'     labs(title = 'Cell Type') +
#'     theme_void() +
#'     coord_equal(),
#'   Results %>%
#'     dplyr::select(BoundaryResult, R, Sinusoid) %>%
#'     tidyr::unnest(BoundaryResult) %>%
#'     dplyr::mutate(BoundaryScore = pmin(BoundaryScore, 0.5)) %>%
#'     ggplot(aes(X, Y, color = BoundaryScore)) +
#'     facet_grid(R ~ Sinusoid) +
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
#'   Results %>%
#'     dplyr::select(BoundaryResult, R, Sinusoid) %>%
#'     tidyr::unnest(BoundaryResult) %>%
#'     ggplot(aes(X, Y, color = factor(SynoraAnnotation))) +
#'     facet_grid(R ~ Sinusoid) +
#'     geom_point(size = 1) +
#'     scale_color_manual(values = c(Boundary = '#4DAF4AFF', Nest = '#377EB8FF', Outside = '#984EA3FF'),
#'                        name = "Synora Annotation") +
#'     labs(title = "Synora Annotation") +
#'     theme_void() +
#'     coord_equal(),
#'   Results %>%
#'     dplyr::select(ShapeMetrics, R, Sinusoid) %>%
#'     dplyr::mutate(ShapeMetrics %>%
#'                     dplyr::as_tibble()) %>%
#'     dplyr::select(-c(ShapeMetrics)) %>%
#'     dplyr::ungroup() %>%
#'     tidyr::pivot_longer(cols = !c(R, Sinusoid), names_to = 'ShapeMetric', values_to = 'Value') %>%
#'     dplyr::nest_by(ShapeMetric) %>%
#'     dplyr::mutate(plot = {data %>%
#'         ggplot(aes(Sinusoid, R, fill = Value)) +
#'         geom_tile() +
#'         geom_text(aes(label = sprintf('%.2f', Value))) +
#'         labs(title = ShapeMetric) +
#'         theme_void() +
#'         scale_fill_gradient2(high = '#CB2314', guide = 'none') +
#'         coord_fixed(ratio = 1)} %>%
#'           list()) %>%
#'     dplyr::pull(plot) %>%
#'     patchwork::wrap_plots()) +
#'   patchwork::plot_layout(guides = "collect")
#' print(p)

GetShapeMetrics <- function(INPUT, X_POSITION, Y_POSITION,
                            CELL_ID_COLUMN, CELL_ID_PREFIX,
                            ANNO_COLUMN = 'SynoraAnnotation',
                            ANNO_OF_BOUNDARY = 'Boundary', ANNO_OF_NEST = 'Nest',
                            SHAPE_METRICS = c('Boundary2NestRatio'),
                            SEPARATE_NESTS = F, DIST_THRESHOLD = Inf, MIN_NEST_SIZE = 10) {
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

  valid_metrics <- c('Boundary2NestRatio', 'Nest2BoundaryDistCoV', 'NestSolidity', 'NestFractalDimension', 'NestClarkEvansIndex')
  invalid_metrics <- setdiff(SHAPE_METRICS, valid_metrics)

  if (length(invalid_metrics) > 0) {
    stop(
      "Invalid SHAPE_METRICS detected: \033[1;4;41m",
      paste(invalid_metrics, collapse = ", "),
      "\033[0m.\nValid SHAPE_METRICS: \033[1;32m",
      paste(valid_metrics, collapse = "\033[0m, \033[1;32m"),
      "\033[0m."
    )
  }

  INPUT <- INPUT %>%
    dplyr::transmute(Cell_ID = !!as.name(CELL_ID_COLUMN),
                     X = !!as.name(X_POSITION),
                     Y = !!as.name(Y_POSITION),
                     Anno = !!as.name(ANNO_COLUMN) %>%
                       forcats::fct_recode(Boundary = ANNO_OF_BOUNDARY, Nest = ANNO_OF_NEST) %>%
                       forcats::fct_drop()) %>%
    dplyr::arrange(Anno)

  calculateSolidity <- function(nest_cells) {
    if (nrow(nest_cells) < 3) return(NA_real_)
    areas <- c(1, 10000) %>%
      purrr::map(~ {
        nest_cells %>%
          as.matrix() %>%
          concaveman::concaveman(concavity = .x) %>%
          list() %>%
          sf::st_polygon() %>%
          sf::st_area()
      })
    areas[[1]] / areas[[2]]
  }

  calculateFractalDimension <- function(nest_cells) {
    if (nrow(nest_cells) < 1) return(NA_real_)
    nest_cells <- nest_cells %>% as.matrix()
    min_x <- min(nest_cells[, 1])
    max_x <- max(nest_cells[, 1])
    min_y <- min(nest_cells[, 2])
    max_y <- max(nest_cells[, 2])
    x_range <- max_x - min_x
    y_range <- max_y - min_y
    data_range <- max(x_range, y_range)
    bmin <- 0.1 * mean(dbscan::kNN(nest_cells, k = 1)$dist) * 2
    bmax <- data_range / 4
    bfact = 1.5
    oset = 10

    s_values <- numeric(0)
    s <- bmax
    while (s >= bmin) {
      s_values <- c(s_values, s)
      s <- s / bfact
    }

    if (length(s_values) == 0) stop("Invalid box size parameters")

    n_s <- numeric(length(s_values))

    for (i in seq_along(s_values)) {
      s <- s_values[i]
      min_count <- Inf

      for (j in 1:oset) {
        offset_x <- runif(1, 0, s)
        offset_y <- runif(1, 0, s)

        grid_x <- floor((nest_cells[, 1] - (min_x - offset_x)) / s)
        grid_y <- floor((nest_cells[, 2] - (min_y - offset_y)) / s)

        cell_ids <- paste(grid_x, grid_y, sep = ",")
        count <- length(unique(cell_ids))

        if (count < min_count) min_count <- count
      }

      n_s[i] <- min_count
    }
    valid <- n_s > 1 & s_values < Inf
    df <- data.frame(
      log_inv_s = log(1/s_values[valid]),
      log_n = log(n_s[valid])
    )
    if (nrow(df) < 3) {
      warning("Not enough valid nest cells for regression")
      return(NA)
    } else {
      result <- 3:nrow(df) %>%
        purrr::set_names() %>%
        purrr::map(function(k) {lm(log_n ~ log_inv_s, data = df[1:k, ])}) %>%
        tibble::enframe(name = 'k', value = 'fit') %>%
        dplyr::rowwise() %>%
        dplyr::mutate(r2 = summary(fit)$r.squared,
                      slope = coef(fit)[2]) %>%
        dplyr::ungroup() %>%
        dplyr::filter(r2 >= 0.99) %>%
        dplyr::slice_max(order_by = slope, n = 1) %>%
        dplyr::pull(slope) %>%
        .[[1]]
      return(result)
    }
  }
  calculateCE <- function(nest_cells) {
    if (nrow(nest_cells) < 1) return(NA_real_)
    nest_cells <- nest_cells %>% as.matrix()
    spatstat.geom::ppp(
      x = nest_cells[,1],
      y = nest_cells[,2],
      window = spatstat.geom::owin(range(nest_cells[,1]), range(nest_cells[,2]))) %>%
      spatstat.explore::clarkevans(correction = 'Donnelly')
  }

  calculateDistCoV <- function(nest_cells, boundary_cells) {
    if (nrow(nest_cells) < 1 || nrow(boundary_cells) < 2) return(NA_real_)
    distances <- dbscan::kNN(
      x = as.matrix(boundary_cells),
      query = as.matrix(nest_cells),
      k = 1
    )$dist[, 1]
    if (mean(distances) == 0) return(NA_real_)
    stats::sd(distances) / mean(distances)
  }

  if (!SEPARATE_NESTS) {
    RESULT <- SHAPE_METRICS %>%
      purrr::set_names() %>%
      purrr::map(function(X) {
        switch(X,
               "Boundary2NestRatio" = {
                 boundary_count <- sum(INPUT$Anno == 'Boundary', na.rm = TRUE)
                 nest_count <- sum(INPUT$Anno == 'Nest', na.rm = TRUE)
                 if (nest_count < MIN_NEST_SIZE || boundary_count == 0) NA_real_ else boundary_count / nest_count
               },
               "NestSolidity" = calculateSolidity(
                 nest_cells = INPUT %>% dplyr::filter(Anno == "Nest") %>% dplyr::select(X, Y)),
               "NestFractalDimension" = calculateFractalDimension(
                 nest_cells = INPUT %>% dplyr::filter(Anno == "Nest") %>% dplyr::select(X, Y)),
               "Nest2BoundaryDistCoV" = calculateDistCoV(
                 nest_cells = INPUT %>% dplyr::filter(Anno == "Nest") %>% dplyr::select(X, Y),
                 boundary_cells = INPUT %>% dplyr::filter(Anno == "Boundary") %>% dplyr::select(X, Y)),
               "NestClarkEvansIndex" = calculateCE(
                 nest_cells = INPUT %>% dplyr::filter(Anno == "Nest") %>% dplyr::select(X, Y)),
               NA_real_
        )
      })
  } else if (SEPARATE_NESTS & (sum(INPUT$Anno == 'Nest', na.rm = TRUE) >= MIN_NEST_SIZE)) {
    SeparatedNests <- .SeparateNests(
      INPUT = INPUT,
      X_POSITION = 'X',
      Y_POSITION = 'Y',
      ANNO_COLUMN = 'Anno',
      CELL_ID_COLUMN = 'Cell_ID',
      DIST_THRESHOLD = DIST_THRESHOLD,
      MIN_NEST_SIZE = MIN_NEST_SIZE,
      ANNO_OF_NEST = 'Nest') %>%
      dplyr::filter(!is.na(Nest_ID)) %>%
      base::split(.$Nest_ID)

    SHAPE_METRICS <- c("NestSize", "NestIDs", SHAPE_METRICS[SHAPE_METRICS != "Boundary2NestRatio"])

    RESULT <- SeparatedNests %>%
      purrr::map(function(INPUT_i) {
        SHAPE_METRICS %>%
          purrr::set_names() %>%
          purrr::map(function(X) {
            switch(X,
                   "NestSize" = sum(INPUT_i$Anno == 'Nest', na.rm = TRUE),
                   "NestIDs" = INPUT_i %>% dplyr::filter(Anno == "Nest") %>% dplyr::pull(Cell_ID),
                   "NestSolidity" = calculateSolidity(
                     nest_cells = INPUT_i %>% dplyr::filter(Anno == "Nest") %>% dplyr::select(X, Y)),
                   "NestFractalDimension" = calculateFractalDimension(
                     nest_cells = INPUT_i %>% dplyr::filter(Anno == "Nest") %>% dplyr::select(X, Y)),
                   "Nest2BoundaryDistCoV" = calculateDistCoV(
                     nest_cells = INPUT_i %>% dplyr::filter(Anno == "Nest") %>% dplyr::select(X, Y),
                     boundary_cells = INPUT %>% dplyr::filter(Anno == "Boundary") %>% dplyr::select(X, Y)),
                   "NestClarkEvansIndex" = calculateCE(
                     nest_cells = INPUT_i %>% dplyr::filter(Anno == "Nest") %>% dplyr::select(X, Y)),
                   NA_real_
            )
          })
      })
  } else {
    warning("Not enough valid nest cells for separation")
    RESULT <- list() %>% setNames(character(0))
  }
  return(RESULT)
}

.SeparateNests <- function(INPUT, X_POSITION, Y_POSITION,
                           ANNO_COLUMN, CELL_ID_COLUMN,
                           ANNO_OF_NEST = 'Nest',
                           DIST_THRESHOLD = Inf,
                           MIN_NEST_SIZE = 10) {
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

  DELDIR <- deldir::deldir(x = INPUT[[X_POSITION]],
                           y = INPUT[[Y_POSITION]],
                           rw = c(min(INPUT[[X_POSITION]]),
                                  max(INPUT[[X_POSITION]]),
                                  min(INPUT[[Y_POSITION]]),
                                  max(INPUT[[Y_POSITION]]))
  )

  EDGES <- DELDIR$delsgs %>%
    dplyr::mutate(dist = sqrt((x1 - x2) ^ 2 + (y1 - y2) ^ 2),
                  from = INPUT[[CELL_ID_COLUMN]][ind1],
                  to = INPUT[[CELL_ID_COLUMN]][ind2],
                  anno_from = INPUT[[ANNO_COLUMN]][ind1],
                  anno_to = INPUT[[ANNO_COLUMN]][ind2]) %>%
    dplyr::filter(anno_from == ANNO_OF_NEST & anno_to == ANNO_OF_NEST) %>%
    dplyr::filter(dist < DIST_THRESHOLD)

  RESULT <- tidygraph::tbl_graph(
    nodes = INPUT,
    edges = EDGES,
    node_key = 'Cell_ID',
    directed = F) %>%
    tidygraph::activate(nodes) %>%
    dplyr::mutate(n_group_id = tidygraph::group_components()) %>%
    dplyr::group_by(n_group_id) %>%
    dplyr::mutate(n_group_size = n()) %>%
    dplyr::ungroup() %>%
    tidygraph::activate(nodes) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(Nest_ID = ifelse(n_group_size >= MIN_NEST_SIZE & !!as.name(ANNO_COLUMN) == ANNO_OF_NEST, paste0('Nest_', n_group_id), NA))


  return(RESULT)
}




