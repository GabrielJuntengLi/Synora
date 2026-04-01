# .ringPerimeter ----
.ringPerimeter <- function(coords) {
  dx <- diff(coords[, 1]); dy <- diff(coords[, 2])
  sum(sqrt(dx^2 + dy^2))
}


# .deldirEdges ----
.deldirEdges <- function(df, CELL_ID_COLUMN, ANNO_COLUMN, ANNO_OF_NEST,
                         DIST_THRESHOLD, x_vec, y_vec) {
  tryCatch({
    n   <- nrow(df)
    DEL <- deldir::deldir(
      x  = x_vec, y  = y_vec,
      rw = c(min(x_vec), max(x_vec), min(y_vec), max(y_vec))
    )
    DEL$delsgs %>%
      dplyr::filter(ind1 <= n, ind2 <= n) %>%
      dplyr::mutate(
        dist      = sqrt((x1 - x2)^2 + (y1 - y2)^2),
        from      = df[[CELL_ID_COLUMN]][ind1],
        to        = df[[CELL_ID_COLUMN]][ind2],
        anno_from = df[[ANNO_COLUMN]][ind1],
        anno_to   = df[[ANNO_COLUMN]][ind2]
      ) %>%
      dplyr::filter(!is.na(from), !is.na(to),
                    anno_from == ANNO_OF_NEST, anno_to == ANNO_OF_NEST,
                    dist < DIST_THRESHOLD) %>%
      dplyr::select(from, to)
  }, error = function(e) {
    warning("[Synora] deldir() failed: ", conditionMessage(e), call. = FALSE)
    NULL
  })
}


# .buildNestGraph ----
.buildNestGraph <- function(INPUT, CELL_ID_COLUMN, ANNO_COLUMN, ANNO_OF_NEST,
                            DIST_THRESHOLD, WINDOW_SIZE, WINDOW_OVERLAP,
                            LARGE_N_THRESHOLD) {
  x_vec   <- INPUT[["X"]]
  y_vec   <- INPUT[["Y"]]
  n_cells <- nrow(INPUT)
  
  if (n_cells <= LARGE_N_THRESHOLD) {
    message("[Synora] Global Delaunay (n = ", n_cells, ").")
    EDGES <- .deldirEdges(INPUT, CELL_ID_COLUMN, ANNO_COLUMN, ANNO_OF_NEST,
                          DIST_THRESHOLD, x_vec, y_vec)
  } else {
    message("[Synora] Windowed Delaunay (n = ", n_cells, ").")
    x_min <- min(x_vec); x_max <- max(x_vec)
    y_min <- min(y_vec); y_max <- max(y_vec)
    x_breaks <- seq(x_min, x_max, by = WINDOW_SIZE)
    if (utils::tail(x_breaks, 1) < x_max) x_breaks <- c(x_breaks, x_max)
    y_breaks <- seq(y_min, y_max, by = WINDOW_SIZE)
    if (utils::tail(y_breaks, 1) < y_max) y_breaks <- c(y_breaks, y_max)
    
    all_tile_edges <- list()
    tile_idx <- 0L
    for (xi in seq_len(length(x_breaks) - 1L)) {
      for (yi in seq_len(length(y_breaks) - 1L)) {
        tile_idx <- tile_idx + 1L
        cx_lo <- x_breaks[xi];   cx_hi <- x_breaks[xi + 1L]
        cy_lo <- y_breaks[yi];   cy_hi <- y_breaks[yi + 1L]
        ex_lo <- cx_lo - WINDOW_OVERLAP; ex_hi <- cx_hi + WINDOW_OVERLAP
        ey_lo <- cy_lo - WINDOW_OVERLAP; ey_hi <- cy_hi + WINDOW_OVERLAP
        in_ext  <- x_vec >= ex_lo & x_vec <= ex_hi & y_vec >= ey_lo & y_vec <= ey_hi
        in_core <- x_vec >= cx_lo & x_vec <  cx_hi & y_vec >= cy_lo & y_vec <  cy_hi
        if (xi == length(x_breaks) - 1L) in_core <- in_core | (x_vec == cx_hi)
        if (yi == length(y_breaks) - 1L) in_core <- in_core | (y_vec == cy_hi)
        tile_df  <- INPUT[in_ext, , drop = FALSE]; row.names(tile_df) <- NULL
        core_ids <- INPUT[[CELL_ID_COLUMN]][in_core]
        if (nrow(tile_df) < 3L || length(core_ids) < 1L) next
        te <- .deldirEdges(tile_df, CELL_ID_COLUMN, ANNO_COLUMN, ANNO_OF_NEST,
                           DIST_THRESHOLD, tile_df[["X"]], tile_df[["Y"]])
        if (!is.null(te) && nrow(te) > 0)
          all_tile_edges[[tile_idx]] <- te %>%
          dplyr::filter(from %in% core_ids | to %in% core_ids)
      }
    }
    EDGES <- dplyr::bind_rows(all_tile_edges)
    if (!is.null(EDGES) && nrow(EDGES) > 0)
      EDGES <- EDGES %>%
      dplyr::mutate(e_lo = pmin(from, to), e_hi = pmax(from, to)) %>%
      dplyr::distinct(e_lo, e_hi, .keep_all = TRUE) %>%
      dplyr::select(-e_lo, -e_hi)
  }
  
  if (is.null(EDGES) || nrow(EDGES) == 0)
    EDGES <- data.frame(from = character(0), to = character(0))
  
  tidygraph::tbl_graph(nodes = INPUT, edges = EDGES,
                       node_key = CELL_ID_COLUMN, directed = FALSE) %>%
    tidygraph::activate(nodes) %>%
    dplyr::mutate(n_group_id = tidygraph::group_components())
}


# .SeparateNests ----
# Separate nests via Delaunay graph components; returns tibble with Nest_ID.
.SeparateNests <- function(INPUT, ANNO_COLUMN, CELL_ID_COLUMN,
                           ANNO_OF_NEST, DIST_THRESHOLD, MIN_NEST_SIZE,
                           WINDOW_SIZE, WINDOW_OVERLAP, LARGE_N_THRESHOLD) {
  .buildNestGraph(INPUT, CELL_ID_COLUMN, ANNO_COLUMN, ANNO_OF_NEST,
                  DIST_THRESHOLD, WINDOW_SIZE, WINDOW_OVERLAP, LARGE_N_THRESHOLD) %>%
    tidygraph::activate(nodes) %>%
    dplyr::group_by(n_group_id) %>%
    dplyr::mutate(n_group_size = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(
      Nest_ID = dplyr::if_else(
        n_group_size >= MIN_NEST_SIZE & !!as.name(ANNO_COLUMN) == ANNO_OF_NEST,
        paste0("Nest_", n_group_id),
        NA_character_
      )
    )
}


# .chaikinSmooth ----
# Smooth a closed polygon ring using Chaikin's corner-cutting algorithm.
# Iteratively replaces each edge with two points at 1/4 and 3/4 along it.
# iterations: more = smoother (2-3 is usually sufficient).
.chaikinSmooth <- function(ring, iterations = 3L) {
  pts <- ring[-nrow(ring), , drop = FALSE]   # drop closing duplicate
  n   <- nrow(pts)
  if (n < 3L) return(ring)
  
  for (iter in seq_len(iterations)) {
    n        <- nrow(pts)
    idx_curr <- seq_len(n)
    idx_next <- c(seq(2L, n), 1L)
    
    q <- 0.75 * pts[idx_curr, , drop = FALSE] +
      0.25 * pts[idx_next, , drop = FALSE]
    r <- 0.25 * pts[idx_curr, , drop = FALSE] +
      0.75 * pts[idx_next, , drop = FALSE]
    
    new_pts                          <- matrix(0, nrow = 2L * n, ncol = 2L)
    new_pts[seq(1L, 2L * n, 2L), ]  <- q
    new_pts[seq(2L, 2L * n, 2L), ]  <- r
    pts <- new_pts
  }
  
  rbind(pts, pts[1L, , drop = FALSE])   # re-close
}


# .buildNestHull ----
# Build a smoothed concave hull for one nest.
#   1. concaveman on nest cells  -> correct shape, jagged
#   2. Chaikin smoothing         -> smooth perimeter for metrics
#   3. Convex hull fallback      -> < 3 unique positions
.buildNestHull <- function(nc_ids, all_pts, anno_nest,
                           boundary_cells, dist_threshold, concavity,
                           smooth_iterations = 5L) {
  
  nc     <- all_pts[all_pts$Cell_ID %in% nc_ids, , drop = FALSE]
  if (nrow(nc) < 3L) return(NULL)
  mat_nc <- unique(as.matrix(nc[, c("X", "Y")]))
  if (nrow(mat_nc) < 3L) return(NULL)
  
  # Case A: concave hull
  ring <- tryCatch(
    concaveman::concaveman(mat_nc, concavity = concavity),
    error = function(e) NULL
  )
  
  # Case B: convex hull fallback
  if (is.null(ring) || nrow(ring) < 3L) {
    ring <- tryCatch({
      mp  <- sf::st_multipoint(mat_nc)
      raw <- sf::st_coordinates(sf::st_convex_hull(mp))
      raw[, c("X", "Y"), drop = FALSE]
    }, error = function(e) NULL)
  }
  
  if (is.null(ring)) return(NULL)
  
  # Chaikin smoothing — removes cell-level jaggedness while preserving
  # overall concave shape; does NOT affect IsClosedNest or cell-based metrics.
  if (smooth_iterations > 0L && nrow(ring) >= 3L) {
    smoothed <- tryCatch(
      .chaikinSmooth(ring, iterations = smooth_iterations),
      error = function(e) ring
    )
    valid <- tryCatch(
      sf::st_is_valid(sf::st_polygon(list(smoothed))),
      error = function(e) FALSE
    )
    if (isTRUE(valid)) ring <- smoothed
  }
  
  ring
}


# .classifyNestClosure ----
# Classify a nest as open (touches tissue margin) or closed.
# A nest is CLOSED when <= MARGIN_TOLERANCE fraction of hull vertices
# lie within MARGIN_DIST of the tissue outer margin line.
.classifyNestClosure <- function(hull_ring, margin_line,
                                 MARGIN_DIST, MARGIN_TOLERANCE) {
  if (is.null(hull_ring) || nrow(hull_ring) < 2)
    return(list(IsClosedNest = FALSE, ClosedFraction = NA_real_))
  dists <- vapply(seq_len(nrow(hull_ring) - 1L), function(k)
    as.numeric(sf::st_distance(sf::st_point(hull_ring[k, ]), margin_line)),
    numeric(1))
  margin_frac <- mean(dists <= MARGIN_DIST)
  list(IsClosedNest   = margin_frac <= MARGIN_TOLERANCE,
       ClosedFraction = 1 - margin_frac)
}


# .calcPerNestMetrics ----
# Compute all requested per-nest shape metrics; returns a named list.
.calcPerNestMetrics <- function(nc_i, hull_ring, is_closed,
                                bc_near, margin_line,
                                requested, BOUNDARY_RADIUS, SEED) {
  out <- list()
  
  # Hull-derived geometry
  ch_poly  <- sf::st_convex_hull(sf::st_multipoint(as.matrix(nc_i[, c("X", "Y")])))
  ch_area  <- as.numeric(sf::st_area(ch_poly))
  ch_raw   <- sf::st_coordinates(ch_poly)
  ch_perim <- .ringPerimeter(ch_raw[, c("X", "Y"), drop = FALSE])
  
  cc_area  <- tryCatch(
    as.numeric(sf::st_area(sf::st_polygon(list(hull_ring)))),
    error = function(e) NA_real_)
  cc_perim <- if (!is.null(hull_ring) && nrow(hull_ring) >= 2)
    .ringPerimeter(hull_ring) else NA_real_
  
  if ("NestSolidity" %in% requested)
    out$NestSolidity <- if (anyNA(c(cc_area, ch_area)) || ch_area == 0)
      NA_real_ else cc_area / ch_area
  
  if ("BoundaryRoughness" %in% requested)
    out$BoundaryRoughness <- if (anyNA(c(cc_perim, ch_perim)) || ch_perim == 0)
      NA_real_ else cc_perim / ch_perim
  
  if ("NestCompactness" %in% requested)
    out$NestCompactness <- if (anyNA(c(cc_area, cc_perim)) || cc_perim == 0)
      NA_real_ else (4 * pi * cc_area) / cc_perim^2
  
  if ("NestElongation" %in% requested) {
    mat <- as.matrix(nc_i[, c("X", "Y")])
    pca <- tryCatch(stats::prcomp(mat, center = TRUE, scale. = FALSE),
                    error = function(e) NULL)
    out$NestElongation <- if (is.null(pca) || length(pca$sdev) < 2 || pca$sdev[1] == 0)
      NA_real_ else pca$sdev[2] / pca$sdev[1]
  }
  
  if ("NestFractalDimension" %in% requested)
    out$NestFractalDimension <- .calcFractalDimension(nc_i, SEED)
  
  # Closed-nest-only metrics
  closed_metrics <- c("BoundaryIntegrity", "NestBoundaryBalance",
                      "NestEnclosure", "NestBoundaryDensity")
  for (m in intersect(requested, closed_metrics)) {
    out[[m]] <- if (!is_closed) NA_real_ else switch(m,
                                                     "BoundaryIntegrity" = {
                                                       if (nrow(bc_near) == 0) 0
                                                       else {
                                                         mat_nc <- as.matrix(nc_i[, c("X", "Y")])
                                                         mat_bc <- as.matrix(bc_near[, c("X", "Y")])
                                                         mean(vapply(seq_len(nrow(mat_nc)), function(i)
                                                           any(sqrt((mat_bc[, 1] - mat_nc[i, 1])^2 +
                                                                      (mat_bc[, 2] - mat_nc[i, 2])^2) <= BOUNDARY_RADIUS),
                                                           logical(1)))
                                                       }
                                                     },
                                                     "NestBoundaryBalance" = {
                                                       if (nrow(bc_near) == 0) return(NA_real_)
                                                       cx     <- mean(nc_i$X); cy <- mean(nc_i$Y)
                                                       angles <- atan2(bc_near$Y - cy, bc_near$X - cx)
                                                       breaks <- seq(-pi, pi, length.out = 9L)
                                                       counts <- tabulate(findInterval(angles, breaks, rightmost.closed = TRUE),
                                                                          nbins = 8L)
                                                       counts <- counts[counts > 0]
                                                       if (length(counts) < 2) 0
                                                       else { p <- counts / sum(counts); -sum(p * log(p)) / log(8) }
                                                     },
                                                     "NestEnclosure" = {
                                                       mat <- as.matrix(nc_i[, c("X", "Y")])
                                                       mean(vapply(seq_len(nrow(mat)), function(k)
                                                         as.numeric(sf::st_distance(sf::st_point(mat[k, ]), margin_line)),
                                                         numeric(1)))
                                                     },
                                                     "NestBoundaryDensity" = {
                                                       if (is.na(cc_perim) || cc_perim == 0) NA_real_
                                                       else nrow(bc_near) / cc_perim
                                                     }
    )
  }
  
  out
}


# .calcFractalDimension ----
# Box-counting fractal dimension for a nest point cloud.
.calcFractalDimension <- function(nc, SEED = 42) {
  if (nrow(nc) < 2) return(NA_real_)
  mat        <- as.matrix(nc[, c("X", "Y")])
  data_range <- max(diff(range(mat[, 1])), diff(range(mat[, 2])))
  if (data_range == 0) return(NA_real_)
  k_nn  <- min(1L, nrow(mat) - 1L)
  bmin  <- 0.1 * mean(dbscan::kNN(mat, k = k_nn)$dist) * 2
  bmax  <- data_range / 4
  if (bmax <= bmin) return(NA_real_)
  s_values <- numeric(0); s <- bmax
  while (s >= bmin) { s_values <- c(s_values, s); s <- s / 1.5 }
  if (length(s_values) == 0) return(NA_real_)
  min_x <- min(mat[, 1]); min_y <- min(mat[, 2])
  set.seed(SEED)
  n_s <- vapply(s_values, function(s) {
    counts <- vapply(seq_len(10), function(j) {
      ox <- stats::runif(1, 0, s); oy <- stats::runif(1, 0, s)
      length(unique(paste(floor((mat[, 1] - (min_x - ox)) / s),
                          floor((mat[, 2] - (min_y - oy)) / s), sep = ",")))
    }, numeric(1))
    min(counts)
  }, numeric(1))
  set.seed(NULL)
  valid <- n_s > 1 & is.finite(s_values)
  if (sum(valid) < 3) return(NA_real_)
  df   <- data.frame(log_inv_s = log(1 / s_values[valid]), log_n = log(n_s[valid]))
  fits <- lapply(3:nrow(df), function(k)
    stats::lm(log_n ~ log_inv_s, data = df[seq_len(k), ]))
  r2     <- vapply(fits, function(f) summary(f)$r.squared, numeric(1))
  slopes <- vapply(fits, function(f) stats::coef(f)[2], numeric(1))
  good   <- which(r2 >= 0.99)
  if (length(good) > 0) return(max(slopes[good]))
  slopes[which.max(r2)]
}


# GetShapeMetrics ----
#' @title GetShapeMetrics
#' @description Compute shape and morphology metrics for tumour nests detected
#'   by \code{GetBoundary}. Supports the same three input modes as
#'   \code{GetBoundary} (list, data frame with \code{SAMPLE_COLUMN}, single
#'   data frame).
#'
#' @param INPUT A data frame or named/unnamed list of data frames containing
#'   \code{SynoraAnnotation} output from \code{GetBoundary}.
#' @param X_POSITION Name of the X coordinate column.
#' @param Y_POSITION Name of the Y coordinate column.
#' @param ANNO_COLUMN Annotation column name. Default \code{"SynoraAnnotation"}.
#' @param CELL_ID_COLUMN (optional) Cell ID column name. Row indices used if
#'   omitted.
#' @param CELL_ID_PREFIX (optional) Prefix for auto-generated cell IDs.
#' @param SAMPLE_COLUMN (optional) Column name identifying samples within a
#'   single data frame. Default \code{NULL}.
#' @param ANNO_OF_BOUNDARY Annotation value(s) for boundary cells.
#'   Default \code{"Boundary"}. Auto-detects stratified variants
#'   (\code{"Boundary_Inner"} / \code{"Boundary_Outer"}) when
#'   \code{"Boundary"} is absent.
#' @param ANNO_OF_NEST Annotation value for nest cells.
#'   Default \code{"Nest"}.
#' @param SHAPE_METRICS Character vector of metrics to compute. Global metrics:
#'   \code{"Boundary2NestRatio"}, \code{"NestFragmentation"},
#'   \code{"NestDispersion"}, \code{"NestFractalDimension"}.
#'   Per-nest metrics (require \code{SEPARATE_NESTS = TRUE}):
#'   \code{"NestSolidity"}, \code{"BoundaryRoughness"},
#'   \code{"NestCompactness"}, \code{"NestElongation"},
#'   \code{"NestFractalDimension"}, \code{"BoundaryIntegrity"},
#'   \code{"NestBoundaryBalance"}, \code{"NestEnclosure"},
#'   \code{"NestBoundaryDensity"}.
#'   Default \code{"Boundary2NestRatio"}.
#' @param SEPARATE_NESTS Logical; individualise nests via Delaunay graph
#'   components. Required for per-nest metrics. Default \code{FALSE}.
#' @param DIST_THRESHOLD Maximum edge length in Delaunay triangulation used
#'   for nest separation. Default \code{Inf}.
#' @param MIN_NEST_SIZE Minimum number of cells for a component to be
#'   considered a nest. Default \code{10}.
#' @param WINDOW_SIZE Tile width/height for windowed Delaunay on large
#'   datasets. Default \code{500}.
#' @param WINDOW_OVERLAP Overlap between adjacent tiles. Default \code{50}.
#' @param LARGE_N_THRESHOLD Cell count above which windowed Delaunay is used.
#'   Default \code{10000}.
#' @param CONCAVITY Concavity parameter passed to \code{concaveman}.
#'   Lower = tighter hull. Default \code{1.5}.
#' @param MARGIN_DIST Distance threshold (same units as coordinates) for
#'   classifying hull vertices as touching the tissue margin.
#'   Default \code{10}.
#' @param MARGIN_TOLERANCE Maximum fraction of hull vertices allowed near the
#'   margin for a nest to be classified as closed. Default \code{0.05}.
#' @param BOUNDARY_RADIUS Search radius for associating boundary cells with a
#'   nest (used in closed-nest-only metrics). Default \code{50}.
#' @param SEED Random seed for fractal dimension estimation. Default \code{42}.
#'
#' @return A list (or named list of lists for multi-sample input) with:
#' \describe{
#'   \item{Global}{One-row data frame of sample-level metrics.}
#'   \item{PerNest}{One-row-per-nest data frame. \code{NULL} when
#'     \code{SEPARATE_NESTS = FALSE}.}
#'   \item{HullMap}{Hull polygon vertices with all per-nest metrics joined.
#'     \code{NULL} when \code{SEPARATE_NESTS = FALSE}.}
#'   \item{CellMap}{Cell-to-nest ID mapping data frame.
#'     \code{NULL} when \code{SEPARATE_NESTS = FALSE}.}
#' }
#'
#' @examples
#' library(Synora)
#' data("DummyData")
#'
#' anno_list <- GetBoundary(
#'   INPUT                = DummyData,
#'   CELL_ID_COLUMN       = "Cell_ID",
#'   X_POSITION           = "X",
#'   Y_POSITION           = "Y",
#'   ANNO_COLUMN          = "CT",
#'   RADIUS               = 20,
#'   NEST_SPECIFICITY     = 0.25,
#'   BOUNDARY_SPECIFICITY = 0.05
#' )
#'
#' # --- List mode — global metrics across samples -----------
#' shape_list <- GetShapeMetrics(
#'   INPUT          = anno_list,
#'   X_POSITION     = "X",
#'   Y_POSITION     = "Y",
#'   CELL_ID_COLUMN = "Cell_ID",
#'   SHAPE_METRICS  = c("Boundary2NestRatio", "NestFragmentation")
#' )
#' # Each element has $Global and $PerNest slots
#' shape_list[["Sample1"]]$Global
#'
#' # --- Single data frame — per-nest metrics ----------------
#' shape_single <- GetShapeMetrics(
#'   INPUT          = anno_list[[1]],
#'   X_POSITION     = "X",
#'   Y_POSITION     = "Y",
#'   CELL_ID_COLUMN = "Cell_ID",
#'   SHAPE_METRICS  = c("NestSolidity", "NestCompactness",
#'                       "BoundaryRoughness"),
#'   SEPARATE_NESTS = TRUE
#' )
#' shape_single$Global
#' shape_single$PerNest[, c("Nest_ID", "NestSize", "NestSolidity",
#'                           "NestCompactness")]
#'
#' # --- Closure classification + hull map -------------------
#' shape_closure <- GetShapeMetrics(
#'   INPUT             = anno_list[[1]],
#'   X_POSITION        = "X",
#'   Y_POSITION        = "Y",
#'   CELL_ID_COLUMN    = "Cell_ID",
#'   SHAPE_METRICS     = "NestEnclosure",
#'   SEPARATE_NESTS    = TRUE,
#'   MARGIN_DIST       = 10,
#'   MARGIN_TOLERANCE  = 0.05
#' )
#' table(shape_closure$PerNest$IsClosedNest)
#' head(shape_closure$HullMap)   # polygon vertices ready for ggplot2
#'
#' @export
GetShapeMetrics <- function(INPUT,
                            X_POSITION,
                            Y_POSITION,
                            ANNO_COLUMN        = "SynoraAnnotation",
                            CELL_ID_COLUMN,
                            CELL_ID_PREFIX,
                            SAMPLE_COLUMN      = NULL,
                            ANNO_OF_BOUNDARY   = "Boundary",
                            ANNO_OF_NEST       = "Nest",
                            SHAPE_METRICS      = "Boundary2NestRatio",
                            SEPARATE_NESTS     = FALSE,
                            DIST_THRESHOLD     = Inf,
                            MIN_NEST_SIZE      = 10,
                            WINDOW_SIZE        = 500,
                            WINDOW_OVERLAP     = 50,
                            LARGE_N_THRESHOLD  = 10000,
                            CONCAVITY          = 1.5,
                            MARGIN_DIST        = 10,
                            MARGIN_TOLERANCE   = 0.05,
                            BOUNDARY_RADIUS    = 50,
                            SEED               = 42) {
  
  ARGS        <- as.list(match.call())[-1]
  ARGS$INPUT  <- NULL
  
  .dispatch <- function(df, nm) {
    tryCatch(
      do.call(.GetShapeMetrics_Single, c(list(INPUT = df), ARGS)),
      error = function(e) stop("Error in sample '", nm, "': ", conditionMessage(e))
    )
  }
  
  if (is.list(INPUT) && !is.data.frame(INPUT)) {
    if (length(INPUT) == 0) stop("INPUT list is empty.")
    if (is.null(names(INPUT)))
      INPUT <- purrr::set_names(INPUT, paste0("Sample_", seq_along(INPUT)))
    message("[Synora] List mode — ", length(INPUT), " sample(s).")
    return(purrr::imap(INPUT, .progress = "[Synora] Processing", .f = .dispatch))
  }
  
  if (!is.data.frame(INPUT)) stop("INPUT must be a data frame or list of data frames.")
  
  if (!is.null(SAMPLE_COLUMN)) {
    if (!SAMPLE_COLUMN %in% names(INPUT))
      stop("SAMPLE_COLUMN '", SAMPLE_COLUMN, "' not found.")
    message("[Synora] Data frame mode with SAMPLE_COLUMN — ",
            dplyr::n_distinct(INPUT[[SAMPLE_COLUMN]]), " sample(s).")
    return(INPUT %>% base::split(.[[SAMPLE_COLUMN]]) %>%
             purrr::imap(.progress = "[Synora] Processing", .f = .dispatch))
  }
  
  message("[Synora] Single data frame mode.")
  .dispatch(INPUT, "single")
}


# .GetShapeMetrics_Single ----
# Internal single-sample worker called by GetShapeMetrics.
# Returns list(Global, PerNest, HullMap, CellMap).
.GetShapeMetrics_Single <- function(INPUT,
                                    X_POSITION,
                                    Y_POSITION,
                                    ANNO_COLUMN        = "SynoraAnnotation",
                                    CELL_ID_COLUMN,
                                    CELL_ID_PREFIX,
                                    ANNO_OF_BOUNDARY   = "Boundary",
                                    ANNO_OF_NEST       = "Nest",
                                    SHAPE_METRICS      = "Boundary2NestRatio",
                                    SEPARATE_NESTS     = FALSE,
                                    DIST_THRESHOLD     = Inf,
                                    MIN_NEST_SIZE      = 10,
                                    WINDOW_SIZE        = 500,
                                    WINDOW_OVERLAP     = 50,
                                    LARGE_N_THRESHOLD  = 10000,
                                    CONCAVITY          = 2,
                                    MARGIN_DIST        = 5,
                                    MARGIN_TOLERANCE   = 0.1,
                                    BOUNDARY_RADIUS    = 50,
                                    SEED               = 42) {
  
  # Validation
  if (!is.data.frame(INPUT)) stop("INPUT must be a data frame.")
  for (col in c(X_POSITION, Y_POSITION, ANNO_COLUMN))
    if (!col %in% names(INPUT)) stop("`", col, "` not found in INPUT.")
  
  GLOBAL_METRICS   <- c("Boundary2NestRatio", "NestFragmentation",
                        "NestDispersion", "NestFractalDimension")
  PER_NEST_METRICS <- c("NestSolidity", "BoundaryRoughness", "NestCompactness",
                        "NestElongation", "NestFractalDimension",
                        "BoundaryIntegrity", "NestBoundaryBalance",
                        "NestEnclosure", "NestBoundaryDensity")
  VALID_METRICS    <- unique(c(GLOBAL_METRICS, PER_NEST_METRICS))
  
  bad <- setdiff(SHAPE_METRICS, VALID_METRICS)
  if (length(bad) > 0)
    stop("Invalid SHAPE_METRICS: [", paste(bad, collapse = ", "), "].")
  
  per_nest_requested <- intersect(SHAPE_METRICS, PER_NEST_METRICS)
  global_requested   <- intersect(SHAPE_METRICS, GLOBAL_METRICS)
  
  if (!SEPARATE_NESTS && length(per_nest_requested) > 0)
    stop("Metric(s) [", paste(per_nest_requested, collapse = ", "),
         "] require SEPARATE_NESTS = TRUE.")
  
  # Cell ID
  if (missing(CELL_ID_COLUMN)) {
    pfx    <- if (missing(CELL_ID_PREFIX)) "" else paste0(CELL_ID_PREFIX, "_")
    INPUT  <- INPUT %>% dplyr::mutate(Cell_ID = paste0(pfx, dplyr::row_number()))
    CELL_ID_COLUMN <- "Cell_ID"
  } else {
    if (!CELL_ID_COLUMN %in% names(INPUT))
      stop("CELL_ID_COLUMN '", CELL_ID_COLUMN, "' not found.")
  }
  
  # Annotation resolution
  anno_vals <- unique(INPUT[[ANNO_COLUMN]][!is.na(INPUT[[ANNO_COLUMN]])])
  STRAT     <- c("Boundary_Inner", "Boundary_Outer")
  if (identical(ANNO_OF_BOUNDARY, "Boundary") && !"Boundary" %in% anno_vals) {
    detected <- intersect(STRAT, anno_vals)
    if (length(detected) == 0)
      stop("ANNO_OF_BOUNDARY 'Boundary' not found. Available: ",
           paste(anno_vals, collapse = ", "))
    message("[Synora] Auto-detected boundary variant(s): [",
            paste(detected, collapse = ", "), "].")
    ANNO_OF_BOUNDARY <- detected
  }
  if (!ANNO_OF_NEST %in% anno_vals)
    stop("ANNO_OF_NEST '", ANNO_OF_NEST, "' not found.")
  
  # Normalise to working columns
  INPUT_NORM <- INPUT %>%
    dplyr::transmute(
      Cell_ID = as.character(!!as.name(CELL_ID_COLUMN)),
      X       = !!as.name(X_POSITION),
      Y       = !!as.name(Y_POSITION),
      Anno    = !!as.name(ANNO_COLUMN)
    )
  
  nest_cells     <- INPUT_NORM %>% dplyr::filter(Anno == ANNO_OF_NEST)
  boundary_cells <- INPUT_NORM %>% dplyr::filter(Anno %in% ANNO_OF_BOUNDARY)
  nest_count     <- nrow(nest_cells)
  
  # Global metrics
  global_row <- list()
  
  if ("Boundary2NestRatio" %in% global_requested)
    global_row$Boundary2NestRatio <-
    if (nest_count == 0 || nrow(boundary_cells) == 0) NA_real_
  else nrow(boundary_cells) / nest_count
  
  if ("NestFragmentation" %in% global_requested || SEPARATE_NESTS) {
    g <- .buildNestGraph(INPUT_NORM, "Cell_ID", "Anno", ANNO_OF_NEST,
                         DIST_THRESHOLD, WINDOW_SIZE, WINDOW_OVERLAP,
                         LARGE_N_THRESHOLD)
    nest_nodes <- g %>% tidygraph::activate(nodes) %>%
      dplyr::as_tibble() %>% dplyr::filter(Anno == ANNO_OF_NEST)
    if ("NestFragmentation" %in% global_requested)
      global_row$NestFragmentation <-
      if (nrow(nest_nodes) == 0) NA_real_
    else dplyr::n_distinct(nest_nodes$n_group_id) / nrow(nest_nodes)
  }
  
  if ("NestDispersion" %in% global_requested) {
    mat_all  <- as.matrix(INPUT_NORM[, c("X", "Y")])
    mat_nest <- as.matrix(nest_cells[, c("X", "Y")])
    A <- as.numeric(sf::st_area(
      sf::st_convex_hull(sf::st_multipoint(mat_all))))
    global_row$NestDispersion <-
      if (nest_count < 2 || A == 0) NA_real_
    else {
      mean_nnd <- mean(dbscan::kNN(mat_nest, k = 1)$dist)
      mean_nnd / (1 / (2 * sqrt(nest_count / A)))
    }
  }
  
  if ("NestFractalDimension" %in% global_requested)
    global_row$NestFractalDimension <- .calcFractalDimension(nest_cells, SEED)
  
  Global_df <- as.data.frame(global_row)
  
  # Early return if no per-nest work needed
  if (!SEPARATE_NESTS)
    return(list(Global = Global_df, PerNest = NULL))
  
  if (nest_count < MIN_NEST_SIZE) {
    warning("[Synora] nest_count (", nest_count, ") < MIN_NEST_SIZE. ",
            "Returning empty PerNest.", call. = FALSE)
    return(list(Global = Global_df, PerNest = data.frame()))
  }
  
  # Separate nests
  SeparatedNests <- .SeparateNests(
    INPUT_NORM, "Anno", "Cell_ID", ANNO_OF_NEST,
    DIST_THRESHOLD, MIN_NEST_SIZE, WINDOW_SIZE, WINDOW_OVERLAP,
    LARGE_N_THRESHOLD
  ) %>%
    dplyr::filter(!is.na(Nest_ID)) %>%
    base::split(.$Nest_ID)
  
  if (length(SeparatedNests) == 0) {
    warning("[Synora] No nests passed MIN_NEST_SIZE.", call. = FALSE)
    return(list(Global = Global_df, PerNest = data.frame()))
  }
  
  # Tissue margin line (outer contour of entire point cloud)
  all_mat     <- as.matrix(INPUT_NORM[, c("X", "Y")])
  margin_ring <- concaveman::concaveman(all_mat, concavity = CONCAVITY)
  margin_line <- sf::st_linestring(margin_ring)
  
  # Per-nest loop
  per_nest_rows <- list()
  hull_map_list <- list()
  
  for (nest_nm in names(SeparatedNests)) {
    nest_i <- SeparatedNests[[nest_nm]]
    nc_i   <- nest_i %>% dplyr::filter(Anno == ANNO_OF_NEST)
    
    hull_ring <- tryCatch(
      .buildNestHull(
        nc_ids         = nc_i$Cell_ID,
        all_pts        = INPUT_NORM,
        anno_nest      = ANNO_OF_NEST,
        boundary_cells = boundary_cells,
        dist_threshold = DIST_THRESHOLD,
        concavity      = CONCAVITY
      ),
      error = function(e) {
        warning("[Synora] Hull failed for ", nest_nm, ": ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    
    if (!is.null(hull_ring) && nrow(hull_ring) >= 3)
      hull_map_list[[nest_nm]] <- data.frame(
        Nest_ID = nest_nm,
        X       = hull_ring[, 1],
        Y       = hull_ring[, 2]
      )
    
    closure   <- .classifyNestClosure(hull_ring, margin_line,
                                      MARGIN_DIST, MARGIN_TOLERANCE)
    is_closed <- closure$IsClosedNest
    
    bc_near <- if (is_closed && nrow(boundary_cells) > 0) {
      mat_nc    <- as.matrix(nc_i[, c("X", "Y")])
      mat_bc    <- as.matrix(boundary_cells[, c("X", "Y")])
      min_dists <- FNN::knnx.dist(mat_nc, mat_bc, k = 1L)[, 1]
      boundary_cells[min_dists <= BOUNDARY_RADIUS, , drop = FALSE]
    } else {
      boundary_cells[0, , drop = FALSE]
    }
    
    metrics <- .calcPerNestMetrics(
      nc_i, hull_ring, is_closed, bc_near, margin_line,
      per_nest_requested, BOUNDARY_RADIUS, SEED
    )
    
    row <- data.frame(
      Nest_ID        = nest_nm,
      NestSize       = nrow(nc_i),
      IsClosedNest   = is_closed,
      ClosedFraction = closure$ClosedFraction,
      stringsAsFactors = FALSE
    )
    for (m in per_nest_requested)
      row[[m]] <- if (!is.null(metrics[[m]])) metrics[[m]] else NA_real_
    
    per_nest_rows[[nest_nm]] <- row
  }
  
  PerNest_df <- dplyr::bind_rows(per_nest_rows)
  
  # HullMap: hull vertices with per-nest metrics joined for easy ggplot2 use
  HullMap_df <- dplyr::bind_rows(hull_map_list) %>%
    dplyr::left_join(
      PerNest_df %>% dplyr::select(-NestSize),
      by = "Nest_ID"
    )
  
  # CellMap: cell-to-nest ID mapping
  CellMap_df <- dplyr::bind_rows(SeparatedNests) %>%
    dplyr::filter(Anno == ANNO_OF_NEST) %>%
    dplyr::select(Cell_ID, Nest_ID)
  
  n_closed <- sum(PerNest_df$IsClosedNest, na.rm = TRUE)
  message("[Synora] ", nrow(PerNest_df), " nest(s); ",
          n_closed, " closed (MARGIN_TOLERANCE=", MARGIN_TOLERANCE * 100,
          "%, MARGIN_DIST=", MARGIN_DIST, ").")
  
  list(
    Global  = Global_df,
    PerNest = PerNest_df,
    HullMap = HullMap_df,
    CellMap = CellMap_df
  )
}