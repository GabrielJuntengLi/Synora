# PlotShapeMetrics ----
#' @title PlotShapeMetrics
#' @description Multi-panel diagnostic plot for \code{GetShapeMetrics} output.
#'   Produces up to three panel types, controlled by \code{PLOT_TYPE}:
#'   \describe{
#'     \item{\code{"spatial"}}{Spatial scatter of all cells coloured by
#'       \code{SynoraAnnotation}, with optional per-nest hull overlays.}
#'     \item{\code{"per_nest"}}{Jitter + mean±SD plot of numeric per-nest
#'       metrics, faceted by metric and split by closure status.}
#'     \item{\code{"summary"}}{Bar chart of mean±SD per metric comparing
#'       all nests vs closed nests only.}
#'     \item{\code{"all"}}{All three panels stacked via \pkg{patchwork}.}
#'   }
#'   Supports multi-sample input: when both \code{SPATIAL_INPUT} and
#'   \code{METRICS_RESULT} are named lists, returns a named list of plots.
#'
#' @param SPATIAL_INPUT A data frame (or named list of data frames) containing
#'   cell coordinates and \code{SynoraAnnotation} output from
#'   \code{GetBoundary}.
#' @param METRICS_RESULT Output from \code{GetShapeMetrics}. A single result
#'   list or a named list of results matching \code{SPATIAL_INPUT}.
#' @param X_POSITION Name of the X coordinate column. Default \code{"X"}.
#' @param Y_POSITION Name of the Y coordinate column. Default \code{"Y"}.
#' @param ANNO_COLUMN Annotation column name. Default \code{"SynoraAnnotation"}.
#' @param CELL_ID_COLUMN Cell ID column name. Default \code{"Cell_ID"}.
#' @param PLOT_TYPE Character vector of panel types to produce. Any combination
#'   of \code{"spatial"}, \code{"per_nest"}, \code{"summary"}, or
#'   \code{"all"}. Default \code{"all"}.
#' @param HIGHLIGHT_CLOSED Logical; colour closed-nest hulls differently from
#'   open-nest hulls in the spatial panel. Default \code{TRUE}.
#' @param HULL_TYPE Hull geometry for the spatial panel: \code{"concave"}
#'   (uses \pkg{concaveman}) or \code{"convex"}. Default \code{"concave"}.
#' @param CONCAVITY Concavity parameter passed to \code{concaveman}.
#'   Default \code{2}.
#' @param COLOR_NEST Point colour for Nest cells. Default \code{"#377EB8FF"}.
#' @param COLOR_BOUNDARY Point colour for Boundary cells.
#'   Default \code{"#4DAF4AFF"}.
#' @param COLOR_OUTSIDE Point colour for Outside cells.
#'   Default \code{"#984EA3FF"}.
#' @param COLOR_CLOSED_HULL Hull outline colour for closed nests.
#'   Default \code{"#E63946"}.
#' @param COLOR_OPEN_HULL Hull outline colour for open nests.
#'   Default \code{"#ADB5BD"}.
#' @param COLOR_CLOSED_BAR Bar colour for the "Closed Nests" group in the
#'   summary panel. Default \code{"#E63946"}.
#' @param COLOR_ALL_BAR Bar colour for the "All Nests" group in the summary
#'   panel. Default \code{"#457B9D"}.
#' @param POINT_SIZE Point size for the spatial scatter. Default \code{0.8}.
#' @param HULL_LINEWIDTH Linewidth for hull outlines. Default \code{0.6}.
#' @param SAMPLE_NAME (optional) String used as the plot title.
#'   Default \code{NULL}.
#'
#' @return A \code{ggplot} or \pkg{patchwork} object (single sample), or a
#'   named list of such objects (multi-sample input).
#'
#' @examples
#' library(Synora)
#' data("DummyData")
#'
#' anno <- GetBoundary(
#'   INPUT = DummyData[[1]], CELL_ID_COLUMN = "Cell_ID",
#'   X_POSITION = "X", Y_POSITION = "Y", ANNO_COLUMN = "CT",
#'   RADIUS = 20, NEST_SPECIFICITY = 0.25, BOUNDARY_SPECIFICITY = 0.05
#' )
#' shapes <- GetShapeMetrics(
#'   INPUT = anno, X_POSITION = "X", Y_POSITION = "Y",
#'   CELL_ID_COLUMN = "Cell_ID", SEPARATE_NESTS = TRUE,
#'   SHAPE_METRICS = c("NestSolidity", "NestCompactness")
#' )
#'
#' # --- All panels (spatial + per-nest + summary) -----------
#' PlotShapeMetrics(SPATIAL_INPUT = anno, METRICS_RESULT = shapes,
#'                  PLOT_TYPE = "all")
#'
#' # --- Spatial panel only ----------------------------------
#' PlotShapeMetrics(SPATIAL_INPUT = anno, METRICS_RESULT = shapes,
#'                  PLOT_TYPE = "spatial", HIGHLIGHT_CLOSED = TRUE,
#'                  HULL_TYPE = "concave", CONCAVITY = 2)
#'
#' # --- Per-nest scatter panel only -------------------------
#' PlotShapeMetrics(SPATIAL_INPUT = anno, METRICS_RESULT = shapes,
#'                  PLOT_TYPE = "per_nest")
#'
#' # --- Custom colour scheme --------------------------------
#' PlotShapeMetrics(
#'   SPATIAL_INPUT     = anno,    METRICS_RESULT = shapes,
#'   PLOT_TYPE         = "all",
#'   COLOR_NEST        = "#1D3557", COLOR_BOUNDARY    = "#457B9D",
#'   COLOR_OUTSIDE     = "#A8DADC", COLOR_CLOSED_HULL = "#E63946",
#'   COLOR_OPEN_HULL   = "#F1FAEE", POINT_SIZE        = 1.0,
#'   HULL_LINEWIDTH    = 0.8,       SAMPLE_NAME       = "Sample 1"
#' )
#'
#' # --- Multi-sample list input -----------------------------
#' anno_list   <- GetBoundary(
#'   INPUT = DummyData, CELL_ID_COLUMN = "Cell_ID",
#'   X_POSITION = "X", Y_POSITION = "Y", ANNO_COLUMN = "CT",
#'   RADIUS = 20, NEST_SPECIFICITY = 0.25, BOUNDARY_SPECIFICITY = 0.05
#' )
#' shapes_list <- GetShapeMetrics(
#'   INPUT = anno_list, X_POSITION = "X", Y_POSITION = "Y",
#'   CELL_ID_COLUMN = "Cell_ID", SEPARATE_NESTS = TRUE,
#'   SHAPE_METRICS = c("NestSolidity", "NestCompactness")
#' )
#' plot_list <- PlotShapeMetrics(
#'   SPATIAL_INPUT  = anno_list, METRICS_RESULT = shapes_list,
#'   PLOT_TYPE      = "summary"
#' )
#' plot_list[["Sample1"]]
#'
#' @export
PlotShapeMetrics <- function(SPATIAL_INPUT,
                             METRICS_RESULT,
                             X_POSITION        = "X",
                             Y_POSITION        = "Y",
                             ANNO_COLUMN       = "SynoraAnnotation",
                             CELL_ID_COLUMN    = "Cell_ID",
                             PLOT_TYPE         = "all",
                             HIGHLIGHT_CLOSED  = TRUE,
                             HULL_TYPE         = "concave",
                             CONCAVITY         = 2,
                             COLOR_NEST        = "#377EB8FF",
                             COLOR_BOUNDARY    = "#4DAF4AFF",
                             COLOR_OUTSIDE     = "#984EA3FF",
                             COLOR_CLOSED_HULL = "#E63946",
                             COLOR_OPEN_HULL   = "#ADB5BD",
                             COLOR_CLOSED_BAR  = "#E63946",
                             COLOR_ALL_BAR     = "#457B9D",
                             POINT_SIZE        = 0.8,
                             HULL_LINEWIDTH    = 0.6,
                             SAMPLE_NAME       = NULL) {
  
  VALID_TYPES <- c("spatial", "per_nest", "summary", "all")
  bad_types   <- setdiff(PLOT_TYPE, VALID_TYPES)
  if (length(bad_types) > 0)
    stop("Invalid PLOT_TYPE: [", paste(bad_types, collapse = ", "),
         "]. Valid: spatial, per_nest, summary, all.")
  if ("all" %in% PLOT_TYPE) PLOT_TYPE <- c("spatial", "per_nest", "summary")
  PLOT_TYPE <- unique(PLOT_TYPE)
  
  # Multi-sample dispatch
  if (is.list(SPATIAL_INPUT) && !is.data.frame(SPATIAL_INPUT)) {
    if (!is.list(METRICS_RESULT) || is.null(names(METRICS_RESULT)))
      stop("When SPATIAL_INPUT is a list, METRICS_RESULT must be a named list.")
    shared <- intersect(names(SPATIAL_INPUT), names(METRICS_RESULT))
    if (length(shared) == 0)
      stop("No matching names between SPATIAL_INPUT and METRICS_RESULT.")
    if (length(shared) < length(names(SPATIAL_INPUT)))
      warning("Skipping unmatched samples: [",
              paste(setdiff(names(SPATIAL_INPUT), shared), collapse = ", "),
              "].", call. = FALSE)
    return(purrr::imap(
      purrr::set_names(shared),
      function(nm, ...) .PlotShapeMetrics_Single(
        SPATIAL_INPUT     = SPATIAL_INPUT[[nm]],
        METRICS_RESULT    = METRICS_RESULT[[nm]],
        X_POSITION        = X_POSITION,      Y_POSITION        = Y_POSITION,
        ANNO_COLUMN       = ANNO_COLUMN,      CELL_ID_COLUMN    = CELL_ID_COLUMN,
        PLOT_TYPE         = PLOT_TYPE,        HIGHLIGHT_CLOSED  = HIGHLIGHT_CLOSED,
        HULL_TYPE         = HULL_TYPE,        CONCAVITY         = CONCAVITY,
        COLOR_NEST        = COLOR_NEST,       COLOR_BOUNDARY    = COLOR_BOUNDARY,
        COLOR_OUTSIDE     = COLOR_OUTSIDE,    COLOR_CLOSED_HULL = COLOR_CLOSED_HULL,
        COLOR_OPEN_HULL   = COLOR_OPEN_HULL,  COLOR_CLOSED_BAR  = COLOR_CLOSED_BAR,
        COLOR_ALL_BAR     = COLOR_ALL_BAR,    POINT_SIZE        = POINT_SIZE,
        HULL_LINEWIDTH    = HULL_LINEWIDTH,   SAMPLE_NAME       = nm
      )
    ))
  }
  
  .PlotShapeMetrics_Single(
    SPATIAL_INPUT     = SPATIAL_INPUT,    METRICS_RESULT    = METRICS_RESULT,
    X_POSITION        = X_POSITION,       Y_POSITION        = Y_POSITION,
    ANNO_COLUMN       = ANNO_COLUMN,       CELL_ID_COLUMN    = CELL_ID_COLUMN,
    PLOT_TYPE         = PLOT_TYPE,         HIGHLIGHT_CLOSED  = HIGHLIGHT_CLOSED,
    HULL_TYPE         = HULL_TYPE,         CONCAVITY         = CONCAVITY,
    COLOR_NEST        = COLOR_NEST,        COLOR_BOUNDARY    = COLOR_BOUNDARY,
    COLOR_OUTSIDE     = COLOR_OUTSIDE,     COLOR_CLOSED_HULL = COLOR_CLOSED_HULL,
    COLOR_OPEN_HULL   = COLOR_OPEN_HULL,   COLOR_CLOSED_BAR  = COLOR_CLOSED_BAR,
    COLOR_ALL_BAR     = COLOR_ALL_BAR,     POINT_SIZE        = POINT_SIZE,
    HULL_LINEWIDTH    = HULL_LINEWIDTH,    SAMPLE_NAME       = SAMPLE_NAME
  )
}


# .PlotShapeMetrics_Single ----
#' @keywords internal
.PlotShapeMetrics_Single <- function(SPATIAL_INPUT,
                                     METRICS_RESULT,
                                     X_POSITION,
                                     Y_POSITION,
                                     ANNO_COLUMN,
                                     CELL_ID_COLUMN,
                                     PLOT_TYPE,
                                     HIGHLIGHT_CLOSED,
                                     HULL_TYPE,
                                     CONCAVITY,
                                     COLOR_NEST,
                                     COLOR_BOUNDARY,
                                     COLOR_OUTSIDE,
                                     COLOR_CLOSED_HULL,
                                     COLOR_OPEN_HULL,
                                     COLOR_CLOSED_BAR,
                                     COLOR_ALL_BAR,
                                     POINT_SIZE,
                                     HULL_LINEWIDTH,
                                     SAMPLE_NAME) {
  
  # Guards
  has_per_nest <- !is.null(METRICS_RESULT$PerNest) &&
    is.data.frame(METRICS_RESULT$PerNest) &&
    nrow(METRICS_RESULT$PerNest) > 0
  has_global   <- !is.null(METRICS_RESULT$Global) &&
    is.data.frame(METRICS_RESULT$Global) &&
    ncol(METRICS_RESULT$Global) > 0
  has_cell_map <- !is.null(METRICS_RESULT$CellMap) &&
    is.data.frame(METRICS_RESULT$CellMap) &&
    nrow(METRICS_RESULT$CellMap) > 0
  
  FIXED_COLS      <- c("Nest_ID", "NestSize", "IsClosedNest", "ClosedFraction")
  numeric_metrics <- if (has_per_nest) {
    candidate <- setdiff(names(METRICS_RESULT$PerNest), FIXED_COLS)
    candidate[vapply(candidate, function(m)
      is.numeric(METRICS_RESULT$PerNest[[m]]), logical(1))]
  } else character(0)
  
  # Annotation colour map
  anno_vals <- unique(as.character(SPATIAL_INPUT[[ANNO_COLUMN]]))
  anno_vals <- anno_vals[!is.na(anno_vals)]
  color_map <- c()
  if ("Boundary"       %in% anno_vals) color_map["Boundary"]       <- COLOR_BOUNDARY
  if ("Boundary_Inner" %in% anno_vals) color_map["Boundary_Inner"] <- COLOR_BOUNDARY
  if ("Boundary_Outer" %in% anno_vals) color_map["Boundary_Outer"] <- scales::alpha(COLOR_BOUNDARY, 0.6)
  if ("Nest"           %in% anno_vals) color_map["Nest"]           <- COLOR_NEST
  if ("Outside"        %in% anno_vals) color_map["Outside"]        <- COLOR_OUTSIDE
  if ("Noise"          %in% anno_vals) color_map["Noise"]          <- "#CCCCCC"
  
  plots <- list()
  
  # Spatial panel
  if ("spatial" %in% PLOT_TYPE) {
    
    p_spatial <- ggplot2::ggplot(
      SPATIAL_INPUT,
      ggplot2::aes(x     = !!as.name(X_POSITION),
                   y     = !!as.name(Y_POSITION),
                   color = as.character(!!as.name(ANNO_COLUMN)))
    ) +
      ggplot2::geom_point(size = POINT_SIZE, na.rm = TRUE) +
      ggplot2::scale_color_manual(values   = color_map,
                                  name     = "Annotation",
                                  na.value = "#E0E0E0") +
      ggplot2::coord_equal() +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.title  = ggplot2::element_text(size = 9, face = "bold"),
        legend.text   = ggplot2::element_text(size = 8),
        plot.title    = ggplot2::element_text(size = 10, face = "bold",
                                              margin = ggplot2::margin(b = 4)),
        plot.subtitle = ggplot2::element_text(size = 8, color = "grey40",
                                              margin = ggplot2::margin(b = 4))
      )
    
    # Per-nest hull overlays
    if (has_per_nest && has_cell_map && HULL_TYPE != "none") {
      
      cell_map              <- METRICS_RESULT$CellMap
      cell_map[["Cell_ID"]] <- as.character(cell_map[["Cell_ID"]])
      if (CELL_ID_COLUMN != "Cell_ID")
        names(cell_map)[names(cell_map) == "Cell_ID"] <- CELL_ID_COLUMN
      
      spatial_mapped <- SPATIAL_INPUT %>%
        dplyr::mutate(!!CELL_ID_COLUMN := as.character(!!as.name(CELL_ID_COLUMN))) %>%
        dplyr::left_join(cell_map, by = CELL_ID_COLUMN)
      
      closure_lut <- stats::setNames(METRICS_RESULT$PerNest$IsClosedNest,
                                     METRICS_RESULT$PerNest$Nest_ID)
      
      hull_df_list <- lapply(METRICS_RESULT$PerNest$Nest_ID, function(nest_nm) {
        nc_i <- spatial_mapped[!is.na(spatial_mapped$Nest_ID) &
                                 spatial_mapped$Nest_ID == nest_nm, ]
        if (nrow(nc_i) < 3) return(NULL)
        mat  <- as.matrix(nc_i[, c(X_POSITION, Y_POSITION)])
        ring <- if (HULL_TYPE == "concave")
          tryCatch(concaveman::concaveman(mat, concavity = CONCAVITY),
                   error = function(e) {
                     raw <- sf::st_coordinates(sf::st_convex_hull(sf::st_multipoint(mat)))
                     raw[, c("X", "Y"), drop = FALSE]
                   })
        else {
          raw <- sf::st_coordinates(sf::st_convex_hull(sf::st_multipoint(mat)))
          raw[, c("X", "Y"), drop = FALSE]
        }
        data.frame(X = ring[, 1], Y = ring[, 2],
                   Nest_ID  = nest_nm,
                   IsClosed = isTRUE(closure_lut[[nest_nm]]))
      })
      
      hull_df <- dplyr::bind_rows(
        hull_df_list[!vapply(hull_df_list, is.null, logical(1))])
      
      if (nrow(hull_df) > 0) {
        # Vectorised outline: two geom_polygon calls (closed / open)
        closed_h <- hull_df %>% dplyr::filter( IsClosed)
        open_h   <- hull_df %>% dplyr::filter(!IsClosed)
        hull_color <- if (HIGHLIGHT_CLOSED)
          list(closed = COLOR_CLOSED_HULL, open = COLOR_OPEN_HULL)
        else
          list(closed = COLOR_OPEN_HULL,   open = COLOR_OPEN_HULL)
        
        if (nrow(closed_h) > 0)
          p_spatial <- p_spatial + ggplot2::geom_polygon(
            data    = closed_h,
            mapping = ggplot2::aes(x = X, y = Y, group = Nest_ID),
            fill = NA, color = hull_color$closed,
            linewidth = HULL_LINEWIDTH, inherit.aes = FALSE
          )
        if (nrow(open_h) > 0)
          p_spatial <- p_spatial + ggplot2::geom_polygon(
            data    = open_h,
            mapping = ggplot2::aes(x = X, y = Y, group = Nest_ID),
            fill = NA, color = hull_color$open,
            linewidth = HULL_LINEWIDTH, inherit.aes = FALSE
          )
        
        if (HIGHLIGHT_CLOSED) {
          n_closed  <- sum(METRICS_RESULT$PerNest$IsClosedNest, na.rm = TRUE)
          n_open    <- nrow(METRICS_RESULT$PerNest) - n_closed
          p_spatial <- p_spatial + ggplot2::labs(
            subtitle = paste0(n_closed, " closed nest(s)  \u2014  ",
                              n_open,   " open nest(s)"))
        }
      }
      
    } else if (has_per_nest && !has_cell_map && HULL_TYPE != "none") {
      warning("[Synora] METRICS_RESULT$CellMap is missing. ",
              "Re-run GetShapeMetrics with SEPARATE_NESTS = TRUE.",
              call. = FALSE)
    }
    
    # Global metrics as caption
    if (has_global) {
      global_str <- paste(
        vapply(names(METRICS_RESULT$Global), function(nm) {
          v <- METRICS_RESULT$Global[[nm]]
          sprintf("%s: %s", nm,
                  if (is.numeric(v) && !is.na(v)) round(v, 4) else "NA")
        }, character(1)),
        collapse = "\n"
      )
      p_spatial <- p_spatial +
        ggplot2::labs(caption = global_str) +
        ggplot2::theme(
          plot.caption = ggplot2::element_text(size = 7, hjust = 0,
                                               color = "grey40",
                                               margin = ggplot2::margin(t = 4))
        )
    }
    
    if (!is.null(SAMPLE_NAME))
      p_spatial <- p_spatial + ggplot2::labs(title = SAMPLE_NAME)
    
    plots$spatial <- p_spatial
  }
  
  # Per-nest panel
  if ("per_nest" %in% PLOT_TYPE && has_per_nest &&
      length(numeric_metrics) > 0) {
    
    per_nest_long <- METRICS_RESULT$PerNest %>%
      dplyr::select(dplyr::all_of(c("Nest_ID", "IsClosedNest",
                                    numeric_metrics))) %>%
      tidyr::pivot_longer(cols      = dplyr::all_of(numeric_metrics),
                          names_to  = "Metric",
                          values_to = "Value") %>%
      dplyr::mutate(
        NestType = ifelse(IsClosedNest, "Closed", "Open"),
        Metric   = factor(Metric, levels = numeric_metrics)
      )
    
    plots$per_nest <- ggplot2::ggplot(
      per_nest_long,
      ggplot2::aes(x = NestType, y = Value, color = NestType)
    ) +
      ggplot2::geom_jitter(width = 0.15, size = 2, alpha = 0.8, na.rm = TRUE) +
      ggplot2::stat_summary(
        fun     = mean,
        fun.min = function(x) mean(x) - stats::sd(x),
        fun.max = function(x) mean(x) + stats::sd(x),
        geom = "pointrange", color = "black", size = 0.5, na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(
        values = c(Closed = COLOR_CLOSED_HULL, Open = COLOR_OPEN_HULL),
        guide  = "none"
      ) +
      ggplot2::facet_wrap(ggplot2::vars(Metric), scales = "free_y",
                          ncol = min(3L, length(numeric_metrics))) +
      ggplot2::labs(
        x     = NULL, y = "Value",
        title = if (!is.null(SAMPLE_NAME))
          paste0(SAMPLE_NAME, " \u2014 Per-nest metrics") else "Per-nest metrics"
      ) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(strip.text       = ggplot2::element_text(face = "bold", size = 9),
                     axis.text.x      = ggplot2::element_text(size = 9),
                     panel.grid.minor = ggplot2::element_blank())
  }
  
  # Summary panel
  if ("summary" %in% PLOT_TYPE && has_per_nest &&
      length(numeric_metrics) > 0) {
    
    .summarise_group <- function(df) {
      df %>%
        dplyr::select(dplyr::all_of(numeric_metrics)) %>%
        tidyr::pivot_longer(dplyr::everything(),
                            names_to = "Metric", values_to = "Value") %>%
        dplyr::group_by(Metric) %>%
        dplyr::summarise(Mean = mean(Value, na.rm = TRUE),
                         SD   = stats::sd(Value, na.rm = TRUE),
                         N    = sum(!is.na(Value)),
                         .groups = "drop")
    }
    
    summary_df <- dplyr::bind_rows(
      .summarise_group(METRICS_RESULT$PerNest) %>%
        dplyr::mutate(Group = "All Nests"),
      .summarise_group(METRICS_RESULT$PerNest %>%
                         dplyr::filter(IsClosedNest)) %>%
        dplyr::mutate(Group = "Closed Nests")
    ) %>%
      dplyr::mutate(
        SD     = dplyr::coalesce(SD, 0),
        Group  = factor(Group, levels = c("All Nests", "Closed Nests")),
        Metric = factor(Metric, levels = numeric_metrics)
      ) %>%
      dplyr::filter(!is.nan(Mean))
    
    if (nrow(summary_df) > 0)
      plots$summary <- ggplot2::ggplot(
        summary_df,
        ggplot2::aes(x = Group, y = Mean, fill = Group)
      ) +
      ggplot2::geom_col(width = 0.55, alpha = 0.9) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = pmax(Mean - SD, 0), ymax = Mean + SD),
        width = 0.2, linewidth = 0.5, color = "grey30"
      ) +
      ggplot2::scale_fill_manual(
        values = c("All Nests"    = COLOR_ALL_BAR,
                   "Closed Nests" = COLOR_CLOSED_BAR),
        guide  = "none"
      ) +
      ggplot2::facet_wrap(ggplot2::vars(Metric), scales = "free_y",
                          ncol = min(3L, dplyr::n_distinct(summary_df$Metric))) +
      ggplot2::labs(
        x     = NULL, y = "Mean \u00b1 SD",
        title = if (!is.null(SAMPLE_NAME))
          paste0(SAMPLE_NAME, " \u2014 Summary") else "Summary"
      ) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(strip.text       = ggplot2::element_text(face = "bold", size = 9),
                     axis.text.x      = ggplot2::element_text(size = 9),
                     panel.grid.minor = ggplot2::element_blank())
  }
  
  # Assemble
  if (length(plots) == 0) {
    warning("[Synora] No plots produced. Check METRICS_RESULT contains PerNest data.",
            call. = FALSE)
    return(invisible(NULL))
  }
  if (length(plots) == 1) return(plots[[1]])
  
  patchwork::wrap_plots(plots, ncol = 1) +
    patchwork::plot_annotation(
      title = SAMPLE_NAME,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 12, face = "bold"))
    )
}