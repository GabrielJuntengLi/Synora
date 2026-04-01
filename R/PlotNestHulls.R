# PlotNestHulls ----
#' @title PlotNestHulls
#' @description Visualise per-nest hull polygons from \code{GetShapeMetrics}
#'   output. Requires \code{SEPARATE_NESTS = TRUE} in \code{GetShapeMetrics}.
#'
#'   Supports multi-sample input: when \code{METRICS_RESULT} is a named list
#'   (from list-mode \code{GetShapeMetrics}), returns a named list of ggplots.
#'
#' @param METRICS_RESULT Output from \code{GetShapeMetrics}. A single result
#'   list (with \code{$HullMap}) or a named list of such results.
#' @param FILL_MODE Fill style for hull polygons. One of:
#'   \describe{
#'     \item{\code{"identity"}}{Each nest receives a distinct categorical colour.}
#'     \item{\code{"metric"}}{Continuous heatmap fill by a numeric metric
#'       (set via \code{FILL_METRIC}).}
#'     \item{\code{"closure"}}{Two-colour fill: closed nests vs open nests.}
#'     \item{\code{"none"}}{No fill; outline only.}
#'   }
#'   Default \code{"identity"}.
#' @param FILL_METRIC Name of the numeric column in \code{HullMap} to use as
#'   the heatmap fill variable. Required when \code{FILL_MODE = "metric"}.
#'   Default \code{NULL}.
#' @param FILL_PALETTE Colour palette name. \code{"viridis"} uses
#'   \code{scale_fill_viridis_c}; any other string is passed to
#'   \code{scale_fill_distiller}. Default \code{"viridis"}.
#' @param FILL_DIRECTION Palette direction (\code{1} or \code{-1}).
#'   Default \code{1}.
#' @param FILL_ALPHA Polygon fill transparency (0–1). Default \code{0.55}.
#' @param CONTOUR_CLOSED Logical; draw a thicker/coloured border around closed
#'   nests and a thinner grey border around open nests. Default \code{TRUE}.
#' @param COLOR_CLOSED Border (and closure-fill) colour for closed nests.
#'   Default \code{"#E63946"}.
#' @param COLOR_OPEN Border (and closure-fill) colour for open nests.
#'   Default \code{"#ADB5BD"}.
#' @param LINEWIDTH_CLOSED Border linewidth for closed nests. Default \code{1.2}.
#' @param LINEWIDTH_OPEN Border linewidth for open nests. Default \code{0.5}.
#' @param BACKGROUND_POINTS (optional) A data frame (or named list of data
#'   frames for multi-sample mode) of cells to plot as a background scatter
#'   layer. Default \code{NULL}.
#' @param BG_ANNO_COLUMN Annotation column in \code{BACKGROUND_POINTS} used
#'   for point colouring. Default \code{"SynoraAnnotation"}.
#' @param BG_X X coordinate column in \code{BACKGROUND_POINTS}.
#'   Default \code{"X"}.
#' @param BG_Y Y coordinate column in \code{BACKGROUND_POINTS}.
#'   Default \code{"Y"}.
#' @param BG_COLOR_MAP Named character vector mapping annotation values to
#'   colours for the background scatter layer.
#' @param BG_POINT_SIZE Point size for background cells. Default \code{0.5}.
#' @param BG_POINT_ALPHA Point transparency for background cells.
#'   Default \code{0.4}.
#' @param LABEL_NESTS Logical; overlay \code{Nest_ID} text labels at hull
#'   centroids. Default \code{FALSE}.
#' @param LABEL_SIZE Text size for nest labels. Default \code{2.5}.
#' @param SAMPLE_NAME (optional) String used as the plot title.
#'   Default \code{NULL}.
#'
#' @return A \code{ggplot} object, or a named list of \code{ggplot} objects
#'   when \code{METRICS_RESULT} is a named list.
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
#'   CELL_ID_COLUMN = "Cell_ID", SEPARATE_NESTS = TRUE
#' )
#'
#' # --- Categorical fill (one colour per nest) ---------------
#' PlotNestHulls(METRICS_RESULT = shapes, FILL_MODE = "identity",
#'               FILL_ALPHA = 0.5)
#'
#' # --- Heatmap fill by a numeric metric --------------------
#' PlotNestHulls(METRICS_RESULT = shapes, FILL_MODE = "metric",
#'               FILL_METRIC = "NestSolidity", FILL_PALETTE = "viridis",
#'               FILL_ALPHA = 0.6)
#'
#' # --- Closure colouring (closed vs open) ------------------
#' PlotNestHulls(METRICS_RESULT = shapes, FILL_MODE = "closure",
#'               COLOR_CLOSED = "#E63946", COLOR_OPEN = "#ADB5BD",
#'               CONTOUR_CLOSED = TRUE, LINEWIDTH_CLOSED = 1.4)
#'
#' # --- Background cell layer + nest ID labels --------------
#' PlotNestHulls(
#'   METRICS_RESULT    = shapes,    FILL_MODE      = "identity",
#'   BACKGROUND_POINTS = anno,      BG_ANNO_COLUMN = "SynoraAnnotation",
#'   BG_X = "X",                    BG_Y           = "Y",
#'   BG_POINT_SIZE     = 0.4,       BG_POINT_ALPHA = 0.35,
#'   LABEL_NESTS       = TRUE,      LABEL_SIZE     = 2.5
#' )
#'
#' # --- Multi-sample list → named list of ggplots -----------
#' anno_list   <- GetBoundary(
#'   INPUT = DummyData, CELL_ID_COLUMN = "Cell_ID",
#'   X_POSITION = "X", Y_POSITION = "Y", ANNO_COLUMN = "CT",
#'   RADIUS = 20, NEST_SPECIFICITY = 0.25, BOUNDARY_SPECIFICITY = 0.05
#' )
#' shapes_list <- GetShapeMetrics(
#'   INPUT = anno_list, X_POSITION = "X", Y_POSITION = "Y",
#'   CELL_ID_COLUMN = "Cell_ID", SEPARATE_NESTS = TRUE
#' )
#' hull_plots <- PlotNestHulls(
#'   METRICS_RESULT    = shapes_list, FILL_MODE = "closure",
#'   BACKGROUND_POINTS = anno_list    # named list matched by sample name
#' )
#' hull_plots[["Sample1"]]
#'
#' @export
PlotNestHulls <- function(METRICS_RESULT,
                          FILL_MODE          = "identity",
                          FILL_METRIC        = NULL,
                          FILL_PALETTE       = "viridis",
                          FILL_DIRECTION     = 1,
                          FILL_ALPHA         = 0.55,
                          CONTOUR_CLOSED     = TRUE,
                          COLOR_CLOSED       = "#E63946",
                          COLOR_OPEN         = "#ADB5BD",
                          LINEWIDTH_CLOSED   = 1.2,
                          LINEWIDTH_OPEN     = 0.5,
                          BACKGROUND_POINTS  = NULL,
                          BG_ANNO_COLUMN     = "SynoraAnnotation",
                          BG_X               = "X",
                          BG_Y               = "Y",
                          BG_COLOR_MAP       = c(Nest       = "#BFDBFE",
                                                 Boundary   = "#BBF7D0",
                                                 Outside    = "#E9D5FF",
                                                 Noise      = "#E5E7EB"),
                          BG_POINT_SIZE      = 0.5,
                          BG_POINT_ALPHA     = 0.4,
                          LABEL_NESTS        = FALSE,
                          LABEL_SIZE         = 2.5,
                          SAMPLE_NAME        = NULL) {
  
  # Multi-sample dispatch
  if (is.list(METRICS_RESULT) && is.null(METRICS_RESULT$HullMap)) {
    bg_list <- if (is.list(BACKGROUND_POINTS) && !is.data.frame(BACKGROUND_POINTS))
      BACKGROUND_POINTS else NULL
    return(purrr::imap(METRICS_RESULT, function(res, nm) {
      PlotNestHulls(
        METRICS_RESULT    = res,
        FILL_MODE         = FILL_MODE,        FILL_METRIC       = FILL_METRIC,
        FILL_PALETTE      = FILL_PALETTE,     FILL_DIRECTION    = FILL_DIRECTION,
        FILL_ALPHA        = FILL_ALPHA,        CONTOUR_CLOSED    = CONTOUR_CLOSED,
        COLOR_CLOSED      = COLOR_CLOSED,      COLOR_OPEN        = COLOR_OPEN,
        LINEWIDTH_CLOSED  = LINEWIDTH_CLOSED,  LINEWIDTH_OPEN    = LINEWIDTH_OPEN,
        BACKGROUND_POINTS = if (!is.null(bg_list)) bg_list[[nm]] else NULL,
        BG_ANNO_COLUMN    = BG_ANNO_COLUMN,    BG_X              = BG_X,
        BG_Y              = BG_Y,              BG_COLOR_MAP      = BG_COLOR_MAP,
        BG_POINT_SIZE     = BG_POINT_SIZE,     BG_POINT_ALPHA    = BG_POINT_ALPHA,
        LABEL_NESTS       = LABEL_NESTS,       LABEL_SIZE        = LABEL_SIZE,
        SAMPLE_NAME       = nm
      )
    }))
  }
  
  # Guards
  hull_df <- METRICS_RESULT$HullMap
  pernest <- METRICS_RESULT$PerNest
  if (is.null(hull_df) || nrow(hull_df) == 0)
    stop("[Synora] HullMap is empty. Re-run GetShapeMetrics with SEPARATE_NESTS = TRUE.")
  
  VALID_MODES <- c("identity", "metric", "closure", "none")
  if (!FILL_MODE %in% VALID_MODES)
    stop("FILL_MODE must be one of: ", paste(VALID_MODES, collapse = ", "))
  if (FILL_MODE == "metric") {
    if (is.null(FILL_METRIC))
      stop("FILL_MODE = 'metric' requires FILL_METRIC to be specified.")
    if (!FILL_METRIC %in% names(hull_df))
      stop("FILL_METRIC '", FILL_METRIC, "' not found in HullMap. ",
           "Available: ", paste(names(hull_df), collapse = ", "))
  }
  
  # Pre-compute per-vertex border aesthetics and centroid labels
  hull_df <- hull_df %>%
    dplyr::mutate(
      .border_color = dplyr::if_else(IsClosedNest, COLOR_CLOSED, COLOR_OPEN),
      .linewidth    = dplyr::if_else(IsClosedNest, LINEWIDTH_CLOSED, LINEWIDTH_OPEN)
    )
  centroids <- hull_df %>%
    dplyr::group_by(Nest_ID) %>%
    dplyr::summarise(cx = mean(X, na.rm = TRUE),
                     cy = mean(Y, na.rm = TRUE), .groups = "drop")
  
  # Base plot
  p <- ggplot2::ggplot()
  
  # Optional background cell scatter
  if (!is.null(BACKGROUND_POINTS) && is.data.frame(BACKGROUND_POINTS)) {
    bg <- BACKGROUND_POINTS
    if (BG_ANNO_COLUMN %in% names(bg)) {
      bg[[BG_ANNO_COLUMN]] <- as.character(bg[[BG_ANNO_COLUMN]])
      p <- p +
        ggplot2::geom_point(
          data = bg,
          mapping = ggplot2::aes(x     = !!as.name(BG_X),
                                 y     = !!as.name(BG_Y),
                                 color = !!as.name(BG_ANNO_COLUMN)),
          size = BG_POINT_SIZE, alpha = BG_POINT_ALPHA, inherit.aes = FALSE
        ) +
        ggplot2::scale_color_manual(values   = BG_COLOR_MAP,
                                    name     = "Cell type",
                                    na.value = "#D1D5DB")
    } else {
      p <- p +
        ggplot2::geom_point(
          data = bg,
          mapping = ggplot2::aes(x = !!as.name(BG_X), y = !!as.name(BG_Y)),
          size = BG_POINT_SIZE, alpha = BG_POINT_ALPHA,
          color = "#9CA3AF", inherit.aes = FALSE
        )
    }
  }
  
  # Hull fill layer
  if (FILL_MODE == "identity") {
    nest_cols        <- scales::hue_pal()(dplyr::n_distinct(hull_df$Nest_ID))
    names(nest_cols) <- unique(hull_df$Nest_ID)
    p <- p +
      ggplot2::geom_polygon(
        data    = hull_df,
        mapping = ggplot2::aes(x = X, y = Y, group = Nest_ID, fill = Nest_ID),
        alpha = FILL_ALPHA, color = NA, inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_manual(
        values = nest_cols, name = "Nest ID",
        guide  = ggplot2::guide_legend(ncol = 2,
                                       override.aes = list(alpha = 0.8))
      )
    
  } else if (FILL_MODE == "metric") {
    p <- p +
      ggplot2::geom_polygon(
        data    = hull_df,
        mapping = ggplot2::aes(x = X, y = Y, group = Nest_ID,
                               fill = !!as.name(FILL_METRIC)),
        alpha = FILL_ALPHA, color = NA, inherit.aes = FALSE
      )
    p <- p + if (FILL_PALETTE == "viridis")
      ggplot2::scale_fill_viridis_c(name      = FILL_METRIC,
                                    direction = FILL_DIRECTION,
                                    na.value  = "#E5E7EB")
    else
      ggplot2::scale_fill_distiller(palette   = FILL_PALETTE,
                                    direction = FILL_DIRECTION,
                                    name      = FILL_METRIC,
                                    na.value  = "#E5E7EB")
    
  } else if (FILL_MODE == "closure") {
    p <- p +
      ggplot2::geom_polygon(
        data    = hull_df,
        mapping = ggplot2::aes(x = X, y = Y, group = Nest_ID,
                               fill = IsClosedNest),
        alpha = FILL_ALPHA, color = NA, inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_manual(
        values = c("TRUE"  = scales::alpha(COLOR_CLOSED, 0.6),
                   "FALSE" = scales::alpha(COLOR_OPEN,   0.4)),
        labels = c("TRUE" = "Closed", "FALSE" = "Open"),
        name   = "Nest type"
      )
    
  } else {
    # FILL_MODE == "none" — placeholder polygon, no fill or stroke
    p <- p +
      ggplot2::geom_polygon(
        data    = hull_df,
        mapping = ggplot2::aes(x = X, y = Y, group = Nest_ID),
        fill = NA, color = NA, inherit.aes = FALSE
      )
  }
  
  # Hull outlines — vectorised: split closed/open, one geom_polygon each
  if (CONTOUR_CLOSED) {
    closed_hulls <- hull_df %>% dplyr::filter( IsClosedNest)
    open_hulls   <- hull_df %>% dplyr::filter(!IsClosedNest)
    if (nrow(closed_hulls) > 0)
      p <- p + ggplot2::geom_polygon(
        data    = closed_hulls,
        mapping = ggplot2::aes(x = X, y = Y, group = Nest_ID),
        fill = NA, color = COLOR_CLOSED, linewidth = LINEWIDTH_CLOSED,
        inherit.aes = FALSE
      )
    if (nrow(open_hulls) > 0)
      p <- p + ggplot2::geom_polygon(
        data    = open_hulls,
        mapping = ggplot2::aes(x = X, y = Y, group = Nest_ID),
        fill = NA, color = COLOR_OPEN, linewidth = LINEWIDTH_OPEN,
        inherit.aes = FALSE
      )
  } else {
    # Variable colour/linewidth per nest — still vectorised via aes()
    p <- p + ggplot2::geom_polygon(
      data    = hull_df,
      mapping = ggplot2::aes(x = X, y = Y, group = Nest_ID,
                             color     = I(.border_color),
                             linewidth = I(.linewidth)),
      fill = NA, inherit.aes = FALSE
    )
  }
  
  # Optional nest ID labels at centroids
  if (LABEL_NESTS)
    p <- p + ggplot2::geom_text(
      data    = centroids,
      mapping = ggplot2::aes(x = cx, y = cy, label = Nest_ID),
      size = LABEL_SIZE, fontface = "bold", color = "black",
      inherit.aes = FALSE
    )
  
  # Subtitle: closed/open counts
  n_closed <- sum(pernest$IsClosedNest, na.rm = TRUE)
  n_open   <- nrow(pernest) - n_closed
  
  p +
    ggplot2::coord_equal() +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::labs(
      title    = SAMPLE_NAME,
      subtitle = paste0(n_closed, " closed nest(s)  \u2014  ", n_open, " open nest(s)")
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 12,
                                            margin = ggplot2::margin(b = 3)),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey40",
                                            margin = ggplot2::margin(b = 4)),
      legend.title  = ggplot2::element_text(size = 9, face = "bold"),
      legend.text   = ggplot2::element_text(size = 8)
    )
}