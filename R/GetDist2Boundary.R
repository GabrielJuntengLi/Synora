# GetDist2Boundary
#' @title GetDist2Boundary
#' @description Calculate the distance of each cell to the nearest boundary
#'   cell(s), with optional signed directionality.
#'
#'   Supports three input modes (mirroring \code{GetBoundary}):
#'   \enumerate{
#'     \item A \strong{list of data frames} — processed per element, returns a
#'       named list of data frames.
#'     \item A \strong{single data frame with \code{SAMPLE_COLUMN}} — split by
#'       sample, processed independently, returns one reassembled data frame.
#'     \item A \strong{single data frame without \code{SAMPLE_COLUMN}} — treated
#'       as one sample, returns one data frame.
#'   }
#'
#' @param INPUT A data frame or a named/unnamed list of data frames.
#' @param X_POSITION Name of the X coordinate column.
#' @param Y_POSITION Name of the Y coordinate column.
#' @param ANNO_COLUMN Name of the annotation column containing boundary labels.
#'   Default \code{"SynoraAnnotation"}.
#' @param CELL_ID_COLUMN (optional) Cell ID column name. Row indices are used
#'   if omitted.
#' @param CELL_ID_PREFIX (optional) Prefix for auto-generated cell IDs.
#' @param SAMPLE_COLUMN (optional) Column name identifying samples within a
#'   single data frame. Default \code{NULL}.
#' @param ANNO_OF_BOUNDARY Annotation value(s) representing boundary cells.
#'   Accepts a character scalar or vector, e.g.
#'   \code{c("Boundary_Inner", "Boundary_Outer")} for stratified output.
#'   Default \code{"Boundary"}. Auto-detection rules:
#'   \itemize{
#'     \item If \code{"Boundary"} is not found but both stratified variants are
#'       present, both are used automatically.
#'     \item If only one stratified variant is present, only that one is used.
#'     \item If an explicit vector is supplied and only some values are absent
#'       in a given sample (common in multi-sample runs), the missing values
#'       are dropped with a warning and processing continues.
#'     \item If \emph{all} supplied values are absent, the function stops with
#'       an error.
#'   }
#' @param K Number of nearest boundary neighbors whose distances are averaged.
#'   Default \code{5}.
#' @param DIRECTION Named list with elements \code{negative} and
#'   \code{positive} specifying which annotation values receive negative /
#'   positive signed distances respectively. Set to \code{NULL} to return
#'   unsigned distances for all cells. Default
#'   \code{list(negative = "Outside", positive = "Nest")}.
#'
#' @return
#' All original INPUT columns are preserved. \code{GetDist2Boundary} appends:
#' \describe{
#'   \item{Distance2Boundary}{Mean Euclidean distance to the K nearest boundary
#'     cells. Signed if \code{DIRECTION} is set. Boundary cells receive
#'     \code{0}. Cells whose annotation is not covered by \code{DIRECTION}
#'     (e.g. Noise) receive \code{NA}. Cells with \code{NA} in
#'     \code{ANNO_COLUMN} receive \code{NA}.}
#'   \item{Boundary_kNN_IDs}{List-column of length K containing the cell IDs
#'     of the nearest boundary neighbors. \code{NULL} for boundary cells and
#'     cells excluded from analysis.}
#' }
#' 
#' @examples
#' library(Synora)
#' data("DummyData")
#'
#' # First annotate boundaries
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
#' # --- List mode (mirrors GetBoundary list output) ---------
#' dist_list <- GetDist2Boundary(
#'   INPUT          = anno_list,
#'   X_POSITION     = "X",
#'   Y_POSITION     = "Y",
#'   CELL_ID_COLUMN = "Cell_ID",
#'   K              = 5
#' )
#' # Each element gains Distance2Boundary and Boundary_kNN_IDs
#' summary(dist_list[["Sample1"]]$Distance2Boundary)
#'
#' # --- Single data frame, signed distances -----------------
#' dist_single <- GetDist2Boundary(
#'   INPUT          = anno_list[[1]],
#'   X_POSITION     = "X",
#'   Y_POSITION     = "Y",
#'   CELL_ID_COLUMN = "Cell_ID",
#'   K              = 5,
#'   DIRECTION      = list(negative = "Outside", positive = "Nest")
#' )
#' hist(dist_single$Distance2Boundary,
#'      main = "Signed Distance to Boundary", xlab = "Distance")
#'
#' # --- Unsigned distances for all cells --------------------
#' dist_unsigned <- GetDist2Boundary(
#'   INPUT          = anno_list[[1]],
#'   X_POSITION     = "X",
#'   Y_POSITION     = "Y",
#'   CELL_ID_COLUMN = "Cell_ID",
#'   K              = 3,
#'   DIRECTION      = NULL          # no sign applied
#' )
#'
#' # --- Stratified boundary (Boundary_Inner / Boundary_Outer)
#' anno_strat <- GetBoundary(
#'   INPUT                = DummyData[[1]],
#'   CELL_ID_COLUMN       = "Cell_ID",
#'   X_POSITION           = "X",
#'   Y_POSITION           = "Y",
#'   ANNO_COLUMN          = "CT",
#'   RADIUS               = 20,
#'   STRATIFY_BOUNDARY    = TRUE,
#'   NEST_SPECIFICITY     = 0.25,
#'   BOUNDARY_SPECIFICITY = 0.05
#' )
#' dist_strat <- GetDist2Boundary(
#'   INPUT              = anno_strat,
#'   X_POSITION         = "X",
#'   Y_POSITION         = "Y",
#'   CELL_ID_COLUMN     = "Cell_ID",
#'   ANNO_OF_BOUNDARY   = c("Boundary_Inner", "Boundary_Outer"),
#'   K                  = 5
#' )
#' table(dist_strat$SynoraAnnotation, is.na(dist_strat$Distance2Boundary))
#' 
#' @export
#' @importFrom magrittr `%>%`
GetDist2Boundary <- function(INPUT,
                             X_POSITION,
                             Y_POSITION,
                             ANNO_COLUMN      = "SynoraAnnotation",
                             CELL_ID_COLUMN,
                             CELL_ID_PREFIX,
                             SAMPLE_COLUMN    = NULL,
                             ANNO_OF_BOUNDARY = "Boundary",
                             K                = 5,
                             DIRECTION        = list(negative = "Outside",
                                                     positive = "Nest")) {
  
  FORWARDED_ARGS <- list(
    X_POSITION       = X_POSITION,
    Y_POSITION       = Y_POSITION,
    ANNO_COLUMN      = ANNO_COLUMN,
    ANNO_OF_BOUNDARY = ANNO_OF_BOUNDARY,
    K                = K,
    DIRECTION        = DIRECTION
  )
  if (!missing(CELL_ID_COLUMN)) FORWARDED_ARGS$CELL_ID_COLUMN <- CELL_ID_COLUMN
  if (!missing(CELL_ID_PREFIX)) FORWARDED_ARGS$CELL_ID_PREFIX <- CELL_ID_PREFIX
  
  if (is.list(INPUT) && !is.data.frame(INPUT)) {
    
    if (length(INPUT) == 0) stop("INPUT list is empty.")
    
    non_df <- which(!vapply(INPUT, is.data.frame, logical(1)))
    if (length(non_df) > 0)
      stop("INPUT list contains non-data-frame element(s) at position(s): \033[1;4;41m",
           paste(non_df, collapse = ", "), "\033[0m")
    
    if (is.null(names(INPUT))) {
      message("[Synora] INPUT is an unnamed list. ",
              "Assigning positional names: Sample_1, Sample_2, ...")
      INPUT <- purrr::set_names(INPUT, paste0("Sample_", seq_along(INPUT)))
    } else if (any(names(INPUT) == "")) {
      blank_idx <- which(names(INPUT) == "")
      names(INPUT)[blank_idx] <- paste0("Sample_", blank_idx)
      message("[Synora] Partially unnamed list. Filled blank name(s) at position(s): ",
              paste(blank_idx, collapse = ", "))
    }
    
    message("[Synora] Mode: List of data frames — ",
            length(INPUT), " sample(s) detected.")
    
    return(purrr::imap(
      .x        = INPUT,
      .progress = "[Synora] Processing",
      .f        = function(df, nm) {
        tryCatch(
          do.call(.GetDist2Boundary_Single, c(list(INPUT = df), FORWARDED_ARGS)),
          error = function(e)
            stop("Error in sample '\033[1;4;41m", nm, "\033[0m': ",
                 conditionMessage(e))
        )
      }
    ))
    
  } else if (is.data.frame(INPUT)) {
    
    if (!is.null(SAMPLE_COLUMN)) {
      
      if (!SAMPLE_COLUMN %in% names(INPUT))
        stop("SAMPLE_COLUMN: `\033[1;4;41m", SAMPLE_COLUMN,
             "\033[0m` not found in INPUT.")
      
      sample_ids <- unique(INPUT[[SAMPLE_COLUMN]])
      n_samples  <- length(sample_ids)
      
      if (n_samples == 1)
        warning("[Synora] SAMPLE_COLUMN \033[1;4;43m", SAMPLE_COLUMN,
                "\033[0m has only one unique value ('", sample_ids,
                "'). Processing as a single sample.", call. = FALSE)
      
      message("[Synora] Mode: Data frame with SAMPLE_COLUMN '", SAMPLE_COLUMN,
              "' — ", n_samples, " sample(s) detected.")
      
      return(
        INPUT %>%
          base::split(.[[SAMPLE_COLUMN]]) %>%
          purrr::imap(
            .progress = "[Synora] Processing",
            .f = function(df, nm) {
              tryCatch(
                do.call(.GetDist2Boundary_Single, c(list(INPUT = df), FORWARDED_ARGS)),
                error = function(e)
                  stop("Error in sample '\033[1;4;41m", nm, "\033[0m': ",
                       conditionMessage(e))
              )
            }
          ) %>%
          dplyr::bind_rows()
      )
      
    } else {
      
      message("[Synora] Mode: Single data frame.")
      return(do.call(.GetDist2Boundary_Single, c(list(INPUT = INPUT), FORWARDED_ARGS)))
      
    }
    
  } else {
    stop("INPUT must be a data frame or a list of data frames. ",
         "Received: \033[1;4;41m", paste(class(INPUT), collapse = ", "), "\033[0m")
  }
}


# .GetDist2Boundary_Single
#' @keywords internal
.GetDist2Boundary_Single <- function(INPUT,
                                     X_POSITION,
                                     Y_POSITION,
                                     ANNO_COLUMN      = "SynoraAnnotation",
                                     CELL_ID_COLUMN,
                                     CELL_ID_PREFIX,
                                     ANNO_OF_BOUNDARY = "Boundary",
                                     K                = 5,
                                     DIRECTION        = list(negative = "Outside",
                                                             positive = "Nest")) {
  
  # 1. Input validation ----
  if (missing(INPUT))        stop("INPUT data frame must be provided.")
  if (!is.data.frame(INPUT)) stop("INPUT must be a data frame.")
  if (missing(X_POSITION) || missing(Y_POSITION))
    stop("Both X_POSITION and Y_POSITION must be provided.")
  if (!X_POSITION %in% names(INPUT))
    stop("X_POSITION: `\033[1;4;41m", X_POSITION, "\033[0m` not found in INPUT.")
  if (!Y_POSITION %in% names(INPUT))
    stop("Y_POSITION: `\033[1;4;41m", Y_POSITION, "\033[0m` not found in INPUT.")
  if (!ANNO_COLUMN %in% names(INPUT))
    stop("ANNO_COLUMN: `\033[1;4;41m", ANNO_COLUMN, "\033[0m` not found in INPUT.")
  
  if (missing(CELL_ID_COLUMN)) {
    INPUT <- if (missing(CELL_ID_PREFIX)) {
      INPUT %>% dplyr::mutate(Cell_ID = dplyr::row_number())
    } else {
      INPUT %>% dplyr::mutate(Cell_ID = paste0(CELL_ID_PREFIX, "_",
                                               dplyr::row_number()))
    }
    CELL_ID_COLUMN <- "Cell_ID"
  } else {
    if (!CELL_ID_COLUMN %in% names(INPUT))
      stop("CELL_ID_COLUMN: `\033[1;4;41m", CELL_ID_COLUMN,
           "\033[0m` not found in INPUT.")
    if (any(duplicated(INPUT[[CELL_ID_COLUMN]])))
      stop("CELL_ID_COLUMN `\033[1;4;41m", CELL_ID_COLUMN,
           "\033[0m` contains duplicate IDs. IDs must be unique.")
  }
  
  out_cols     <- c("Distance2Boundary", "Boundary_kNN_IDs")
  existing_out <- intersect(out_cols, names(INPUT))
  if (length(existing_out) > 0) {
    warning("Overwriting existing column(s): \033[1;4;43m[",
            paste(existing_out, collapse = ", "),
            "]\033[0m. Check if GetDist2Boundary has already been run on this data.",
            call. = FALSE)
    INPUT <- INPUT %>% dplyr::select(-dplyr::all_of(existing_out))
  }
  
  if (!is.numeric(K) || length(K) != 1 || K <= 0 || K != as.integer(K))
    stop("K must be a positive integer. Got: \033[1;4;41m", K, "\033[0m.")
  
  # 2. ANNO_OF_BOUNDARY resolution ----
  anno_vals        <- unique(INPUT[[ANNO_COLUMN]][!is.na(INPUT[[ANNO_COLUMN]])])
  STRATIFIED_VARIANTS <- c("Boundary_Inner", "Boundary_Outer")
  
  if (identical(ANNO_OF_BOUNDARY, "Boundary") && !"Boundary" %in% anno_vals) {
    detected <- intersect(STRATIFIED_VARIANTS, anno_vals)
    if (length(detected) == 0)
      stop("ANNO_OF_BOUNDARY: `\033[1;4;41mBoundary\033[0m` not found in ",
           "ANNO_COLUMN '", ANNO_COLUMN, "' and no stratified variants ",
           "(Boundary_Inner / Boundary_Outer) detected either.\n",
           "  Available values: ", paste(anno_vals, collapse = ", "))
    if (length(detected) == 1) {
      message("[Synora] 'Boundary' not found. Using stratified variant '",
              detected, "' as boundary reference.")
    } else {
      message("[Synora] 'Boundary' not found. Using stratified variants [",
              paste(detected, collapse = ", "), "] as boundary reference.")
    }
    ANNO_OF_BOUNDARY <- detected
    
  } else {
    missing_vals <- setdiff(ANNO_OF_BOUNDARY, anno_vals)
    
    if (length(missing_vals) == length(ANNO_OF_BOUNDARY)) {
      stop("ANNO_OF_BOUNDARY value(s) \033[1;4;41m[",
           paste(missing_vals, collapse = ", "),
           "]\033[0m not found in ANNO_COLUMN '", ANNO_COLUMN, "'.\n",
           "  Available values: ", paste(anno_vals, collapse = ", "))
      
    } else if (length(missing_vals) > 0) {
      warning("[Synora] ANNO_OF_BOUNDARY value(s) \033[1;4;43m[",
              paste(missing_vals, collapse = ", "),
              "]\033[0m not found in this sample and will be ignored. ",
              "Using: [", paste(setdiff(ANNO_OF_BOUNDARY, missing_vals),
                                collapse = ", "), "].",
              call. = FALSE)
      ANNO_OF_BOUNDARY <- intersect(ANNO_OF_BOUNDARY, anno_vals)
    }
  }
  
  # 3. DIRECTION validation ----
  if (!is.null(DIRECTION)) {
    if (!is.list(DIRECTION) || !all(c("negative", "positive") %in% names(DIRECTION)))
      stop("DIRECTION must be a named list with elements 'negative' and 'positive', ",
           "or NULL.")
    if (identical(DIRECTION$negative, DIRECTION$positive))
      stop("DIRECTION$negative and DIRECTION$positive cannot be the same value.")
    missing_dir <- setdiff(c(DIRECTION$negative, DIRECTION$positive), anno_vals)
    if (length(missing_dir) > 0)
      warning("DIRECTION value(s) \033[1;4;43m[",
              paste(missing_dir, collapse = ", "),
              "]\033[0m not found in ANNO_COLUMN '", ANNO_COLUMN, "'. ",
              "Affected cells will receive NA for Distance2Boundary.",
              call. = FALSE)
  }
  
  # 4. Separate boundary / non-boundary / NA-anno ----
  na_mask  <- is.na(INPUT[[ANNO_COLUMN]])
  n_na     <- sum(na_mask)
  
  if (n_na > 0)
    warning("[Synora] ANNO_COLUMN `", ANNO_COLUMN, "` contains ",
            n_na, " NA value(s). These cells will receive NA for Distance2Boundary.",
            call. = FALSE)
  
  INPUT_VALID    <- INPUT[!na_mask, ]
  INPUT_NA       <- INPUT[na_mask,  ]
  
  is_boundary    <- INPUT_VALID[[ANNO_COLUMN]] %in% ANNO_OF_BOUNDARY
  BOUNDARY_CELLS <- INPUT_VALID[is_boundary,  ]
  QUERY_CELLS    <- INPUT_VALID[!is_boundary, ]
  
  n_boundary <- nrow(BOUNDARY_CELLS)
  
  # 5. Insufficient boundary cells ----
  if (n_boundary < 2) {
    warning("[Synora] Fewer than 2 boundary cells found (n = ", n_boundary, "). ",
            "Distance2Boundary will be NA for all cells.",
            call. = FALSE)
    return(
      INPUT %>%
        dplyr::mutate(
          Distance2Boundary = NA_real_,
          Boundary_kNN_IDs  = list(NULL)
        )
    )
  }
  
  k_eff <- min(K, n_boundary)
  if (k_eff < K)
    warning("[Synora] K = ", K, " but only ", n_boundary,
            " boundary cells available. Using K = ", k_eff, ".",
            call. = FALSE)
  
  # 6. kNN: query = non-boundary, reference = boundary ----
  BOUNDARY_COORDS <- BOUNDARY_CELLS %>%
    dplyr::select(dplyr::all_of(c(X_POSITION, Y_POSITION))) %>%
    as.data.frame()
  
  QUERY_COORDS <- QUERY_CELLS %>%
    dplyr::select(dplyr::all_of(c(X_POSITION, Y_POSITION))) %>%
    as.data.frame()
  
  knn_result <- dbscan::kNN(
    x     = BOUNDARY_COORDS,
    query = QUERY_COORDS,
    k     = k_eff
  )
  
  boundary_ids <- BOUNDARY_CELLS[[CELL_ID_COLUMN]]
  knn_cell_ids <- matrix(boundary_ids[knn_result$id],
                         nrow = nrow(knn_result$id),
                         ncol = ncol(knn_result$id))
  
  # 7. Assemble results ----
  QUERY_RESULT <- QUERY_CELLS %>%
    dplyr::mutate(
      Distance2Boundary = rowMeans(knn_result$dist),
      Boundary_kNN_IDs  = unname(as.list(as.data.frame(t(knn_cell_ids))))
    )
  
  BOUNDARY_RESULT <- BOUNDARY_CELLS %>%
    dplyr::mutate(
      Distance2Boundary = 0,
      Boundary_kNN_IDs  = list(NULL)
    )
  
  RESULT <- dplyr::bind_rows(QUERY_RESULT, BOUNDARY_RESULT) %>%
    { if (n_na > 0)
      dplyr::bind_rows(.,
                       INPUT_NA %>% dplyr::mutate(Distance2Boundary = NA_real_,
                                                  Boundary_kNN_IDs  = list(NULL)))
      else . } %>%
    dplyr::arrange(match(!!as.name(CELL_ID_COLUMN), INPUT[[CELL_ID_COLUMN]]))
  
  # 8. Apply directionality ----
  if (!is.null(DIRECTION)) {
    
    uncovered <- !is.na(RESULT[[ANNO_COLUMN]]) &
      !(RESULT[[ANNO_COLUMN]] %in% ANNO_OF_BOUNDARY) &
      RESULT[[ANNO_COLUMN]] != DIRECTION$negative &
      RESULT[[ANNO_COLUMN]] != DIRECTION$positive
    
    if (any(uncovered)) {
      uncovered_vals <- unique(RESULT[[ANNO_COLUMN]][uncovered])
      warning("[Synora] ", sum(uncovered), " cell(s) with annotation(s) [",
              paste(uncovered_vals, collapse = ", "),
              "] are not covered by DIRECTION and not boundary. ",
              "Their Distance2Boundary will be set to NA.",
              call. = FALSE)
    }
    
    RESULT <- RESULT %>%
      dplyr::mutate(
        Distance2Boundary = dplyr::case_when(
          .data[[ANNO_COLUMN]] %in% ANNO_OF_BOUNDARY  ~  Distance2Boundary,
          .data[[ANNO_COLUMN]] == DIRECTION$positive  ~  Distance2Boundary,
          .data[[ANNO_COLUMN]] == DIRECTION$negative  ~ -Distance2Boundary,
          TRUE                                        ~  NA_real_
        )
      )
  }
  
  return(RESULT)
}