# GetBoundary ----
#' @title GetBoundary
#' @description Detect and annotate boundary cells in spatial omics data.
#'   Supports three input modes:
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
#' @param X_POSITION Name of X coordinate column.
#' @param Y_POSITION Name of Y coordinate column.
#' @param ANNO_COLUMN Name of annotation column. 1 = cell of interest (e.g.
#'   tumor), 0 = other. May be binary or continuous.
#' @param CELL_ID_COLUMN (optional) Cell ID column name. Row indices used if
#'   omitted. Default NULL.
#' @param CELL_ID_PREFIX (optional) Prefix for auto-generated cell IDs.
#'   Default NULL.
#' @param SAMPLE_COLUMN (optional) Column name identifying samples within a
#'   single data frame. When provided, the data frame is split by this column
#'   and each sample is processed independently. Default NULL.
#' @param ANNO_RANGE Annotation range. Default \code{c(0, 1)}. Use
#'   \code{"auto"} to derive from data.
#' @param ANNO_MIDPOINT Annotation midpoint. Numeric or \code{"auto"} (Otsu
#'   thresholding). Default \code{0.5}.
#' @param NEIGHBOR_METHOD Neighborhood method: \code{"radius"},
#'   \code{"knn"}, or \code{"hybrid"}. Default \code{"radius"}.
#' @param RADIUS Neighborhood radius (required for \code{"radius"} /
#'   \code{"hybrid"}).
#' @param KNN_K Number of nearest neighbors (required for \code{"knn"} /
#'   \code{"hybrid"}).
#' @param DENOISE Logical; iterative DBSCAN noise removal. Default \code{TRUE}.
#' @param NEST_MIN_SIZE Minimum cluster size for denoising. Default \code{5}.
#' @param NEST_SPECIFICITY Nest membership threshold (0-1). Default \code{0.25}.
#' @param BOUNDARY_SPECIFICITY Boundary score threshold (0-1). Default
#'   \code{0.05}.
#' @param STRATIFY_BOUNDARY Logical; sub-classify Boundary cells into
#'   \code{Boundary_Inner} (nest-facing) and \code{Boundary_Outer}
#'   (stroma-facing). Default \code{TRUE}.
#' @param VERBOSE Logical; print progress messages. Default \code{FALSE}.
#'
#' @return
#' \describe{
#'   \item{List input}{A named list of data frames, one per input element.}
#'   \item{Data frame + \code{SAMPLE_COLUMN}}{A single data frame with all
#'     samples reassembled, retaining the \code{SAMPLE_COLUMN}.}
#'   \item{Data frame (no \code{SAMPLE_COLUMN})}{A single data frame.}
#' }
#' Each output data frame contains all original columns plus:
#' \describe{
#'   \item{Nb_Count}{Number of neighbors.}
#'   \item{Anno_Midpoint}{Threshold used for Nest/Outside separation.}
#'   \item{Mixedness}{Local heterogeneity score.}
#'   \item{Orientedness}{Synora directional asymmetry metric.}
#'   \item{BoundaryScore}{Mixedness x Orientedness.}
#'   \item{SynoraAnnotation}{Factor annotation (see \code{STRATIFY_BOUNDARY}).
#'     Cells with missing \code{ANNO_COLUMN} values receive \code{NA}.}
#' }
#' 
#' @examples
#' library(Synora)
#' data("DummyData")
#'
#' # --- Single data frame (no SAMPLE_COLUMN) ----------------
#' result_single <- GetBoundary(
#'   INPUT                = DummyData[[1]],
#'   CELL_ID_COLUMN       = "Cell_ID",
#'   X_POSITION           = "X",
#'   Y_POSITION           = "Y",
#'   ANNO_COLUMN          = "CT",
#'   RADIUS               = 20,
#'   NEST_SPECIFICITY     = 0.25,
#'   BOUNDARY_SPECIFICITY = 0.05
#' )
#' head(result_single[, c("Cell_ID", "SynoraAnnotation", "BoundaryScore")])
#'
#' # --- List of data frames (batch mode) --------------------
#' result_list <- GetBoundary(
#'   INPUT                = DummyData,          # named list
#'   CELL_ID_COLUMN       = "Cell_ID",
#'   X_POSITION           = "X",
#'   Y_POSITION           = "Y",
#'   ANNO_COLUMN          = "CT",
#'   RADIUS               = 20,
#'   NEST_SPECIFICITY     = 0.25,
#'   BOUNDARY_SPECIFICITY = 0.05
#' )
#' # result_list is a named list; access individual samples:
#' table(result_list[["Sample1"]]$SynoraAnnotation)
#'
#' # --- Single data frame WITH SAMPLE_COLUMN ----------------
#' combined_df <- do.call(rbind, lapply(
#'   names(DummyData),
#'   function(nm) { d <- DummyData[[nm]]; d$Sample <- nm; d }
#' ))
#' result_combined <- GetBoundary(
#'   INPUT                = combined_df,
#'   SAMPLE_COLUMN        = "Sample",
#'   CELL_ID_COLUMN       = "Cell_ID",
#'   X_POSITION           = "X",
#'   Y_POSITION           = "Y",
#'   ANNO_COLUMN          = "CT",
#'   RADIUS               = 20,
#'   NEST_SPECIFICITY     = 0.25,
#'   BOUNDARY_SPECIFICITY = 0.05
#' )
#' # Returns one reassembled data frame retaining the Sample column
#' table(result_combined$Sample, result_combined$SynoraAnnotation)
#'
#' # --- KNN neighborhood + Otsu midpoint (continuous anno) --
#' result_knn <- GetBoundary(
#'   INPUT                = DummyData[[1]],
#'   CELL_ID_COLUMN       = "Cell_ID",
#'   X_POSITION           = "X",
#'   Y_POSITION           = "Y",
#'   ANNO_COLUMN          = "CT",
#'   NEIGHBOR_METHOD      = "knn",
#'   KNN_K                = 10,
#'   ANNO_RANGE           = "auto",
#'   ANNO_MIDPOINT        = "auto",   # Otsu thresholding
#'   NEST_SPECIFICITY     = 0.25,
#'   BOUNDARY_SPECIFICITY = 0.05,
#'   STRATIFY_BOUNDARY    = TRUE,
#'   VERBOSE              = TRUE
#' )
#' table(result_knn$SynoraAnnotation)
#'
#' @export
#' @importFrom magrittr `%>%`
GetBoundary <- function(INPUT, X_POSITION, Y_POSITION,
                        ANNO_COLUMN, CELL_ID_COLUMN, CELL_ID_PREFIX,
                        SAMPLE_COLUMN        = NULL,
                        ANNO_RANGE           = c(0, 1),
                        ANNO_MIDPOINT        = 0.5,
                        NEIGHBOR_METHOD      = c("radius", "knn", "hybrid"),
                        RADIUS, KNN_K,
                        DENOISE              = TRUE,
                        NEST_MIN_SIZE        = 5,
                        NEST_SPECIFICITY     = 0.25,
                        BOUNDARY_SPECIFICITY = 0.05,
                        STRATIFY_BOUNDARY    = TRUE,
                        VERBOSE              = FALSE) {
  
  FORWARDED_ARGS <- list(
    X_POSITION           = X_POSITION,
    Y_POSITION           = Y_POSITION,
    ANNO_COLUMN          = ANNO_COLUMN,
    ANNO_RANGE           = ANNO_RANGE,
    ANNO_MIDPOINT        = ANNO_MIDPOINT,
    NEIGHBOR_METHOD      = match.arg(NEIGHBOR_METHOD),
    DENOISE              = DENOISE,
    NEST_MIN_SIZE        = NEST_MIN_SIZE,
    NEST_SPECIFICITY     = NEST_SPECIFICITY,
    BOUNDARY_SPECIFICITY = BOUNDARY_SPECIFICITY,
    STRATIFY_BOUNDARY    = STRATIFY_BOUNDARY,
    VERBOSE              = VERBOSE
  )
  if (!missing(CELL_ID_COLUMN)) FORWARDED_ARGS$CELL_ID_COLUMN <- CELL_ID_COLUMN
  if (!missing(CELL_ID_PREFIX)) FORWARDED_ARGS$CELL_ID_PREFIX <- CELL_ID_PREFIX
  if (!missing(RADIUS))         FORWARDED_ARGS$RADIUS         <- RADIUS
  if (!missing(KNN_K))          FORWARDED_ARGS$KNN_K          <- KNN_K
  
  if (is.list(INPUT) && !is.data.frame(INPUT)) {
    # Mode 1: List of data frames ----
    
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
    
    RESULT_LIST <- purrr::imap(
      .x        = INPUT,
      .progress = "[Synora] Processing",
      .f        = function(df, nm) {
        tryCatch(
          do.call(.GetBoundary_Single, c(list(INPUT = df), FORWARDED_ARGS)),
          error = function(e)
            stop("Error in sample '\033[1;4;41m", nm, "\033[0m': ",
                 conditionMessage(e))
        )
      }
    )
    return(RESULT_LIST)
    
  } else if (is.data.frame(INPUT)) {
    
    if (!is.null(SAMPLE_COLUMN)) {
    # Mode 2: SAMPLE_COLUMN provided — split, process, reassemble ----
      
      if (!SAMPLE_COLUMN %in% names(INPUT))
        stop("SAMPLE_COLUMN: `\033[1;4;41m", SAMPLE_COLUMN,
             "\033[0m` not found in INPUT.")
      
      sample_ids <- unique(INPUT[[SAMPLE_COLUMN]])
      n_samples  <- length(sample_ids)
      
      if (n_samples == 1)
        warning("[Synora] SAMPLE_COLUMN \033[1;4;43m", SAMPLE_COLUMN,
                "\033[0m has only one unique value ('", sample_ids,
                "'). Processing as a single sample.")
      
      message("[Synora] Mode: Data frame with SAMPLE_COLUMN '", SAMPLE_COLUMN,
              "' — ", n_samples, " sample(s) detected.")
      
      RESULT_DF <- INPUT %>%
        base::split(.[[SAMPLE_COLUMN]]) %>%
        purrr::imap(
          .progress = "[Synora] Processing",
          .f = function(df, nm) {
            tryCatch(
              do.call(.GetBoundary_Single, c(list(INPUT = df), FORWARDED_ARGS)),
              error = function(e)
                stop("Error in sample '\033[1;4;41m", nm, "\033[0m': ",
                     conditionMessage(e))
            )
          }
        ) %>%
        dplyr::bind_rows()
      
      return(RESULT_DF)
      
    } else {
      # Mode 3: Plain single data frame ----

      message("[Synora] Mode: Single data frame.")
      return(do.call(.GetBoundary_Single, c(list(INPUT = INPUT), FORWARDED_ARGS)))
      
    }
    
    # Invalid input type-----
  } else {
    stop("INPUT must be a data frame or a list of data frames. ",
         "Received: \033[1;4;41m", paste(class(INPUT), collapse = ", "), "\033[0m")
  }
}

# .GetBoundary_Single
#' @keywords internal
.GetBoundary_Single <- function(INPUT, X_POSITION, Y_POSITION,
                                ANNO_COLUMN, CELL_ID_COLUMN, CELL_ID_PREFIX,
                                ANNO_RANGE           = c(0, 1),
                                ANNO_MIDPOINT        = 0.5,
                                NEIGHBOR_METHOD      = c("radius", "knn", "hybrid"),
                                RADIUS               = NULL,
                                KNN_K                = NULL,
                                DENOISE              = TRUE,
                                NEST_MIN_SIZE        = 5,
                                NEST_SPECIFICITY     = 0.25,
                                BOUNDARY_SPECIFICITY = 0.05,
                                STRATIFY_BOUNDARY    = TRUE,
                                VERBOSE              = FALSE) {
  
  # 1. Input validation ----
  if (missing(INPUT))        stop("INPUT data frame must be provided.")
  if (!is.data.frame(INPUT)) stop("INPUT must be a data frame.")
  if (missing(X_POSITION) || missing(Y_POSITION))
    stop("Both X_POSITION and Y_POSITION must be provided.")
  if (!X_POSITION %in% names(INPUT))
    stop("X_POSITION: `\033[1;4;41m", X_POSITION, "\033[0m` not found in INPUT.")
  if (!Y_POSITION %in% names(INPUT))
    stop("Y_POSITION: `\033[1;4;41m", Y_POSITION, "\033[0m` not found in INPUT.")
  
  if (missing(CELL_ID_COLUMN)) {
    if (VERBOSE) message("  Creating Cell_ID...")
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
  
  if (any(duplicated(INPUT[, c(X_POSITION, Y_POSITION)])))
    stop("Duplicate (X, Y) coordinates detected. ",
         "Please remove overlapping cells first.")
  
  # Drop any pre-existing Synora output columns to avoid join conflicts
  synora_cols     <- c("Nb_Count", "Anno_Midpoint", "Mixedness",
                       "Orientedness", "BoundaryScore", "SynoraAnnotation")
  existing_synora <- intersect(synora_cols, names(INPUT))
  if (length(existing_synora) > 0) {
    warning("Overwriting existing Synora column(s): \033[1;4;43m[",
            paste(existing_synora, collapse = ", "),
            "]\033[0m. Check if GetBoundary has already been run on this data.",
            call. = FALSE)
    INPUT <- INPUT %>% dplyr::select(-dplyr::all_of(existing_synora))
  }
  
  if (missing(ANNO_COLUMN))
    stop("ANNO_COLUMN must be provided.")
  if (!ANNO_COLUMN %in% names(INPUT))
    stop("ANNO_COLUMN: `\033[1;4;41m", ANNO_COLUMN, "\033[0m` not found in INPUT.")
  if (!is.numeric(INPUT[[ANNO_COLUMN]]))
    stop("ANNO_COLUMN `\033[1;4;41m", ANNO_COLUMN,
         "\033[0m` must be numeric. Convert to numeric first.")

  NEIGHBOR_METHOD <- match.arg(NEIGHBOR_METHOD)
  if (NEIGHBOR_METHOD %in% c("radius", "hybrid") && is.null(RADIUS))
    stop("RADIUS is required for NEIGHBOR_METHOD \033[1;4;41m'",
         NEIGHBOR_METHOD, "'\033[0m.")
  if (NEIGHBOR_METHOD %in% c("radius", "hybrid") &&
      (!is.numeric(RADIUS) || RADIUS <= 0))
    stop("RADIUS must be a positive numeric. Got: \033[1;4;41m",
         RADIUS, "\033[0m.")
  if (NEIGHBOR_METHOD %in% c("knn", "hybrid") && is.null(KNN_K))
    stop("KNN_K is required for NEIGHBOR_METHOD \033[1;4;41m'",
         NEIGHBOR_METHOD, "'\033[0m.")
  if (NEIGHBOR_METHOD %in% c("knn", "hybrid") &&
      (!is.numeric(KNN_K) || KNN_K <= 0 || KNN_K != as.integer(KNN_K)))
    stop("KNN_K must be a positive integer. Got: \033[1;4;41m",
         KNN_K, "\033[0m.")
  
  na_mask  <- !is.finite(INPUT[[ANNO_COLUMN]])
  n_na     <- sum(na_mask)
  
  if (n_na > 0) {
    na_ids     <- INPUT[[CELL_ID_COLUMN]][na_mask]
    id_display <- if (length(na_ids) <= 10) {
      paste(na_ids, collapse = ", ")
    } else {
      paste0(paste(head(na_ids, 10), collapse = ", "),
             " ... [", length(na_ids) - 10, " more]")
    }
    warning(
      "\033[1;43m[Synora] ANNO_COLUMN `", ANNO_COLUMN, "` contains ",
      "\033[1;4;43m", n_na, " non-finite value(s)\033[0;1;43m ",
      "(NA / NaN / Inf).\033[0m\n",
      "  Affected cell ID(s): \033[1;4;43m", id_display, "\033[0m\n",
      "  These cells will be \033[1;43mexcluded from analysis\033[0m ",
      "and returned with NA in all Synora output columns.",
      call. = FALSE
    )
  }
  
  # 2. Preparation ----
  INPUT_FULL  <- INPUT %>% 
    dplyr::arrange(!!as.name(CELL_ID_COLUMN))
  INPUT_CLEAN <- INPUT_FULL %>% 
    dplyr::filter(is.finite(!!as.name(ANNO_COLUMN)))
  
  if (identical(ANNO_RANGE, "auto")) {
    ANNO_RANGE <- range(INPUT_CLEAN[[ANNO_COLUMN]], na.rm = TRUE)
    if (diff(ANNO_RANGE) == 0) ANNO_RANGE <- ANNO_RANGE + c(-1e-6, 1e-6)
  } else {
    if (!is.numeric(ANNO_RANGE) || length(ANNO_RANGE) != 2 ||
        ANNO_RANGE[1] >= ANNO_RANGE[2])
      stop("ANNO_RANGE must be \033[1;4;41m'auto'\033[0m or a numeric vector ",
           "of length 2 with min < max.")
  }
  
  if (!(is.numeric(ANNO_MIDPOINT) && length(ANNO_MIDPOINT) == 1) &&
      !identical(ANNO_MIDPOINT, "auto"))
    stop("ANNO_MIDPOINT must be a single numeric or \033[1;4;41m'auto'\033[0m.")
  
  if (identical(ANNO_MIDPOINT, "auto")) {
    ANNO_MIDPOINT <- .OtsuThreshold(INPUT_CLEAN[[ANNO_COLUMN]])
    if (!dplyr::between(ANNO_MIDPOINT, ANNO_RANGE[1], ANNO_RANGE[2])) {
      warning("\033[1;43m[Synora] Otsu midpoint (",
              round(ANNO_MIDPOINT, 3),
              ") is outside ANNO_RANGE. Falling back to range midpoint.\033[0m",
              call. = FALSE)
      ANNO_MIDPOINT <- mean(ANNO_RANGE)
    }
    if (VERBOSE)
      message("  Otsu ANNO_MIDPOINT: ", round(ANNO_MIDPOINT, 3))
  }
  
  ANNO_SCALED   <- paste0(ANNO_COLUMN, "_scaled")
  ANNO_DENOISED <- paste0(ANNO_COLUMN, "_denoised")
  
  INPUT_CLEAN <- INPUT_CLEAN %>%
    dplyr::mutate(!!as.name(ANNO_SCALED) := ifelse(
      !!as.name(ANNO_COLUMN) >= ANNO_MIDPOINT,
      (!!as.name(ANNO_COLUMN) - ANNO_MIDPOINT) / (ANNO_RANGE[2] - ANNO_MIDPOINT),
      (!!as.name(ANNO_COLUMN) - ANNO_MIDPOINT) / (ANNO_MIDPOINT - ANNO_RANGE[1])
    ))
  
  NN <- .GetNeighbors(
    COORDS          = INPUT_CLEAN %>%
      dplyr::select(dplyr::all_of(c(X_POSITION, Y_POSITION))),
    NEIGHBOR_METHOD = NEIGHBOR_METHOD,
    RADIUS          = RADIUS,
    KNN_K           = KNN_K
  )
  
  # 3. Core logic ----
  if (DENOISE) {
    
    RESULT_1 <- INPUT_CLEAN %>%
      .GetMO(CELL_ID_COLUMN, X_POSITION, Y_POSITION, ANNO_SCALED, NN, FALSE) %>%
      dplyr::mutate(Nest = (0.5 * Mean_Anno + 0.5) >= NEST_SPECIFICITY)
    
    NOISE_COUNT <- c(); i <- 1
    repeat {
      RESULT_1 <- RESULT_1 %>%
        dplyr::nest_by(Nest, NeighborhoodRadius) %>%
        dplyr::mutate(
          !!as.name(paste0("Noise_", i)) :=
            data %>%
            dplyr::select(dplyr::all_of(c(X_POSITION, Y_POSITION))) %>%
            dbscan::dbscan(eps = NeighborhoodRadius, minPts = NEST_MIN_SIZE) %>%
            .$cluster %>% { . == 0 } %>% list()
        ) %>%
        tidyr::unnest(cols = c(data, !!as.name(paste0("Noise_", i)))) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(!!as.name(CELL_ID_COLUMN)) %>%
        dplyr::mutate(Nest = ifelse(!!as.name(paste0("Noise_", i)), !Nest, Nest))
      
      NOISE_COUNT <- c(NOISE_COUNT, sum(RESULT_1[[paste0("Noise_", i)]]))
      converged   <- if (length(NOISE_COUNT) == 1) NOISE_COUNT == 0 else
        diff(tail(NOISE_COUNT, 2)) == 0
      
      if (converged) {
        RESULT_1 <- RESULT_1 %>%
          dplyr::mutate(!!as.name(ANNO_DENOISED) := dplyr::case_when(
            !(!!as.name(paste0("Noise_", i))) &  Nest ~  1,
            !(!!as.name(paste0("Noise_", i))) & !Nest ~ -1,
            TRUE ~ NA_real_
          ))
        break
      }
      i <- i + 1
    }
    
    RESULT <- RESULT_1 %>%
      .GetMO(CELL_ID_COLUMN, X_POSITION, Y_POSITION, ANNO_DENOISED, NN, TRUE) %>%
      dplyr::mutate(
        BoundaryScore    = Mixedness * Orientedness,
        SynoraAnnotation = dplyr::case_when(
          BoundaryScore >= BOUNDARY_SPECIFICITY                                    ~ "Boundary",
          BoundaryScore <  BOUNDARY_SPECIFICITY & !!as.name(ANNO_DENOISED) ==  1   ~ "Nest",
          BoundaryScore <  BOUNDARY_SPECIFICITY & !!as.name(ANNO_DENOISED) == -1   ~ "Outside",
          TRUE                                                                     ~ "Noise"
        ) %>% factor(levels = c("Boundary", "Nest", "Outside", "Noise"))
      ) %>%
      dplyr::select(
        dplyr::all_of(c(CELL_ID_COLUMN, X_POSITION, Y_POSITION)),
        Nb_Count, Mixedness, Orientedness, BoundaryScore, SynoraAnnotation
      )
    
  } else {
    
    RESULT <- INPUT_CLEAN %>%
      .GetMO(CELL_ID_COLUMN, X_POSITION, Y_POSITION, ANNO_SCALED, NN, TRUE) %>%
      dplyr::mutate(
        BoundaryScore    = Mixedness * Orientedness,
        SynoraAnnotation = dplyr::case_when(
          BoundaryScore >= BOUNDARY_SPECIFICITY        ~ "Boundary",
          (0.5 * Mean_Anno + 0.5) >= NEST_SPECIFICITY ~ "Nest",
          TRUE                                         ~ "Outside"
        ) %>% factor(levels = c("Boundary", "Nest", "Outside", "Noise"))
      ) %>%
      dplyr::select(
        dplyr::all_of(c(CELL_ID_COLUMN, X_POSITION, Y_POSITION)),
        Nb_Count, Mixedness, Orientedness, BoundaryScore, SynoraAnnotation
      )
  }
  
  # 4. Final assembly ----
  JOIN_KEY <- c(CELL_ID_COLUMN, X_POSITION, Y_POSITION)
  
  RESULT_CLEAN <- INPUT_CLEAN %>%
    dplyr::select(-dplyr::all_of(ANNO_SCALED)) %>% 
    dplyr::left_join(RESULT, by = JOIN_KEY) %>%
    dplyr::mutate(
      Anno_Midpoint = ANNO_MIDPOINT,
      Nb_Count      = Nb_Count,
      Mixedness     = ifelse(is.na(Mixedness),     0, Mixedness),
      Orientedness  = ifelse(is.na(Orientedness),  0, Orientedness),
      BoundaryScore = ifelse(is.na(BoundaryScore), 0, BoundaryScore)
    )
  
  if (STRATIFY_BOUNDARY) {
    RESULT_CLEAN <- .StratifyBoundary(RESULT = RESULT_CLEAN, NN = NN)
  }
  
  if (n_na > 0) {
    final_levels <- levels(RESULT_CLEAN$SynoraAnnotation)
    RESULT_FINAL <- INPUT_FULL %>%
      dplyr::left_join(
        RESULT_CLEAN %>%
          dplyr::select(dplyr::all_of(c(JOIN_KEY,
                                        "Nb_Count", "Anno_Midpoint",
                                        "Mixedness", "Orientedness",
                                        "BoundaryScore", "SynoraAnnotation"))),
        by = JOIN_KEY
      ) %>%
      dplyr::mutate(
        SynoraAnnotation = factor(SynoraAnnotation, levels = final_levels)
      )
  } else {
    RESULT_FINAL <- RESULT_CLEAN
  }
  
  return(RESULT_FINAL)
}


# Internal Helpers ----
.OtsuThreshold <- function(data, n_bins = 256) {
  tryCatch({
    data <- data[is.finite(data)]
    if (length(unique(data)) < 2) {
      message("Insufficient variation for Otsu thresholding, returning median")
      return(median(data))
    }
    breaks   <- seq(min(data), max(data), length.out = n_bins + 1)
    h        <- graphics::hist(data, breaks = breaks, plot = FALSE)
    bin_mids <- h$mids
    p        <- h$counts / sum(h$counts)
    
    best_thresh <- bin_mids[1]
    best_var    <- -Inf
    for (i in seq_len(length(p) - 1)) {
      w0 <- sum(p[seq_len(i)])
      w1 <- sum(p[(i + 1):length(p)])
      if (w0 == 0 || w1 == 0) next
      mu0 <- sum(p[seq_len(i)]        * bin_mids[seq_len(i)])        / w0
      mu1 <- sum(p[(i + 1):length(p)] * bin_mids[(i + 1):length(p)]) / w1
      var_between <- w0 * w1 * (mu0 - mu1)^2
      if (var_between > best_var) {
        best_var    <- var_between
        best_thresh <- (bin_mids[i] + bin_mids[i + 1]) / 2
      }
    }
    return(best_thresh)
  }, error = function(e) {
    message("Otsu fallback to median: ", e$message)
    median(data)
  })
}

.GetNeighbors <- function(COORDS,
                          NEIGHBOR_METHOD = c("radius", "knn", "hybrid"),
                          RADIUS = NULL, KNN_K = NULL) {
  NEIGHBOR_METHOD <- match.arg(NEIGHBOR_METHOD)
  COORDS <- as.data.frame(COORDS)
  
  if (NEIGHBOR_METHOD == "radius") {
    nn     <- dbscan::frNN(COORDS, eps = RADIUS, sort = FALSE)
    id     <- nn$id;  dist <- nn$dist;  radius <- RADIUS
    
  } else if (NEIGHBOR_METHOD == "knn") {
    nn     <- dbscan::kNN(COORDS, k = KNN_K, sort = FALSE)
    radius <- stats::quantile(nn$dist[, KNN_K], 0.75)
    id     <- split(nn$id,   seq(nrow(nn$id)))
    dist   <- split(nn$dist, seq(nrow(nn$dist)))
    
  } else if (NEIGHBOR_METHOD == "hybrid") {
    nn         <- dbscan::kNN(COORDS, k = KNN_K, sort = FALSE)
    valid_mask <- nn$dist <= RADIUS
    id   <- lapply(seq_len(nrow(nn$id)),   \(i) nn$id[i,   valid_mask[i, ]])
    dist <- lapply(seq_len(nrow(nn$dist)), \(i) nn$dist[i, valid_mask[i, ]])
    radius <- RADIUS
  }
  list(id = id, dist = dist, radius = radius)
}

.GetMO <- function(INPUT, CELL_ID_COLUMN, X_POSITION, Y_POSITION,
                   ANNO_COLUMN, NN, ORIENTEDNESS = TRUE) {
  INPUT_SLIM <- INPUT %>%
    dplyr::transmute(
      Cell_ID = !!as.name(CELL_ID_COLUMN),
      X       = !!as.name(X_POSITION),
      Y       = !!as.name(Y_POSITION),
      Anno    = !!as.name(ANNO_COLUMN)
    )
  NeighborhoodRadius <- NN$radius
  
  if (ORIENTEDNESS) {
    RESULT <- INPUT_SLIM %>%
      dplyr::mutate(ID = NN$id, Dist = NN$dist) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        X_u     = list((.$X[ID] - X) / Dist),
        Y_u     = list((.$Y[ID] - Y) / Dist),
        Nb_Anno = list(.$Anno[ID])
      ) %>%
      dplyr::mutate(
        Nb_Count  = length(Nb_Anno), Mean_Anno = mean(Nb_Anno),
        X_u_SumBg = sum(X_u),        Y_u_SumBg = sum(Y_u),
        X_u_Sum   = sum(Nb_Anno * X_u), Y_u_Sum = sum(Nb_Anno * Y_u)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        Mixedness          = 1 - (Mean_Anno)^2,
        Orientedness       = (sqrt(X_u_Sum^2 + Y_u_Sum^2) -
                                sqrt(X_u_SumBg^2 + Y_u_SumBg^2)) / Nb_Count,
        Orientedness       = ifelse(Orientedness < 0, 0, Orientedness),
        NeighborhoodRadius = NeighborhoodRadius
      ) %>%
      dplyr::select(Cell_ID, X, Y, Anno, Nb_Count,
                    NeighborhoodRadius, Mean_Anno, Mixedness, Orientedness)
  } else {
    RESULT <- INPUT_SLIM %>%
      dplyr::mutate(ID = NN$id) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Nb_Anno = list(.$Anno[ID])) %>%
      dplyr::mutate(Nb_Count = length(Nb_Anno), Mean_Anno = mean(Nb_Anno)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        Mixedness          = 1 - (Mean_Anno)^2,
        NeighborhoodRadius = NeighborhoodRadius
      ) %>%
      dplyr::select(Cell_ID, X, Y, Anno, Nb_Count,
                    NeighborhoodRadius, Mean_Anno, Mixedness)
  }
  RESULT %>%
    dplyr::rename(
      !!as.name(CELL_ID_COLUMN) := Cell_ID,
      !!as.name(X_POSITION)     := X,
      !!as.name(Y_POSITION)     := Y,
      !!as.name(ANNO_COLUMN)    := Anno
    )
}

.StratifyBoundary <- function(RESULT, NN) {
  anno_vec     <- as.character(RESULT$SynoraAnnotation)
  boundary_idx <- which(anno_vec == "Boundary")
  new_levels   <- c("Boundary_Inner", "Boundary_Outer", "Nest", "Outside", "Noise")
  
  if (length(boundary_idx) == 0) {
    return(RESULT %>%
             dplyr::mutate(SynoraAnnotation = factor(anno_vec, levels = new_levels)))
  }
  
  new_labels <- vapply(boundary_idx, function(i) {
    nb_anno   <- anno_vec[NN$id[[i]]]
    n_nest    <- sum(nb_anno == "Nest",    na.rm = TRUE)
    n_outside <- sum(nb_anno == "Outside", na.rm = TRUE)
    if (n_nest >= n_outside) "Boundary_Inner" else "Boundary_Outer"
  }, character(1))
  
  anno_vec[boundary_idx] <- new_labels
  RESULT %>%
    dplyr::mutate(SynoraAnnotation = factor(anno_vec, levels = new_levels))
}

