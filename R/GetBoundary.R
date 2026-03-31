# GetBoundary ----
#' @title GetBoundary
#' @description Detect and annotate boundary cells in spatial omics data.
#' @param INPUT Input data frame containing cell coordinates and annotations. Each row represents a single cell.
#' @param X_POSITION Name of X coordinate column.
#' @param Y_POSITION Name of Y coordinate column.
#' @param ANNO_COLUMN Name of annotation column. 1 indicates cells of interest (e.g., tumor) and 0 indicates others. 
#'   Can be binary or continuous.
#' @param CELL_ID_COLUMN (optional) Cell ID column. If not provided, row indices will be used. Default NULL.
#' @param CELL_ID_PREFIX (optional) Prefix for auto-generated cell IDs. Default NULL.
#' @param ANNO_RANGE Annotation range. Default c(0, 1). Use "auto" to use min/max of data.
#' @param ANNO_MIDPOINT Annotation midpoint. Numeric value or "auto" to use Otsu's thresholding method. Default 0.5.
#' @param NEIGHBOR_METHOD Neighborhood detection: "radius" (fixed-radius), "knn" (k-nearest), or "hybrid" (kNN within radius). Default "radius".
#' @param RADIUS Neighborhood radius (required for "radius" and "hybrid").
#' @param KNN_K Number of nearest neighbors (required for "knn" and "hybrid").
#' @param DENOISE Logical; whether to perform iterative noise removal via DBSCAN. Default TRUE.
#' @param NEST_MIN_SIZE Minimum nest size for denoising. Default 5.
#' @param NEST_SPECIFICITY Threshold for nest membership (0 to 1). Default 0.25.
#' @param BOUNDARY_SPECIFICITY Threshold for boundary detection (0 to 1). Default 0.05.
#' @param STRATIFY_BOUNDARY Logical; whether to sub-classify Boundary cells into 'Boundary_Inner' (nest-facing) and 'Boundary_Outer' (stroma-facing). Default TRUE.
#' @param VERBOSE Logical; whether to print progress messages. Default FALSE.
#'
#' @return A data frame with original columns plus:
#' \describe{
#'   \item{Nb_Count}{Number of neighbors found for the cell.}
#'   \item{Anno_Midpoint}{The threshold used to separate Nest from Outside.}
#'   \item{Mixedness}{Local heterogeneity score.}
#'   \item{Orientedness}{Synora metric for local cell-type orientation.}
#'   \item{BoundaryScore}{Product of Mixedness and Orientedness.}
#'   \item{SynoraAnnotation}{Factor: \code{Boundary_Inner},
#'     \code{Boundary_Outer}, \code{Nest}, \code{Outside}, \code{Noise}
#'     (when \code{STRATIFY_BOUNDARY = TRUE}); or \code{Boundary},
#'     \code{Nest}, \code{Outside}, \code{Noise} otherwise.}
#' }
#' @export
#' @importFrom magrittr `%>%`
GetBoundary <- function(INPUT, X_POSITION, Y_POSITION,
                        ANNO_COLUMN, CELL_ID_COLUMN, CELL_ID_PREFIX,
                        ANNO_RANGE = c(0, 1), ANNO_MIDPOINT = 0.5,
                        NEIGHBOR_METHOD = c("radius", "knn", "hybrid"),
                        RADIUS, KNN_K,
                        DENOISE             = TRUE,
                        NEST_MIN_SIZE       = 5,
                        NEST_SPECIFICITY    = 0.25,
                        BOUNDARY_SPECIFICITY = 0.05,
                        STRATIFY_BOUNDARY   = TRUE,
                        VERBOSE             = FALSE) {
  
  # 1. Input validation and data cleaning ----
  if (missing(INPUT))        stop("INPUT data frame must be provided")
  if (!is.data.frame(INPUT)) stop("INPUT must be a data frame")
  if (missing(X_POSITION) || missing(Y_POSITION))
    stop("Both X_POSITION and Y_POSITION column names must be provided")
  if (!X_POSITION %in% names(INPUT))
    stop("X_POSITION: `", X_POSITION, "` not found in INPUT")
  if (!Y_POSITION %in% names(INPUT))
    stop("Y_POSITION: `", Y_POSITION, "` not found in INPUT")
  
  # Handle Cell_ID creation or validation
  if (missing(CELL_ID_COLUMN)) {
    if (VERBOSE) message("Creating Cell_ID...")
    INPUT <- if (missing(CELL_ID_PREFIX)) {
      INPUT %>% dplyr::mutate(Cell_ID = dplyr::row_number())
    } else {
      INPUT %>% dplyr::mutate(Cell_ID = paste0(CELL_ID_PREFIX, "_", dplyr::row_number()))
    }
    CELL_ID_COLUMN <- "Cell_ID"
  } else {
    if (!CELL_ID_COLUMN %in% names(INPUT)) stop("Specified CELL_ID_COLUMN not found in INPUT")
    if (any(duplicated(INPUT[[CELL_ID_COLUMN]]))) stop("CELL_ID_COLUMN contains duplicate IDs. IDs must be unique.")
  }
  
  # Check for duplicate coordinates
  if (any(duplicated(INPUT[, c(X_POSITION, Y_POSITION)]))) {
    stop("Duplicate coordinates (X, Y) detected. Please remove overlapping cells first.")
  }
  
  # Clean existing Synora columns to prevent join conflicts
  synora_cols <- c("Nb_Count", "Anno_Midpoint", "Mixedness", "Orientedness", "BoundaryScore", "SynoraAnnotation")
  existing_synora <- intersect(synora_cols, names(INPUT))
  if (length(existing_synora) > 0) {
    warning("Existing Synora result columns detected: [", paste(existing_synora, collapse = ", "), 
            "]. Overwriting these columns. Check if you have already run GetBoundary on this dataset.")
    INPUT <- INPUT %>% dplyr::select(-dplyr::all_of(existing_synora))
  }
  
  if (missing(ANNO_COLUMN)) stop("ANNO_COLUMN must be provided")
  if (!ANNO_COLUMN %in% names(INPUT)) stop("ANNO_COLUMN: `", ANNO_COLUMN, "` not found in INPUT")
  if (!is.numeric(INPUT[[ANNO_COLUMN]])) stop("ANNO_COLUMN must be numeric.")
  if (any(!is.finite(INPUT[[ANNO_COLUMN]]))) stop("ANNO_COLUMN contains NA, NaN, or Inf values.")
  
  NEIGHBOR_METHOD <- match.arg(NEIGHBOR_METHOD)
  if (NEIGHBOR_METHOD %in% c("radius", "hybrid") && missing(RADIUS)) 
    stop("RADIUS required")
  if (NEIGHBOR_METHOD %in% c("radius", "hybrid") && (!is.numeric(RADIUS) || RADIUS <= 0))
    stop("RADIUS must be a positive numeric")
  
  if (NEIGHBOR_METHOD %in% c("knn", "hybrid") && missing(KNN_K)) 
    stop("KNN_K required")
  if (NEIGHBOR_METHOD %in% c("knn", "hybrid") &&
      (!is.numeric(KNN_K) || KNN_K <= 0 || KNN_K != as.integer(KNN_K)))
    stop("KNN_K must be a positive integer")
  
  # 2. Preparation ----
  INPUT <- INPUT %>% dplyr::arrange(!!as.name(CELL_ID_COLUMN))
  
  if (identical(ANNO_RANGE, "auto")) {
    ANNO_RANGE <- range(INPUT[[ANNO_COLUMN]], na.rm = TRUE)
    if (diff(ANNO_RANGE) == 0) ANNO_RANGE <- ANNO_RANGE + c(-1e-6, 1e-6)
  }
  
  if (!(is.numeric(ANNO_MIDPOINT) && length(ANNO_MIDPOINT) == 1) &&
      !identical(ANNO_MIDPOINT, "auto"))
    stop("ANNO_MIDPOINT must be either a single numeric or 'auto'")
  
  if (identical(ANNO_MIDPOINT, "auto")) {
    ANNO_MIDPOINT <- .OtsuThreshold(INPUT[[ANNO_COLUMN]])
    if (VERBOSE) message("Otsu ANNO_MIDPOINT: ", round(ANNO_MIDPOINT, 3))
  }
  
  ANNO_SCALED <- paste0(ANNO_COLUMN, "_scaled")
  INPUT <- INPUT %>%
    dplyr::mutate(!!as.name(ANNO_SCALED) := ifelse(
      !!as.name(ANNO_COLUMN) >= ANNO_MIDPOINT,
      (!!as.name(ANNO_COLUMN) - ANNO_MIDPOINT) / (ANNO_RANGE[2] - ANNO_MIDPOINT),
      (!!as.name(ANNO_COLUMN) - ANNO_MIDPOINT) / (ANNO_MIDPOINT - ANNO_RANGE[1])
    ))
  
  NN <- .GetNeighbors(
    COORDS          = INPUT %>% dplyr::select(dplyr::all_of(c(X_POSITION, Y_POSITION))),
    NEIGHBOR_METHOD = NEIGHBOR_METHOD,
    RADIUS          = if (NEIGHBOR_METHOD %in% c("radius", "hybrid")) RADIUS else NULL,
    KNN_K           = if (NEIGHBOR_METHOD %in% c("knn",    "hybrid")) KNN_K  else NULL
  )
  
  ANNO_DENOISED <- paste0(ANNO_COLUMN, "_denoised")
  
  # 3. Core Logic ----
  if (DENOISE) {
    RESULT_1 <- INPUT %>%
      .GetMO(CELL_ID_COLUMN, X_POSITION, Y_POSITION, ANNO_SCALED, NN, FALSE) %>%
      dplyr::mutate(Nest = (0.5 * Mean_Anno + 0.5) >= NEST_SPECIFICITY)
    
    NOISE_COUNT <- c(); i <- 1
    repeat {
      RESULT_1 <- RESULT_1 %>%
        dplyr::nest_by(Nest, NeighborhoodRadius) %>%
        dplyr::mutate(!!as.name(paste0("Noise_", i)) := data %>%
                        dplyr::select(dplyr::all_of(c(X_POSITION, Y_POSITION))) %>%
                        dbscan::dbscan(eps = NeighborhoodRadius, minPts = NEST_MIN_SIZE) %>%
                        .$cluster %>% { . == 0 } %>% list()) %>%
        tidyr::unnest(cols = c(data, !!as.name(paste0("Noise_", i)))) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(!!as.name(CELL_ID_COLUMN)) %>% 
        dplyr::mutate(Nest = ifelse(!!as.name(paste0("Noise_", i)), !Nest, Nest))
      
      NOISE_COUNT <- c(NOISE_COUNT, sum(RESULT_1[[paste0("Noise_", i)]]))
      if (ifelse(length(NOISE_COUNT) == 1, NOISE_COUNT == 0, diff(tail(NOISE_COUNT, 2)) == 0)) {
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
        BoundaryScore = Mixedness * Orientedness,
        SynoraAnnotation = dplyr::case_when(
          BoundaryScore >= BOUNDARY_SPECIFICITY                                  ~ "Boundary",
          BoundaryScore <  BOUNDARY_SPECIFICITY & !!as.name(ANNO_DENOISED) ==  1 ~ "Nest",
          BoundaryScore <  BOUNDARY_SPECIFICITY & !!as.name(ANNO_DENOISED) == -1 ~ "Outside",
          TRUE                                                                   ~ "Noise"
        ) %>% factor(levels = c("Boundary", "Nest", "Outside", "Noise"))
      )
  } else {
    
    RESULT <- INPUT %>%
      .GetMO(CELL_ID_COLUMN, X_POSITION, Y_POSITION, ANNO_SCALED, NN, TRUE) %>%
      dplyr::mutate(
        BoundaryScore = Mixedness * Orientedness,
        SynoraAnnotation = dplyr::case_when(
          BoundaryScore >= BOUNDARY_SPECIFICITY        ~ "Boundary",
          (0.5 * Mean_Anno + 0.5) >= NEST_SPECIFICITY ~ "Nest",
          TRUE                                         ~ "Outside"
        ) %>% factor(levels = c("Boundary", "Nest", "Outside", "Noise"))
      )
    
  }
  
  # 4. Final Assembly ----
  RESULT <- dplyr::left_join(INPUT, RESULT, by = c(CELL_ID_COLUMN, X_POSITION, Y_POSITION)) %>%
    dplyr::transmute(
      !!as.name(CELL_ID_COLUMN),
      !!as.name(X_POSITION),
      !!as.name(Y_POSITION),
      !!as.name(ANNO_COLUMN),
      Nb_Count,
      Anno_Midpoint = ANNO_MIDPOINT,
      Mixedness     = ifelse(is.na(Mixedness),    0, Mixedness),
      Orientedness  = ifelse(is.na(Orientedness), 0, Orientedness),
      BoundaryScore = ifelse(is.na(BoundaryScore), 0, BoundaryScore),
      SynoraAnnotation
    )
  
  if (STRATIFY_BOUNDARY) {
    RESULT <- .StratifyBoundary(RESULT = RESULT, NN = NN)
  }
  
  return(RESULT)
}


# Internal Helpers ----
.OtsuThreshold <- function(data, n_bins = 256) {
  tryCatch({
    data <- data[is.finite(data)]
    if (length(unique(data)) < 2) {
      message("Insufficient variation for Otsu thresholding, returning median")
      return(median(data))
    }
    
    breaks    <- seq(min(data), max(data), length.out = n_bins + 1)
    h         <- graphics::hist(data, breaks = breaks, plot = FALSE)
    bin_mids  <- h$mids
    p         <- h$counts / sum(h$counts)
    
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
                          RADIUS = NULL,
                          KNN_K  = NULL) {
  NEIGHBOR_METHOD <- match.arg(NEIGHBOR_METHOD)
  COORDS <- as.data.frame(COORDS)
  
  if (NEIGHBOR_METHOD == "radius") {
    nn     <- dbscan::frNN(COORDS, eps = RADIUS, sort = FALSE)
    id     <- nn$id
    dist   <- nn$dist
    radius <- RADIUS
    
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
  
  return(list(id = id, dist = dist, radius = radius))
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
      dplyr::mutate(ID   = NN$id,
                    Dist = NN$dist) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        X_u      = list((.$X[ID] - X) / Dist),
        Y_u      = list((.$Y[ID] - Y) / Dist),
        Nb_Anno  = list(.$Anno[ID])
      ) %>%
      dplyr::mutate(
        Nb_Count    = length(Nb_Anno),
        Mean_Anno   = mean(Nb_Anno),
        X_u_SumBg   = sum(X_u),
        Y_u_SumBg   = sum(Y_u),
        X_u_Sum     = sum(Nb_Anno * X_u),
        Y_u_Sum     = sum(Nb_Anno * Y_u)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        Mixedness    = 1 - (Mean_Anno)^2,
        Orientedness = ((sqrt(X_u_Sum^2 + Y_u_Sum^2) -
                           sqrt(X_u_SumBg^2 + Y_u_SumBg^2)) / Nb_Count),
        Orientedness = ifelse(Orientedness < 0, 0, Orientedness),
        NeighborhoodRadius = NeighborhoodRadius
      ) %>%
      dplyr::select(Cell_ID, X, Y, Anno, Nb_Count,
                    NeighborhoodRadius, Mean_Anno, Mixedness, Orientedness)
    
  } else {
    RESULT <- INPUT_SLIM %>%
      dplyr::mutate(ID = NN$id) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Nb_Anno = list(.$Anno[ID])) %>%
      dplyr::mutate(
        Nb_Count  = length(Nb_Anno),
        Mean_Anno = mean(Nb_Anno)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        Mixedness          = 1 - (Mean_Anno)^2,
        NeighborhoodRadius = NeighborhoodRadius
      ) %>%
      dplyr::select(Cell_ID, X, Y, Anno, Nb_Count,
                    NeighborhoodRadius, Mean_Anno, Mixedness)
  }
  
  RESULT <- RESULT %>%
    dplyr::rename(
      !!as.name(CELL_ID_COLUMN) := Cell_ID,
      !!as.name(X_POSITION)     := X,
      !!as.name(Y_POSITION)     := Y,
      !!as.name(ANNO_COLUMN)    := Anno
    )
  
  return(RESULT)
}

.StratifyBoundary <- function(RESULT, NN) {
  anno_vec     <- as.character(RESULT$SynoraAnnotation)
  boundary_idx <- which(anno_vec == "Boundary")
  new_levels   <- c("Boundary_Inner", "Boundary_Outer", "Nest", "Outside", "Noise")
  
  if (length(boundary_idx) == 0) {
    RESULT <- RESULT %>%
      dplyr::mutate(SynoraAnnotation = factor(anno_vec, levels = new_levels))
    return(RESULT)
  }
  
  new_labels <- vapply(boundary_idx, function(i) {
    nb_anno   <- anno_vec[NN$id[[i]]]
    n_nest    <- sum(nb_anno == "Nest",    na.rm = TRUE)
    n_outside <- sum(nb_anno == "Outside", na.rm = TRUE)
    if (n_nest >= n_outside) "Boundary_Inner" else "Boundary_Outer"
  }, character(1))
  
  anno_vec[boundary_idx] <- new_labels
  RESULT <- RESULT %>%
    dplyr::mutate(SynoraAnnotation = factor(anno_vec, levels = new_levels))
  return(RESULT)
}

