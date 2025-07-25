#' @import GenomicRanges
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#' @importFrom data.table set
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table setkeyv
#' @importFrom data.table setnames
#' @importFrom data.table transpose
#' @import rsvd
#' @importFrom magrittr %>%
#' @importFrom gUtils gr2dt
#' @importFrom gUtils dt2gr
#' @importFrom gUtils gr.nochr
#' @importFrom gUtils %Q%
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @importFrom stats dist
#' @importFrom stats median
#' @importFrom stats na.omit
#' @importFrom MASS ginv
#' @importFrom utils globalVariables
#' @importFrom parallel mclapply
#' @import pbmcapply


globalVariables(c("target_resolution", "sortSeqlevels", "field", ".", "..ix", "L", "L1", "V1", "black_list_pct", "blacklisted", "decomposed_cov", "germline.status", "log.reads", "mclapply", "median.chr", "normal_cov", "foreground", "input.read.counts", "foreground.log", "reads.corrected", "background", "file.available", "mt", "background.log", ".N", ".SD", ":=", "median.idx", ".GRP", "reads.corrected.org", "%>%", "signal", "signal.org"))


##############################
## thresh
##############################
#' @name thresh
#'
#' @title This is second phase in preparing the Panel of Normals
#' @description Internal shrinkage function for values in S vector:
#' y = sgn(x)max(|x| - mu, 0)
#' @keywords internal
#' @param x numeric vector of length == number of markers in the sample. Vector to which shrinkage is to be applied.
#' @param mu integer. Shrinkage parameter to use (default = lambda2 = 1/(sqrt(no. of markers in input vector)))
#' @return numeric vector. Regularized vector y
#' @author Aditya Deshpande

thresh <- function(x, mu) {
  y <- pmax(x - mu, 0, na.rm = TRUE)
  y <- y + pmin(x + mu, 0, na.rm = TRUE)
  return(y)
}


##############################
## apg_project
##############################
#' @name apg_progect
#'
#' @title Accelerated Proximal Gradient based projection method
#' @description project new sample into the burnin space
#' solving projection by Accelerated Proximal Gradient
#'
#' @keywords internal
#' @param m.vec numeric vector of length m. Vector of GC corrected coverage data of sample in question
#' @param U (m  markers x n samples) numeric matrix. The basis of low rank subspace. The dimensions are same as burnin matrix
#' @param lambda1 integer. Tuning parameter (default = 1/(sqrt(no. of markers in input vector))
#' @param lambda2 integer. Tuning parameter (default = 1/(sqrt(no. of markers in input vector))
#' @return list with s and v vectors
#' @author Aditya Deshpande

apg_project <- function(m.vec, U, lambda1, lambda2) {
  q <- dim(U)[1]
  p <- dim(U)[2]
  v <- matrix(0, p, 1)
  s <- matrix(0, q, 1)
  I <- diag(p)
  converged <- FALSE
  k <- 0
  maxiter <- 200
  UUt <- (ginv(crossprod(U) + lambda1 * I) %*% t(U))
  while (converged == FALSE) {
    k <- k + 1
    v.old <- v
    v <- UUt %*% (m.vec - s)
    s.old <- s
    s <- thresh(m.vec - (U %*% v), lambda2)
    e <- max(norm(v - v.old, "F"), norm(s - s.old, "F")) / q
    if (e < 1e-6 || k > maxiter) {
      converged <- TRUE
    }
  }
  return(list(v, s))
}



##############################
## update_cols
##############################
#' @name update_cols
#'
#' @title subspace update based on projection of sample in question
#' @description Internal function that updates basis of low rank subspace after
#' projecting sample in question by warm restart
#'
#' @keywords internal
#' @param U (m  markers x n samples) numeric matrix. The basis of low rank subspace. The dimensions are same as burnin matrix
#' @param lambda1 integer. Tuning parameter (default = 1/(sqrt(no. of markers in input vector))
#' @param A, B m x n numeric matrices. Matrices to update subspace basis
#'
#' @return U (m  markers x n samples) numeric matrix. updated basis U of subspace
#' @author Aditya Deshpande

update_cols <- function(U, A, B, lambda1) {
  r <- dim(U)[2]
  A <- A + (lambda1 * diag(r))
  for (i in 1:r) {
    bi <- B[, i]
    ui <- U[, i]
    ai <- A[, i]
    numerator <- ((bi - U %*% ai) / A[i, i]) + ui
    U[, i] <- numerator / max(norm(numerator, "2"), 1)
  }
  return(U)
}


##############################
## wash_cycle
##############################
#' @name wash_cycle
#'
#' @title function to perform online rPCA descomposition
#' @description Function that calls the online rPCA on sample under question and updates subspace
#' basis and does decomposition
#'
#'
#' @param m.vec numeric vector of length m. Vector of GC corrected coverage data of sample in question
#'
#' @param lambda1, integer. Tuning parameter (default = 1/(sqrt(no. of markers in input vector))
#'
#' @param lambda2, integer. Tuning parameter (default = 1/(sqrt(no. of markers in input vector))
#'
#' @param L.burnin (m  markers x n samples) numeric matrix. L matrix of panel of normals after batch rPCA decomposition
#'
#' @param S.burnin (m  markers x n samples) numeric matrix. S matrix of panel of normals after batch rPCA decomposition
#'
#' @param r integer. Estimated rank of panel of normals after batch rPCA decomposition
#'
#' @param N (m  markers x n samples) numeric matrix. Matrix of GC corrected coverage data of panel of normals for batch rPCA decomposition
#'
#' @param U.hat (m  markers x n samples) numeric matrix. Right singular matrix of L.burnin
#'
#' @param V.hat (m  markers x n samples) numeric matrix. Left singular matrix of L.burnin
#'
#' @param sigma.hat numeric vector of length r. Singular values of L.burnin
#'
#' @param verbose boolean (default == TRUE). Outputs progress.
#'
#' @return S and L numeric vectors each of length m. Vectors for sample in question
#' @author Aditya Deshpande

wash_cycle <- function(m.vec, L.burnin, S.burnin, r, N, U.hat, V.hat, sigma.hat, lambda1 = NA, lambda2 = NA, verbose = TRUE) {
  if (verbose == TRUE) {
    message("Using the detergent provided to start washing")
  }

  m <- dim(L.burnin)[1]
  n <- dim(L.burnin)[2]

  if (is.na(lambda1)) {
    lambda1 <- 1 / (sqrt(m))
    lambda2 <- 1 / (sqrt(m))

    if (verbose == TRUE) {
      message("lambdas calculated")
    }
  }

  U <- U.hat %*% sqrt(diag(sigma.hat))
  ## A = matrix(0, r, r)
  ## B = matrix(0, m, r)

  if (verbose == TRUE) {
    message("calculating A and B")
  }

  # for (i in 1:dim(V.hat)[2]){
  #    A = A + outer(V.hat[, i], V.hat[, i])
  #    B = B + outer((L.burnin[, i] - S.burnin[, i]), V.hat[, i])
  # }

  if (verbose == TRUE) {
    message("calculating v and s")
  }

  projection <- apg_project(m.vec, U, lambda1, lambda2)
  vi <- as.numeric(projection[[1]])
  si <- as.numeric(projection[[2]])
  vi.subtract <- V.hat[, 1]
  V.hat <- cbind(V.hat, vi)
  # A = A + outer(vi, vi) - outer(vi.subtract, vi.subtract)
  # B = B + outer((as.numeric(m.vec) - si), vi) - outer((L.burnin[, 1] - S.burnin[, 1]), vi.subtract)

  if (verbose == TRUE) {
    message("Calculating b")
  }

  ## U = update_cols(U, A, B, lambda1)
  L.vec <- U %*% vi

  return(list(L.vec, si))
}


##############################
## prep_cov
##############################
#' @name prep_cov
#' @title function to prepare covearge data for decomposition and perform normalization
#' @description prepares the GC corrected coverage data for decomposition
#' @keywords internal
#' @param m.vec GRanges object. GRanges containing GC corrected reads as a column in metadata
#' @param blacklist, boolean (default == FALSE). Whether to exclude off-target markers in case of Exomes or targeted sequqnecing. If set to TRUE, needs a GRange marking if each marker is set to be excluded or not.
#' @param burnin.samples.path, character (default = NA). Path to balcklist markers file
#' @param center, character (default == FALSE). Transform the data to centered?
#' @param centering, enum of "mean" or "median" ONLY (default = "mean"). Paradigm of centering.
#' @return vector of length m with processed coverage data
#' @author Aditya Deshpande / Johnathan Rafailov

prep_cov <- function(m.vec = m.vec, PAR.file = NA, use.blacklist = FALSE, blacklist = NA, center = FALSE, centering = "mean", build = "hg19") {
  if (is.na(PAR.file)) {
    if (build == "hg38") {
          message("PAR file not provided, using hg38 default. If this is not the correct build, please provide a GRanges object delineating for corresponding build")
          par.path <- system.file("extdata", "PAR_hg38.rds", package = "dryclean")
        } else {
          message("PAR file not provided, using hg19 default. If this is not the correct build, please provide a GRanges object delineating for corresponding build")
          par.path <- system.file("extdata", "PAR_hg19.rds", package = "dryclean")
        }
        par.gr <- readRDS(par.path)
  } else {
        par.gr <- readRDS(PAR.file)
      }

  mcols(m.vec)$signal[which(is.infinite(mcols(m.vec)$signal))] <- NA
  mcols(m.vec)$og.signal <- mcols(m.vec)$signal
  m.vec = gr2dt(m.vec)
  m.vec[, median.idx := .GRP, by = seqnames]
  m.vec$mt <- suppressWarnings(gr.match(dt2gr(m.vec), par.gr))
  m.vec[, median.idx := ifelse(is.na(mt), median.idx, mt + 24)]
  m.vec[, median.chr := median(og.signal, na.rm = T), by = median.idx]

  ## We calculate mean as well to capture Y chr distribution as will become clear below
  m.vec[, mean.chr := mean(og.signal, na.rm = T), by = median.idx]
  
  if (center == T) {
    if (centering == "mean") {
      m.vec$center.all <- mean(m.vec$og.signal, na.rm = T)
    } else {
      m.vec$center.all <- median(m.vec$og.signal, na.rm = T)
    }
  }

  ## TODO:
  ## This threshold is empirical, to distinguish Y chr in male and female
  m.vec[, signal := ifelse(!(seqnames %in% c("X", "Y")), og.signal / center.all,
                       ifelse(mean.chr < 5, 0, og.signal / center.all)
                       )]

  m.vec[, signal := ifelse(is.na(median.chr), 1, signal)]
  m.vec <- m.vec[, .(seqnames, start, end, og.signal, signal, center.all, median.chr)]
  m.vec[is.na(signal), signal := median.chr]
  min.cov <- min(m.vec[signal > 0]$signal, na.rm = T)
  m.vec[is.infinite(signal), signal := min.cov]
  m.vec[signal < 0, signal := min.cov]
  m.vec[signal == 0, signal := min.cov]
  m.vec[, signal := log(signal)]
  
  return(dt2gr(m.vec))
}


collapse_cov <- function(cov.gr, bin.size = target_resolution, this.field = field) {
  BINSIZE.ROUGH <- bin.size
  cov.gr <- cov.gr[, this.field]
  cov.gr <- gr2dt(cov.gr)
  setnames(cov.gr, this.field, "signal")
  cov_og <- cov.gr
  cov.gr <- cov.gr[!is.infinite(signal), .(signal = median(signal, na.rm = TRUE)),
    by = .(seqnames, start = floor(start / BINSIZE.ROUGH) * BINSIZE.ROUGH + 1)
  ]
  cov.gr[, end := (start + BINSIZE.ROUGH) - 1]
  setnames(cov.gr, "signal", this.field)
  cov.gr <- dt2gr(cov.gr)
  return(cov.gr)
}

generate_template <- function(cov, wgs = wgs, target_resolution = target_resolution, this.field = field, nochr = TRUE, all.chr = c(as.character(1:22), "X", "Y")) {
  if (wgs) {
    inferred_resolution <- median(width(cov), na.rm = T)
    if (target_resolution > inferred_resolution) {
      cov <- collapse_cov(cov, this.field = this.field, bin.size = target_resolution)
    }
  }
  cov <- sortSeqlevels(cov)
  cov <- sort(cov)
  if (nochr) {
    cov <- gUtils::gr.nochr(cov)
  }
  cov <- cov %Q% (seqnames %in% all.chr)
  template <- cov[, c()]
  return(template)
}

#' @title validate_pon_paths
#' @details Checks that the provided Panel of Normals (PON) file paths are valid and accessible. This function validates that each file path supplied in the input is a well-formed path and that the file exists. It raises an informative error if any path does not meet the required criteria for subsequent processing. Returns as list of rebinned coverages. The function performs necessary checks for file existence and proper format. Ensure that the provided paths correspond to files that are accessible in the current environment.
#' @name validate_pon_paths
#' @param pon_paths A character vector of file paths to PON files.
#' @param field
#' @param target_resolution
#' @param cores
#' @param all.chr
#' @return A character vector of validated PON file paths.
#' @examples
#' \dontrun{
#'   # Validate PON file paths
#'   valid_paths <- validate_pon_paths(c("/path/to/pon1.txt", "/path/to/pon2.txt"))
#' }
validate_pon_paths <- function(pon_paths, field, target_resolution, cores, all.chr = c(1:22, "X", "Y"), verbose = TRUE) {

  if (verbose) {
    message("Validating PON paths...")
  }

  skip.files <- !file.exists(pon_paths)
  

  if(length(which(skip.files)) > 0) {
    message(sprintf("Number of files that do not exist: %d", sum(skip.files)))
    message(paste0("Skipping files (do not exist): ", paste(pon_paths[skip.files], collapse = ", ")))
  }

  result <- mclapply(pon_paths[!skip.files], function(path) {
    cov.dt <- readRDS(path) %>% gr.nochr() %>% gr2dt()
    bin.size <- median(cov.dt$width, na.rm = TRUE)
    n.bins.chr <- cov.dt[seqnames %in% all.chr, .N, by = seqnames]
    setnames(n.bins.chr, "N", paste0("normal_", match(path, pon_paths)))
    list(bin.size = bin.size, n.bins.chr = n.bins.chr)
  }, mc.cores = cores)

  bin.sizes <- sapply(result, function(x) x$bin.size)
  names(bin.sizes) <- paste0("normal_", seq_along(bin.sizes))

  rebin <- bin.sizes != target_resolution

  if (rebin){
    message("Rebinning")
  }
  
  file.reads <- pbmcapply::pbmcmapply(function(nm, path, rebin){
    cov.gr <- readRDS(path) %>% gr.nochr()
    og.resolution <- median(width(cov.gr), na.rm = TRUE)
    if (rebin) {
      cov.gr <- collapse_cov(cov.gr, this.field = field, bin.size = target_resolution)
    }
    return(cov.gr) 
  }, names(rebin), pon_paths[!skip.files], rebin, mc.cores = cores)

  return(file.reads)

}
