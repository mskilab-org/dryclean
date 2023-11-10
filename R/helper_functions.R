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
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @importFrom stats dist
#' @importFrom stats median
#' @importFrom stats na.omit
#' @importFrom MASS ginv
#' @importFrom utils globalVariables
#' @importFrom parallel mclapply
#' @import gUtils
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

thresh <- function(x, mu){
  y = pmax(x - mu, 0, na.rm = TRUE)
  y = y + pmin(x + mu, 0, na.rm = TRUE)
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

apg_project <- function(m.vec, U, lambda1, lambda2){
  q = dim(U)[1]
  p = dim(U)[2]
  v = matrix(0, p, 1)
  s = matrix(0, q, 1)
  I = diag(p)
  converged = FALSE
  k = 0
  maxiter = 200
  UUt = (ginv(crossprod(U) + lambda1*I) %*% t(U))
  while (converged == FALSE){
    k =  k + 1
    v.old = v
    v = UUt %*% (m.vec - s)
    s.old = s
    s = thresh(m.vec - (U %*% v), lambda2)
    e = max(norm(v - v.old, "F"), norm(s - s.old, "F"))/q
    if (e < 1e-6 || k > maxiter){
      converged = TRUE
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

update_cols <- function(U, A, B, lambda1){
  r = dim(U)[2]
  A = A + (lambda1 * diag(r))
  for (i in 1:r){
    bi = B[, i]
    ui = U[, i]
    ai = A[, i]
    numerator = ((bi - U %*% ai)/ A[i, i]) + ui
    U[, i] = numerator / max(norm(numerator, "2"), 1)
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

wash_cycle <- function(m.vec, L.burnin, S.burnin , r, N, U.hat, V.hat, sigma.hat, lambda1 = NA, lambda2 = NA, verbose = TRUE){
  
  if (verbose == TRUE){
    message("Using the detergent provided to start washing")
  }
  
  L.burnin = L.burnin
  r = r
  S.burnin = S.burnin
  U.hat = U.hat
  V.hat = V.hat
  sigma.hat = sigma.hat
  m = dim(L.burnin)[1]
  n = dim(L.burnin)[2]
  
  if (is.na(lambda1)){
    
    lambda1 = 1/(sqrt(m))
    lambda2 = 1/(sqrt(m))
    
    if (verbose == TRUE){
      message("lambdas calculated")
    }
  }
  
  U = U.hat %*% sqrt(diag(sigma.hat))
  ## A = matrix(0, r, r)
  ## B = matrix(0, m, r)
  
  if (verbose == TRUE){
    message("calculating A and B")
  }
  
  #for (i in 1:dim(V.hat)[2]){
  #    A = A + outer(V.hat[, i], V.hat[, i])
  #    B = B + outer((L.burnin[, i] - S.burnin[, i]), V.hat[, i])
  #}
  
  if (verbose == TRUE){
    message("calculating v and s")
  }
  
  projection = apg_project(m.vec, U, lambda1, lambda2)
  vi = as.numeric(projection[[1]])
  si = as.numeric(projection[[2]])
  vi.subtract = V.hat[, 1]
  V.hat = cbind(V.hat, vi)
  #A = A + outer(vi, vi) - outer(vi.subtract, vi.subtract)
  #B = B + outer((as.numeric(m.vec) - si), vi) - outer((L.burnin[, 1] - S.burnin[, 1]), vi.subtract)
  
  if (verbose == TRUE){
    message("Calculating b")
  }
  
  ## U = update_cols(U, A, B, lambda1)
  L.vec = U %*% vi
  
  return(list(L.vec, si))
}


##############################
## prep_cov
##############################
#' @name prep_cov
#'
#' @title function to prepare covearge data for decomposition 
#' @description prepares the GC corrected coverage data for decomposition
#' 
#' @keywords internal
#' 
#' @param m.vec GRanges object. GRanges containing GC corrected reads as a cloumn in metadata
#' 
#' @param blacklist, boolean (default == FALSE). Whether to exclude off-target markers in case of Exomes or targeted sequqnecing. If set to TRUE, needs a GRange marking if each marker is set to be excluded or not.
#'
#' @param burnin.samples.path, character (default = NA). Path to balcklist markers file
#'
#' @return vector of length m with processed coverage data
#' 
#' @author Aditya Deshpande

prep_cov <- function(m.vec = m.vec, blacklist = FALSE, blacklist_path = NA){
  m.vec = gr2dt(m.vec)
  m.vec = m.vec[, .(seqnames, signal)]
  m.vec[, median.chr := median(.SD$signal, na.rm = T), by = seqnames]
  m.vec[is.na(signal), signal := median.chr]
  min.cov = min(m.vec[signal > 0]$signal, na.rm = T)
  m.vec[is.infinite(signal), signal := min.cov]
  m.vec[signal == 0, signal := min.cov]
  m.vec[signal < 0, signal := min.cov]
  
  if (blacklist){
    blacklist.pon =  readRDS(blacklist_path)
    m.vec$blacklisted = blacklist.pon$blacklisted
    m.vec[blacklisted == TRUE, signal := NA]
    m.vec = na.omit(m.vec)
  }
  
  m.vec[, signal := log(signal)]
  
  return(m.vec)
}


# collapse_cov <- function(cov.gr, bin.size = target_resolution, this.field = field){
#   BINSIZE.ROUGH = bin.size
#   cov.gr = cov.gr[, this.field]
#   cov.gr = gr2dt(cov.gr)
#   setnames(cov.gr, this.field, "signal")
#   cov.gr = cov.gr[!is.infinite(signal), .(signal = median(signal, na.rm = TRUE)),
#                   by = .(seqnames, start = floor(start/BINSIZE.ROUGH)*BINSIZE.ROUGH+1)]
#   cov.gr[, end := (start + BINSIZE.ROUGH) - 1]
#   setnames(cov.gr, "signal", this.field)
#   cov.gr = dt2gr(cov.gr)
#   return(cov.gr)
# }

collapse_cov <- function(cov.gr, bin.size = target_resolution, this.field = field){
  BINSIZE.ROUGH = bin.size
  cov.gr = cov.gr[, this.field]
  cov.gr = gr2dt(cov.gr)
  setnames(cov.gr, this.field, "signal")
  cov_og = cov.gr
  cov.gr = cov.gr[!is.infinite(signal), .(signal = median(signal, na.rm = TRUE)),
                  by = .(seqnames, start = floor(start/BINSIZE.ROUGH)*BINSIZE.ROUGH+1)]  
  cov.gr[, end := (start + BINSIZE.ROUGH) - 1]
  s <- cov_og %>% group_by(seqnames) %>% summarize(end = max(end)) %>% data.table() %>% filter(end %% BINSIZE.ROUGH == 0) %>% select(seqnames)  
  for (chr in s$seqnames){
    cov_fixed <- cov.gr %>% filter(seqnames == chr) %>% select(seqnames, start, end, signal)
    cov_fixed <- cov_fixed[nrow(cov_fixed)]
    cov_fixed <- cov_fixed %>%
      mutate(start = end + 1, end = end + BINSIZE.ROUGH, signal = NA)
    cov.gr <- rbind(cov.gr,cov_fixed)
    cov.gr <- cov.gr %>% arrange(seqnames, start)
  }
  setnames(cov.gr, "signal", this.field)
  cov.gr = dt2gr(cov.gr)
  return(cov.gr)
}



generate_template <- function(cov, wgs = wgs, target_resolution = target_resolution, this.field = field, nochr = TRUE, all.chr = c(as.character(1:22), "X")){
  if (wgs){
    inferred_resolution = median(width(cov), na.rm = T)
    if (target_resolution > inferred_resolution){
      cov = collapse_cov(cov, this.field = this.field, bin.size = target_resolution)
    }
  }
  cov = sortSeqlevels(cov)
  cov = sort(cov)
  if (nochr) {
    cov = gUtils::gr.nochr(cov)
  } 
  cov = cov %Q% (seqnames %in% all.chr)
  template = cov[, c()]
  return(template)
}



