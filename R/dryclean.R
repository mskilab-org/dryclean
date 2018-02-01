#' @import GenomicRanges
#' @import data.table
#' @import rsvd
#' @import gUtils
#' @import Matrix
#' @import stats
#' @import MASS


##############################
## thresh
##############################
#' @name thresh
#'
#' @title Internal shrinkage function for values in S vector
#' @description Internal shrinkage function for values in S vector:
#' y = sgn(x)max(|x| - mu, 0)
#' @keywords internal
#' @param x numeric vector of length == number of markers in the sample. Vector to which shrinkage is to be applied.
#' @param mu integer. Shrinkage parameter to use (default = lambda2 = 1/(sqrt(no. of markers in input vector)))
#' @return numeric vector. Regularized vector y
#' @author Aditya Deshpande

thresh = function(x, mu){
    y = pmax(x - mu, 0, na.rm = TRUE)
    y = y + pmin(x + mu, 0, na.rm = TRUE)
    return(y)}


##############################
## apg_project
##############################
#' @name apg_progect
#'
#' @title Accelerated Proximal Gradient based projection method
#' @description project new sample into the burnin space
#' solving projection by Accelerated Proximal Gradient
#' min_{v, s} 0.5*|m-Uv-s|_2^2 + 0.5*lambda1*|v|^2 + lambda2*|s|_1
#' @keywords internal
#' @param m.vec numeric vector of length m. Vector of GC corrected coverage data of sample in question
#' @param U (m  markers x n samples) numeric matrix. The basis of low rank subspace. The dimensions are same as burnin matrix
#' @param lambda1 integer. Tuning parameter (default = 1/(sqrt(no. of markers in input vector))
#' @param lambda2 integer. Tuning parameter (default = 1/(sqrt(no. of markers in input vector))
#' @return list with s and v vectors 
#' @author Aditya Deshpande

apg_project = function(m.vec, U, lambda1, lambda2){
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
        ### message("on ", k)
        v.old = v
        v = UUt %*% (m.vec - s)
        s.old = s
        s = thresh(m.vec - (U %*% v), lambda2)
        ### print(norm(v - v.old, "2"))
        ### print(norm(s - s.old, "2"))
        e = max(norm(v - v.old, "2"), norm(s - s.old, "2"))/q
        if (e < 1e-6 || k > maxiter){
            converged = TRUE
        }
    }
    return(list(v, s))}

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

update_cols = function(U, A, B, lambda1){
    r = dim(U)[2]
    A = A + (lambda1 * diag(r))
    for (i in 1:r){
        bi = B[, i]
        ui = U[, i]
        ai = A[, i]
        numerator = ((bi - U %*% ai)/ A[i, i]) + ui
        U[, i] = numerator / max(norm(numerator, "2"), 1)
    }
    return(U)}


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
#' @param lambda1, lambda2 integer. Tuning parameter (default = 1/(sqrt(no. of markers in input vector))
#' @param L.burnin (m  markers x n samples) numeric matrix. L matrix of panel of normals after batch rPCA decomposition 
#' @param S.burnin (m  markers x n samples) numeric matrix. S matrix of panel of normals after batch rPCA decomposition 
#' @param r integer. Estimated rank of panel of normals after batch rPCA decomposition
#' @param N (m  markers x n samples) numeric matrix. Matrix of GC corrected coverage data of panel of normals for batch rPCA decomposition
#' @param U.hat (m  markers x n samples) numeric matrix. Right singular matrix of L.burnin
#' @param V.hat (m  markers x n samples) numeric matrix. Left singular matrix of L.burnin
#' @param sigma.hat numeric vector of length r. Singular values of L.burnin 
#' @param decomp boolean (default = FALSE). If TRUE, carry batch rPCA decomposition on N matrix
#' @export
#' @return S and L numeric vectors each of length m. Vectors for sample in question
#' @author Aditya Deshpande

wash_cycle = function(m.vec, L.burnin, S.burnin , r, N, U.hat, V.hat, sigma.hat, lambda1 = NA, lambda2 = NA, decomp = FALSE, verbose = TRUE){
    if (decomp == TRUE){
        if (verbose == TRUE){message("Starting decomposition, this will take a while!")}
        rpca.N = rrpca(N, trace = TRUE)
        L.burnin = rpca.N$L
        r = rpca.N$k
        S.burnin = rpca.N$S
    }
    else{
        if (verbose == TRUE){message("Using default set of normals a.k.a solvent")}
        L.burnin = L.burnin
        r = r
        S.burnin = S.burnin
        U.hat = U.hat
        V.hat = V.hat
        sigma.hat = sigma.hat
    }
    m = dim(L.burnin)[1]
    n = dim(L.burnin)[2]
    if (is.na(lambda1)){
        lambda1 = 1/(sqrt(m))
        lambda2 = 1/(sqrt(m))
        if (verbose == TRUE){message("lambdas calculated")}
    }
    ## initialization
    if (verbose == TRUE){message("And thus we begin the 'wash cycle'")}
    U = U.hat %*% sqrt(diag(sigma.hat))
    A = matrix(0, r, r)
    B = matrix(0, m, r)
    if (verbose == TRUE){message("calculating A and B")}

    for (i in 1:dim(V.hat)[2]){
        A = A + outer(V.hat[, i], V.hat[, i])
        B = B + outer((L.burnin[, i] - S.burnin[, i]), V.hat[, i])
    }

    if (verbose == TRUE){message("calculating v and s")}
    projection = apg_project(m.vec, U, lambda1, lambda2)
    vi = as.numeric(projection[[1]])
    si = as.numeric(projection[[2]])
    vi.subtract = V.hat[, 1]
    V.hat = cbind(V.hat, vi)
    A = A + outer(vi, vi) - outer(vi.subtract, vi.subtract)
    B = B + outer((as.numeric(m.vec) - si), vi) - outer((L.burnin[, 1] - S.burnin[, 1]), vi.subtract)
    if (verbose == TRUE){message("Updating subspace")}
    U = update_cols(U, A, B, lambda1)
    L.vec = U %*% vi
    return(list(L.vec, si))}


##############################
## prep_cov
##############################
#' @name prep_cov
#'
#' @title function to prepare covearge data for decomposition 
#' @description preapares the GC corrected coverage data for decomposition
#'
#' 
#' @keywords internal
#' @param m.vec GRanges object. GRanges containing GC corrected reads as a cloumn in metadata
#' @return vector of length m with processed coverage data
#' @author Aditya Deshpande

prep_cov = function(m.vec = m.vec){
    m.vec = gr2dt(m.vec)
    ##m.vec = as.data.table(m.vec)
    m.vec = m.vec[, .(seqnames, reads.corrected)]
    m.vec[, reads.corrected := log(reads.corrected)]
    invisible(lapply(names(m.vec),function(.name) set(m.vec, which(is.infinite(m.vec[[.name]])), j = .name,value = 0)))
    invisible(lapply(names(m.vec),function(.name) set(m.vec, which(is.na(m.vec[[.name]])), j = .name,value = 0)))
    invisible(lapply(names(m.vec),function(.name) set(m.vec, which(is.nan(m.vec[[.name]])), j = .name,value = 0)))
    return(m.vec)}


##############################
## start_wash_cycle
##############################
#' @name start_wash_cycle
#'
#' @title function begins the online rPCA process and parallelizes the chromosmes for speed.  
#' @description function begins the online rPCA process and parallelizes the chromosmes for speed.  
#' It is the wrapper that takes in GRanges and outputs GRanges with decomposition 
#' 
#' @param cov GRanges object containig the GC corrected cov data outputed from fragCounter. Needs metadata with header "reads.corrected"
#' @param mc.cores interger. Number of cores to use for parallelization
#' @param burnin.samples.path string. Path to burnin samples for each chromosome
#' @param whole_genome boolean (default = TRUE). For processing chromosome or whole genome
#' @param chr string. if a single chromosome is to be processed, name of chromosome. Requires whole_genome = FALSE
#' @export
#' @return GRange object with decomposition
#' @author Aditya Deshpande

start_wash_cycle = function(cov, mc.cores = 1, burnin.samples.path = NA, verbose = TRUE, chr = NA, whole_genome = TRUE){
    if (whole_genome == TRUE){
        chr.nm = c(as.character(1:22), "X", "Y")}
    else {chr.nm = chr}
    ##print(chr.nm)
    chr.list = mclapply(chr.nm, function(x){
        if(verbose == TRUE){
            message("Loading burnins")}
        if (is.na(burnin.samples.path)){
            stop('Need burnin files to procced')}
        rpca.1 = readRDS(paste0(burnin.samples.path, "/rpca.burnin.chr", x, ".rds"))
        if(verbose == TRUE){
            message(paste0("Let's begin, this is chr", x))}
        m.vec = prep.cov(cov)
        m.vec = m.vec[seqnames == x]
        m.vec = as.matrix(m.vec$reads.corrected)
        L.burnin = rpca.1$L
        S.burnin = rpca.1$S
        r = rpca.1$k
        U.hat = rpca.1$U.hat
        V.hat = rpca.1$V.hat
        sigma.hat = rpca.1$sigma.hat
        decomposed = dryclean::wash_cycle(m.vec = m.vec, L.burnin = L.burnin, S.burnin = S.burnin, r = r, U.hat = U.hat, V.hat = V.hat, sigma.hat = sigma.hat)
        if (verbose == TRUE){
            message("combining matrices with gRanges")}
        cov = cov[seqnames == x]
        cov = cbind(decomposed[[2]], cov)
        colnames(cov)[1] = 'S'
        cov[reads.corrected == 0, S := 0]
        cov[is.na(reads.corrected), S := NA]
        cov[, S1 := exp(S)]
        cov = cbind(decomposed[[1]], cov)
        colnames(cov)[1] = 'L'
        cov[reads.corrected == 0, L := 0]
        cov[is.na(reads.corrected), L := NA]
        cov[, L1 := exp(L) ]
        cov[, log.reads := log(reads.corrected)]
        cov[is.infinite(log.reads), log.reads := NA]
        return(cov)}, mc.cores = mc.cores, mc.preschedule = FALSE)

    results = rbindlist(chr.list, fill = T)
    results = dt2gr(results)
    return(results)}


    
