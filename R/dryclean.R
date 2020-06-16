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
#' @importFrom gUtils gr2dt
#' @importFrom gUtils dt2gr
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @importFrom stats dist
#' @importFrom stats median
#' @importFrom stats na.omit
#' @importFrom MASS ginv
#' @importFrom utils globalVariables
#' @import gUtils   

globalVariables(c(".", "..ix", "L", "L1", "V1", "black_list_pct", "blacklisted", "decomposed_cov", "germline.status", "log.reads", "mclapply", "median.chr", "normal_cov", "foreground", "input.read.counts", "foreground.log", "reads.corrected", "background", "background.log", ".N", ".SD", ":=", "median.idx", ".GRP", "reads.corrected.org", "%>%", "signal", "signal.org"))




##############################
## prepare_detergent
##############################
#' @name prepare_detergent
#'
#' @title This is the frist stage and involves preparing the Panel of Normals (PON)
#' that will be used for online decomposition
#' 
#' @description This function takes in gRanges outputs from fragCounter and extracts GC corrected
#' read count data and carries rPCA decomposition on the matrix thus created. The normal samples used to form the PON can be selected randomly or by clustering the genomic bacground or all samples can be used.
#' 
#' @export
#' @param normal.table.path character path to data.table containing two columns "sample" and "normal_cov". See manual for details.
#' 
#' @param use.all boolean (default == TRUE). If all normal samples are to be used for creating PON.
#' 
#' @param choose.randomly boolean (default == FALSE). If a random subset of normal samples are to be used for creating PON.
#' 
#' @param choose.by.clustering boolean (default == FALSE). Clusters normal samples based on the genomic background and takes a random sample from within the clusters.
#' 
#' @param number.of.samples interger (default == 50). If choose.by.clustering == TRUE, this is the number of clusters at which to cut tree.
#'
#' @param tolerance numeric (default == 0.0001). Tolerance for error for batch rPCA. We suggest keeping this value.
#' 
#' @param save.pon boolean (default == FALSE). If PON needs to be saved.
#' 
#' @param path.to.save character (default == NA). Path to save the PON created if save.pon == TRUE.
#' 
#' @param num.cores interger (default == 1). Number of cores to use for parallelization.
#'
#' @param verbose boolean (default == TRUE). Outputs progress.
#'
#' @param is.human boolean (default == TRUE). Organism type.
#' 
#' @param build genome build to define PAR region in chromosome X.
#'
#' @param field character (default == "reads.corrected"). Field to use for processing.
#' 
#' @return \code{prepare_detergent} returns a list containing the following components:
#' 
#'    \item{L}{  array_like; \cr
#'              low-rank component; \eqn{(m, n)} dimensional array.
#'    }
#'    \item{S}{  array_like; \cr
#'              sparse component; \eqn{(m, n)} dimensional array.
#'    }
#'    \item{k}{  numeric; \cr
#'              estimated rank of the matrix; used in the online implementation.
#'    }
#'    \item{U.hat}{  array_like; \cr
#'              left singular vectors of L; \eqn{(m, k)} dimensional array.
#'    }
#'    \item{V.hat}{  array_like; \cr
#'              right singular vectors of L; \eqn{(n, k)} dimensional array.
#'    }
#'    \item{sigma.hat}{  array_like; \cr
#'              singular values; vector of length \eqn{(k)}.
#'    }
#' 
#' @author Aditya Deshpande


prepare_detergent <- function(normal.table.path = NA, use.all = TRUE, choose.randomly = FALSE, choose.by.clustering = FALSE, number.of.samples = 50, save.pon = FALSE, path.to.save = NA, verbose = TRUE, num.cores = 1, tolerance = 0.0001, is.human = TRUE, build = "hg19", field = "reads.corrected"){
    
    if (verbose){
        message("Starting the preparation of Panel of Normal samples a.k.a detergent")
    }

    if (is.na(normal.table.path)){
        stop("Need a table with paths to normal samples to create a PON")
    }

    if (is.na(path.to.save) & save.pon == TRUE){
        stop("Need a path to save decomposed PON")
    }

    normal.table = readRDS(normal.table.path)
    setkeyv(normal.table, "sample")

    num.samp = nrow(normal.table)

    if (verbose){
        message(paste0(num.samp, " samples available"))
    }
        
    if (use.all & choose.randomly | use.all & choose.by.clustering | choose.randomly & choose.by.clustering | use.all & choose.randomly & choose.by.clustering){
        stop("only one of use.all, choose.randomly, choose.by.clustering can be set to TRUE. Rectify and restart")
    }

    if (choose.randomly){
        if (verbose){ 
            message(paste0("Selecting ", number.of.samples, " normal samples randomly"))
        }
        set.seed(12)
        samp.final = sample(1:num.samp, number.of.samples)
        samp.final = normal.table[samp.final]
        setkey(samp.final, "sample")
    }

    if (choose.by.clustering){

        if (verbose){ 
            message("Starting the clustering")
        }
        
        mat.small = mclapply(normal.table[, sample], function(nm){
            this.cov = tryCatch(readRDS(normal.table[nm, normal_cov]), error = function(e) NULL)
            if (!is.null(this.cov)){
                this.cov = this.cov[, field] %>% gr2dt() %>% setnames(., field, "signal")
                reads = this.cov[seqnames == "22", .(seqnames, signal)]
                reads[, median.chr := median(.SD$signal, na.rm = T), by = seqnames]
                reads[is.na(signal), signal := median.chr]
                min.cov = min(reads[signal > 0]$signal, na.rm = T)
                reads[signal == 0, signal := min.cov]
                reads[signal < 0, signal := min.cov]
                reads = log(reads[, .(signal)])
                reads = transpose(reads)
                reads = cbind(reads, nm)
            } else {reads = data.table(NA)}
            return(reads)}, mc.cores = num.cores)
        
        gc()
        
        mat.sub = rbindlist(mat.small, fill = T)
        mat.sub = na.omit(mat.sub)
        ix = ncol(mat.sub)
        samp.names = mat.sub[, ..ix]
        mat.sub.t = transpose(mat.sub[, 1:(ncol(mat.sub) - 1)])
        rm(mat.sub)
        gc()

        if (verbose){ 
            message("Starting decomposition on a small section of genome")
        }

        mat.sub.t = as.matrix(mat.sub.t)
        gc()

        rpca.mat = rrpca.mod(mat.sub.t, tol = tolerance, trace = F)

        rm(mat.sub.t)
        gc()
        
        if (verbose){ 
            message("Starting clustering")
        }

        l.mat = rpca.mat$L
        
        rm(rpca.mat)
        gc()
        
        l.mat = t(l.mat)
        rownames(l.mat) = samp.names$nm
        clust.out = hclust(dist(l.mat))
        
        memb = cutree(clust.out, k = number.of.samples)
        memb = setDT(as.data.frame(data.matrix(memb)), keep.rownames = T)
        
        set.seed(12)
        samp.selected = memb[, .SD[sample(1:.N, 1)], by = V1]$rn
        samp.final = normal.table[samp.selected]
        setkey(samp.final, "sample")
    }

    if (use.all){
        if (verbose){ 
            message("Using all samples")
        }

        samp.final = normal.table
        setkey(samp.final, "sample")
    }

    
    message("Balancing pre-decomposition")

    if (is.human){
        if (build == "hg19"){
            par = 2700000
        } else if (build == "hg38"){
            par = 2782000
        } else {
            stop("provide either hg19 or hg38 build")
        }
    }


    samp.final[, file.available := file.exists(normal_cov)]

    message("Checking for existence of files")
    
    samp.final = samp.final[file.available == TRUE]

    mat.n = mclapply(samp.final[, sample], function(nm){
        this.cov = tryCatch(readRDS(samp.final[nm, normal_cov]), error = function(e) NULL)
        if (!is.null(this.cov)){
            all.chr = c(as.character(1:22), "X")
            this.cov = this.cov %Q% (seqnames %in% all.chr)
            this.cov = this.cov[, field] %>% gr2dt() %>% setnames(., field, "signal")
            setnames(this.cov, "signal", "signal.org")
            this.cov[, median.idx := .GRP, by = seqnames]
            if (is.human){
                this.cov[, median.idx := ifelse(seqnames == "X" & start < par, 24, median.idx)]
            }
            this.cov[, median.chr := median(signal.org, na.rm = T), by = median.idx]
            this.cov[, signal := signal.org/median.chr]
            message(nm)
            reads = this.cov[, .(seqnames, signal, median.chr)]
            ##reads[, median.chr := median(.SD$signal, na.rm = T), by = seqnames]
            reads[is.na(signal), signal := median.chr]
            min.cov = min(reads[signal > 0]$signal, na.rm = T)
            reads[is.infinite(signal), signal := min.cov]
            reads[signal == 0, signal := min.cov]
            reads[signal < 0, signal := min.cov]
            reads[, signal := log(signal)]
            reads = reads[, .(signal)]
            if (!any(is.infinite(reads$signal))){
                reads = transpose(reads)
                return(reads)
            } 
        } 
    }, mc.cores = num.cores)

    gc()

    mat.bind = rbindlist(mat.n, fill = T)
    mat.bind = na.omit(mat.bind)
    mat.bind.t = transpose(mat.bind)

    rm(mat.bind)
    gc()

    if (verbose){ 
        message("Starting decomposition")
    }

    mat.bind.t = as.matrix(mat.bind.t)
    gc()

    detergent = rrpca.mod(mat.bind.t, trace = F, tol = tolerance)

    rm(mat.bind.t)
    gc()
    
    rsvd.L.burnin = rsvd(detergent$L, k = detergent$k)
    detergent$U.hat = rsvd.L.burnin$u
    detergent$V.hat = t(rsvd.L.burnin$v)
    detergent$sigma.hat = rsvd.L.burnin$d

    if (verbose){ 
        message("Finished making the PON or detergent and saving it to the path provided")
    }

    if (save.pon){
        saveRDS(detergent, paste0(path.to.save, "/detergent.rds"))
    }
    
    return(detergent)
}


##############################
## identify_germline
##############################
#' @name identify_germline
#'
#' @title This is first phase in preparing the Panel of Normals (PON)
#' that will be used for online decomposition
#' 
#' @description This function takes in gRanges outputs from fragCounter and extracts GC corrected
#' read count data and carries rPCA decomposition on the matrix thus created
#' 
#' @export
#' @param normal.table.path character path to data.table containing two columns "sample", "normal_cov" and additional column "decomposed_cov" that contains drycleaned normal sample outputs. See manual for details.
#' 
#' @param signal.thresh numeric (default == 0.5). This is the threshold to be used to identify an amplification (markers with signal intensity > 0.5) or deletions (markers with signal intensity < -0.5) in log space from dryclean outputs.
#'
#' @param pct.thresh numeric (default == 0.98). Proportion of samples in which a given marker is free of germline event.
#'
#' @param save.grm boolean (default == FALSE). If the germline list needs to be saved.     
#' 
#' @param path.to.save charater (default == NA). Path to save the germline list created if save.grm == TRUE.
#' 
#' @param num.cores interger (default == 1). Number of cores to use for parallelization.
#'
#' @param verbose boolean (default == TRUE). Outputs progress.
#' 
#' @return A GRange object with a metadata field annotating germline markers.
#' 
#' @author Aditya Deshpande

identify_germline <- function(normal.table.path = NA, signal.thresh = 0.5, pct.thresh = 0.98, verbose = TRUE, save.grm = FALSE, path.to.save = NA, num.cores = 1){

    if (verbose){
        message("Starting the preparation of Panel of Normal samples a.k.a detergent")
    }

    if (is.na(normal.table.path)){
        stop("Need a table with paths to decomposed normal samples to identify germline events")
    }
    
    if (is.na(path.to.save) & save.grm == TRUE){
        stop("Need a path to save identified germline events")
    }

    normal.table = readRDS(normal.table.path)
    setkeyv(normal.table, "sample")
    
    mat.pon = mclapply(normal.table[, sample], function(nm){
        this.cov = tryCatch(readRDS(normal.table[nm, decomposed_cov]), error = function(e) NULL)
        if (!is.null(this.cov)){
            this.cov = gr2dt(this.cov)
            reads = this.cov[, .(foreground.log)] 
            reads = transpose(reads)
        } else {reads = data.table(NA)}
        return(reads)}, mc.cores = num.cores)
    
    gc()
    
    mat.bind = rbindlist(mat.pon, fill = T)
    ## mat.bind = na.omit(mat.bind)
    mat.bind.t = transpose(mat.bind)
    
    rm(mat.bind)
    gc()
    
    for(col in names(mat.bind.t)) set(mat.bind.t, i = which(abs(mat.bind.t[[col]]) > signal.thresh), j = col, value = NA) 
    for(col in names(mat.bind.t)) set(mat.bind.t, i = which(!is.na(mat.bind.t[[col]])), j = col, value = 1)
    for(col in names(mat.bind.t)) set(mat.bind.t, i = which(is.na(mat.bind.t[[col]])), j = col, value = 0)

    gc()
    
    mat.bind.t[, black_list_pct := rowSums(.SD)/dim(mat.bind.t)[2]]
    
    mat.bind.t[, germline.status := ifelse(black_list_pct > pct.thresh, FALSE, TRUE)]

    if (nrow(mat.bind.t[germline.status == FALSE]) < 0.5 * nrow(mat.bind.t)){
        warning("More than 50% markers classified as germline, consider adjusting thresholds.")
    }
    
    template = readRDS(normal.table[1, normal_cov])
    values(template) <- NULL
    template$germline.status <- mat.bind.t$germline.status

    
    rm(mat.bind.t)
    gc()

    if (verbose){ 
        message("Finished identifying germline markers based pn thresholds provided and saving it to the path provided")
    }

    if (save.grm){
        saveRDS(template, paste0(path.to.save, "/germline.markers.rds"))
    }

    return(template)

}


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
#' @export
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
    A = matrix(0, r, r)
    B = matrix(0, m, r)
    
    if (verbose == TRUE){
        message("calculating A and B")
    }

    for (i in 1:dim(V.hat)[2]){
        A = A + outer(V.hat[, i], V.hat[, i])
        B = B + outer((L.burnin[, i] - S.burnin[, i]), V.hat[, i])
    }

    if (verbose == TRUE){
        message("calculating v and s")
    }

    projection = apg_project(m.vec, U, lambda1, lambda2)
    vi = as.numeric(projection[[1]])
    si = as.numeric(projection[[2]])
    vi.subtract = V.hat[, 1]
    V.hat = cbind(V.hat, vi)
    A = A + outer(vi, vi) - outer(vi.subtract, vi.subtract)
    B = B + outer((as.numeric(m.vec) - si), vi) - outer((L.burnin[, 1] - S.burnin[, 1]), vi.subtract)

    if (verbose == TRUE){
        message("Updating subspace")
    }

    U = update_cols(U, A, B, lambda1)
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

prep_cov <- function(m.vec = m.vec, blacklist = FALSE, burnin.samples.path = NA){
    m.vec = gr2dt(m.vec)
    m.vec = m.vec[, .(seqnames, signal)]
    m.vec[, median.chr := median(.SD$signal, na.rm = T), by = seqnames]
    m.vec[is.na(signal), signal := median.chr]
    min.cov = min(m.vec[signal > 0]$signal, na.rm = T)
    m.vec[is.infinite(signal), signal := min.cov]
    m.vec[signal == 0, signal := min.cov]
    m.vec[signal < 0, signal := min.cov]

    if (blacklist){
        blacklist.pon =  readRDS(paste0(burnin.samples.path, "/blacklist.rds"))
        m.vec$blacklisted = blacklist.pon$blacklisted
        m.vec[blacklisted == TRUE, signal := NA]
        m.vec = na.omit(m.vec)
    }

    m.vec[, signal := log(signal)]

    return(m.vec)
}


##############################
## start_wash_cycle
##############################
#' @name start_wash_cycle
#'
#' @title function begins the online rPCA process. Use this function if you performed batch rPCA on samples as whole withiout dividing into chromosomes.For exomes and whole genomes where number of normal samples are small (<=100).  
#' @description function begins the online rPCA process.  
#' It is the wrapper that takes in GRanges and outputs GRanges with decomposition 
#' 
#' @param cov GRanges object containig the GC corrected cov data outputed from fragCounter. Needs metadata with header "reads.corrected".
#' 
#' @param mc.cores interger (default == 1). Number of cores to use for parallelization.
#' 
#' @param detergent.pon.path string. Path to pon/detergent genrated using normal samples.
#' 
#' @param whole_genome boolean (default = TRUE). For this function always set this parameter to TRUE.
#' 
#' @param use.blacklist boolean (default = FALSE). Whether to exclude off-target markers in case of Exomes or targeted sequqnecing. If set to TRUE, needs a GRange marking if each marker is set to be excluded or not.
#'
#' @param germline.filter boolean (default == TRUE). If germline markers need to be removed from decomposition.
#'
#' @param germline.file character (default == NA). Path to file with germline markers.
#' 
#' @param verbose boolean (default == TRUE). Outputs progress.
#' 
#' @param is.human boolean (default == TRUE). Organism type.
#' 
#' @param field character (default == "reads.corrected"). Field to use for processing.
#' 
#' @param chr integer (default == NA). Depricated. Can be used to decompose a single chromosome.
#' 
#' @export
#' 
#' @return GRange object with decomposition
#' @author Aditya Deshpande



start_wash_cycle <- function(cov, mc.cores = 1, detergent.pon.path = NA, verbose = TRUE, whole_genome = TRUE, use.blacklist = FALSE, chr = NA, germline.filter = FALSE, germline.file = NA, field = "reads.corrected", is.human = TRUE){

    if(verbose == TRUE){
        message("Loading PON a.k.a detergent from path provided")
    }

    if (is.na(detergent.pon.path)){
        stop('Need pon/detergent file to procced. Use prepare_detergent() command')
    }

    rpca.1 = readRDS(detergent.pon.path)
    
    if(verbose == TRUE){
        message(paste0("Let's begin, this is whole exome/genome"))
    }

    if (germline.filter & is.na(germline.file)){
        stop("If germiline.filter is set to TRUE, provide path to germline marker file")
    }

    all.chr = c(as.character(1:22), "X")
    
    if (is.human){
        cov = cov %Q% (seqnames %in% all.chr)
    }
    
    cov = cov[, field] %>% gr2dt() %>% setnames(., field, "signal") %>% dt2gr()
    m.vec = prep_cov(cov, blacklist = use.blacklist,
                     burnin.samples.path = detergent.pon.path)
    
    m.vec = as.matrix(m.vec$signal)
    L.burnin = rpca.1$L
    S.burnin = rpca.1$S
    r = rpca.1$k
    U.hat = rpca.1$U.hat
    V.hat = rpca.1$V.hat
    sigma.hat = rpca.1$sigma.hat

    if(verbose == TRUE){
        message("Initializing wash cycle")
    }

    decomposed = wash_cycle(m.vec = m.vec, L.burnin = L.burnin,
                            S.burnin = S.burnin, r = r, U.hat = U.hat,
                            V.hat = V.hat, sigma.hat = sigma.hat)
    
    if (verbose == TRUE){
        message("Combining matrices with gRanges")
    }

    cov = gr2dt(cov)
    cov[, median.chr := median(.SD$signal, na.rm = T), by = seqnames]

    if (use.blacklist){
        cov[is.na(signal), signal := median.chr]
        cov[is.infinite(signal), signal := median.chr]
        blacklist.pon =  readRDS(paste0(detergent.pon.path, "/blacklist.rds"))
        cov$blacklisted = blacklist.pon$blacklisted
        cov[blacklisted == TRUE, signal := NA]
        cov = na.omit(cov)
    }
    
    setnames(cov, "signal", "input.read.counts")
    cov = cbind(decomposed[[2]], cov)
    colnames(cov)[1] = 'foreground.log'
    cov[is.na(input.read.counts), foreground.log := NA]
    cov[, foreground := exp(foreground.log)]
    cov[input.read.counts == 0, foreground := 0]
    cov[is.na(input.read.counts), foreground := NA]
    cov = cbind(decomposed[[1]], cov)
    colnames(cov)[1] = 'background.log'
    cov[is.na(input.read.counts), background.log := NA]
    cov[, background := exp(background.log) ]
    cov[input.read.counts == 0, background := 0]
    cov[is.na(input.read.counts), background := NA]
    cov[, log.reads := log(input.read.counts)]
    cov[is.infinite(log.reads), log.reads := NA]

    if (germline.filter){
        germ.file = readRDS(germline.file)
        cov$germline.status = germ.file$germline.status
        cov[germline.status == TRUE, foreground := NA]
        cov[germline.status == TRUE, foreground.log := NA]
        cov = na.omit(cov)
    }

    cov = dt2gr(cov)
    return(cov)
}

message("Giddy up!")

    
