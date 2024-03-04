#' @importFrom data.table setkeyv
#' @import GenomicRanges
#' @import rsvd
#' @importFrom MASS ginv
#' @import gUtils
#' @import DNAcopy
#' @importFrom dplyr mutate filter


##############################
## dryclean R6 object
##############################
#' @name dryclean
#' @title dryclean
#' @description dryclean R6 class storing all methods and values necessary for "drycleaning"
#' @details Add more details 
#' 
#' @export
#' 
#' @author Aditya Deshpande <asd3002@med.cornell.edu>, Sebastian Brylka <sebastian.brylka@nyulangone.org> 


dryclean <- R6::R6Class("dryclean",
    private = list(
    history = NULL,
    pon = NULL,
    dt_mismatch = NULL
    ),
  
  public = list(
    
    #' @method initialize() initialize()
    #' @description Initialize dryclean object. Authors: Aditya Deshpande, Sebastian Brylka
    #' 
    #' @param pon PON object
    #' 

    
    initialize = function(pon) {
      
      private$history <- data.table(action = character(), date = character())
      private$history <- rbindlist(list(private$history, data.table(action = "Created dryclean object", date = as.character(Sys.time()))))
      
      private$pon = pon

    },
    
    #' @method clean() clean()
    #' @description Function begins the online rPCA process. Use this function if you performed batch rPCA on samples as whole without dividing into chromosomes. For exomes and whole genomes where number of normal samples are small (<=100). Author: Aditya Deshpande 
    #' @details The function begins the online rPCA process. It is the wrapper that takes in GRanges and outputs GRanges with decomposition 
    #' 
    #' @param cov character path to the granges coverage file to be drycleaned
    #'
    #' @param center boolean (default == TRUE) whether to center the coverage before drycleaning
    #' 
    #' @param cbs boolean (default == FALSE) whether to perform cbs on the drycleaned coverage
    #' 
    #' @param cnsignif integer (default = 1e-5) the significance levels for the tests in cbs to accept change-points
    #'
    #' @param mc.cores interger (default == 1) number of cores to use for parallelization
    #' 
    #' @param use.blacklist boolean (default = FALSE) whether to exclude off-target markers in case of Exomes or targeted sequencing. If set to TRUE, needs a GRange marking if each marker is set to be excluded or not or will use a default mask
    #'
    #' @param blacklist_path character (default = NA) if use.blacklist == TRUE, path a GRanges object marking if each marker is set to be excluded or not
    #'
    #' @param germline.filter boolean (default == FALSE) if germline markers need to be removed from decomposition
    #'
    #' @param verbose boolean (default == TRUE) outputs progress
    #' 
    #' @param field character (default == "reads.corrected") field to use for processing
    #' 
    #' @param testing boolean (default = FALSE) DO NOT CHANGE
    #' 
    #' @return Drycleaned coverage in GRanges format
    
    clean = function(cov, center = TRUE, centering = "mean", cbs = FALSE, cnsignif = 1e-5, mc.cores = 1, verbose = TRUE, use.blacklist = FALSE, blacklist_path = NA, germline.filter = FALSE, field = "reads.corrected", testing = FALSE){
      
      message("Loading coverage")
      private$history <- rbindlist(list(private$history, data.table(action = paste("Loaded coverage from", cov), date = as.character(Sys.time()))))
      cov = readRDS(cov)
      
      is.chr = FALSE
      if(any(grepl("chr", as.character(seqnames(cov))))){
        cov = gr.sub(cov)
        is.chr = TRUE
      }
      
      cov <- cov %>% gr2dt() %>% filter(seqnames != 'Y') %>% dt2gr()
      #cov <- cov %>% gr2dt() %>% filter(seqnames != 'X') %>% dt2gr()
      
      if(verbose == TRUE){
        message("Loading PON a.k.a detergent")
      }
      
      
      private$history <- rbindlist(list(private$history, data.table(action = "Loaded PON", date = as.character(Sys.time()))))
        
      tumor.binsize = median(gr2dt(cov)$width)
      pon.binsize = median(gr2dt(private$pon$get_template())$width)
      
      if (tumor.binsize != pon.binsize & testing == FALSE){
        message(paste0("WARNING: Input tumor bin size = ", tumor.binsize,"bp. PON bin size = ", pon.binsize,"bp. Rebinning tumor to bin size of PON..."))
        private$history <- rbindlist(list(private$history, data.table(action = paste("Rebinning tumor to", pon.binsize, "bp bin size"), date = as.character(Sys.time()))))
        suppressWarnings({ 
          cov = gr.val(query = private$pon$get_template(), cov, val = field)
          })
      }
      
      pon.length <- private$pon$get_template() %>%
        gr2dt() %>%
        dplyr::filter(seqnames != "Y") %>%
        dt2gr() %>%
        length()
      
      if (length(cov) != pon.length & testing == FALSE) {
        dt_mismatch = data.table(chr = c(), coverage = c(), pon = c())                                                                            
        for(chr in c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "\
X")){                                                                                                                                     
          dt_mismatch <- rbind(dt_mismatch,                                                                                                       
                               data.table(                                                                                                        
                                 chr = chr,                                                                                                       
                                 coverage = cov %>% gr2dt() %>% filter(seqnames == chr) %>% nrow(),                                               
                                 pon = private$pon$get_template() %>% gr2dt() %>% filter(seqnames == chr) %>% nrow()                               
                               ))                                                                                                                 
        }                                                                                                                                         
        dt_mismatch <- dt_mismatch %>% filter(coverage != pon)  
        private$dt_mismatch = dt_mismatch
        if(testing == FALSE){
          message("WARNING: Number of bins of coverage and PON does not match. Use get_mismatch() function to see mismatched chromosomes\nAligning coverage to the PON")
          suppressWarnings({ 
            cov = gr.val(query = private$pon$get_template(), cov, val = field)
          })
        }
      }
      
      if(verbose == TRUE){
        message(paste0("Let's begin, this is whole exome/genome"))
      }
      
      private$history <- rbindlist(list(private$history, data.table(action = paste("Started drycleaning the coverage file"), date = as.character(Sys.time()))))
      
      if (germline.filter & is.null(private$pon$get_inf_germ())){
        stop("If germiline.filter is set to TRUE, pon must have a inf_germ element, see prepare_detergent for details")
      }
      
      all.chr = c(as.character(1:22), "X")
      
      local.all.chr = all.chr
      cov = cov %Q% (seqnames %in% local.all.chr)
      cov = cov[, field] %>% gr2dt() %>% setnames(., field, "signal")
      cov = cov %>% dt2gr()
      
      if(center == TRUE){
        
        
        message("Median-centering the sample")
        private$history <- rbindlist(list(private$history, data.table(action = paste("Median-normalization of coverage"), date = as.character(Sys.time()))))
        mcols(cov)[which(is.na(mcols(cov)[, "signal"])), "signal"] = 0
        mcols(cov)[which(is.infinite(mcols(cov)[, "signal"])), "signal"] = NA
        if(centering == "mean"){
          values(cov)[, "signal"] = values(cov)[, "signal"] / mean(values(cov)[, "signal"], na.rm = TRUE)
        }
        if(centering == "median"){
          values(cov)[, "signal"] = values(cov)[, "signal"] / median(values(cov)[, "signal"], na.rm = TRUE)
        }
      }      
      
      
      cov = sortSeqlevels(cov)
      cov = sort(cov)
      
      if(use.blacklist == TRUE ){
        if((is.na(blacklist_path) | blacklist_path == "NA")){
          blacklist_path = system.file("extdata", "blacklist_A.rds", package = 'dryclean')
          message(paste0("Applying the default mask to the coverage"))
          private$history <- rbindlist(list(private$history, data.table(action = paste("Applying the defualt mask to coverage"), date = as.character(Sys.time()))))
        }else{
          blacklist_path = blacklist_path
          message(paste0("Applying the provided mask to the coverage"))
          private$history <- rbindlist(list(private$history, data.table(action = paste("Applying the provided mask to coverage"), date = as.character(Sys.time()))))
        }
        
        suppressWarnings({ 
          blacklist_pon <- gUtils::gr.val(query = private$pon$get_template(),
                                          target = readRDS(blacklist_path),                                                                                               
                                          val = "blacklisted",                                                                                                       
                                          na.rm = TRUE)$blacklisted 
          blacklist_pon <- ifelse(blacklist_pon == 1, TRUE, FALSE)
          
          
          blacklist_cov <- gUtils::gr.val(query = cov,
                                          target = readRDS(blacklist_path),                                                                                               
                                          val = "blacklisted",                                                                                                       
                                          na.rm = TRUE)$blacklisted                                                                                                  
          blacklist_cov <- ifelse(blacklist_cov == 1, TRUE, FALSE)
        })
         
      }
      
      
      m.vec = prep_cov(cov, use.blacklist = use.blacklist, blacklist = blacklist_cov)
      
      
      m.vec = as.matrix(m.vec$signal)
      L.burnin = private$pon$get_L()
      S.burnin = private$pon$get_S()
      r = private$pon$get_k()
      U.hat = private$pon$get_U_hat()
      V.hat = private$pon$get_V_hat()
      sigma.hat = private$pon$get_sigma_hat()
      
      
      if(verbose == TRUE){
        message("Initializing wash cycle")
      }
      
      if (use.blacklist){
        blacklisted <- !blacklist_pon
        L.burnin = L.burnin[blacklisted,]
        S.burnin = S.burnin[blacklisted,]
        U.hat = U.hat[blacklisted, ]
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
        cov_template <- cov
        cov[is.na(signal), signal := median.chr]
        cov[is.infinite(signal), signal := median.chr]
        cov = cov[!blacklist_cov,]
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
        germ.file = private$pon$get_inf_germ()
        cov$germline.status = germ.file$germline.status
        cov[germline.status == TRUE, foreground := NA]
        cov[germline.status == TRUE, foreground.log := NA]
        cov = na.omit(cov)
      }
      
      cov = dt2gr(cov)

      if(use.blacklist){
        cov_template = dt2gr(cov_template)
        suppressWarnings({ 
          cov <- gUtils::gr.val(query = cov_template,
                                target = cov,                                                                                               
                                val =c("background.log",
                                        "foreground.log", 
                                        "input.read.counts", 
                                        "median.chr",
                                        "foreground",
                                        "background",
                                        "log.reads"),                                                                                                       
                                na.rm = TRUE)
          })

      }
      
      if (is.chr){
        cov = gr.chr(cov)
      }
      
      private$history <- rbindlist(list(private$history, data.table(action = paste("Finished drycleaning the coverage file"), date = as.character(Sys.time()))))
    
      if(cbs == TRUE){
        
        message("Starting CBS on the drycleaned sample")
        
        tcov = cov
        
        tcov$ratio = values(tcov)[, "foreground"]
        ss.n = NULL
        new.sl = seqlengths(tcov)
        
        if (any(is.na(new.sl)))
        {
          tmp.sl = data.table(sn = as.character(seqnames(tcov)),
                              end = end(tcov))[ , max(end, na.rm = T), by = sn][ ,  structure(V1, names = sn)]
          new.sl[is.na(new.sl)] = tmp.sl[is.na(new.sl)]
          new.sl = new.sl[!is.na(new.sl)]
        }
        
        
        ix = which(!is.na(tcov$ratio))
        cat('sending ', length(ix), ' segments\n')
        suppressWarnings({
          cna = DNAcopy::CNA(log(tcov$ratio[ix]), as.character(seqnames(tcov))[ix], start(tcov)[ix], data.type = 'logratio')
        })
        gc()
        cat('finished making cna\n')
        seg = DNAcopy::segment(DNAcopy::smooth.CNA(cna), alpha = cnsignif, verbose = T) ## 1e-5!!! TODO URGENT
        cat('finished segmenting\n')
        utils::capture.output({seg_dt = print(seg); setDT(seg_dt)}, type = "output", file = "/dev/null") #### KH
        out = gUtils::seg2gr(seg_dt[!(is.na(seg.mean) | is.na(loc.start) | is.na(loc.end))], new.sl) ## remove seqlengths that have not been segmented #### KH
        out = gUtils::gr.fix(out, new.sl, drop = T)
        cat(length(out), ' segments produced\n')
        names(out) = NULL
        
        gc()
        
        private$history <- rbindlist(list(private$history, data.table(action = paste("Applied CBS correction to the drycleaned coverage file"), date = as.character(Sys.time()))))
        
        #return(out)
        saveRDS(out, "cbs_output.rds")
        
        private$history <- rbindlist(list(private$history, data.table(action = paste("Saved CBS output in current directory as cbs_output.rds"), date = as.character(Sys.time()))))
        
      }
      
      return(cov)
    },

     
    #' @method get_history() get_history()
    #' @description Function returns the history of the dryclean object
    #'
    #' @return Prints the history of the dryclean object as data table 
    get_history = function(){
      for (i in 1:nrow(private$history)){
        cat(paste0(private$history$date[i], "\t", private$history$action[i], "\n"))
      }
    },
    
    #' @method get_mismatch() get_mismatch()
    #' @description Function returns the data table with mismatching in seqlengths between pon and coverage
    #'
    #' @return Returns the data table with mismatching chromosomes of pon and coverage
    get_mismatch = function(){
      return(private$dt_mismatch)
    }

  )
)

##############################
## pon R6 object
##############################
#' @name pon
#' @title pon
#' @description pon R6 class storing PON (panel of normals) necessary for "drycleaning"
#' @details Add more details 
#' 
#' @export
#' 
#' @author Aditya Deshpande <asd3002@med.cornell.edu>, Sebastian Brylka <sebastian.brylka@nyulangone.org> 

pon <- R6::R6Class("pon",
                   private = list(
                     L = NULL,
                     S = NULL,
                     err = NULL,
                     k = NULL,
                     U.hat = NULL,
                     V.hat = NULL,
                     sigma.hat = NULL,
                     inf.germ = NULL,
                     template = NULL,
                     
                     prepare_pon = function(save_pon, use.all, choose.randomly, choose.by.clustering, number.of.samples, verbose, num.cores, tolerance, is.human, build, field, PAR.file, balance, infer.germline, signal.thresh, pct.thresh, wgs, target_resolution, nochr, all.chr){
                       
                       normal.table = private$normal.table

                       if(save_pon == TRUE){
                         message(paste0("WARNING: New PON will be generated and saved at ",private$pon.path))
                         Sys.sleep(3)
                         message("\nGiving you some time to think...\n")
                         Sys.sleep(5)
                         
                         path.to.save = private$pon.path
                       }
                       
                       private$history <- rbindlist(list(private$history, data.table(action = paste("Started PON preparation"), date = as.character(Sys.time()))))
                       if (verbose){
                         message("Starting the preparation of Panel of Normal samples a.k.a detergent")
                       }
                       
                       #normal.table = readRDS(normal.table.path)
                       setkeyv(normal.table, "sample")
                       
                       num.samp = nrow(normal.table)
                       
                       if (verbose){
                         message(paste0(num.samp, " samples available"))
                       }
                       
                       if (use.all & choose.randomly | use.all & choose.by.clustering | choose.randomly & choose.by.clustering | use.all & choose.randomly & choose.by.clustering){
                         stop("only one of use.all, choose.randomly, choose.by.clustering can be set to TRUE. Rectify and restart")
                       }
                       
                       if (nochr) {
                         template = generate_template(cov = gUtils::gr.nochr(readRDS(normal.table[1]$normal_cov)), wgs = wgs, target_resolution = target_resolution, this.field = field, nochr = nochr, all.chr = all.chr)
                       } else {
                         template = generate_template(cov = readRDS(normal.table[1]$normal_cov), wgs = wgs, target_resolution = target_resolution, this.field = field, nochr = nochr, all.chr = all.chr)
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
                             if (nochr) {
                               this.cov = gUtils::gr.nochr(this.cov)
                             }
                             ## this.cov = standardize_coverage(gr.nochr(this.cov), template = template, wgs = wgs, target_resolution = target_resolution, this.field = field)
                             this.cov = this.cov[, field] %>% gr2dt() %>% setnames(., field, "signal")
                             ## reads = this.cov[seqnames == "22", .(seqnames, signal)]
                             reads = this.cov[seqnames == seqnames[1], .(seqnames, signal)]
                             reads[, median.chr := median(.SD$signal, na.rm = T), by = seqnames]
                             reads[is.na(signal), signal := median.chr]
                             min.cov = min(reads[signal > 0]$signal, na.rm = T)
                             reads[signal == 0, signal := min.cov]
                             reads[signal < 0, signal := min.cov]
                             reads = log(reads[, .(signal)])
                             reads = transpose(reads)
                             reads = cbind(reads, nm)
                           } ## else {reads = data.table(NA)}
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
                       
                       
                       if(!is.null(PAR.file)){
                         if(PAR.file == "NA"){PAR.file = NA}
                       }
                       if(is.null(PAR.file)){PAR.file = NA}
                       
                       
                       if (is.na(PAR.file)){
                         if (build == "hg38") {
                           message("PAR file not provided, using hg38 default. If this is not the correct build, please provide a GRanges object delineating for corresponding build")
                           par.path = system.file("extdata", "PAR_hg38.rds", package = 'dryclean')
                         } else {
                           message("PAR file not provided, using hg19 default. If this is not the correct build, please provide a GRanges object delineating for corresponding build")
                           par.path = system.file("extdata", "PAR_hg19.rds", package = 'dryclean')
                         }
                         par.gr = readRDS(par.path)
                       } else {par.gr = readRDS(PAR.file)}
                       
                       message("PAR read")
                       
                       #samp.final[, file.available := file.exists(normal_cov)]
                       samp.final <- samp.final %>%
                         mutate(file.available = file.exists(normal_cov))
                       
                       
                       message("Checking for existence of files")
                       
                       samp.final = samp.final[file.available == TRUE]
                       
                       message(paste0(nrow(samp.final), " files present"))
                       
                       
                       mat.n = pbmcapply::pbmclapply(samp.final[, sample], function(nm, all.chr){
                         this.cov = tryCatch(readRDS(samp.final[sample == nm, normal_cov]), error = function(e) NULL)
                         if (!is.null(this.cov)){
                           ## this.cov = standardize_coverage(cov = gr.nochr(this.cov), template = template, wgs = wgs, target_resolution = target_resolution, this.field = field)
                           if (nochr) {
                             this.cov = gr.nochr(this.cov)
                           }
                           this.cov = this.cov %Q% (seqnames %in% all.chr)
                           this.cov = sortSeqlevels(this.cov)
                           this.cov = sort(this.cov)
                           this.cov = this.cov[, field] %>% gr2dt() %>% setnames(., field, "signal.org")
                           if (balance){
                             this.cov[, median.idx := .GRP, by = seqnames]
                             this.cov$mt = suppressWarnings(gr.match(dt2gr(this.cov), par.gr))
                             this.cov[, median.idx := ifelse(is.na(mt), median.idx, mt+24)]
                             this.cov[, median.chr := median(signal.org, na.rm = T), by = median.idx]
                             this.cov[, signal := ifelse(seqnames != "X", signal.org,
                                                         ifelse(median.chr == 0, 1, signal.org/median.chr))]
                           } else {
                             this.cov[, signal := signal.org]
                             this.cov[, median.chr := median(signal.org, na.rm = T)]
                           }
                           reads = this.cov[, .(seqnames, signal, median.chr)]
                           reads[is.na(signal), signal := median.chr]
                           min.cov = min(reads[signal > 0]$signal, na.rm = T)
                           reads[is.infinite(signal), signal := min.cov]
                           reads[signal == 0, signal := min.cov]
                           reads[signal < 0, signal := min.cov]
                           reads[, signal := log(signal)]
                           reads = reads[, .(signal)]
                           if (!any(is.infinite(reads$signal))){
                             #reads = transpose(reads)
                             return(reads$signal)
                           }
                         }
                       }, all.chr, mc.cores = num.cores, ignore.interactive = T)

                       if (verbose){
                         message("Starting decomposition")
                       }
                       

                       mat.bind.t = matrix(unlist(mat.n), ncol = length(mat.n))
                       #print(nrow(mat.bind.t))
                       #print(ncol(mat.bind.t))
                       rm(mat.n)
                       gc()

                       detergent = rrpca.mod(mat.bind.t, trace = F, tol = tolerance)
                       
                       rm(mat.bind.t)
                       gc()
                       
                       rsvd.L.burnin = rsvd(detergent$L, k = detergent$k)
                       detergent$U.hat = rsvd.L.burnin$u
                       detergent$V.hat = t(rsvd.L.burnin$v)
                       detergent$sigma.hat = rsvd.L.burnin$d
                       ##browser()
                       ##this.template = readRDS(samp.final[1]$normal_cov)
                       ##this.template = sortSeqlevels(this.template)
                       ##this.template = sort(this.template)
                       ##
                       detergent$template = template
                       
                       if (infer.germline){
                         this.s = as.data.table(detergent$S)
                         gc()
                         for(col in names(this.s)) set(this.s, i = which(abs(this.s[[col]]) > signal.thresh), j = col, value = NA)
                         for(col in names(this.s)) set(this.s, i = which(!is.na(this.s[[col]])), j = col, value = 1)
                         for(col in names(this.s)) set(this.s, i = which(is.na(this.s[[col]])), j = col, value = 0)
                         gc()
                         this.germ = this.s[, .(black_list_pct = rowSums(.SD)/dim(this.s)[2])]
                         this.germ[, germline.status := ifelse(black_list_pct > pct.thresh, FALSE, TRUE)]
                         this.germ = dt2gr(cbind(gr2dt(detergent$template), this.germ))
                         detergent$inf_germ = this.germ
                       }
                       
                       private$L = detergent$L
                       private$S = detergent$S
                       private$err = detergent$err
                       private$k = detergent$k
                       private$U.hat = detergent$U.hat
                       private$V.hat = detergent$V.hat
                       private$inf.germ = detergent$inf_germ
                       private$sigma.hat = detergent$sigma.hat
                       
                       private$template = detergent$template
                       private$history <- rbindlist(list(private$history, data.table(action = paste("Created new PON from the normal samples"), date = as.character(Sys.time()))))
                      
                      if (verbose){
                         message("Finished making the PON")
                       }
                       
                       if(save_pon == TRUE){
                         saveRDS(detergent, paste0(path.to.save))
                         private$history <- rbindlist(list(private$history, data.table(action = paste("Saved new PON at", private$pon.path), date = as.character(Sys.time()))))
                         if (verbose){
                           message("Finished saving the PON to the provided path")
                         }
                         
                        }
                       
                       rm(detergent)
                       
                     },
                     
                     history = NULL,
                     normal.table = NULL,
                     pon.path = NULL
                  
                   ),
                   
                   public = list(
                     
                     
                     #' @method initialize() initialize()
                     #' @description Initialize PON object. Authors: Aditya Deshpande, Sebastian Brylka
                     #' 
                     #' @param normal_vector character (default == c()) vector of paths to normal samples
                     #' 
                     #' @param pon_path character (default == NULL) path to PON/detergent
                     #' 
                     #' @param create_new_pon boolean (default == FALSE) whether to create a new PON from normal samples
                     #' 
                     #' @param save_pon boolean (default == FALSE) if create_new_pon == TRUE, whether to save pon to path given by pon_path
                     #' 
                     #' @param field character (default == "reads.corrected") field to use for processing
                     #'
                     #' @param use.all boolean (default == TRUE) if all normal samples are to be used for creating PON
                     #' 
                     #' @param choose.randomly boolean (default == FALSE) if a random subset of normal samples are to be used for creating PON
                     #' 
                     #' @param choose.by.clustering boolean (default == FALSE) clusters normal samples based on the genomic background and takes a random sample from within the clusters
                     #' 
                     #' @param number.of.samples interger (default == 50) if choose.by.clustering == TRUE, this is the number of clusters at which to cut tree
                     #'
                     #' @param tolerance numeric (default == 0.0001) tolerance for error for batch rPCA. We suggest keeping this value
                     #' 
                     #' @param num.cores integer (default == 1) number of cores to use for parallelization
                     #'
                     #' @param verbose boolean (default == TRUE) outputs progress
                     #'
                     #' @param is.human boolean (default == TRUE) organism type
                     #' 
                     #' @param build character (default == "hg19") genome build to define PAR region in chromosome X
                     #'
                     #' @param PAR.file character (default == NULL) this is a GRanges with the boundaries of PAR region in X chr
                     #'
                     #' @param balance boolean (default == TRUE) experimental variable to take into consideration 1 copy of X chr in male sample
                     #'
                     #' @param infer.germline boolean (default = FALSE) if use the L matrix to infer germline events
                     #'
                     #' @param signal.thresh numeric (default == 0.5) this is the threshold to be used to identify an amplification (markers with signal intensity > 0.5) or deletions (markers with signal intensity < -0.5) in log space from dryclean outputs
                     #'
                     #' @param pct.thresh numeric (default == 0.98) proportion of samples in which a given marker is free of germline event
                     #'
                     #' @param wgs boolean (default == TRUE) whether whole genome is being used
                     #'
                     #' @param target_resolution numeric (default == 1e3) resolution at which to conduct analyses
                     #'
                     #' @param nochr logical (default = TRUE) remove chr prefix
                     #'
                     #' @param all.chr list (default = c(as.character(1:22), "X")) list of chromosomes
                     
                     initialize = function(normal_vector = c(), pon_path = NULL, create_new_pon = FALSE, save_pon = FALSE, field = "reads.corrected", use.all = TRUE, choose.randomly = FALSE, choose.by.clustering = FALSE, number.of.samples = 50, verbose = TRUE, num.cores = 1, tolerance = 0.0001, is.human = TRUE, build = "hg19", PAR.file = NULL, balance = TRUE, infer.germline = FALSE, signal.thresh = 0.3, pct.thresh = 0.80, wgs = TRUE, target_resolution = 1000, nochr = TRUE, all.chr = c(as.character(1:22), "X")) {
                       
                       message("Loading PON...")
                       
                       private$history <- data.table(action = character(), date = character())
                       private$history <- rbindlist(list(private$history, data.table(action = "Created pon object", date = as.character(Sys.time()))))
                       
                       if(is.null(pon_path) & create_new_pon == FALSE){
                         stop("ERROR: Provide pon_path or set create_new_pon = TRUE")
                       }
                       
                       if(!is.null(pon_path) & create_new_pon == FALSE){
                         if(!file.exists(pon_path)){
                           stop("ERROR: PON data file not found or is invalid.")
                         }
                         pon = readRDS(pon_path)
                         private$L = pon$L
                         private$S = pon$S
                         private$err = pon$err
                         private$k = pon$k
                         private$U.hat = pon$U.hat
                         private$V.hat = pon$V.hat
                         private$sigma.hat = pon$sigma.hat
                         private$template = pon$template
                         rm(pon)
                       }
                       
                       if(create_new_pon == TRUE){
                         
                         if(length(normal_vector) == 0){
                           stop("ERROR: Invalid input vector of paths to normal samples")
                         }
                         
                         if (!file.exists(normal_vector[1])){
                           stop("ERROR: Invalid input vector of paths to normal samples")
                         }
                         
                         if(save_pon == TRUE){
                           if(is.null(pon_path)){
                             stop("ERROR: If you want to save the pon, provide pon_path")
                           }
                           private$pon.path = pon_path
                         }
                         
                         private$normal.table <- data.table(normal_cov = normal_vector) %>% mutate(sample = paste0("sample_",row_number())) %>% setcolorder(c("sample","normal_cov"))
                         private$history <- rbindlist(list(private$history, data.table(action = paste("Loaded normal vector"), date = as.character(Sys.time()))))
                         
                         private$prepare_pon(save_pon = save_pon, use.all = use.all, choose.randomly = choose.randomly, choose.by.clustering = choose.by.clustering, number.of.samples = number.of.samples, verbose = verbose, num.cores = num.cores, tolerance = tolerance, is.human = is.human, build = build, field = field, PAR.file = PAR.file, balance = balance, infer.germline = infer.germline, signal.thresh = signal.thresh, pct.thresh = pct.thresh, wgs = wgs, target_resolution = target_resolution, nochr = nochr, all.chr = all.chr)
                         
                       }
                         
                       message("PON loaded")
                       
                     },
                     
                     #' @method get_history() get_history()
                     #' @description Function returns the history of the pon object
                     #'
                     #' @return Prints the history of the pon object as data table 
                     get_history = function(){
                       for (i in 1:nrow(private$history)){
                         cat(paste0(private$history$date[i], "\t", private$history$action[i], "\n"))
                       }
                     },
                     
                     #' @method get_template() get_template()
                     #' @description Function returns the template of the pon object
                     #'
                     #' @return Template of the pon object 
                     get_template = function(){
                       return(private$template)
                     },
                     
                     #' @method get_L() get_L()
                     #' @description Function returns the L matrix of the pon object
                     #'
                     #' @return L matrix of the pon object 
                     get_L = function(){
                       return(private$L)
                     },
                     
                     #' @method get_S() get_S()
                     #' @description Function returns the S matrix of the pon object
                     #'
                     #' @return S matrix of the pon object 
                     get_S = function(){
                       return(private$S)
                     },
                     
                     #' @method get_err() get_err()
                     #' @description Function returns the err matrix of the pon object
                     #'
                     #' @return err matrix of the pon object 
                     get_err = function(){
                       return(private$err)
                     },
                     
                     #' @method get_k() get_k()
                     #' @description Function returns the k matrix of the pon object
                     #'
                     #' @return k matrix of the pon object 
                     get_k = function(){
                       return(private$k)
                     },
                     
                     #' @method get_U_hat() get_U_hat()
                     #' @description Function returns the U_hat matrix of the pon object
                     #'
                     #' @return U_hat matrix of the pon object 
                     get_U_hat = function(){
                       return(private$U.hat)
                     },
                     
                     #' @method get_V_hat() get_V_hat()
                     #' @description Function returns the V_hat matrix of the pon object
                     #'
                     #' @return V_hat matrix of the pon object 
                     get_V_hat = function(){
                       return(private$V.hat)
                     },
                     
                     #' @method get_sigma_hat() get_sigma_hat()
                     #' @description Function returns the sigma_hat matrix of the pon object
                     #'
                     #' @return sigma_hat matrix of the pon object 
                     get_sigma_hat = function(){
                       return(private$sigma.hat)
                     },
                     
                     #' @method get_inf_germ() get_inf_germ()
                     #' @description Function returns the inf_germ matrix of the pon object
                     #'
                     #' @return inf_germ matrix of the pon object 
                     get_inf_germ = function(){
                       return(private$inf.germ)
                     },
                     
                     #' @method get_seqlengths() get_seqlengths()
                     #' @description Function returns the seqlengths of the template of pon object
                     #'
                     #' @return seqlengths of template of the pon object 
                     get_seqlengths = function(){
                       return(seqlengths(private$template))
                     }
                   )
)
