% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dryclean.R
\name{pon}
\alias{pon}
\title{pon}
\description{
pon R6 class storing PON (panel of normals) necessary for "drycleaning"
}
\details{
Add more details
}
\author{
Aditya Deshpande <asd3002@med.cornell.edu>, Sebastian Brylka <sebastian.brylka@nyulangone.org>
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-pon-new}{\code{pon$new()}}
\item \href{#method-pon-get_history}{\code{pon$get_history()}}
\item \href{#method-pon-get_template}{\code{pon$get_template()}}
\item \href{#method-pon-get_L}{\code{pon$get_L()}}
\item \href{#method-pon-get_S}{\code{pon$get_S()}}
\item \href{#method-pon-get_err}{\code{pon$get_err()}}
\item \href{#method-pon-get_k}{\code{pon$get_k()}}
\item \href{#method-pon-get_U_hat}{\code{pon$get_U_hat()}}
\item \href{#method-pon-get_V_hat}{\code{pon$get_V_hat()}}
\item \href{#method-pon-get_sigma_hat}{\code{pon$get_sigma_hat()}}
\item \href{#method-pon-get_inf_germ}{\code{pon$get_inf_germ()}}
\item \href{#method-pon-get_seqlengths}{\code{pon$get_seqlengths()}}
\item \href{#method-pon-save_pon}{\code{pon$save_pon()}}
\item \href{#method-pon-clone}{\code{pon$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-pon-new"></a>}}
\if{latex}{\out{\hypertarget{method-pon-new}{}}}
\subsection{Method \code{new()}}{
Initialize PON object. Authors: Aditya Deshpande, Sebastian Brylka
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{pon$new(
  normal_vector = c(),
  pon_path = NULL,
  create_new_pon = FALSE,
  save_pon = FALSE,
  field = "reads.corrected",
  use.all = TRUE,
  choose.randomly = FALSE,
  choose.by.clustering = FALSE,
  number.of.samples = 50,
  verbose = TRUE,
  num.cores = 1,
  tolerance = 1e-04,
  is.human = TRUE,
  build = "hg19",
  PAR.file = NULL,
  balance = TRUE,
  infer.germline = FALSE,
  signal.thresh = 0.3,
  pct.thresh = 0.8,
  wgs = TRUE,
  target_resolution = 1000,
  nochr = TRUE,
  all.chr = c(as.character(1:22), "X", "Y")
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{normal_vector}}{character (default == c()) vector of paths to normal samples}

\item{\code{pon_path}}{character (default == NULL) path to PON/detergent}

\item{\code{create_new_pon}}{boolean (default == FALSE) whether to create a new PON from normal samples}

\item{\code{save_pon}}{boolean (default == FALSE) if create_new_pon == TRUE, whether to save pon to path given by pon_path}

\item{\code{field}}{character (default == "reads.corrected") field to use for processing}

\item{\code{use.all}}{boolean (default == TRUE) if all normal samples are to be used for creating PON}

\item{\code{choose.randomly}}{boolean (default == FALSE) if a random subset of normal samples are to be used for creating PON}

\item{\code{choose.by.clustering}}{boolean (default == FALSE) clusters normal samples based on the genomic background and takes a random sample from within the clusters}

\item{\code{number.of.samples}}{interger (default == 50) if choose.by.clustering == TRUE, this is the number of clusters at which to cut tree}

\item{\code{verbose}}{boolean (default == TRUE) outputs progress}

\item{\code{num.cores}}{integer (default == 1) number of cores to use for parallelization}

\item{\code{tolerance}}{numeric (default == 0.0001) tolerance for error for batch rPCA. We suggest keeping this value}

\item{\code{is.human}}{boolean (default == TRUE) organism type}

\item{\code{build}}{character (default == "hg19") genome build to define PAR region in chromosome X}

\item{\code{PAR.file}}{character (default == NULL) this is a GRanges with the boundaries of PAR region in X chr}

\item{\code{balance}}{boolean (default == TRUE) experimental variable to take into consideration 1 copy of X chr in male sample}

\item{\code{infer.germline}}{boolean (default = FALSE) if use the L matrix to infer germline events}

\item{\code{signal.thresh}}{numeric (default == 0.5) this is the threshold to be used to identify an amplification (markers with signal intensity > 0.5) or deletions (markers with signal intensity < -0.5) in log space from dryclean outputs}

\item{\code{pct.thresh}}{numeric (default == 0.98) proportion of samples in which a given marker is free of germline event}

\item{\code{wgs}}{boolean (default == TRUE) whether whole genome is being used}

\item{\code{target_resolution}}{numeric (default == 1e3) resolution at which to conduct analyses}

\item{\code{nochr}}{logical (default = TRUE) remove chr prefix}

\item{\code{all.chr}}{list (default = c(as.character(1:22), "X")) list of chromosomes}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-pon-get_history"></a>}}
\if{latex}{\out{\hypertarget{method-pon-get_history}{}}}
\subsection{Method \code{get_history()}}{
Function returns the history of the pon object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{pon$get_history()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
Prints the history of the pon object as data table
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-pon-get_template"></a>}}
\if{latex}{\out{\hypertarget{method-pon-get_template}{}}}
\subsection{Method \code{get_template()}}{
Function returns the template of the pon object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{pon$get_template()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
Template of the pon object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-pon-get_L"></a>}}
\if{latex}{\out{\hypertarget{method-pon-get_L}{}}}
\subsection{Method \code{get_L()}}{
Function returns the L matrix of the pon object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{pon$get_L()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
L matrix of the pon object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-pon-get_S"></a>}}
\if{latex}{\out{\hypertarget{method-pon-get_S}{}}}
\subsection{Method \code{get_S()}}{
Function returns the S matrix of the pon object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{pon$get_S()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
S matrix of the pon object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-pon-get_err"></a>}}
\if{latex}{\out{\hypertarget{method-pon-get_err}{}}}
\subsection{Method \code{get_err()}}{
Function returns the err matrix of the pon object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{pon$get_err()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
err matrix of the pon object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-pon-get_k"></a>}}
\if{latex}{\out{\hypertarget{method-pon-get_k}{}}}
\subsection{Method \code{get_k()}}{
Function returns the k matrix of the pon object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{pon$get_k()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
k matrix of the pon object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-pon-get_U_hat"></a>}}
\if{latex}{\out{\hypertarget{method-pon-get_U_hat}{}}}
\subsection{Method \code{get_U_hat()}}{
Function returns the U_hat matrix of the pon object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{pon$get_U_hat()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
U_hat matrix of the pon object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-pon-get_V_hat"></a>}}
\if{latex}{\out{\hypertarget{method-pon-get_V_hat}{}}}
\subsection{Method \code{get_V_hat()}}{
Function returns the V_hat matrix of the pon object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{pon$get_V_hat()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
V_hat matrix of the pon object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-pon-get_sigma_hat"></a>}}
\if{latex}{\out{\hypertarget{method-pon-get_sigma_hat}{}}}
\subsection{Method \code{get_sigma_hat()}}{
Function returns the sigma_hat matrix of the pon object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{pon$get_sigma_hat()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
sigma_hat matrix of the pon object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-pon-get_inf_germ"></a>}}
\if{latex}{\out{\hypertarget{method-pon-get_inf_germ}{}}}
\subsection{Method \code{get_inf_germ()}}{
Function returns the inf_germ matrix of the pon object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{pon$get_inf_germ()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
inf_germ matrix of the pon object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-pon-get_seqlengths"></a>}}
\if{latex}{\out{\hypertarget{method-pon-get_seqlengths}{}}}
\subsection{Method \code{get_seqlengths()}}{
Function returns the seqlengths of the template of pon object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{pon$get_seqlengths()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
seqlengths of template of the pon object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-pon-save_pon"></a>}}
\if{latex}{\out{\hypertarget{method-pon-save_pon}{}}}
\subsection{Method \code{save_pon()}}{
Function to save_pon if not saved
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{pon$save_pon(path_to_save = "./detergent.rds")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path_to_save}}{character (default == "./") path to save the pon object}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-pon-clone"></a>}}
\if{latex}{\out{\hypertarget{method-pon-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{pon$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
