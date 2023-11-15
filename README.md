[![Build Status](https://app.travis-ci.com/mskilab/dryclean.svg?branch=master)](https://app.travis-ci.com/mskilab/dryclean)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/dryclean.svg)](https://codecov.io/github/mskilab/dryclean?branch=master)

### Dockerized installation

To make our life easier, we have created a `Docker` container with the latest stable release of Dryclean and its dependencies. This can be found [here](https://hub.docker.com/r/mskilab/dryclean/tags/). The latest updated version is `0.0.2`, so make sure to select the correct `tag`.

---
title: dryclean tutorial
---

![](inst/extdata/DNAhanger.png)

### <font color=black> Robust PCA based method to de-noise genomic coverage data.</font>

## <font color=black> Installations </font>

Install devtools from CRAN


```R
install.packages('devtools')
```

Set this to allow dependencies that throw warnings to be installed.


```R
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)
```

Install dependent packages and latest Bioconductor (if you haven't already)


```R
source('https://bioconductor.org/biocLite.R')
biocLite('GenomicRanges')
```

Install mskilab R dependencies (gUtils)


```R
devtools::install_github('mskilab/gUtils')
```

Install dryclean


```R
devtools::install_github('mskilab-org/dryclean')
```

(after installing R package) Add dryclean directory to PATH and test the executable 

```{bash}
$ export PATH=${PATH}:$(Rscript -e 'cat(paste0(installed.packages()["dryclean", "LibPath"], "/dryclean/extdata/"))')
$ drcln -h ## to see the help message
```

## <font color=black> Tutorial </font>

Dryclean is a robust principal component analysis (rPCA) based method. It uses a panel of normal (PON) samples to learn the landscape of both biological and technical noise in read depth data. Dryclean then uses this landscape to significantly reduce noise and artifacts in the signal for tumor samples. The input to the algorithm is a [GenomicsRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) object containing read depth. You can use read counts from your favorite tool (there are many fast tools out there, for example:  [megadepth](https://bioconductor.org/packages/release/bioc/html/megadepth.html)). Using uncorrected read counts as input for Dryclean works well based on our experience, but if you wish, you can use the GC and mappability corrected read depth data from fragCounter, which can be found at: https://github.com/mskilab/fragCounter.

###  <font color=black> 1. Creating Panel of Normal aka detergent </font>

There are 2 options for instantiating the PON object: 

Option 1: Load an existing PON from a path.

To load an existing PON from the path into the <code>pon_object</code>, run:

```R
pon_object = pon$new(pon_path = "~/git/dryclean/inst/extdata/detergent.rds")
```
Option 2: Create a new PON from normal samples.

To create a new PON, the vector with paths to the normal samples is needed.

Following is an example of such a vector


```R
normal_vector_example = readRDS("~/git/dryclean/inst/extdata/normal_vector.rds")
normal_vector_example

[1] "/git/dryclean/extdata/samp1.rds"
[2] "/git/dryclean/extdata/samp2.rds"
[3] "/git/dryclean/extdata/samp3.rds"
```

To make a new PON, you need to instantiate a PON object and set <code>create_new_pon = TRUE</code>. 

```R
pon_object = pon$new(
    create_new_pon = TRUE, 
    normal_vector = normal_vector_example
    )
```

The parameters that could be used in PON generation:

<table style="border: 1px solid black; border-collapse: collapse;">
  <tbody>
    <tr>
      <th style="border: 1px solid black; padding: 5px;">Parameter</th>
      <th style="border: 1px solid black; padding: 5px;">Default value</th>
      <th style="border: 1px solid black; padding: 5px;">Description</th>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">create_new_pon</td>
      <td style="border: 1px solid black; padding: 5px;">FALSE</td>
      <td style="border: 1px solid black; padding: 5px;">Whether to create a new PON from normal samples</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">pon_path</td>
      <td style="border: 1px solid black; padding: 5px;">NULL</td>
      <td style="border: 1px solid black; padding: 5px;">If create_new_pon==FALSE, the path to the existing PON; If create_new_pon==TRUE and save_pon == TRUE, the path to save the new PON</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">normal_vector</td>
      <td style="border: 1px solid black; padding: 5px;">c()</td>
      <td style="border: 1px solid black; padding: 5px;">Vector of paths to normal samples</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">save_pon</td>
      <td style="border: 1px solid black; padding: 5px;">FALSE</td>
      <td style="border: 1px solid black; padding: 5px;">If create_new_pon==TRUE, whether to save pon to the path given by pon_path</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">field</td>
      <td style="border: 1px solid black; padding: 5px;">"reads.corrected"</td>
      <td style="border: 1px solid black; padding: 5px;">Field name in GRanges metadata of normal samples to use for PON generation</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">use.all</td>
      <td style="border: 1px solid black; padding: 5px;">TRUE</td>
      <td style="border: 1px solid black; padding: 5px;">Whether all normal samples are to be used for creating PON</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">choose.randomly</td>
      <td style="border: 1px solid black; padding: 5px;">FALSE</td>
      <td style="border: 1px solid black; padding: 5px;">If use.all==FALSE, whether a random subset of normal samples are to be used for creating PON</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">choose.by.clustering</td>
      <td style="border: 1px solid black; padding: 5px;">FALSE</td>
      <td style="border: 1px solid black; padding: 5px;">Whether to cluster normal samples based on the genomic background and take a random sample from within the clusters</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">number.of.samples</td>
      <td style="border: 1px solid black; padding: 5px;">50</td>
      <td style="border: 1px solid black; padding: 5px;">If choose.by.clustering==TRUE or choose.randomly==TRUE, the number of clusters/samples to use</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">tolerance</td>
      <td style="border: 1px solid black; padding: 5px;">0.0001</td>
      <td style="border: 1px solid black; padding: 5px;">Tolerance for error for batch rPCA; we suggest keeping this value</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">num.cores</td>
      <td style="border: 1px solid black; padding: 5px;">1</td>
      <td style="border: 1px solid black; padding: 5px;">Number of cores to use for parallelization</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">verbose</td>
      <td style="border: 1px solid black; padding: 5px;">TRUE</td>
      <td style="border: 1px solid black; padding: 5px;">Whether to output progress</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">is.human</td>
      <td style="border: 1px solid black; padding: 5px;">TRUE</td>
      <td style="border: 1px solid black; padding: 5px;">Organism type</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">build</td>
      <td style="border: 1px solid black; padding: 5px;">"hg19"</td>
      <td style="border: 1px solid black; padding: 5px;">Genome build to define PAR region in chromosome X</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">PAR.file</td>
      <td style="border: 1px solid black; padding: 5px;">NULL</td>
      <td style="border: 1px solid black; padding: 5px;">GRanges with the boundaries of PAR region in X chr</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">balance</td>
      <td style="border: 1px solid black; padding: 5px;">TRUE</td>
      <td style="border: 1px solid black; padding: 5px;">Experimental variable to take into consideration 1 copy of X chr in male sample</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">infer.germline</td>
      <td style="border: 1px solid black; padding: 5px;">FALSE</td>
      <td style="border: 1px solid black; padding: 5px;">Whether to use the L matrix to infer germline events</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">signal.thresh</td>
      <td style="border: 1px solid black; padding: 5px;">0.5</td>
      <td style="border: 1px solid black; padding: 5px;">The threshold to be used to identify an amplification (markers with signal intensity > 0.5) or deletions (markers with signal intensity < -0.5) in log space from dryclean outputs</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">pct.thresh</td>
      <td style="border: 1px solid black; padding: 5px;">0.98</td>
      <td style="border: 1px solid black; padding: 5px;">Proportion of samples in which a given marker is free of germline event</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">wgs</td>
      <td style="border: 1px solid black; padding: 5px;">TRUE</td>
      <td style="border: 1px solid black; padding: 5px;">Whether whole genome is being used</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">target_resolution</td>
      <td style="border: 1px solid black; padding: 5px;">1,000</td>
      <td style="border: 1px solid black; padding: 5px;">Desired bin size of the PON</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">nochr</td>
      <td style="border: 1px solid black; padding: 5px;">TRUE</td>
      <td style="border: 1px solid black; padding: 5px;">Whether to remove chr prefix</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">all.chr</td>
      <td style="border: 1px solid black; padding: 5px;">c(as.character(1:22), "X")</td>
      <td style="border: 1px solid black; padding: 5px;">List of chromosomes</td>
    </tr>
  </tbody>
</table>

<br>

The pon_object contains the following methods:

<table>
<tbody>
<tr><td>1. <code>get_L()</code> - returns L, the low ranked matrix of all the PONs calculated by batch robust PCA method</td></tr>
<tr><td>2. <code>get_S()</code> - returns S, the sparse matrix of all the PONs calculated by batch robust PCA method</td></tr>
<tr><td>3. <code>get_k()</code> - returns k, the estimated rank of a matrix where coverage values from each normal sample forms a column</td></tr>
<tr><td>4. <code>get_U_hat()</code> - returns U.hat, svd decompsed left sigular matrix of L required for online implentation of rPCA</td></tr>
<tr><td>5. <code>get_V_hat()</code> - returns V.hat, svd decompsed right sigular matrix of L required for online implentation of rPCA</td></tr>
<tr><td>6. <code>get_sigma_hat()</code> - returns sigma.hat, svd decompsed first k sigular values of L required for online implentation of rPCA</td></tr>
<tr><td>7. <code>get_inf_germ()</code> - returns inf.germ, the inferred germline obtained from the normal samples</td></tr>
<tr><td>8. <code>get_seqlengths()</code> - returns seqlengths of each chromosome of the PON objects</td></tr>
<tr><td>9. <code>get_history()</code> - returns the history of actions on the pon object with timestamps</td></tr>
</tbody>
</table>


###  <font color=black> 2. Normalizing the coverage aka drycleaning </font>

Following is a dummy example. The data directory has a dummy coverage gRanges object with "reads.corrected" field.


```R
coverage_file = readRDS("~/git/dryclean/inst/extdata/dummy_coverage.rds")
coverage_file
```

```R
GRanges object with 50 ranges and 1 metadata column:
       seqnames    ranges strand | reads.corrected
          <Rle> <IRanges>  <Rle> |       <numeric>
   [1]       22       1-3      * |        2.869742
   [2]       22       3-5      * |        2.221168
   [3]       22       5-7      * |        3.576461
   [4]       22       7-9      * |        3.289552
   [5]       22      9-11      * |        0.013421
   ...      ...       ...    ... .             ...
  [46]       22     91-93      * |         1.89621
  [47]       22     93-95      * |         4.16527
  [48]       22     95-97      * |         1.22947
  [49]       22     97-99      * |         2.79558
  [50]       22    99-101      * |         1.88191
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

In order to run dryclean, instantiate a dryclean class object first with previously created pon object.


```R
dryclean_object <- dryclean$new(pon = pon_object)

```

After initializing the dryclean object, use the <code>clean</code> function to normalize the coverage with the path to the coverage data as the required <code>cov</code> parameter. For the sake of example, we set the parameter <code>testing=TRUE</code> but typically, you would leave it at its default value.


```R
dryclean_object$clean(cov = "~/git/dryclean/inst/extdata/dummy_coverage.rds", testing = TRUE)

```

```R

GRanges object with 50 ranges and 7 metadata columns:
       seqnames    ranges strand | background.log foreground.log input.read.counts median.chr foreground background log.reads
          <Rle> <IRanges>  <Rle> |      <numeric>      <numeric>         <numeric>  <numeric>  <numeric>  <numeric> <numeric>
   [1]       22       1-3      * |    -0.00161363      0.0130496        1.16515870    1.04791 1.01313509   0.998388  0.152857
   [2]       22       3-5      * |    -0.00515368      0.0000000        0.90182764    1.04791 1.00000000   0.994860 -0.103332
   [3]       22       5-7      * |    -0.00162024      0.2332078        1.45209733    1.04791 1.26264386   0.998381  0.373009
   [4]       22       7-9      * |    -0.00176592      0.1497312        1.33560805    1.04791 1.16152201   0.998236  0.289387
   [5]       22      9-11      * |    -0.00147886     -5.0694027        0.00544911    1.04791 0.00628617   0.998522 -5.212303
   ...      ...       ...    ... .            ...            ...               ...        ...        ...        ...       ...
  [46]       22     91-93      * |    -0.00161919      -0.118465          0.769892    1.04791   0.888283   0.998382 -0.261505
  [47]       22     93-95      * |    -0.00515368       0.389148          1.691161    1.04791   1.475722   0.994860  0.525415
  [48]       22     95-97      * |    -0.00515368      -0.548208          0.499183    1.04791   0.577985   0.994860 -0.694783
  [49]       22     97-99      * |    -0.00176489       0.000000          1.135049    1.04791   1.000000   0.998237  0.126676
  [50]       22    99-101      * |    -0.00176960      -0.125884          0.764086    1.04791   0.881717   0.998232 -0.269075
  -------
  seqinfo: 1 sequence from an unspecified genome

```

The output has following metadata fields: 

<table>
<tbody>
<tr><td>1. background.log: This is the L low ranked vector after decomposition and represent the background noise separated by dryclean in the log space</td></tr>
<tr><td>2. foreground.log: The S vector with the inferred copy number signal separated by dryclean, that forms foreground, in the log space</td></tr>
<tr><td>3. input.read.counts: This is the mean-normalized count input in linear space</td></tr>
<tr><td>4. median.chr: median chromosome signal</td></tr>
<tr><td>5. foreground: Foreground signal, that forms SCNAs (S vector) in read count/ratio space</td></tr>
<tr><td>6. background: This is the L low ranked vector after decomposition and represent the background noise separated by dryclean in read count/ratio space </td></tr>
<tr><td>7. log.reads: log of the mean-normalized count</td></tr>
</tbody>
</table>

The parameters that can be used in clean() function:

<table style="border: 1px solid black; border-collapse: collapse;">
  <tbody>
    <tr>
      <th style="border: 1px solid black; padding: 5px;">Parameter</th>
      <th style="border: 1px solid black; padding: 5px;">Default value</th>
      <th style="border: 1px solid black; padding: 5px;">Description</th>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">cov</td>
      <td style="border: 1px solid black; padding: 5px;">REQUIRED</td>
      <td style="border: 1px solid black; padding: 5px;">Path to the granges coverage file to be normalized</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">field</td>
      <td style="border: 1px solid black; padding: 5px;">"reads.corrected"</td>
      <td style="border: 1px solid black; padding: 5px;">Field name in GRanges metadata of coverage to use for drycleaning</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">centered</td>
      <td style="border: 1px solid black; padding: 5px;">FALSE</td>
      <td style="border: 1px solid black; padding: 5px;">Whether a coverage has already been mean-normalized</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">cbs</td>
      <td style="border: 1px solid black; padding: 5px;">FALSE</td>
      <td style="border: 1px solid black; padding: 5px;">Whether to perform cbs on the drycleaned coverage</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">cnsignif</td>
      <td style="border: 1px solid black; padding: 5px;">1e-5</td>
      <td style="border: 1px solid black; padding: 5px;">The significance levels for the tests in cbs to accept change-points</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">mc.cores</td>
      <td style="border: 1px solid black; padding: 5px;">1</td>
      <td style="border: 1px solid black; padding: 5px;">Number of cores to use for parallelization</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">use.blacklist</td>
      <td style="border: 1px solid black; padding: 5px;">FALSE</td>
      <td style="border: 1px solid black; padding: 5px;">Whether to exclude off-target markers in case of Exomes or targeted sequencing; If set to TRUE, needs a GRange marking if each marker is set to be excluded or not</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">blacklist_path</td>
      <td style="border: 1px solid black; padding: 5px;">NA</td>
      <td style="border: 1px solid black; padding: 5px;">If use.blacklist == TRUE, path a GRanges object marking if each marker is set to be excluded or not</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">germline.filter</td>
      <td style="border: 1px solid black; padding: 5px;">FALSE</td>
      <td style="border: 1px solid black; padding: 5px;">Whether germline markers need to be removed from decomposition</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">germline.file</td>
      <td style="border: 1px solid black; padding: 5px;">NA</td>
      <td style="border: 1px solid black; padding: 5px;">Path to file with germline markers</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">verbose</td>
      <td style="border: 1px solid black; padding: 5px;">TRUE</td>
      <td style="border: 1px solid black; padding: 5px;">Outputs progress</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">is.human</td>
      <td style="border: 1px solid black; padding: 5px;">TRUE</td>
      <td style="border: 1px solid black; padding: 5px;">Organism type</td>
    </tr>
    <tr>
      <td style="border: 1px solid black; padding: 5px;">all.chr</td>
      <td style="border: 1px solid black; padding: 5px;">c(as.character(1:22), "X")</td>
      <td style="border: 1px solid black; padding: 5px;">List of chromosomes</td>
    </tr>
  </tbody>
</table>

<br>

For 'dryclean' to work correctly, the lengths of each sequence on each chromosome in the coverage and PON (Panel of Normal) data must match. If you attempt to normalize the coverage with PON data of different sequence lengths, you will encounter an error. In the event of such an error, you can utilize the <code>get_mismatch()</code> method to obtain a data table of all chromosomes with mismatched lengths. Additionally, you can use the <code>get_history()</code> method to review all actions performed on the object with timestamps.



###  <font color=black> 3. Running `dryclean` on tumor sample from command line</font>

Dryclean CLI offers two modes of operation: PON generation and normalization using an existing PON.

Mode 1 (Default): Coverage Normalization with an existing PON. To select this mode, set <code>--mode 'coverage'</code>. In this mode, Dryclean employs an existing PON specified by <code>--pon</code> to normalize the GRanges coverage provided with <code>--input</code>. The normalized coverage is saved as GRanges in the <code>--outdir</code> directory (default = './').

Example (Note: <code>--testing TRUE</code> only for example purposes; typically, you would use the default value):


```R
./drcln --input inst/extdata/samp1.rds --pon inst/extdata/detergent.rds --testing TRUE
```

```R

▓█████▄   ██▀███  ██   ██▓ ▄████▄   ██▓    ▓█████ ▄▄▄       ███▄    █
 ██▀ ██▌ ▓██   ██  ██  ██  ██▀ ▀█  ▓██▒    ▓█   ▀ ████▄     ██ ▀█   █
░██   █▌ ▓██ ░▄█    ██ ██  ▓█    ▄  ██░    ░███   ██  ▀█▄   ██  ▀█ ██▒
░▓█▄   ▌ ▒██▀▀█▄   ░ ▐██▓ ▒▓▓▄ ▄██▒ ██░    ░▓█  ▄ ██▄▄▄▄█   ██▒  ▐▌██▒
░▒████▓  ░██▓  ██  ░ ██▒    ▓███▀ ░░█████ ▒█████▒ █     █▒ ██░   ▓██░
 ▒ ▓  ▒  ░  ▓ ░▒▓░  ██    ░ ░▒ ▒  ░░ ▒░▓  ░░░ ▒░ ░▒▒   ▓▒█░░ ▒░   ▒ ▒
 ░ ▒  ▒    ░▒ ░  ░  ░░▒░   ░  ▒   ░ ░ ▒  ░ ░ ░  ░ ▒   ▒▒ ░░ ░░   ░ ▒░
 ░ ░  ░    ░░   ░   ░  ░░  ░          ░ ░  ░    ░    ░   ▒      ░   ░ ░
   ░        ░     ░ ░     ░ ░          ░  ░   ░  ░     ░  ░     ░   ░
 ░               ░ ░     ░       ░    ░     ░     ░      ░     ░ 


(Let's dryclean the genomes!)

ℹ Loading dryclean
Loading PON...
PON loaded
Loading coverage
Loading PON a.k.a detergent
Let's begin, this is whole exome/genome
Centering the sample
Initializing wash cycle
Using the detergent provided to start washing
lambdas calculated
calculating A and B
calculating v and s
Calculating b
Combining matrices with gRanges
Giddy Up!

```

Mode 2: PON Generation. To select this mode, set <code>--mode 'pon'</code>. In this mode, a new Panel of Normals (PON) is generated using a vector of normal samples saved as .rds, specified with <code>--normal_vector</code> flag. The newly created PON is then saved in the <code>--outdir</code> directory (default = './'). 

Example: 

```R
./drcln --mode "pon" --normal_vector inst/extdata/normal_vector.rds

```

```R

▓█████▄   ██▀███  ██   ██▓ ▄████▄   ██▓    ▓█████ ▄▄▄       ███▄    █
 ██▀ ██▌ ▓██   ██  ██  ██  ██▀ ▀█  ▓██▒    ▓█   ▀ ████▄     ██ ▀█   █
░██   █▌ ▓██ ░▄█    ██ ██  ▓█    ▄  ██░    ░███   ██  ▀█▄   ██  ▀█ ██▒
░▓█▄   ▌ ▒██▀▀█▄   ░ ▐██▓ ▒▓▓▄ ▄██▒ ██░    ░▓█  ▄ ██▄▄▄▄█   ██▒  ▐▌██▒
░▒████▓  ░██▓  ██  ░ ██▒    ▓███▀ ░░█████ ▒█████▒ █     █▒ ██░   ▓██░
 ▒ ▓  ▒  ░  ▓ ░▒▓░  ██    ░ ░▒ ▒  ░░ ▒░▓  ░░░ ▒░ ░▒▒   ▓▒█░░ ▒░   ▒ ▒
 ░ ▒  ▒    ░▒ ░  ░  ░░▒░   ░  ▒   ░ ░ ▒  ░ ░ ░  ░ ▒   ▒▒ ░░ ░░   ░ ▒░
 ░ ░  ░    ░░   ░   ░  ░░  ░          ░ ░  ░    ░    ░   ▒      ░   ░ ░
   ░        ░     ░ ░     ░ ░          ░  ░   ░  ░     ░  ░     ░   ░
 ░               ░ ░     ░       ░    ░     ░     ░      ░     ░ 


(Let's dryclean the genomes!)

ℹ Loading dryclean
Loading PON...
WARNING: New PON will be generated and saved at ./pon.rds

Giving you some time to think...

Starting the preparation of Panel of Normal samples a.k.a detergent
3 samples available
Using all samples
PAR file not provided, using hg19 default. If this is not the correct build, please provide a GRanges object delineating for corresponding build
PAR read
Checking for existence of files
3 files present
Starting decomposition
This is version 2
Finished making the PON
Finished saving the PON to the provided path
PON loaded
Giddy Up!

```

All CLI options:

```R
./drcln -h

```

```R

Options:
	--mode=MODE
		Mode of operation: 'pon' or 'coverage'. Set to 'pon' for PON generation and 'coverage' for normalizing a sample using existing PON

	-p PON, --pon=PON
		path to the existing Panel Of Normal (PON) saved as .rds

	-i INPUT, --input=INPUT
		path to the coverage file in GRanges format saved as .rds

	-t CENTERED, --centered=CENTERED
		are the samples centered

	-s CBS, --cbs=CBS
		whether to perform cbs on the drycleaned coverage

	-n CNSIGNIF, --cnsignif=CNSIGNIF
		the significance levels for the test to accept change-points in cbs

	-c CORES, --cores=CORES
		number of cores to use

	-w WHOLEGENOME, --wholeGenome=WHOLEGENOME
		whether whole genome is being used

	-b BLACKLIST, --blacklist=BLACKLIST
		whether there are blacklisted makers

	-l BLACKLIST_PATH, --blacklist_path=BLACKLIST_PATH
		if --blacklist == TRUE, path to a GRanges object marking if each marker is set to be excluded or not

	-g GERMLINE.FILTER, --germline.filter=GERMLINE.FILTER
		if PON based germline filter is to be used for removing some common germline events, if set to TRUE, give path to germline annotated file

	-f GERMLINE.FILE, --germline.file=GERMLINE.FILE
		path to file annotated with germline calls, if germline.filter == TRUE

	-m HUMAN, --human=HUMAN
		whther the samples under consideration are human

	-F FIELD, --field=FIELD
		field name in GRanges metadata to use for drycleaning

	-C ALL.CHR, --all.chr=ALL.CHR
		list of chromosomes to dryclean

	-B BUILD, --build=BUILD
		hg19/hg38 build for human samples

	-T TESTING, --testing=TESTING
		DO NOT CHANGE

	--normal_vector=NORMAL_VECTOR
		if mode = 'pon', path to a vector containing normal coverages in GRanges format saved as .rds

	--field_pon=FIELD_PON
		field name in GRanges metadata of normal samples to use for pon generation

	-o OUTDIR, --outdir=OUTDIR
		output directory

	-h, --help
		show this help message and exit
```


## <font color=black> Panel of Normal for 1kb WGS (hg19) </font>

The Panel of Normal samples (PON) of 395 TCGA WGS normal samples was created using hierarchical clustering approach described above and filtered for CNPs.

The file is 17G in size. 

WGS 1 kb PON: https://mskilab.s3.amazonaws.com/hg19/WGS/detergent.rds
