[![Build Status](https://travis-ci.org/mskilab/dryclean.svg?branch=master)](https://travis-ci.org/mskilab/dryclean)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/dryclean.svg)](https://codecov.io/github/mskilab/dryclean?branch=master)

# <font color=black> dryclean </font>

### <font color=black> Robust PCA based method to de-noise genomic coverage data.</font>

## <font color=black> Installations </font>

Install devtools from CRAN


```R
install.packages('devtools')
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
devtools::install_github('mskilab/dryclean')
```

## <font color=black> Tutorial </font>

This tool is for read depth normalization robustly and independent of paired normal sample. There are two parts to it, first part requires batch rPCA decomposition of a Panel of Normals we call "solvent" (keeping true to the dryclean analogy!). We will make a standard PON available from TCGA data for whole exomes and whole genomes soon. You can also create your own solvents but this can be expensive in terms of time and memory for whole genomes. We highly suggest using our TCGA PONs for WGS. The second part is online implemetation of rPCA that decompses the tumor sample in to $S$ and $L$ where the former matrix captures copy number variation and latter captures thechnical noise and background copy number based on the PON.

Important note: A prerequisite for running dryclean is GC correction by our method called fragCounter which will become a part of dryclean in near future update.

###  <font color=black> 1. Creating solvent </font>
library(dryclean)setwd("~/git/dryclean/data")
If you were to create your own PON for dryclean to use, you will need to run fragCounter on your samples. The output from fragCounter is used to create solvent. You should have all your normals in one directory. The third thing you need is a data.table with two columns:
1. "paris" column contains the sample name you will use to index the sample
2. "normal_cov" is a column with paths to the normal samples to be used

Following is an example of such a table


```R
normal_table_example = readRDS("~/git/dryclean/data/normal_table.rds")
normal_table_example
```


<table>
<thead><tr><th scope=col>pair</th><th scope=col>normal_cov</th></tr></thead>
<tbody>
	<tr><td>samp1                        </td><td>~/git/dryclean/data/samp1.rds</td></tr>
	<tr><td>samp2                        </td><td>~/git/dryclean/data/samp2.rds</td></tr>
	<tr><td>samp3                        </td><td>~/git/dryclean/data/samp3.rds</td></tr>
</tbody>
</table>



Once you have all the essential elements, run phase 1 of creating solvents.
In this pahse, each GRange object is read in and pertinent columns extracted. A mtrix is created with each sample as a column and rPCA is carried out. We use randomized version

Note: The `rrPCA()` is purposefully set to verbose mode as it displays the "target rank" of matrix. This is not a standard output at this point and needs to be noted for second phase. This is crucial for getting pahse2 to work


```R
solvent = prepare_solvent_phase1(normal_table_path = "~/git/dryclean/data/normal_table.rds", mc.cores = 1)
```

    Starting the prep for first phase requiring randomized rPCA
    This process may take time depending on dimensions of input


    
     Iteration: 1  predicted rank = 1  target rank k = 2  Fro. error = 0.418420723344439
     Iteration: 2  predicted rank = 1  target rank k = 2  Fro. error = 0.0606708071355258
     Iteration: 3  predicted rank = 1  target rank k = 2  Fro. error = 0.0159516219506272
     Iteration: 4  predicted rank = 1  target rank k = 2  Fro. error = 0.00837253226261099
     Iteration: 5  predicted rank = 1  target rank k = 2  Fro. error = 0.00687830281345693
     Iteration: 6  predicted rank = 1  target rank k = 2  Fro. error = 0.0053365739356178
     Iteration: 7  predicted rank = 1  target rank k = 2  Fro. error = 0.0036876138150812
     Iteration: 8  predicted rank = 1  target rank k = 2  Fro. error = 0.00190013440514129
     Iteration: 9  predicted rank = 1  target rank k = 2  Fro. error = 0.000783951852238052
     Iteration: 10  predicted rank = 1  target rank k = 2  Fro. error = 0.000340787584414362
     Iteration: 11  predicted rank = 1  target rank k = 2  Fro. error = 0.000274445726843516
     Iteration: 12  predicted rank = 1  target rank k = 2  Fro. error = 0.000263312013672022
     Iteration: 13  predicted rank = 1  target rank k = 2  Fro. error = 0.000197490734067693
     Iteration: 14  predicted rank = 1  target rank k = 2  Fro. error = 0.000117371650209188
     Iteration: 15  predicted rank = 1  target rank k = 2  Fro. error = 6.12925948234308e-05
     Iteration: 16  predicted rank = 1  target rank k = 2  Fro. error = 3.88266083109873e-05
     Iteration: 17  predicted rank = 1  target rank k = 2  Fro. error = 3.18871504840978e-05
     Iteration: 18  predicted rank = 1  target rank k = 2  Fro. error = 2.41617367212416e-05
     Iteration: 19  predicted rank = 1  target rank k = 2  Fro. error = 1.51291492408584e-05
     Iteration: 20  predicted rank = 1  target rank k = 2  Fro. error = 7.95531579597329e-06


```R
names(solvent)
head(solvent$L)
```


<ol class=list-inline>
	<li>'L'</li>
	<li>'S'</li>
	<li>'err'</li>
</ol>




<table>
<tbody>
	<tr><td>0.1186491</td><td>0.1293625</td><td>0.1181923</td></tr>
	<tr><td>0.3428964</td><td>0.3738583</td><td>0.3415761</td></tr>
	<tr><td>0.1189428</td><td>0.1296828</td><td>0.1184849</td></tr>
	<tr><td>0.1230901</td><td>0.1342045</td><td>0.1226161</td></tr>
	<tr><td>0.1094483</td><td>0.1193309</td><td>0.1090269</td></tr>
	<tr><td>0.1233053</td><td>0.1344392</td><td>0.1228305</td></tr>
</tbody>
</table>



Phase2 consists of carrying SVD on subspace matrix $L$ which is necessary to carry on the online uodate of subspace later on. We use randomizec version of SVD. 

As noted earlier, at this stage, user must note the final `k` printed from phase1 to be used in this pahse


```R
solvent = prepare_solvent_phase2(solvent = solvent, k = 2) #Note that k needs to be set and this is printed in the following output
```

    Starting the prep for second pahse involving svd
    Starting randomizec SVD on L matrix

names(solvent)'L' 'S' 'err' 'k' 'U.hat' 'V.hat' 'sigma.hat'
Now we are ready for tumor decomposition

###  <font color=black> 2. Running `dryclean` on tumor sample </font>

Following is a dummy example. The data diretory has a dummy coverage gRanges object which requires "reads.corrected" field 


```R
coverage_file = readRDS("~/git/dryclean/data/dummy_coverage.rds")
coverage_file
```


    GRanges object with 100 ranges and 1 metadata column:
            seqnames    ranges strand |    reads.corrected
               <Rle> <IRanges>  <Rle> |          <numeric>
        [1]        1   [1, 10]      * |   2.86974197952077
        [2]        1   [1, 10]      * |   2.22116750082932
        [3]        1   [1, 10]      * |   3.57646101620048
        [4]        1   [1, 10]      * |   3.28955231001601
        [5]        1   [1, 10]      * | 0.0134209531825036
        ...      ...       ...    ... .                ...
       [96]        1   [1, 10]      * |   3.57150336261839
       [97]        1   [1, 10]      * |   3.80716656334698
       [98]        1   [1, 10]      * | 0.0819604389835149
       [99]        1   [1, 10]      * |   3.99150879820809
      [100]        1   [1, 10]      * |   2.46343192528002
      -------
      seqinfo: 1 sequence from an unspecified genome; no seqlengths


The solvent is a list with the following elements: 
1. $L$: This is the $L$ low ranked matrix of all the PONs calculated by batch robust PCA mathod
2. $S$: This is the $S$ sparse matrix of all the PONs calculated by batch robust PCA mathod
3. $k$: This is estimated rank of a matrix where coverage values from each normal sample forms a column
4. $U.hat$: svd decompsed left sigular matrix of $L$ required for online implentation of rPCA
5. $V.hat$: svd decompsed right sigular matrix of $L$ required for online implentation of rPCA
6. $sigma.hat$: svd decompsed first $k$ sigular values of $L$ required for online implentation of rPCA

solvents can be created or one can used precomputed ones provided by us


```R
solvent = readRDS("~/git/dryclean/data/rpca.burnin.chr1.rds")
names(solvent)
```


<ol class=list-inline>
	<li>'L'</li>
	<li>'S'</li>
	<li>'k'</li>
	<li>'U.hat'</li>
	<li>'V.hat'</li>
	<li>'sigma.hat'</li>
</ol>



In order to run dryclean, simply invoke the following function


```R
cov_out = start_wash_cycle(cov = coverage_file, burnin.samples.path = "~/git/dryclean/data", whole_genome = FALSE, chr = "1")
```


```R
head(cov_out)
```


    GRanges object with 6 ranges and 6 metadata columns:
          seqnames    ranges strand |                    L                 S
             <Rle> <IRanges>  <Rle> |            <numeric>         <numeric>
      [1]        1   [1, 10]      * | -0.00330337586216281 0.520359271078844
      [2]        1   [1, 10]      * |  -0.0109178359911416  1.43994339613236
      [3]        1   [1, 10]      * |   0.0139477201074875  1.29967267859486
      [4]        1   [1, 10]      * |  -0.0175408901210239  1.06422653415043
      [5]        1   [1, 10]      * | -0.00953444080993428 -4.55653813644987
      [6]        1   [1, 10]      * | 0.000743513472141894  1.27843817478798
             reads.corrected                 S1                L1         log.reads
                   <numeric>          <numeric>         <numeric>         <numeric>
      [1]   2.86974197952077   1.68263206215469 0.996702074280938  1.05422212312404
      [2]   2.22116750082932   4.22045691605059 0.989141547271502 0.798032958921047
      [3]   3.57646101620048   3.66809582481963   1.0140454433659  1.27437376852101
      [4]   3.28955231001601   2.89859615165932 0.982612055717719  1.19075147953513
      [5] 0.0134209531825036 0.0104983399276166 0.990510867858899 -4.31093812294871
      [6]   3.07809003512375   3.59102678732255   1.0007437899468  1.12430928616619
      -------
      seqinfo: 25 sequences from an unspecified genome

