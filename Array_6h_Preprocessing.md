dHPC gene expression 6 hours post-LPS: Data pre-processing and normalization
================
G. Lewandowski
February 17, 2017

When working with Affymetrix GeneChip Gene ST arrays there are three steps necessary before indentification of differentially expressed genes (DEGs): - Background correction - Normalization - Summarization of probesets or transcripts

Using the `oligo` package in Bioconductor the three steps can be easily done sequentially.

The background-corrected/normalized expression data is then filtered to contain only the main probes (by excluding the Affy-controls) and to exclude the probes that are either not expressed or expressed at very low levels in at least 5 of the 6 samples.

### Load the required libraries

`oligo` package for preprocessing of lpsrawdata
`pd.ragene.1.0.st.v1` and `ragene10stprobeset.db` are used for annotation of the probesets on the Affymetrix Rat Gene ST v. 01 genechip.

``` r
suppressPackageStartupMessages({
    library(oligo)
    library(pd.ragene.1.0.st.v1)
    library(ragene10stprobeset.db)
    library(RColorBrewer)
    library(pander)
        })
```

### Load **lpsrawdata.6h** `GeneFeatureSet` object generated during the QC step.

``` r
load(file = file.path("processedData", "lpsrawdata_6h.Rdata")) 
```

### Preprocessing.

A classic and powerful method for preprocessing Affymetrix gene expression arrays is the RMA method (Robust Multichip Average). The RMA alogorithm applies the largest relative adjustments to the lowest probe intensities and preserves the order of probe intensities relative to the distribution of probe intensities. Applying the `rma()` function from the `oligo` package results in
(1) RMA background correction,
(2) quantile normalization and
(3) probe summarization into probesets via the median-polish method.

When applied to a GeneFeatureSet object, `rma()` can produce summaries at different levels: probeset (as defined in the PGF) and 'core genes' (as defined in the core.mps file). Affymetrix made available meta-probeset files (MPS) that define "new probesets", which allow summarization to the gene-level. For gene arrays, there's "core" MPS. We will summarize to the gene-level.

The background-corrected/normalized expression data is then filtered to contain only the main probes (by excluding the Affy-controls).

``` r
# RMA background correction, quantile normaization, 
# summarization to the core gene level

filepath <- file.path("processedData", "RMA_6h_lpsdata.Rdata")
if (file.exists(filepath)){
        load(filepath)
} else {
        RMA.lpsdata.6h <- rma(lpsrawdata.6h, background=TRUE, normalize=TRUE, target="core") # 29214 features
        
        # get the 'main' type of probesets from the Gene ST arrays
        # using the getMainProbes function from the affycoretools package
        RMA.lpsdata.6h <- suppressMessages(affycoretools::getMainProbes(RMA.lpsdata.6h))
                # 27342 features
        
        save(RMA.lpsdata.6h, file=filepath)
}   
n <- length(featureNames(RMA.lpsdata.6h))
```

There are 27342 main probes in the RMA background-corrected/normalized expression data set (hereafter referred to as `RMA.lpsdata.6h`).

#### Intensity distribution of RMA-background corrected/normalized arrays.

Histograms and boxplots of array intensities:

``` r
usr.col <- brewer.pal(6, "Dark2")
set.seed(1234)
par(mfrow = c(1,2))

hist(RMA.lpsdata.6h, transfo=identity,  
     lty = rep(1, 6), 
     lwd = rep(1.75, 6), 
     col = usr.col, 
     main = "RMA Main Probes: Histograms")
legend("topright", sampleNames(RMA.lpsdata.6h), 
       lty = rep(1, 6), lwd = rep(2, 6), col = usr.col, cex = 0.75)

set.seed(1234)
boxplot(RMA.lpsdata.6h, transfo=identity,  
        col = usr.col, 
        las =3, cex.axis = 0.75, 
        names = sampleNames(RMA.lpsdata.6h), 
        main = "RMA Main Probes: Distributionss")
```

<img src="Array_6h_Preprocessing_files/figure-markdown_github/intensity distribution_01-1.png" style="display: block; margin: auto;" />

### Gene Filtering

Probes with low expression (intensity) in 5 of 6 samples (regardless of sample group) are excluded from the expression data set. The reason for this is because the low expression probes do not add any information to the analysis, but instead add noise that can obscure real information during the selection of differentially expressed genes. Importantly, low variance genes **ARE NOT** filtered out, so that variance shrinkage and moderated t-statistic calculations are not impacted.

Gene filtering procedure:
1. determine the summary statistics for the probe expression distribution.
2. set the filter at the first quantile (25%).
3. filter the expression data to include only genes with intensity values greater than the filter value on 2 of 6 arrays.

The following table shows the minimun, first quantile, median, mean, third quantile and maximum expression values for the `RMA_lpsdata.6h` data set.

``` r
pander(df1, caption = "RMA_lpsdata.6h: Summary Statistics", 
       caption.placement = "top")
```

<table style="width:74%;">
<caption>RMA_lpsdata.6h: Summary Statistics</caption>
<colgroup>
<col width="8%" />
<col width="18%" />
<col width="12%" />
<col width="9%" />
<col width="18%" />
<col width="6%" />
</colgroup>
<thead>
<tr class="header">
<th align="center">Min</th>
<th align="center">1st Quant.</th>
<th align="center">Median</th>
<th align="center">Mean</th>
<th align="center">3rd Quant.</th>
<th align="center">Max</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center">1.72</td>
<td align="center">4.84</td>
<td align="center">6.39</td>
<td align="center">6.43</td>
<td align="center">7.99</td>
<td align="center">14.04</td>
</tr>
</tbody>
</table>

The similarity between the mean and median indicates a normal distribution of expression values. A filter value of 4.84 (first quantile) can be used to exclude low expression probes.

``` r
# filter ExpressionSet for probes with expression > first quantile 
# on at least 2 of the 6 arrays using the `genefilter` library.

filepath <- file.path("processedData", "RMA.lpsdata.6h_filt.Rdata")

if(file.exists(filepath)) {
    load(filepath)
} else {
    lowexprfilter <- firstquant
    f1 <- genefilter::kOverA(2, lowexprfilter)
    i1 <- genefilter::genefilter(RMA.lpsdata.6h, f1)
    RMA.lpsdata.6h_filt <- RMA.lpsdata.6h[i1,]
    # 20891 features
    
    save(RMA.lpsdata.6h_filt, file = filepath)
}
```

Using the low expression filter 6451 probes have been excluded from the expression set, leaving 20891 probes for further processing.

The following table shows the minimun, first quantile, median, mean, third quantile and maximum expression values for the new expression data set.

``` r
pander(df2, caption = "Filtered RMA_lpsdata.6h: Summary Statistics", 
       caption.placement = "top")
```

<table style="width:74%;">
<caption>Filtered RMA_lpsdata.6h: Summary Statistics</caption>
<colgroup>
<col width="8%" />
<col width="18%" />
<col width="12%" />
<col width="9%" />
<col width="18%" />
<col width="6%" />
</colgroup>
<thead>
<tr class="header">
<th align="center">Min</th>
<th align="center">1st Quant.</th>
<th align="center">Median</th>
<th align="center">Mean</th>
<th align="center">3rd Quant.</th>
<th align="center">Max</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center">2.85</td>
<td align="center">5.94</td>
<td align="center">7.14</td>
<td align="center">7.28</td>
<td align="center">8.39</td>
<td align="center">14.04</td>
</tr>
</tbody>
</table>

#### Inspection of the lowly expressed-filtered, RMA-processed array data.

Histograms and boxplots of array intensities:

``` r
set.seed(2345)
par(mfrow = c(1,2))

hist(RMA.lpsdata.6h_filt, transfo=identity,  
     lty = rep(1, 6), lwd = rep(1.75, 6), 
     col = usr.col, 
     main = "Filtered RMA Probes: Histograms")
legend("topright", 
       sampleNames(RMA.lpsdata.6h_filt), 
       lty = rep(1, 6), lwd = rep(2, 6), 
       col = usr.col, cex = 0.75)

set.seed(2345)
boxplot(RMA.lpsdata.6h_filt, transfo=identity,  
        col = usr.col, las =3, cex.axis = 0.75, 
        names = sampleNames(RMA.lpsdata.6h_filt), 
        main = "Filtered RMA Probes: Distributions")
```

<img src="Array_6h_Preprocessing_files/figure-markdown_github/intensity distribution_02-1.png" style="display: block; margin: auto;" />

#### MA plots:

``` r
par(mfrow = c(1,1))
yl = c(-4,2)
# MA plot for the saline control group at 6 hours (rats 1,3,5)
oligo::MAplot(RMA.lpsdata.6h_filt[, RMA.lpsdata.6h_filt$group == "sal.6h"],
              transfo=identity, pairs=TRUE, 
              main="Saline 6 hours", ylim=yl)
```

<img src="Array_6h_Preprocessing_files/figure-markdown_github/core genes within groups MA plots-1.png" style="display: block; margin: auto;" />

``` r

# MA plot for the lps treated group at 6 hours (rats 2,4,24)
MAplot(RMA.lpsdata.6h_filt[, RMA.lpsdata.6h_filt$group == "lps.6h"],
       pairs=TRUE, main ="LPS 6 hours", ylim=yl)
```

<img src="Array_6h_Preprocessing_files/figure-markdown_github/core genes within groups MA plots-2.png" style="display: block; margin: auto;" />

``` r

# MA plot for 6 h lps vs reference 6 h saline
RMA.lpsdata.6h_filt$group <- relevel(as.factor(RMA.lpsdata.6h_filt$group),
                                     ref="sal.6h")
groups <- (RMA.lpsdata.6h_filt$group)
suppressMessages(MAplot(RMA.lpsdata.6h_filt, 
                        groups=groups, which=2, refSamples=1))
## Warning in KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, :
## Binning grid too coarse for current (small) bandwidth: consider increasing
## 'gridsize'
```

<img src="Array_6h_Preprocessing_files/figure-markdown_github/core genes within groups MA plots-3.png" style="display: block; margin: auto;" />

#### Relationships between samples using heirarchical clustering:

``` r
exprDat <- exprs(RMA.lpsdata.6h_filt)
distance <- dist(t(exprDat), method = "maximum")
clusters <- hclust(distance)
plot(clusters)
abline(h = 4.5, col = "red")
```

![](Array_6h_Preprocessing_files/figure-markdown_github/heirarchical%20clustering-1.png)

#### Probe level model

Looking for probe effects in the normalized data; PLM fits a model and outputs a number of statistics that are used for quality diagnostics.

``` r
# probe level model on raw lps data: 
# background corrected and normalized (RMA)
# summarization target = 'core', background correction = T, normalization = T

filepath <- file.path("processedData", "plm.filt.RMA.6h.lps.Rdata")

# load the plm for plm_rma_6h if it already exists, else, generate the plm
if (file.exists(filepath)){
    load(filepath)
    
} else {
    plm.rma.6h <- fitProbeLevelModel(lpsrawdata.6h, target = "core", 
                                     background = TRUE, normalize = TRUE)
    
    save(plm.rma.6h, file = filepath)
}
rm(lpsrawdata.6h)
```

##### Relative log expression and Normalized Unscaled Standard Error:

**Relative Log Expression (RLE)** plot is constructed from the PLM model. For each array the RLE plot is a box plot of the relative expression for each gene- which is defined as the difference of the PLM estimated gene expression an the median expression.

**Normalized Unscaled Standard Error (NUSE)** for each gene is defined as the PLM fit standard error minus the median standard error across all the arrays. The median standard error for any probeset should be 1 across all the arrays. A NUSE value significantly higher from other samples likely indicates a lower quality chip (array).

``` r
par(mfrow = c(1,2))
# box plots of RLE values for the normalized data
RLE(plm.rma.6h, type = "plot",las =3, cex.axis = 0.75, names = sampleNames(RMA.lpsdata.6h_filt), main = "RLE")

# box plots of NUSE values for the normalized data, 
# summarized to the core genes level
NUSE(plm.rma.6h, type = "plot",las =3, cex.axis = 0.75, names = sampleNames(RMA.lpsdata.6h_filt), main = "NUSE")
```

![](Array_6h_Preprocessing_files/figure-markdown_github/plm.rma.6h%20RLE%20Nuse-1.png)

#### R session information:

    ## R version 3.3.2 (2016-10-31)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 14393)
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252 
    ## [2] LC_CTYPE=English_United States.1252   
    ## [3] LC_MONETARY=English_United States.1252
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.1252    
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] pander_0.6.0                RColorBrewer_1.1-2         
    ##  [3] ragene10stprobeset.db_8.5.0 org.Rn.eg.db_3.4.0         
    ##  [5] AnnotationDbi_1.36.2        pd.ragene.1.0.st.v1_3.14.1 
    ##  [7] DBI_0.5-1                   RSQLite_1.1-2              
    ##  [9] oligo_1.38.0                Biostrings_2.42.1          
    ## [11] XVector_0.14.0              IRanges_2.8.1              
    ## [13] S4Vectors_0.12.1            Biobase_2.34.0             
    ## [15] oligoClasses_1.36.0         BiocGenerics_0.20.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.9                BiocInstaller_1.24.0      
    ##  [3] GenomeInfoDb_1.10.3        bitops_1.0-6              
    ##  [5] iterators_1.0.8            tools_3.3.2               
    ##  [7] zlibbioc_1.20.0            digest_0.6.12             
    ##  [9] bit_1.1-12                 evaluate_0.10             
    ## [11] memoise_1.0.0              preprocessCore_1.36.0     
    ## [13] lattice_0.20-34            ff_2.2-13                 
    ## [15] Matrix_1.2-8               foreach_1.4.3             
    ## [17] yaml_2.1.14                stringr_1.1.0             
    ## [19] knitr_1.15.1               affxparser_1.46.0         
    ## [21] rprojroot_1.2              grid_3.3.2                
    ## [23] rmarkdown_1.3              magrittr_1.5              
    ## [25] backports_1.0.5            codetools_0.2-15          
    ## [27] htmltools_0.3.5            GenomicRanges_1.26.2      
    ## [29] splines_3.3.2              SummarizedExperiment_1.4.0
    ## [31] KernSmooth_2.23-15         stringi_1.1.2             
    ## [33] RCurl_1.95-4.8             affyio_1.44.0

    ##           used  (Mb) gc trigger  (Mb) max used (Mb)
    ## Ncells 3718048 198.6    8273852 441.9  6252777  334
    ## Vcells 4219837  32.2   58996814 450.2 73660482  562
