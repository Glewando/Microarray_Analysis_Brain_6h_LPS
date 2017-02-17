dHPC gene expression 6 hours post-LPS: Array QC
================
G. Lewandowski
February 10, 2017

Generation of a GeneFeatureSet container for the **6 hours** LPS-Microarry data and sample information.
-------------------------------------------------------------------------------------------------------

### GeneFeatureSet container generation using the `oligo` package.

Set up GeneFeatureSet container for the LPS\_Microarray data. Read in sample information and use the `oligo` package to read in the raw intensity data (CEL files).

### Load CEL files using the `oligo` package and add the sample info as phenoData to the `GeneFeatureSet` object.

``` r
# load the oligo package
suppressPackageStartupMessages(library(oligo))

# list the CEL files
# the CEL file names will be used as rownames in the phenoData and as sampleNames
celfiles <- list.celfiles(file.path("LPS_Data_Files"))

# set up sample info data frame by 
# (1) reading in the sample_info.txt file (returns a data.frame object)
# (2) change the rownames to the CEL filenames (sample names),
# (3) remove the first (redundant) column

sampleInfo <- read.delim(file = file.path("LPS_Data_Files", "sample_info.txt"), 
                         header = TRUE, stringsAsFactors = FALSE)
# order the sampleInfo according to CEL file name (column x)
sampleInfo <- sampleInfo[order(sampleInfo$X),]
# change rownames from 1-12 to sample names
rownames(sampleInfo) <- substr(celfiles, 1, 5)
# change column 1 name from X to filename
colnames(sampleInfo)[1] <- "filename" 

# set up metadata with row names matching the column names of the sampleInfo df
# because these are single channel arrays, use ALL for channel value
rns <- colnames(sampleInfo)
metadata <- data.frame(labelDescription = rns, channel = factor('ALL'))

# load the CEL files; returns a "GeneFeatureSet" object
lpsrawdata <- read.celfiles(list.celfiles(file.path("LPS_Data_Files"), 
                                          full.name = TRUE)) 

# add the sample info as phenoData to lpsdataraw
# set up the pheno data for lpsdataraw
pd <- new('AnnotatedDataFrame', data = sampleInfo, varMetadata = metadata)

# set the phenoData for lpsdataraw
phenoData(lpsrawdata) <- pd

# set the sample names in the lpsrawdata container to the sample names from CEL file names
sampleNames(lpsrawdata) <- rownames(sampleInfo)

# Subset the data to obtain the 6 hour treatment samples:
lpsrawdata.6h <- lpsrawdata[, lpsrawdata$group %in% c("sal.6h", "lps.6h")]
rm(pd, rns, sampleInfo, metadata)
```

### Save the lpsrawdata GeneFeatureSet object to an external file.

``` r
save(lpsrawdata, file = file.path("processedData", "lpsrawdata.Rdata"))
save(lpsrawdata.6h, file = file.path("processedData", "lpsrawdata_6h.Rdata"))
```

Quality control analysis of the experimental arrays
---------------------------------------------------

### Intensity distribution across the arrays:

#### Intensity (density) Histograms

``` r
library(RColorBrewer)
par(mfrow = c(1,2))
usr.col <- brewer.pal(12, "Paired")
set.seed(1234)
oligo::hist(lpsrawdata.6h, which = "all", lty = rep(1, 6), 
            col = usr.col, main = "Intensity Histograms")

legend("topright", sampleNames(lpsrawdata.6h), lty = rep(1, 6), 
       col = usr.col, cex = 0.6)
```

<img src="Array_6h_QC_files/figure-markdown_github/intensity histograms-1.png" style="display: block; margin: auto;" />

The intensity distributions for all arrays are similiar.

#### Box plots:

``` r
set.seed(1234)
oligo::boxplot(lpsrawdata.6h, which = "all", 
               col = usr.col, las =3, cex.axis = 0.75, 
               names = sampleNames(lpsrawdata.6h), main = "Intensity Distribution")
```

<img src="Array_6h_QC_files/figure-markdown_github/intensity distribution boxplots-1.png" style="display: block; margin: auto;" />

The median intensity and widths of the intensity distributions are similar for all 6 arrays.

### MA plots

Plotting log-ratios, **M**, versus average log-intensities, **A** for two arrays (samples) is a strategy to visualize the relationship between two sample. For the single channel Affy Gene ST arrays, the comparison is done with pairs of samples. If two samples are from the same experimental group, then the arrays should behave in a similar manner. Meaning that the plotted data points should appear symmetrically about a horizontal line through zero. Any deviation from the zero horizontal line represents different responses of the two samples. MA plots of the raw (non-normalized) data can be used to determine if there are any array problems. MA plots are done for pairs in each of the four groups. The groups are: "sal.6h" "lps.6h".

Normalization efforts are then used to get the data points on the MA plot to be fitted symmetrically around the zero horizontal line.

``` r
# relevel groups so that sal.6h is the reference group
lpsrawdata.6h$group <- factor(lpsrawdata.6h$group)
lpsrawdata.6h$group <- relevel(lpsrawdata.6h$group, ref = "sal.6h")
groups <- (lpsrawdata.6h$group)

xl = c(2.8, 4)
yl = c(-1, 1)

# MA plot for the saline control group at 6 hours (rats 1,3,5)
MAplot(lpsrawdata.6h[, groups == "sal.6h"], main = "saline",
              pairs=TRUE, ylim = yl, xlim = xl)
```

    ## Loading required package: RSQLite

    ## Loading required package: DBI

<img src="Array_6h_QC_files/figure-markdown_github/MA plots-1.png" style="display: block; margin: auto;" />

``` r
# MA plot for the lps treated group at 6 hours (rats 2,4,24)
MAplot(lpsrawdata.6h[, groups == "lps.6h"], pairs=TRUE, 
       ylim = yl, xlim = xl, main = "LPS")
```

<img src="Array_6h_QC_files/figure-markdown_github/MA plots-2.png" style="display: block; margin: auto;" />

``` r
# MA plot for 6 h lps vs reference 6 h saline
MAplot(lpsrawdata.6h, groups=groups, which=2, 
       refSamples=1, show.statistics=FALSE)
```

<img src="Array_6h_QC_files/figure-markdown_github/MA plots-3.png" style="display: block; margin: auto;" />

### Probe level models

Looking for probe effects in the raw data; PLM fits a model and outputs a number of statistics that are used for quality diagnostics.

``` r
# Raw data no preprocessing
# summarization target = 'core', background correction = F, normalization = F

# load the plm for lpsrawdata.6h if it already exists, else, generate the plm for the raw data
plm_file <- file.path("processedData", "plm_lpsrawdata_6h.RData")

if (file.exists(plm_file)){
    load(plm_file)
} else {
    plm_fit <- fitProbeLevelModel(lpsrawdata.6h, target = "core", 
                                  background = FALSE, normalize = FALSE)
    save(plm_fit, file = plm_file)
}
```

#### Relative log expression and Normalized Unscaled Standard Error:

**Relative Log Expression (RLE)** plot is constructed from the PLM model. For each array the RLE plot is a box plot of the relative expression for each gene- which is defined as the difference of the PLM estimated gene expression an the median expression.

**Normalized Unscaled Standard Error (NUSE)** for each gene is defined as the PLM fit standard error minus the median standard error across all the arrays. The median standard error for any probeset should be 1 across all the arrays. A NUSE value significantly higher from other samples likely indicates a lower quality chip (array).

``` r
par(mfrow = c(1,2))

# box plots of RLE values from the raw data
RLE(plm_fit, type = "plot",las =3, cex.axis = 0.75, 
    names = sampleNames(lpsrawdata.6h), main = "RLE (raw data)")

# box plots of NUSE values from the raw data
NUSE(plm_fit, type = "plot",las =3, cex.axis = 0.75, 
     names = sampleNames(lpsrawdata.6h), main = "NUSE (raw data)")
```

![](Array_6h_QC_files/figure-markdown_github/RLE%20Nuse%20boxplots-1.png)

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
    ##  [1] pd.ragene.1.0.st.v1_3.14.1 DBI_0.5-1                 
    ##  [3] RSQLite_1.1-2              RColorBrewer_1.1-2        
    ##  [5] oligo_1.38.0               Biostrings_2.42.1         
    ##  [7] XVector_0.14.0             IRanges_2.8.1             
    ##  [9] S4Vectors_0.12.1           Biobase_2.34.0            
    ## [11] oligoClasses_1.36.0        BiocGenerics_0.20.0       
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

    ##            used  (Mb) gc trigger  (Mb) max used  (Mb)
    ## Ncells  5252119 280.5    8273852 441.9  6861544 366.5
    ## Vcells 31833080 242.9   94225919 718.9 82397194 628.7
