---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- 如果要生成github主页上的README.md, 需要将此文件复制到R包的主目录下,然后设置for_github参数为TRUE,然后knit,运行完成后删除主目录下的README.Rmd文件 -->






# SeuratExplorer




<!-- badges: start -->
[![](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![](https://www.r-pkg.org/badges/version/SeuratExplorer)](https://cran.r-project.org/package=SeuratExplorer)
[![](https://img.shields.io/badge/devel%20version-0.1.2-rossellhayes.svg)](https://github.com/fentouxungui/SeuratExplorer)
[![](http://cranlogs.r-pkg.org/badges/grand-total/SeuratExplorer)](https://cran.r-project.org/package=SeuratExplorer)
[![](https://img.shields.io/github/languages/code-size/fentouxungui/SeuratExplorer.svg)](https://github.com/fentouxungui/SeuratExplorer)
[![AskDeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/fentouxungui/SeuratExplorer)
[![AskZreadAI]()](https://zread.ai/fentouxungui/SeuratExplorer)
<!-- badges: end -->

> An ``Shiny`` App for Exploring scRNA-seq Data Processed in ``Seurat``

A simple, one-command package which runs an interactive dashboard capable of common visualizations for single cell RNA-seq. ``SeuratExplorer`` requires a processed ``Seurat`` object, which is saved as ``rds`` or ``qs2`` file.

## Why build this R package

> Currently, there is still no good tools for visualising the analysis results from ``Seurat``, when the bioinformatics analyst hands over the results to the user, if the user does not have any R language foundation, it is still difficult to retrieve the results and re-analysis on their own, and this R package is designed to help such users to visualize and explore the anaysis results. The only thing to do for such users is to configure R and Rstudio on their own computers, and then install ``SeuratExplorer``, without any other operations, an optional way is to upload the ``Seurat object`` file to a server which has been deployed with ``shinyserver`` and ``SeuratExplorer``.

> Essentially, what ``SeuratExplorer`` done is just to perform visual operations for command line tools from ``Seurat`` or other packages.

## Installation

Install the latest version from github - ***Recommended***:


``` r
if(!require(devtools)){install.packages("devtools")}
install_github("fentouxungui/SeuratExplorer", dependencies = TRUE)
```

Or install from CRAN:


``` r
# install none-CRAN dependency
if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
BiocManager::install(c("ComplexHeatmap", "MAST", "limma", "DESeq2"))
if(!require(devtools)){install.packages("devtools")}
devtools::install_github("immunogenomics/presto")

install.packages("SeuratExplorer")
```

## Run app on local
 

``` r
library(SeuratExplorer)
launchSeuratExplorer()
```

## Deploy on server

You can deploy this app on a shiny server, which allows people to view their data on a webpage by uploading the data to server.

A live demo: Upload an Rds or qs2 file, with file size no more than 20GB, to [Demo Site](http://www.nibs.ac.cn:666/SeuratExplorer/). You can download a mini demo data from [github](https://github.com/fentouxungui/SeuratExplorerServer/blob/main/inst/extdata/source-data/fly/Rds-file/G101_PC20res04.rds).


``` r
# app.R
library(SeuratExplorer)
launchSeuratExplorer()
```

## Assay option

> [Seurat Assay](https://github.com/satijalab/seurat/wiki/Assay)

> The Assay class stores single cell data.For typical scRNA-seq experiments, a Seurat object will have a single Assay ("RNA"). This assay will also store multiple 'transformations' of the data, including raw counts (@counts slot), normalized data (@data slot), and scaled data for dimensional reduction (@scale.data slot).


SeuratExplorer allows for assay switching, thereby multiple data types can be supported, including:

- scRNA-seq, usually the default 'RNA' assay.

- scATAC-seq, usually named with "ATAC", and "ACTIVITY" based on former.

- Xenium data, usually named with 'Xenium'.

- Visium HD data, usually named with 'Visium'.

- CITE-seq data,for the data of antibody-derived tags, usually named with 'ADT'.

etc.

Besides:

- SCT assay for using SCT normalization method

- cellbender assay by using cellbender output

- lsi assay from LSI weight reduction

etc.

**Slots**

- counts: Stores unnormalized data such as raw counts or TPMs

- data: Normalized data matrix

- scale.data:	Scaled data matrix

## Introduction

### Load data

- support ``Seurat`` object saved as ``.rds`` or ``.qs2`` file.

- support data processed by ``Seurat`` V5 and older versions. it may takes a while to update ``Seurat`` object when loading data.

<div class="figure">
<img src="/Users/detanw/Library/R/arm64/4.5/library/SeuratExplorer/extdata/www/upload-data.png" alt="plot of chunk unnamed-chunk-8" width="100%" />
<p class="caption">plot of chunk unnamed-chunk-8</p>
</div>

### Cell Metadata

- support download cell metadata in ``csv`` format, which can be used for further analysis.

<div class="figure">
<img src="/Users/detanw/Library/R/arm64/4.5/library/SeuratExplorer/extdata/www/cell-metadata.jpg" alt="plot of chunk unnamed-chunk-9" width="100%" />
<p class="caption">plot of chunk unnamed-chunk-9</p>
</div>

### Dimensional Reduction Plot

- support options for **Dimension Reductions**

- support options for **Cluster Resolution**

- support **split** plots

- support highlight selected clusters

- support adjust the height/width ratio of the plot

- support options for showing **cluster label**

- support adjust label size

- support adjust point size

- support download plot in pdf format, what you see is what you get

**Example plots:**

<div class="figure">
<img src="/Users/detanw/Library/R/arm64/4.5/library/SeuratExplorer/extdata/www/Dimplot-splited.png" alt="plot of chunk unnamed-chunk-10" width="80%" />
<p class="caption">plot of chunk unnamed-chunk-10</p>
</div>

### Feature Plot

- support display multiple genes simultaneous, genes names are case-insensitive. Tips: paste multiple genes from excel

- support options for **Dimension Reductions**

- support **split** plots

- support change colors for the lowest expression and highest expression

- support adjust the height/width ratio of the plot

- support adjust point size

- support download plot in pdf format, what you see is what you get

- support switch Assays which contain any one of the slots: counts, data, scale.data

**Example plots:**

<div class="figure">
<img src="/Users/detanw/Library/R/arm64/4.5/library/SeuratExplorer/extdata/www/Featureplot-splited.png" alt="plot of chunk unnamed-chunk-11" width="50%" />
<p class="caption">plot of chunk unnamed-chunk-11</p>
</div>

### Violin Plot

- support display multiple genes simultaneous, genes names are case-insensitive. Tips: paste multiple genes from excel

- support options for **Cluster Resolution**

- support **split** plots

- support **stack** and **flip** plot, and color mapping selection.

- support adjust point size and transparency

- support adjust font size on x and y axis

- support adjust the height/width ratio of the plot

- support download plot in pdf format, what you see is what you get

- support switch Assays which contain any one of the slots: counts, data, scale.data

**Example plots:**

<div class="figure">
<img src="/Users/detanw/Library/R/arm64/4.5/library/SeuratExplorer/extdata/www/ViolinPlot-splited-Stack.png" alt="plot of chunk unnamed-chunk-12" width="50%" />
<p class="caption">plot of chunk unnamed-chunk-12</p>
</div>

### Dot Plot

- support display multiple genes simultaneous, genes names are case-insensitive. Tips: paste multiple genes from excel

- support options for **Cluster Resolution** and subset clusters

- support **split** plots

- support cluster clusters

- support rotate axis and flip coordinate

- support adjust point size and transparency

- support adjust font size on x and y axis

- support adjust the height/width ratio of the plot

- support download plot in pdf format, what you see is what you get

- support switch Assays which contain slot: data

**Example plots:**

<div class="figure">
<img src="/Users/detanw/Library/R/arm64/4.5/library/SeuratExplorer/extdata/www/DotPlot-Splited.png" alt="plot of chunk unnamed-chunk-13" width="50%" />
<p class="caption">plot of chunk unnamed-chunk-13</p>
</div>

### Heatmap for cell level expression

- support display multiple genes simultaneous, genes names are case-insensitive. Tips: paste multiple genes from excel

- support options for **Cluster Resolution** and reorder clusters

- support adjust font size and rotation angle of cluster label, and flip coordinate

- support adjust the height of group bar

- support adjust the gap size between groups

- support adjust the font size of gene names

- support adjust the height/width ratio of the plot

- support download plot in pdf format, what you see is what you get

- support Assay switch

**Example plots:**

<div class="figure">
<img src="/Users/detanw/Library/R/arm64/4.5/library/SeuratExplorer/extdata/www/Heatmap-CellLevel.png" alt="plot of chunk unnamed-chunk-14" width="100%" />
<p class="caption">plot of chunk unnamed-chunk-14</p>
</div>

### Heatmap for group averaged expression

- support display multiple genes simultaneous, genes names are case-insensitive. Tips: paste multiple genes from excel

- support options for **Cluster Resolution** and reorder clusters

- support adjust font size and rotation angle of cluster label

- support adjust the font size of gene names

- support adjust the height/width ratio of the plot

- support download plot in pdf format, what you see is what you get

- support switch Assays which contain any one of the slots: data, scale.data

**Example plots:**

<div class="figure">
<img src="/Users/detanw/Library/R/arm64/4.5/library/SeuratExplorer/extdata/www/Heatmap-GroupLevel-2.png" alt="plot of chunk unnamed-chunk-15" width="50%" />
<p class="caption">plot of chunk unnamed-chunk-15</p>
</div>

### Ridge Plot

- support display multiple genes simultaneous, genes names are case-insensitive. Tips: paste multiple genes from excel

- support options for **Cluster Resolution** and reorder clusters

- support adjust column numbers

- support stack plot and color mapping

- support adjust font size on x and y axis

- support adjust the height/width ratio of the plot

- support download plot in pdf format, what you see is what you get

- support switch Assays which contain any one of the slots: counts, data, scale.data

**Example plots:**

<div class="figure">
<img src="/Users/detanw/Library/R/arm64/4.5/library/SeuratExplorer/extdata/www/RidgePlot.png" alt="plot of chunk unnamed-chunk-16" width="50%" />
<p class="caption">plot of chunk unnamed-chunk-16</p>
</div>

### Plot Cell Percentage

- support facet

- support adjust the height/width ratio of the plot

- support download plot in pdf format, what you see is what you get

**Example plots:**

<div class="figure">
<img src="/Users/detanw/Library/R/arm64/4.5/library/SeuratExplorer/extdata/www/CellRatio-Splited.png" alt="plot of chunk unnamed-chunk-17" width="50%" />
<p class="caption">plot of chunk unnamed-chunk-17</p>
</div>

### Find Cluster Markers and DEGs Analysis

This usually takes longer, please wait patiently.Please save the results before start a new analysis, the old results will be overwritten by the new results, the results can be downloaded as ``csv`` format.

#### Support two ways

- support find markers for all clusters

- support calculate DEGs for self-defined two groups, you can subset cells before calculate DEGs between two groups, default use all cells of two groups.

You can modify part calculation parameters before a analysis.

- support switch Assays which contain any one of the slots: counts, data

**Screen shots:**

<div class="figure">
<img src="/Users/detanw/Library/R/arm64/4.5/library/SeuratExplorer/extdata/www/DEGs-2.png" alt="plot of chunk unnamed-chunk-18" width="50%" />
<p class="caption">plot of chunk unnamed-chunk-18</p>
</div>

#### Output description

<div class="figure">
<img src="/Users/detanw/Library/R/arm64/4.5/library/SeuratExplorer/extdata/www/DEGs-4.jpg" alt="plot of chunk unnamed-chunk-19" width="100%" />
<p class="caption">plot of chunk unnamed-chunk-19</p>
</div>

> [FindMarkers(object, ...)](https://satijalab.org/seurat/reference/findmarkers)
>
> A data.frame with a ranked list of putative markers as rows, and associated statistics as columns (p-values, ROC score, etc., depending on the test used (test.use)). The following columns are always present:
> 
> avg_logFC: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group
> 
> pct.1: The percentage of cells where the gene is detected in the first group
> 
> pct.2: The percentage of cells where the gene is detected in the second group
> 
> p_val_adj: Adjusted p-value, based on bonferroni correction using all genes in the dataset

### Top Expressed Features

Highly expressed genes can reflect the main functions of cells, there two ways to do this. the first - ``Find Top Genes by Cell`` could find gene only high express in a few cells, while the second - ``Find Top Genes by Accumulated UMI counts`` is biased to find the highly expressed genes in most cells by accumulated UMI counts.

- support Assay switch

#### 1. Find Top Genes by Cell

#### How?

Step1: for each cell, find genes that has high UMI percentage, for example, if a cell has 10000 UMIs, and the ``UMI percentage cutoff`` is set to 0.01, then all genes that has more than 10000 * 0.01 = 100 UMIs is thought to be the highly expressed genes for this cell.

Step2: summary those genes for each cluster, firstly get all highly expressed genes in a cluster, some genes may has less cells, then for each gene, count cells in which this genes is highly expressed, and also calculate the mean and median UMI percentage in those highly expressed cells.


<div class="figure">
<img src="/Users/detanw/Library/R/arm64/4.5/library/SeuratExplorer/extdata/www/Find-Top-Genes-by-Cell.jpg" alt="plot of chunk unnamed-chunk-20" width="80%" />
<p class="caption">plot of chunk unnamed-chunk-20</p>
</div>

#### Output description

- ``celltype``: the cluster name which is define by ``Choose A Cluster Resolution``

- ``total.cells``: total cell in this cluster

- ``Gene``: this Gene is highly expressed in at least 1 cell in this cluster

- ``total.pos.cells``: how many cells express this gene

- ``total.UMI.pct``: (all UMIs of this gene)/(total UMIs of this cluster)

- ``cut.Cells``:  how many cells highly express this gene

- ``cut.pct.mean``: in those highly expressed cells, the mean expression percentage

- ``cut.pct.median``: in those highly expressed cells, the median expression percentage


#### 2. Find Top Genes by Mean UMI counts

for each cluster, calculate the ``top n`` highly expressed genes by Mean UMI counts. if a cluster has less than 3 cells, this cluster will be escaped.

- support switch Assays which contain slot: counts

<div class="figure">
<img src="/Users/detanw/Library/R/arm64/4.5/library/SeuratExplorer/extdata/www/Find-Top-Genes-by-Mean-UMI-counts.jpg" alt="plot of chunk unnamed-chunk-21" width="80%" />
<p class="caption">plot of chunk unnamed-chunk-21</p>
</div>

#### Output description

- ``CellType``: the cluster name which is define by ``Choose A Cluster Resolution``

- ``total.cells``: total cell in this cluster

- ``Gene``: the ``top n`` highly expressed genes 

- ``total.pos.cells``: how many cells express this gene

- ``MeanUMICounts``: (total accumulated UMI counts) / (total cells of this cluster)

- ``PCT``:  (total accumulated UMI counts of the gene) / (total UMIs of cluster cells)

### Feature Summary

Summary interested features by cluster, such as the positive cell percentage and mean/median expression level.

- support switch Assays which contain slot: data

<div class="figure">
<img src="/Users/detanw/Library/R/arm64/4.5/library/SeuratExplorer/extdata/www/gene-short-summary.jpg" alt="plot of chunk unnamed-chunk-22" width="80%" />
<p class="caption">plot of chunk unnamed-chunk-22</p>
</div>

#### Output description

- ``celltype``: the cluster name which is define by ``Choose A Cluster Resolution``

- ``TotalCells``: total cell in this cluster

- ``Gene``: the input genes 

- ``PCT``: the percentage of how many cells express this gene in this cluster

- ``Expr.mean``: the mean normalized expression in this cluster

- ``Expr.median``:  the median normalized expression in this cluster

### Feature Correlation Analysis

Can calculate the correlation value of gene pairs within cells from a cluster, support pearson & spearman methods.

- support switch Assays which contain slot: data

#### 3 ways to do

- ``Find Top Correlated Gene Pairs``: to find top 1000 correlated gene pairs

- ``Find Correlated Genes for A Gene``: to find the most correlated genes for input genes

- ``Calculate Correlation for A Gene List``: to calculate the correlation value for each pair of the input genes

<div class="figure">
<img src="/Users/detanw/Library/R/arm64/4.5/library/SeuratExplorer/extdata/www/featurecorrelation.jpg" alt="plot of chunk unnamed-chunk-23" width="100%" />
<p class="caption">plot of chunk unnamed-chunk-23</p>
</div>

#### Output description

<div class="figure">
<img src="/Users/detanw/Library/R/arm64/4.5/library/SeuratExplorer/extdata/www/feature-correlation-output.jpg" alt="plot of chunk unnamed-chunk-24" width="40%" />
<p class="caption">plot of chunk unnamed-chunk-24</p>
</div>

- ``GeneA``: the first gene in a Gene pair

- ``GeneB``:  the second gene in a Gene pair

- ``correlation``: the correlation value

if nothing return, this is because the input genes has very low expression level, very low expressed genes will be removed before analysis.

### Search Features

all features(genes) extracted from the row names of assay, can be used search features.

- support switch Assays which contain any one of the slots: counts, data, scale.data

### Metadata of Cells

The metadata of all cells extracted from the meta.data slot of Seurat object, ehich contains descriptive information for each cell, such as quality control metrics, cell type classifications, batch information, and experimental conditions. This metadata is crucial for organizing, filtering, integrating, and visualizing single-cell RNA-seq data.


### Structure of Seurat Object

> The Seurat object is an S4 class in R designed to store and manage single-cell expression data and associated analyses. It is a highly structured and self-contained object, allowing for the integration of various data modalities and analytical results.

> Key Slots and their Contents:

> assays: This is a list containing one or more Assay objects. Each Assay object represents a specific type of expression data.
> Each Assay object itself contains slots like counts (raw data), data (normalized data), scale.data (scaled data), and meta.features (feature-level metadata).

> meta.data:
> A data frame storing cell-level metadata. This includes information such as the number of features detected per cell (nFeature_RNA), original identity classes (orig.ident), and can be extended with additional information (e.g., cell type annotations, sample information).

> active.assay:
> A character string indicating the name of the currently active or default assay for analysis.

> active.ident:
> Stores the active cluster identity for each cell, typically resulting from clustering analyses.

> reductions:
> A list of DimReduc objects, each representing a dimensionality reduction technique applied to the data (e.g., PCA, UMAP, tSNE). These objects store the lower-dimensional embeddings of the cells.

> graphs:
> A list of Graph objects, typically storing nearest-neighbor graphs used in clustering and other analyses.
images:
> For spatial transcriptomics data, this slot stores Image objects containing spatial image data and information linking spots to their physical locations.

> project.name:
> A character string holding the name of the project.

> misc:
> A list for storing miscellaneous information not fitting into other specific slots.

## Key related packages

- [satijalab/seurat](https://github.com/satijalab/seurat): Seurat is an R toolkit for single cell genomics, developed and maintained by the Satija Lab at NYGC.

- [Hla-Lab/SeuratExplorer](https://github.com/rwcrocker/SeuratExplorer/): An interactive R shiny application for exploring scRNAseq data processed in Seurat.

- [junjunlab/scRNAtoolVis](https://github.com/junjunlab/scRNAtoolVis): Some useful function to make your scRNA-seq plot more beautiful.

- [rstudio/shiny-server](https://github.com/rstudio/shiny-server): Shiny Server is a server program that makes Shiny applications available over the web.

## Session Info


```
#> R version 4.5.1 (2025-06-13)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Tahoe 26.0
#> 
#> Matrix products: default
#> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: Europe/Zurich
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] SeuratExplorer_0.1.2 shiny_1.11.1        
#> 
#> loaded via a namespace (and not attached):
#>   [1] RcppAnnoy_0.0.22         splines_4.5.1            later_1.4.4              tibble_3.3.0             polyclip_1.10-7         
#>   [6] fastDummies_1.7.5        lifecycle_1.0.4          shinyjqui_0.4.1          rprojroot_2.1.1          globals_0.18.0          
#>  [11] processx_3.8.6           lattice_0.22-7           MASS_7.3-65              crosstalk_1.2.2          magrittr_2.0.4          
#>  [16] plotly_4.11.0            sass_0.4.10              jquerylib_0.1.4          yaml_2.3.10              remotes_2.5.0           
#>  [21] shinyBS_0.61.1           httpuv_1.6.16            Seurat_5.3.0             sctransform_0.4.2        spam_2.11-1             
#>  [26] sp_2.2-0                 sessioninfo_1.2.3        pkgbuild_1.4.8           spatstat.sparse_3.1-0    reticulate_1.43.0       
#>  [31] cowplot_1.2.0            pbapply_1.7-4            RColorBrewer_1.1-3       abind_1.4-8              pkgload_1.4.0           
#>  [36] Rtsne_0.17               purrr_1.1.0              circlize_0.4.16          ggrepel_0.9.6            irlba_2.3.5.1           
#>  [41] listenv_0.9.1            gitcreds_0.1.2           spatstat.utils_3.2-0     testthat_3.2.3           goftest_1.2-3           
#>  [46] RSpectra_0.16-2          spatstat.random_3.4-2    fitdistrplus_1.2-4       parallelly_1.45.1        commonmark_2.0.0        
#>  [51] codetools_0.2-20         DT_0.34.0                tidyselect_1.2.1         shape_1.4.6.1            farver_2.1.2            
#>  [56] shinyWidgets_0.9.0       matrixStats_1.5.0        spatstat.explore_3.5-3   jsonlite_2.0.0           ellipsis_0.3.2          
#>  [61] progressr_0.16.0         ggridges_0.5.7           ggalluvial_0.12.5        survival_3.8-3           tools_4.5.1             
#>  [66] ica_1.0-3                Rcpp_1.1.0               glue_1.8.0               gridExtra_2.3            qs2_0.1.5               
#>  [71] xfun_0.53                usethis_3.2.1            dplyr_1.1.4              shinydashboard_0.7.3     withr_3.0.2             
#>  [76] BiocManager_1.30.26      fastmap_1.2.0            litedown_0.7             callr_3.7.6              digest_0.6.37           
#>  [81] R6_2.6.1                 mime_0.13                colorspace_2.1-2         scattermore_1.2          tensor_1.5.1            
#>  [86] markdown_2.0             spatstat.data_3.1-8      tidyr_1.3.1              generics_0.1.4           renv_1.1.5              
#>  [91] data.table_1.17.8        httr_1.4.7               htmlwidgets_1.6.4        uwot_0.2.3               pkgconfig_2.0.3         
#>  [96] gtable_0.3.6             lmtest_0.9-40            S7_0.2.0                 brio_1.1.5               htmltools_0.5.8.1       
#> [101] profvis_0.4.0            dotCall64_1.2            SeuratObject_5.2.0       scales_1.4.0             png_0.1-8               
#> [106] spatstat.univar_3.1-4    knitr_1.50               rstudioapi_0.17.1        reshape2_1.4.4           nlme_3.1-168            
#> [111] curl_7.0.0               shinydashboardPlus_2.0.6 cachem_1.1.0             zoo_1.8-14               GlobalOptions_0.1.2     
#> [116] stringr_1.5.2            KernSmooth_2.23-26       parallel_4.5.1           miniUI_0.1.2             shinycssloaders_1.1.0   
#> [121] desc_1.4.3               pillar_1.11.1            grid_4.5.1               vctrs_0.6.5              RANN_2.6.2              
#> [126] urlchecker_1.0.1         promises_1.3.3           stringfish_0.17.0        xtable_1.8-4             cluster_2.1.8.1         
#> [131] evaluate_1.0.5           cli_3.6.5                compiler_4.5.1           rlang_1.1.6              future.apply_1.20.0     
#> [136] ps_1.9.1                 plyr_1.8.9               fs_1.6.6                 stringi_1.8.7            viridisLite_0.4.2       
#> [141] deldir_2.0-4             lazyeval_0.2.2           devtools_2.4.5           spatstat.geom_3.6-0      colourpicker_1.3.0      
#> [146] Matrix_1.7-3             RcppHNSW_0.6.0           patchwork_1.3.2          future_1.67.0            ggplot2_4.0.0           
#> [151] ROCR_1.0-11              fontawesome_0.5.3        igraph_2.1.4             memoise_2.0.1            RcppParallel_5.1.11-1   
#> [156] bslib_0.9.0
```

## 中文介绍

[微信公众号： 分析力工厂](https://mp.weixin.qq.com/s/lpvI9OnyN95amOeVGmeyMQ)
