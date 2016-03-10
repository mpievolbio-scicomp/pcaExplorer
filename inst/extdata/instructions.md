This information is also contained in the `pcaExplorer` package vignette.

## Introduction

`pcaExplorer` is a Bioconductor package containing a Shiny application for
analyzing expression data in different conditions and experimental factors. 

It is a general-purpose interactive companion tool for RNA-seq analysis, which 
guides the user in exploring the Principal Components of the data under inspection.

`pcaExplorer` provides tools and functionality to detect outlier samples, genes
that show particular patterns, and additionally provides a functional interpretation of 
the principal components for further quality assessment and hypothesis generation
on the input data. 

Moreover, a novel visualization approach is presented to simultaneously assess 
the effect of more than one experimental factor on the expression levels.

Thanks to its interactive/reactive design, it is designed to become a practical
companion to any RNA-seq dataset analysis, making exploratory data analysis 
accessible also to the bench biologist, while providing additional insight also
for the experienced data analyst.


## Launching the application

The `pcaExplorer` app can be launched in different modes:

- `pcaExplorer(dds = dds, rld = rld)`, where `dds` is a `DESeqDataSet` object and `rld` is a `DESeqTransform`
object, which were created during an existing session for the analysis of an RNA-seq
dataset with the `DESeq2` package

- `pcaExplorer(dds = dds)`, where `dds` is a `DESeqDataSet` object. The `rld` object is automatically 
computed upon launch.

- `pcaExplorer(countmatrix = countmatrix, coldata = coldata)`, where `countmatrix` is a count matrix, generated
after assigning reads to features such as genes via tools such as `HTSeq-count` or `featureCounts`, and `coldata`
is a data frame containing the experimental covariates of the experiments, such as condition, tissue, cell line,
run batch and so on.

- `pcaExplorer()`, and then subsequently uploading the count matrix and the covariates data frame through the 
user interface. These files need to be formatted as tab separated files, which is a common format for storing
such count values.

Additional parameters and objects that can be provided to the main `pcaExplorer` function are:

- `pca2go`, which is an object created by the `pca2go` function, which scans the genes with high loadings in 
each principal component and each direction, and looks for functions (such as GO Biological Processes) that 
are enriched above the background. The offline `pca2go` function is based on the routines and algorithms of 
the `topGO` package, but as an alternative, this object can be computed live during the execution of the app
exploiting the `goana` function, provided by the `limma` package. Although this likely provides more general
(and probably less informative) functions, it is a good compromise for obtaining a further data interpretation.

- `annotation`, a data frame object, with `row.names` as gene identifiers (e.g. ENSEMBL ids) identical to the 
row names of the count matrix or `dds` object, and an extra column `gene_name`, containing e.g. HGNC-based 
gene symbols. This can be used for making information extraction easier, as ENSEMBL ids (a usual choice when
assigning reads to features) do not provide an immediate readout for which gene they refer to. This can be
either passed as a parameter when launching the app, or also uploaded as a tab separated text file.



## The controls sidebar

Most of the input controls are located in the sidebar, some are as well in the individual tabs of the app.
By changing one or more of the input parameters, the user can get a fine control on what is displayed.

### Data upload

These file input controls are available when no `dds` or `countmatrix` + `coldata` are provided. Additionally,
it is possible to upload the `annotation` data frame.

### App settings

Here are the parameters that set input values for most of the tabs. By hovering over with the mouse,
the user can receive additional information on how to set the parameter, powered by the `shinyBS` package.

- x-axis PC - Select the principal component to display on the x axis
- y-axis PC - Select the principal component to display on the y axis
- Group/color by - Select the group of samples to stratify the analysis. Can also assume multiple values.
- Nr of (most variant) genes - Number of genes to select for computing the principal components. The top n genes are
selected ranked by their variance inter-samples
- Alpha - Color transparency for the plots. Can assume values from 0 (transparent) to 1 (opaque)
- Labels size - Size of the labels for the samples in the principal components plots
- Points size - Size of the points to be plotted in the principal components plots
- Variable name size - Size of the labels for the genes PCA - correspond to the samples names
- Scaling factor - Scale value for resizing the arrow corresponding to the variables in the PCA for the genes. It
should be used for mere visualization purposes
- Color palette - Select the color palette to be used in the principal components plots. The number of colors 
is selected automatically according to the number of samples and to the levels of the factors of interest
and their interactions

### Plot export settings        

Width and height for the figures to export are input here in cm.

Additional controls available in the single tabs are also assisted by tooltips that show on hovering the mouse.
Normally they are tightly related to the plot/output they are placed nearby.


           
## The app panels

The `pcaExplorer` app is structured in different panels, each focused on a different aspect of the 
data exploration. 

Most of the panels work extensively with click-based and brush-based interactions, to gain additional
depth in the explorations, for example by zooming, subsetting, selecting. This is possible thanks to the 
recent developments in the `shiny` package/framework.

The available panels are the described in the following subsections.

### About

Contains general information on `pcaExplorer`.

### Instructions

This is where you most likely are reading this text (otherwise in the package vignette).

### Samples View

This panel displays the PCA projections of sample expression profiles onto any pair of components,
a scree plot, a zoomed PCA plot, a plot of the genes with top and bottom loadings, and a PCA plot where
it is possible to remove samples deemed to be outliers in the analysis.

### Genes View

This panel displays the PCA projections of genes abundances onto any pair of components, with samples
as biplot variables, to identify interesting groups of genes. Zooming is also possible, and clicking on single
genes, a boxplot is returned, grouped by the factors of interest. A static and an interactive heatmap are 
provided, including the subset of selected genes. These are also reported in `datatable` objects.

### GeneFinder

The user can search and display the expression values of a gene of interest, either by ID or gene
name, as provided in the `annotation`. A handy panel for quick screening of shortlisted genes, again grouped by
the factors of interest.
### PCA2GO

This panel shows the functional annotation of the principal components, with GO functions enriched in the 
genes with high loadings on the selected principal components. It allows for the live computing of the object,
that can otherwise provided as a parameter when launching the app. The panel displays a PCA plot for the 
samples, surrounded on each side by the tables with the functions enriched in each component and direction.

### Multifactor Exploration

This panel allows for the multifactor exploration of datasets with 2 or more experimental factors. The user has to select 
first the two factors and the levels for each. Then, it is possible to combine samples from Factor1-Level1 in the selected
order by clicking on each sample name. The matrix to build requires an equal number of samples for each level of Factor 1.
A typical case for choosing factors 1 and 2 is for example when different conditions and tissues are present

Once constructed, a plot is returned that tries to represent simultaneously the effect of the two factors on the data.
Each gene is represented by a dot-line-dot structure, with the color that is indicating the tissue (factor 2) where the gene 
is mostly expressed. Each gene has two dots, one for each condition level (factor 1), and the position of the points is dictated
by the scores of the principal components calculated on the matrix object. The line connecting the dots is darker when the 
tissue where the gene is mostly expressed varies throughout the conditions. 

This representation is under active development, and it is promising for identifying interesting sets or clusters of genes
according to their behavior on the Principal Components subspaces. Zooming and exporting of the underlying genes is also
allowed by brushing on the main plot.


## `pcaExplorer` on published datasets

We can run `pcaExplorer` for demonstration purpose on published datasets that are available as SummarizedExperiment
in an experiment Bioconductor packages.

We will use the `airway` dataset, which can be installed with this command

```
source("https://bioconductor.org/biocLite.R")
biocLite("airway")
```

This package provides a RangedSummarizedExperiment object of read counts in genes for an RNA-Seq experiment 
on four human airway smooth muscle cell lines treated with dexamethasone. More details such as gene models and 
count quantifications can be found in the `airway` package vignette. 

To run `pcaExplorer` on this dataset, the following commands are required

```
library(airway)

data(airway)

dds_airway <- DESeqDataSet(airway,design=~dex+cell)
dds_airway
rld_airway <- rlogTransformation(dds_airway)
rld_airway
pcaExplorer(dds = dds_airway,
            rlt = rld_airway)
```
The `annotation` for this dataset can be built by exploiting the `org.Hs.eg.db` package

```
library(org.Hs.eg.db)
genenames_airway <- mapIds(org.Hs.eg.db,keys = rownames(dds_airway),column = "SYMBOL",keytype="ENSEMBL")
annotation_airway <- data.frame(gene_name = genenames_airway,
                                rownames = rownames(dds_airway),
                                stringsAsFactors = FALSE)
head(annotation_airway)                                
# launch the app now with the annotation
pcaExplorer(dds = dds_airway,
            rlt = rld_airway,
            annotation = annotation_airway)

```

If desired, alternatives can be used. See the well written annotation workflow available at the Bioconductor site (https://bioconductor.org/help/workflows/annotation/annotation/).

## `pcaExplorer` on synthetic datasets

For testing and demonstration purposes, a function is also available to generate synthetic datasets whose counts
are generated based on two or more experimental factors.

This can be called with the command

```
dds_multifac <- makeExampleDESeqDataSet_multifac(betaSD_condition = 3,betaSD_tissue = 1)
```

See all the available parameters by typing `?makeExampleDESeqDataSet_multifac`. Credits are given to the initial
implementation by Mike Love in the `DESeq2` package.

The following steps run the app with the synthetic dataset

```
dds_multifac <- makeExampleDESeqDataSet_multifac(betaSD_condition = 1,betaSD_tissue = 3)
dds_multifac
rld_multifac <- rlogTransformation(dds_multifac)
rld_multifac
## checking how the samples cluster on the PCA plot
pcaplot(rld_multifac,intgroup = c("condition","tissue"))
## launching the app
pcaExplorer(dds = dds_multifac,
            rlt = rld_multifac)
```

When such a dataset is provided, 


## Functions exported by the package for standalone usage

export(correlatePCs)
export(genespca)
export(hi_loadings)
export(limmaquickpca2go)
export(makeExampleDESeqDataSet_multifac)
export(pca2go)
export(pcaExplorer)
export(pcaplot)
export(pcascree)
export(plotPCcorrs)
export(topGOtable)




