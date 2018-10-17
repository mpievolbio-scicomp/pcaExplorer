# About pcaExplorer

`pcaExplorer` is a Bioconductor package containing a Shiny application for
analyzing expression data in different conditions and experimental factors. 

`pcaExplorer` guides the user in exploring the Principal Components of the data, 
providing tools and functionality to detect outlier samples, genes that show 
particular patterns, and additionally provides a functional interpretation of 
the principal components for further quality assessment and hypothesis generation
on the input data. 

Thanks to its interactive/reactive design, it is designed to become a practical
companion to any RNA-seq dataset analysis, making exploratory data analysis 
accessible also to the bench biologist, while providing additional insight also
for the experienced data analyst.

Moreover, `pcaExplorer` supports reproducible research with state saving and automated 
report generation. 

`pcaExplorer` was developed in the Bioinformatics Division led by Harald Binder 
at the IMBEI (Institut für Medizinische Biometrie, Epidemiologie und Informatik) 
in the University Medical Center of the Johannes Gutenberg University Mainz.

## Developers

<a href="mailto:mailto:marinif@uni-mainz.de?subject=[pcaExplorer_feedback]" class="btn btn-primary">Federico Marini</a>

## Code

All code for `pcaExplorer` is available on 
<a href="https://github.com/federicomarini/pcaExplorer" target="_blank">GitHub</a>.

# Citation info

If you use `pcaExplorer` for your analysis, please cite it as here below:

```r
citation("pcaExplorer")
```

```
To cite package ‘pcaExplorer’ in publications use:

  Federico Marini (2018). pcaExplorer: Interactive Visualization of RNA-seq Data Using
  a Principal Components Approach. R package version 2.6.0.
  https://github.com/federicomarini/pcaExplorer

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {pcaExplorer: Interactive Visualization of RNA-seq Data Using a Principal Components Approach},
    author = {Federico Marini},
    year = {2018},
    note = {R package version 2.6.0},
    url = {https://github.com/federicomarini/pcaExplorer},
  }
```
