# pcaExplorer 2.8.0

## New features

* Added a `NEWS.md` file to track changes to the package
* PCA plots now are correctly generated with fixed coordinates

# pcaExplorer 2.6.0

## New features

* Automatically computing size factors where required
* Added progress indication when compiling the report

## Bug fixes

* Fixed after changes in threejs package
* Edited dropdown menu to remove unused green badge
* Menus start expanded on the side, again
* `theme_bw` applied when needed, corrected previous behavior

## Other notes

* Updated citation infos
* Slight difference in handling validate/need errors

# pcaExplorer 2.2.0

## New features

* Added Demo data, loadable via demo button

## Bug fixes

* Plots work now without cutting out points when zooming in

## Other notes

* Saved reactive values are now exported to dedicate environments (instead of assigning to global)

# pcaExplorer 1.99.0

## Other notes

* Reflecting the major feature added, will trigger a major version number bump. Welcome soon, pcaExplorer 2.0.0!

# pcaExplorer 1.1.5

## New features

* Automated report generation - template available + editor in the app tab for advance user customization
* Support for state saving, in the global environment as well as with binary data
* All plots generated can be now exported with the dedicated button 
* Added confidence ellipse for PCA plot
* Added 3d pca plot
* Added functions to automatically retrieve the annotation in format ready to use for the app
* Added profile explorer function, for plotting behaviour across all samples for subset of genes
* Added distribution plots
* Added pairwise correlation plot
* Added table to enhance readability of the gene finder plot, also by annotating sample names

## Bug fixes

* Minor typos fixed in the tabs
* Added option row.names to read.delim for allowing row names when uploading the data

## Other notes

* Added extra info in the about section
* Instructions and vignette rewritten to reflect new design of the app

# pcaExplorer 1.1.3

## Bug fixes

* Remove y axis limits to gene boxplots
* Fixed: correct labels and colors assignements for genespca

# pcaExplorer 1.0.0

## Other notes

* Released in Bioconductor 3.3

# pcaExplorer 0.99.1

## Other notes

* Changed format of the NEWS file

# pcaExplorer 0.99.0

## Other notes

* Ready for submission to Bioconductor

# pcaExplorer 0.9.0

## Other notes

* Added TravisCI integration for both branches
* Added appveyor integration - plus badges in the README.md
* Code cleanup
* Added screenshots for the vignette
* Removed some lengthy tests

# pcaExplorer 0.8.0

## New features

* Selection of identifier type available in pca2go

## Bug fixes

* Couple of layout fixes

## Other notes

* MIT license
* Added TravisCI integration
* Added codecov integration
* Enhanced documentation

# pcaExplorer 0.7.0

## New features

* Vignette full draft done

# pcaExplorer 0.6.4

## Other notes

* Updated NEWS file

# pcaExplorer 0.6.3

## New features

* About and Instructions done by now
* Added some missing details on the documentations

# pcaExplorer 0.6.2

## Other notes

* Corrected wordings for (cor)relations of principal components with covariates
* Added a couple of checks if correct objects are provided

# pcaExplorer 0.6.1

## New features

* Added function to remove selected samples suspected to be deemed as outliers, in order to see the effect of clustering on the good ones

# pcaExplorer 0.6.0

## Other notes

* Documentation completed
* Examples fully working, cleaned up further a little more.

# pcaExplorer 0.5.0

## Other notes

* Further steps in direction of R CMD check

# pcaExplorer 0.4.0

## New features

* Added pca2go live functionality

# pcaExplorer 0.3.0

## New features

* Added color palette to choose, and dependent on the samples and factors available/selected

# pcaExplorer 0.2.0

## New features

* Multifactorial exploration completed and adaptable to each dataset

# pcaExplorer 0.1.0

## New features

* Restyling and (re)packaging mostly completed
