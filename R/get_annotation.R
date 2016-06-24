

#' Title
#'
#' @param dds
#' @param biomart_dataset
#' @param idtype
#'
#' @return
#' @export
#'
#' @examples
get_annotation <- function(dds, biomart_dataset, idtype){
  if(is.null(biomart_dataset))
    stop("Select a species to generate the corresponding annotation.
To obtain a list, type mart = useMart('ensembl'), followed by listDatasets(mart).")

  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  host="www.ensembl.org",
                  dataset=biomart_dataset)
  anns <- getBM(attributes = c(idtype, "external_gene_name", "description"),
                filters = idtype,
                values = rownames(dds),
                mart = mart)

  # keep and match with the ones that are actually there
  anns2 <- anns[match(rownames(dds), anns[, 1]), ]
  rownames(anns2) <- rownames(dds)
  # rename the columns rsp. add row names to be consistent with other function
  colnames(anns2) <- c("gene_id","gene_name","description")

  return(anns2)
}





#' Title
#'
#' @param dds
#' @param orgdb_species
#' @param idtype
#'
#' @return
#' @export
#'
#' @examples
get_annotation_orgdb <- function(dds, orgdb_species, idtype){
  if(is.null(orgdb_species))
    stop("Select a species to generate the corresponding annotation")

  orgdbpkgs <- data.frame(
    pkg = c("org.Ag.eg.db","org.At.tair.db","org.Bt.eg.db","org.Ce.eg.db","org.Cf.eg.db","org.Dm.eg.db","org.Dr.eg.db","org.EcK12.eg.db",
            "org.EcSakai.eg.db","org.Gg.eg.db","org.Hs.eg.db","org.Hs.ipi.db","org.Mm.eg.db","org.Mmu.eg.db","org.Pf.plasmo.db",
            "org.Pt.eg.db","org.Rn.eg.db","org.Sc.sgd.db","org.Sco.eg.db","org.Ss.eg.db","org.Tgondii.eg.db","org.Xl.eg.db"),
    descr = c("Anopheles","Arabidopsis","Bovine","Worm","Canine","Fly","Zebrafish","E coli strain K12","E coli strain Sakai","Chicken",
              "Human","org.Hs.ipi.db","Mouse","Rhesus","Malaria","Chimp","Rat","Yeast","Streptomyces coelicolor","Pig","Toxoplasma gondii",
              "Xenopus"),
    stringsAsFactors = FALSE
  )

  if(!(orgdb_species %in% orgdbpkgs$pkg)){
    message("The orgDB package is most likely not existent in Bioconductor")
    message("It should be one of",orgdbpkgs$pkg)
  }

  if(!require(orgdb_species,character.only=TRUE))
    stop("The package ",orgdb_species, " is not installed/available. Try installing it with biocLite('",orgdb_species,"')")


  if(!(idtype %in% keytypes(eval(parse(text=orgdb_species))))){
    stop("The key you provided is not listed as key for the annotation package. Please try one of ",
         paste(keytypes(eval(parse(text=orgdb_species))),collapse=","))
  }

  pkg <- eval(parse(text=orgdb_species))

  if(idtype=="SYMBOL")
    warning("You probably do not need to convert symbol to symbol") # the performance would somehow be affected

  anns_vec <- mapIds(pkg, keys=rownames(dds), column=c("SYMBOL"),
                     keytype=idtype)

  anns <- data.frame(
    gene_id = rownames(dds),
    gene_name = anns_vec,
    stringsAsFactors = FALSE,
    row.names = rownames(dds)
  )

  return(anns)
}

