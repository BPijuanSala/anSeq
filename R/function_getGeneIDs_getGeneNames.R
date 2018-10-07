########################################################################################
## Title: function_getGeneIDs_getGeneNames.R
## Author: Blanca Pijuan-Sala
## Description: Functions to obtain gene IDs or gene Names.
## Date: 25 December 2016
## **anSeq package**
########################################################################################


##----------------------------------------
## getGeneID
##----------------------------------------

#' @title Obtain geneIDs from geneNames
#' @description It takes a vector of geneNames and looks for the geneIDs. 
#' @param geneNames vector.
#' @param Table By default, it uses mm10 geneID - geneName association from Ensembl.
#' The table should have the columns named as 'Gene.ID' and 'Associated.Gene.Name'
#' @return  List with table of associations and vector of geneIDs.
#' @author Blanca Pijuan-Sala.
#' @rdname getGeneID
#' @export

getGeneID <- function(geneNames,Table=geneTable){
  geneTab.select <- Table[Table[,"Associated.Gene.Name"] %in% geneNames,]
  
  geneTab.select.order <- geneTab.select[order(match(geneTab.select[,"Associated.Gene.Name"], geneNames)),]
  geneIDs <- geneTab.select.order[,"Gene.ID"]
  
  
  return(list(
    table = geneTab.select.order,
    geneIDs = as.character(geneIDs)
  ))
}



##----------------------------------------
## getGeneName
##----------------------------------------

#' @title Obtain geneNames from geneIDs
#' @description It takes a vector of geneIDs and looks for the geneNames
#' @param geneIDs vector.
#' @param Table By default, it uses mm10 geneID - geneName association from Ensembl.
#' The table should have the columns named as 'Gene.ID' and 'Associated.Gene.Name'
#' @return  List with table of associations and vector of geneNames.
#' @author Blanca Pijuan-Sala.
#' @rdname getGeneName
#' @export

getGeneName <- function(geneIDs,Table=geneTable){
  row.names(Table) <- Table[,"Gene.ID"]
  geneTab.select <- Table[geneIDs,]
  geneNames <- as.character(geneTab.select[,"Associated.Gene.Name"])
  
  return(list(
    table = geneTab.select,
    geneNames = as.character(geneNames)
  ))
}

