########################################################################################
## Title: function_findSurfaceMarkers.R
## Author: Blanca Pijuan-Sala
## Description: R function to find genes that have been described as surface proteins.
## Date: 6 June 2018
## **anSeq package**
########################################################################################


#' @title Find markers within groups
#' @description It finds markers of particular groups within the dataset by comparing them to 
#' the other groups. 
#' @param data vector of gene Names.
#' @param Reference Vector with names being the ENTREZ gene symbols and elements being the 
#' confidence of being a surface protein. Default: MouseSurface. You can set it to human by doing
#' HumanSurface
#' @author Blanca Pijuan-Sala.
#' @export
#' @rdname findSurfaceMarkers

findSurfaceMarkers<-function(data,Reference=MouseSurface){
  Ref <- Reference
  Reference <- as.character(Reference)
  names(Reference)<-names(Ref)
  sm <- Reference[names(Reference)%in%data]
  return(sm)
}
