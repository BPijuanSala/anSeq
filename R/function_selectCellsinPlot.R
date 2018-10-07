########################################################################################
## Title: function_selectCellsinPlot.R
## Author: Blanca Pijuan Sala
## Description: Select cells within a plot.
## Date: 26 February 2017
## **anSeq package**
########################################################################################


#' @title Select cells from a plot
#' @description It allows you to select cells from a plot
#'  and it returns a vector of cell names or indices. 
#'  Select your points by left clicking, and then right clicking to close the 
#'  polygon.
#' @param data matrix or dataframe with x and y in columns 1 and 2, respectively.
#' Ideally, the rownames should correspond to the cell names.
#' @param colors row vector with the colors for the data. The indices between
#' the data and the colors should match. Default: All points will be black.
#' @param pch Type of point in plot. Default: 20.
#' @param ... Other additional features for plot.
#' @author Blanca Pijuan Sala.
#' @return Vector of cells.
#' @seealso  \code{\link[gatepoints]{fhs}}.
#' @export
#' @rdname selectCellsinPlot
#' @import grDevices
#' @importFrom gatepoints fhs

selectCellsinPlot <- function(data,colors=NULL,pch=20,...){
  if (is.null(colors)){
    colors <- rep("black",nrow(data))
    
  }
  dim = ncol(data)
  
  if (dim==2){
    
    X11()
    
    plot(data[,1],data[,2],col=colors,pch=pch,...)
    
    selectedPoints <- gatepoints::fhs(data)
    
  }  else {
    cat("Error: The number of dimensions should be 2.")
  }
  
}