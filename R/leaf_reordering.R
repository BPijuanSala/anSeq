########################################################################################
## Title: leaf_reordering.R
## Author: Blanca Pijuan-Sala
## Description: It reorders the leaves of the hierarchical clustering so that they 
##              depending on their similarities
## Date: 15/11/2016
## **anSeq package**
########################################################################################

##----------------------------------------
## findHVG
##----------------------------------------

#' @title Reorder leaves of hierarchical clustering
#' @description It reorders the leaves of the hierarchical clustering so that they 
#'              depending on their similarities. 
#' @param data Matrix of counts with genes as rows and cells as columns
#' @param clust.method Clustering method. Default: ward.D2.
#' @param cor.method method of correlation to compute clustering. Default: spearman.
#' @param reorder What would you like to reorder? Options: "column", "row", "both". Default: "column".
#' @return  Resulting list from hierarchical clustering.
#' @author Blanca Pijuan-Sala.
#' @export
#' @rdname leaf_reordering
#' @import cba
#' @export
leaf_reordering <- function(data, clust.method = "ward.D2",reorder="column", cor.method="spearman"){
  #reorder can be column or row.
  
  if (reorder == "column"){
    
    mat1 <- as.matrix(data)
    
    cor.mat1<-cor(mat1,method=cor.method)
    dissim1<-(1-cor.mat1)/2
    d1<-as.dist(dissim1)
    
    #d <- dist(mat1)
    
    hc <- hclust(d1, method = clust.method, members=NULL)
    
    co <- order.optimal(d1, hc$merge)
    
    ho <- hc
    
    ho$merge <- co$merge
    
    ho$order <- co$order
    
    ho$labels
    return(ho)
    
  }
  if (reorder=="row"){
    
    mat2 <- t(mat1)
    
    
    cor.mat2<-cor(mat2,method=cor.method)
    dissim2<-(1-cor.mat2)/2
    d2<-as.dist(dissim2)
    
    #d2 <- dist(mat2)
    
    hc2 <- hclust(d2, method = clust.method, members=NULL)
    
    co2 <- order.optimal(d2, hc2$merge)
    
    ho2 <- hc2
    
    ho2$merge <- co2$merge
    
    ho2$order <- co2$order
    return(ho2)
  }
  if (reorder == "both"){
    
    mat1 <- as.matrix(data)
    
    cor.mat1<-cor(mat1,method=cor.method)
    dissim1<-(1-cor.mat1)/2
    d1<-as.dist(dissim1)
    
    #d <- dist(mat1)
    
    hc <- hclust(d1, method = clust.method, members=NULL)
    
    co <- order.optimal(d1, hc$merge)
    
    ho <- hc
    
    ho$merge <- co$merge
    
    ho$order <- co$order
    
    ho$labels
    
    
    mat2 <- t(mat1)
    
    
    cor.mat2<-cor(mat2,method=cor.method)
    dissim2<-(1-cor.mat2)/2
    d2<-as.dist(dissim2)
    
    #d2 <- dist(mat2)
    
    hc2 <- hclust(d2, method = clust.method, members=NULL)
    
    co2 <- order.optimal(d2, hc2$merge)
    
    ho2 <- hc2
    
    ho2$merge <- co2$merge
    
    ho2$order <- co2$order
    
    return (list(cell.clust = ho,
                 gene.clust = ho2))
    
  }
}