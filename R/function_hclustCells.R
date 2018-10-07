########################################################################################
## Title: hclustCells.R
## Author: Blanca Pijuan Sala
## Description: Cluster cells using leaf reordering.
## Date: 26 February 2017
## **anSeq package**
########################################################################################

#' @title Cluster cells
#' @description Cluster cells using the spearman rank correlation
#' as a similarity metric and (1-correlation)/2 as a distance metric.
#' @param data Dataframe or matrix containing the counts (genes x cells). The 
#' algorithm will cluster on these.
#' @param clust.method Clustering method to use with hclust. Default: WarD2.
#' For more options, see \code{\link[stats]{hclust}}.
#' @param leafreorder Specify whether you want to obtain an hclust object
#' with reordered gene ("row") or cell leaves ("column"). If you want
#' to obtain both as a list, specify "both". Default: column (it will 
#' reorder the cells). 
#' @return hclust object
#' @seealso \code{\link[stats]{hclust}}.
#' @references \code{\link[stats]{hclust}}, \code{\link[cba]{order.optimal}}
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname hclustCells
#' @importFrom cba order.optimal 
#' @importFrom stats hclust reorder cor as.dist

hclustCells <- function(data, clust.method = "ward.D2",leafreorder="column"){
  #reorder can be column or row.
  
  if (reorder == "column"){
    
    mat1 <- as.matrix(data)
    
    cor.mat1<-cor(mat1,method="spearman")
    dissim1<-(1-cor.mat1)/2
    d1<-as.dist(dissim1)
    
    #d <- dist(mat1)
    
    hc <- stats::hclust(d1, method = clust.method, members=NULL)
    
    co <- cba::order.optimal(d1, hc$merge)
    
    ho <- hc
    
    ho$merge <- co$merge
    
    ho$order <- co$order
    
    ho$labels
    return(ho)
    
  }
  if (reorder=="row"){
    
    mat2 <- t(mat1)
    
    
    cor.mat2<-cor(mat2,method="spearman")
    dissim2<-(1-cor.mat2)/2
    d2<-as.dist(dissim2)
    
    #d2 <- dist(mat2)
    
    hc2 <- stats::hclust(d2, method = clust.method, members=NULL)
    
    co2 <- cba::order.optimal(d2, hc2$merge)
    
    ho2 <- hc2
    
    ho2$merge <- co2$merge
    
    ho2$order <- co2$order
    return(ho2)
  }
  if (reorder == "both"){
    
    mat1 <- as.matrix(data)
    
    cor.mat1<-cor(mat1,method="spearman")
    dissim1<-(1-cor.mat1)/2
    d1<-as.dist(dissim1)
    
    #d <- dist(mat1)
    
    hc <- stats::hclust(d1, method = clust.method, members=NULL)
    
    co <- cba::order.optimal(d1, hc$merge)
    
    ho <- hc
    
    ho$merge <- co$merge
    
    ho$order <- co$order
    
    ho$labels
    
    
    mat2 <- t(mat1)
    
    
    cor.mat2<-cor(mat2,method="spearman")
    dissim2<-(1-cor.mat2)/2
    d2<-as.dist(dissim2)
    
    #d2 <- dist(mat2)
    
    hc2 <- stats::hclust(d2, method = clust.method, members=NULL)
    
    co2 <- cba::order.optimal(d2, hc2$merge)
    
    ho2 <- hc2
    
    ho2$merge <- co2$merge
    
    ho2$order <- co2$order
    
    return (list(cell.clust = ho,
                 gene.clust = ho2))
    
  }
}