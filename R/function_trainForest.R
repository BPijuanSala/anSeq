########################################################################################
## Title: trainForest.R
## Author: Blanca Pijuan Sala
## Date: 26 February 2017
## **bioSeq package**
########################################################################################

#' @title Train randomForest on a dataset.
#' @description Cluster cells using the spearman rank correlation
#' as a similarity metric and (1-correlation)/2 as a distance metric.
#' @param data Dataframe or matrix containing the counts (genes x cells).
#' It is recommended to use only the HVG for this.
#' @param groups Vector specifying which group each cell belongs to.
#' Cells should be as names and group ID as vector content.
#' @param numTrees Number of trees to grow in randomForest. Default: 1000.
#' Cells should be as names and group ID as vector content.
#' @param rankData Logical. Rank data before training based on counts? 
#' Default: TRUE
#' @return list with randomForest object (RForest) and a vector of the final
#' genes used for training.
#' @seealso \code{\link[randomForest]{randomForest}},
#' \code{\link[anSeq]{allocateCellsForest}}
#' @references \code{\link[randomForest]{randomForest}}.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname trainForest
#' @importFrom randomForest randomForest importance

trainForest <- function(data,groups,numTrees=1000,rankData=T){
  #Grow forest #######
  dataOriginal <- data
  cells<-colnames(data)
  if (rankData==T){
    cat("Ranking data...\n")
    data.ranks<-apply(data, 2, function(x) rank(x, ties.method = "average"))

  } else {
    data.ranks <- dataOriginal
  }

  clust<-groups[cells]
  names(clust)<-cells
  cat("Training preliminary forest...\n")

  set.seed(0)
  RForest.test=randomForest::randomForest(t(data.ranks), factor(clust), ntree = numTrees)

  #extract top 25% genes and re-run the random forest

  imp<-randomForest::importance(RForest.test)
  g<-row.names(imp)[order(imp[,1], decreasing=T)][1:(ceiling(length(imp)/4))]#Predict data


  most.import.high.var<-g

  data=dataOriginal[g,cells]
  if (rankData==T){
    cat("Ranking data...\n")

    data.ranks<-apply(data, 2, function(x) rank(x, ties.method = "average"))

  } else {
    data.ranks=data
  }

  clust<-groups[cells]
  names(clust)<-cells

  cat("Growing forest...\n")

  set.seed(0)#1
  RForest=randomForest::randomForest(t(data.ranks), factor(clust), ntree=numTrees)

  return(list(
    RForest = RForest,
    genesTrained = most.import.high.var
  ))


}
