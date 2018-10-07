########################################################################################
## Title: allocateCellsForest.R
## Author: Blanca Pijuan-Sala
## Description: Use randomForest classifier to allocate cells.
## Date: 26 February 2017
## **anSeqSeq package**
########################################################################################

#' @title Allocate cells based on trained randomForest.
#' @description For each cell, it allocates it to its most likely cluster
#' it belongs to from the trained random Forest.
#' @param newData Dataframe or matrix containing the counts (genes x cells).
#' These should be from the cells you want to predict the cluster of.
#' @param RForestObj randomForest object already trained. To train a dataset, you
#' can use \code{\link[anSeq]{trainForest}}.
#' @param genesTrained vector of the genes that were used to train
#' the randomForest.
#' @return dataframe with predicted allocations and probabilities.
#' @seealso \code{\link[randomForest]{randomForest}},
#' \code{\link[anSeq]{trainForest}}.
#' @references \code{\link[randomForest]{randomForest}}.
#' @author Blanca Pijuan Sala.
#' @export
#' @rdname allocateCellsForest
#' @importFrom stats predict
allocateCellsForest <- function(newData,RForestObj, genesTrained){
  
  cells.new<-colnames(newData)
  most.import.high.var <- genesTrained
  RForest <- RForestObj
  
  test=newData[most.import.high.var,cells.new]
  test<-apply(test, 2, function(x) rank(x, ties.method = "average"))
  
  set.seed(0)
  prob.RF<-apply(test, 2, function(x) max(predict(RForest, x, type="prob") ))
  alloc.RF<-apply(test, 2, function(x) (predict(RForest, x) ))
  
  
  alloc.prob.RF<-data.frame(prediction=alloc.RF, probability=prob.RF)
  return (alloc.prob.RF)
  
  
}