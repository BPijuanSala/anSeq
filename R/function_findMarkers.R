########################################################################################
## Title: function_findMarkers.R
## Author: Blanca Pijuan-Sala
## Description: R function to find genes expressed in a particular cluster compared to the others.
## Date: 4 January 2017 - updated 06/06/2018
## **anSeq package**
########################################################################################


#' @title Find markers within groups
#' @description It finds markers of particular groups within the dataset by comparing them to 
#' the other groups. 
#' @param counts Dataframe or matrix containing the counts (genes x cells).
#' @param genes vector of genes to use for analysis. Default: rownames(counts). It is important that
#' this list has been processed (i.e. not containing genes which all counts are 0). You can set this 
#' to the highly variable genes of the dataset. You can compute them using \code{findHVG}.
#' @param groups Vector specifying what group each cell belongs to. e.g. 
#' c("group1","group2","group1"...). This vector should be labelled (names()) by the cell name.
#' @param topGeneThres Threshold of upregulated genes to consider. 
#' Only used for 'strict' option. Default: 50.
#' @param logFCThres Threshold of log Fold Change to consider. 
#' Only used for 'strict' option. Default = 0.
#' @param groupTarget Set this to one label to perform DE between that group
#' i.e. group target, vs the rest of the cells.
#' @param qvalue FDR q-value to consider for threshold when testing. Default: 0.05.
#' @param mode "strict" or "general". "strict" will look for upregulated
#' genes by (1) doing pairwise comparisons between each of the groups and then (2) only
#' giving those genes that are upregulated in all the comparisons from group A
#' to all others e.g. if gene 1 is upregulated only when group A is compared to B
#' but not when A is compared to C, then it will not consider it as a marker gene. Instead, 
#' if the total number of groups is A, B, C (3), and if gene 1 is upregulated in A when
#'  comparing  A to B and A to C, then you will retrieve it."general" will make pairwise
#'  comparisons of group A to all the others as a group and output both upregulated and
#'  downregulated genes in a list of matrices. Default: "general".

#' @return  list of markers for each group.

#' @author Blanca Pijuan Sala.
#' @export
#' @rdname findMarkers
#' @import edgeR
#' @importFrom utils combn
#' @importFrom stats p.adjust
findMarkers <- function(counts,genes=NULL,groups,
                        topGeneThres=1000, 
                        logFCThres = 0, 
                        qvalue=0.05, mode="general",groupTarget=NULL){
  
  if (is.null(genes)){
    genes <- rownames(counts)
  }
  #Define groups
  uniq.groups <- unique(groups)
  
  #Take those genes to take into account from the counts table.
  counts <- counts[genes,]

  if (mode == "general"){
    
    #Create an empty matrix where you will add the upregulated genes.
    genesDE <- list()
    Up.genes <-matrix(, nrow=topGeneThres , ncol = length(uniq.groups))
    Down.genes <-matrix(, nrow=topGeneThres , ncol = length(uniq.groups))
    if (is.null(groupTarget)){
      
      for (i in c(1:length(uniq.groups))){
        
        
        groupB <- as.character(uniq.groups[i])
        cat (paste("\n\n\nComparing",groupB,"to others...\n"))
        
        groupA <- "other"
        groups.bias <- as.character(groups)
        names(groups.bias)= names(groups)
        
        
        
        groups.bias[names(groups.bias)[which(groups.bias != groupB)]] = groupA
        groups.bias <- as.factor(groups.bias)
        #Create a DGEList
        dgList <- edgeR::DGEList(counts=as.matrix(counts), genes = rownames(counts), 
                                 group=groups.bias)
        
        
        #Estimate dispersions needed for edgeR.
        
        
        cat("Estimating dispersions...")
        dgList <- edgeR::estimateCommonDisp(dgList, verbose=TRUE)
        dgList <- edgeR::estimateTagwiseDisp(dgList)
        
        
        res <- edgeR::exactTest(dgList,pair = c(groupA , groupB))
        fdr <- stats::p.adjust(res$table$PValue, method="BH")
        res2 <- cbind(res$table,fdr)
        res2.sort <- res2[order(res2$fdr,decreasing = FALSE),]
        
        #res2.sort.top <- res2.sort[which(res2.sort$fdr < qvalue),]
        #up <- rownames(res2.sort.top[which(res2.sort.top$logFC<(-logFCThres)),])[1:topGeneThres]
        
        #down <- rownames(res2.sort.top[which(res2.sort.top$logFC>(logFCThres)),])[1:topGeneThres]
        
        genesDE[[groupB]] <- res2.sort[res2.sort$fdr<qvalue,]
      }
      
      #olnames(Up.genes) = colnames(Down.genes) = uniq.groups
      
      return(genesDE)
    } else {
      
      groupB <- groupTarget
      cat (paste("\n\n\nComparing",groupB,"to others...\n"))
      groupA <- "other"
      groups.bias <- as.character(groups)
      names(groups.bias)= names(groups)
      
      
      groups.bias[names(groups.bias)[which(groups.bias != groupB)]] = groupA
      groups.bias <- as.factor(groups.bias)
      #Create a DGEList
      dgList <- edgeR::DGEList(counts=as.matrix(counts), genes = rownames(counts), 
                               group=groups.bias)
      
      
      #Estimate dispersions needed for edgeR.
      
      
      cat("Estimating dispersions...")
      dgList <- edgeR::estimateCommonDisp(dgList, verbose=TRUE)
      dgList <- edgeR::estimateTagwiseDisp(dgList)
      
      
      res <- edgeR::exactTest(dgList,pair = c(groupA , groupB))
      fdr <- stats::p.adjust(res$table$PValue, method="BH")
      res2 <- cbind(res$table,fdr)
      res2.sort <- res2[order(res2$fdr,decreasing = FALSE),]
      
      res2.sort.top <- res2.sort[which(res2.sort$fdr < qvalue),]
      
    

      return(res2.sort.top)
      
    }
    cat("DONE!\n")
    
  } else if (mode=="strict"){
    #Create a DGEList
    dgList <- edgeR::DGEList(counts=as.matrix(counts), genes = rownames(counts), 
                      group=groups)
    
    
    #Estimate dispersions needed for edgeR.
    
    
    cat("\n\nEstimating dispersions...\n")
    dgList <- edgeR::estimateCommonDisp(dgList, verbose=TRUE)
    dgList <- edgeR::estimateTagwiseDisp(dgList)
    
    
    #Create table for pairwise comparisons (pairs will be unique)
    combin <- t(utils::combn(unique(c(uniq.groups,uniq.groups)) , 2 ))
    
    #Reorder columns
    combin.inv <- combin[,c(2,1)]
    
    #Take all pairs in one direction and in the other (e.g. row with 6 4 and 4 6)
    combin2 <- rbind(combin,combin.inv)
    
    #Create an empty matrix where you will add the upregulated genes.
    Up.genes <-matrix(, nrow=topGeneThres , ncol = nrow(combin2))
    
    
    for (i in c(1:nrow(combin))){
      
      #for each unique combination, define group1 and 2.
      group1 <- uniq.groups[combin[i,1]]
      group2 <- uniq.groups[combin[i,2]]
      
      res <- edgeR::exactTest(dgList,pair = c( group1 , group2))
      fdr <- stats::p.adjust(res$table$PValue, method="BH")
      res2 <- cbind(res$table,fdr)
      res2.sort <- res2[order(res2$fdr,decreasing = FALSE),]
      
      res2.sort.top <- res2.sort[which(res2.sort$fdr < qvalue),]
      up1 <- rownames(res2.sort.top[which(res2.sort.top$logFC<(-logFCThres)),])[1:topGeneThres]
      
      idx1 <- which(((combin2[,1]==combin[i,1])&(combin2[,2]==combin[i,2]))==TRUE)
      Up.genes[,idx1] <- up1
      
      up2 <- rownames(res2.sort.top[which(res2.sort.top$logFC>(logFCThres)),])[1:topGeneThres]
      idx2 <- which(((combin2[,1]==combin[i,2])&(combin2[,2]==combin[i,1]))==TRUE)
      Up.genes[,idx2] <- up2
      
      
    }
    
    up.genes.merged <- list()
    for (t in 1:length(uniq.groups)){
      num <- length(which(combin2[,1]==t))
      tab <- table(as.vector(Up.genes[,which(combin2[,1]==t)]))
      up.genes.merged[[t]] <- names(tab)[which(tab == num)]
      if (is.null(names(tab)[which(tab == num)])){
        up.genes.merged[[t]] = c("none")
      }
    }
    names(up.genes.merged)<-uniq.groups
    return(list(
      up.genes <- up.genes.merged
    ))
    
  cat("DONE!\n")
  }
  
    
}
  
    
  



