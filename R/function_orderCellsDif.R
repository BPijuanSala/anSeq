########################################################################################
## Title: function_orderCellsDiff.R
## Author: Blanca Pijuan-Sala
## Description: Order cells along a trajectory using principal curves
## Date: 04 February 2018
## **anSeq package**
########################################################################################


#' @title Cell ordering along a trajectory
#' @description It orders the cells along a trajectory using principal curves. 
#' For this, it first takes the fb$trall, which has the cells stored in the 
#' order they have been classified, and takes those cells belonging to the 
#' tEnd trajectory, keeping the order. With this order, the algorithm calculates
#' how many cells were assigned to differentiated clusters -specified in DifClusts, 
#' with reference to clustCells - and also how many cells were assigned to 
#' undifferentiated clusters - specified in OriginClusts with reference to 
#' clustCells. Next, it calculates the ratio between them 
#' - #OriginCells/#Target Cells - and sets a cutoff in the iteration where
#' there the maximum ratio is found - i.e. this would be when the classification
#' has reached the area of the undifferentiated state for that target trajectory.
#' With this filtered subset of cells belonging to this trajectory, the algorithm
#' then calculates a principal curve fitting these cells. The values used for this
#' are the ones passed through the 'manifold' parameter. With this principal curve,
#' the order of the cells in respect to it is calculated.
#' @param fb Output object resulting from the fateBias function from FateID. 
#' #Note: the original FateID package had some issues so I have fixed the function.
#' @param tEnd Trajectory to take to order the cells. This should correspond to 
#' one of the rownames of fb$probs.
#' @param coorPlot Matrix or dataframe with 2 columns, corresponding to
#' the coordinates of the cells in a 2D plot - e.g. PCA, tSNE, FR...-
#  The rownames should be labelled with the cellnames.
#' @param clustCells Vector specifying the cluster each cell belongs to.
#' This vector should have the cell names as names. 
#' @param OriginClusts Clusters to consider as origin. These will be used to 
#' define where the randomForest classification should no longer be considered.
#' @param DifClusts Clusters to consider as differentiated. These will be used to 
#' define where the randomForest classification should no longer be considered.
#' @param windowSize window size to calculate the ratio between the number of
#' cells belonging to origin clusters and the number of cells belonging to 
#' the target clusters. Default: 50
#' @param stepSize Step size to consider when going through the cells ordered
#' by the iteration number where they have been classified. This is used when 
#' calculating the ratio between the number of cells belonging to origin 
#' clusters and the number of cells belonging to the target clusters.Default: 20
#' @param cutEndMode Specify the method with which the cut threshold 
#' for filtering out the last classified cells should be if endIdx=NULL.
#' Options: "fit": This will fit a curve with loess and the max point 
#' of this fit will be chosen. "firstMax": This will take the value 
#' where the first maximum is achieved. "lastMax": This will take the value 
#' where the last maximum is achieved. Default: fit.
#' @param manifold Matrix or dataframe with cells in rows and features in columns
#' that is used to calculate the principal curve. This should at least include
#' the cells belonging to the trajectory of interest. This could be a PCA 
#' matrix, diffusion map, tSNE, counts...
#' @param strengthen Logical. If TRUE, cells belonging to the start and to 
#' the end of the filtered trajectory - according to the order in fb$trall -
#' will be repeated n times in the manifold matrix to make sure the 
#' principal curve goes through them. Default: TRUE
#' @param strenCells Numerical. This will be used if strengthen=TRUE.
#' Number of cells to consider as 'start' and 'end' cells. These cells will
#' be the ones repeated in the manifold if strengthen=TRUE. Default: 50
#' @param strenPower Numerical. This will be used if strengthen=TRUE. Number 
#' of times to repeat the values of Start and End cells. Default: 10
#' @param colIntensity Color key to plot the ordering of the trajectory.
#' This will be a gradient. Default: "grey","gold","orange","red","blue"
#' @param endIdx If one has the end cell along fb$trall already defined, you
#' can use this parameter to specify its index. Default: NULL. It will calculate
#' it based on the #Origin/#Target ratio.
#' @note It starts from a FateID object obtained from the function fateBias.
#' @author Blanca Pijuan-Sala.
#' @export
#' @rdname orderCellsDif
#' @importFrom princurve principal.curve
#' @importFrom graphics lines
#' @importFrom stats loess lm predict

##Example
##Heart1.obj <- orderCellsDif(fb,tEnd="theart1",coorPlot=fr[,1:2],clustCells=clustColLouv,
#                            OriginClusts=c("green","darkgreen","lightcyan","salmon", "greenyellow","grey60"),
#                            DifClusts=c("brown","turquoise","black","red","cyan","magenta"),
#                            windowSize=50,stepSize=20,manifold=countsPCA$x[,1:90],
#                            strengthen=TRUE,strenCells=50,strenPower=10,
#                            colIntensity=c("grey","gold","orange","red","blue"),
#                            endIdx=NULL)

##tEnd<-"theart1"
##coorPlot <- fr[,1:2]
##clustCells<-clustColLouv
##OriginClusts=c("green","darkgreen","lightcyan","salmon", "greenyellow","grey60")
##DifClusts=c("brown","turquoise","black","red","cyan","magenta")
##windowSize=50
##stepSize=20
##manifold<-countsPCA$x[,1:90]
##strengthen=TRUE
##strenCells=50
##strenPower=10
##colIntensity=c("grey","gold","orange","red","blue")
##endIdx=NULL

orderCellsDif <- function(fb,tEnd,coorPlot,clustCells,OriginClusts,
                          DifClusts,windowSize=50,stepSize=20,cutEndMode="fit",
                          manifold,strengthen=TRUE,strenCells=50,strenPower=10,
                          colIntensity=c("grey","gold","orange","red","blue"),
                          endIdx=NULL
                          ){
  plotsOutput <- list()
  
  par(mfrow=c(2,3))
  likelihoodLineages <- t(fb$probs)[,rownames(coorPlot)]
  
  lin <- tEnd
  cellsLin <- fb$trall[fb$trall%in%fb$tr[[lin]]]
  values <- sort(seq(1,length(cellsLin),1),decreasing=TRUE)
  names(values)<-cellsLin
  ValuesAll <- rep(0,nrow(coorPlot))
  names(ValuesAll)<- rownames(coorPlot)
  ValuesAll[names(values)] <- values
  
  cellsLinSort <-sort(ValuesAll[ValuesAll>0],decreasing=TRUE)
  
  cat("Plotting cells from trajectory before filtering...\n")
  plot(coorPlot[,1],coorPlot[,2],pch=20,col="gray",
       xlab="Coordinate 1",ylab="Coordinate 2",main="Cells from trajectory")
  points(x=coorPlot[names(cellsLinSort),1],coorPlot[names(cellsLinSort),2],col="blue",pch=20)
  
  
  plotsOutput[["plot2Dbefore"]] <- function(){
    plot(coorPlot[,1],coorPlot[,2],pch=20,col="gray",
         xlab="Coordinate 1",ylab="Coordinate 2",main="Cells from trajectory")
    points(x=coorPlot[names(cellsLinSort),1],coorPlot[names(cellsLinSort),2],col="blue",pch=20)
    
  }
  d <- clustCells[names(cellsLinSort)]
  
  cat("Calculating ratio origin cells vs target cells throughout iterations...\n")
  countTimes <- 0
  window <- windowSize
  step <- stepSize
  nBin <- round(length(d)/step)
  ratioEpiDif <- c()
  for (i in 1:nBin){
    start <- countTimes+1
    end <- countTimes+window
    countTimes <- countTimes+step
    if (length(d)<end){
      end <- length(d)
      subset <- d[start:length(d)]
    } else {
      subset <- d[start:end]
      
    }
    epiFrac <- length(which(subset %in% OriginClusts ))
    difFrac <- length(which(subset %in% DifClusts))
    rat <- round(epiFrac/(difFrac+0.1),3)
    names(rat)<-paste0("end",end)
    ratioEpiDif <-c(ratioEpiDif,rat)
    
    
  }
  cat("Plotting ratio origin vs target...\n")

  plot(x=c(1:length(ratioEpiDif)),y=ratioEpiDif,pch=20,
       ylab="#Ori/#Target",xlab="iteration",main="Ratio OriginCells vs TargetCells")
  lines(x=c(1:length(ratioEpiDif)),y=ratioEpiDif,col="blue")
  fit1 <- lm(ratioEpiDif ~ c(1:length(ratioEpiDif)))
  
  plotsOutput[["ratioOriginvsTar"]] <- function(){
    plot(x=c(1:length(ratioEpiDif)),y=ratioEpiDif,pch=20,
         ylab="#Ori/#Target",xlab="iteration",main="Ratio OriginCells vs TargetCells")
    lines(x=c(1:length(ratioEpiDif)),y=ratioEpiDif,col="blue")
  }
    
  
  lo <- loess(ratioEpiDif~c(1:length(ratioEpiDif)))
  cat("Plotting fitted line on ratio origin vs target...\n")
  
  plot(x=c(1:length(ratioEpiDif)),y=ratioEpiDif,pch=20,
       ylab="#Ori/#Target",xlab="iteration",main="Fitting Line - Ratio OriginCells vs TargetCells")
  lines(predict(lo), col='red', lwd=2)
  
  plotsOutput[["ratioOriginvsTarFittedLine"]] <- function(){
    plot(x=c(1:length(ratioEpiDif)),y=ratioEpiDif,pch=20,
         ylab="#Ori/#Target",xlab="iteration",main="Fitting Line - Ratio OriginCells vs TargetCells")
    lines(predict(lo), col='red', lwd=2)
  }
  
  
  idx <- which(predict(lo) == max(predict(lo)))
  
  
  if(is.null(endIdx)){
    if(cutEndMode == "fit"){
      cat("Taking threshold as maximum of fit...\n")
      endCut <- as.numeric(gsub("end","",names(ratioEpiDif)[idx]))
      
      plotEnd <- idx
      lines(x=rep(as.numeric(plotEnd),2),y=c(min(ratioEpiDif),max(ratioEpiDif)),col="red")
      
      } else if (cutEndMode == "firstMax"){
      cat("Taking threshold as first maximum value...\n")
      endCut <- as.numeric(gsub("end","",names(ratioEpiDif)[which(ratioEpiDif==max(ratioEpiDif))])[1])
      
      plotEnd <- which(ratioEpiDif==max(ratioEpiDif))[1]
      lines(x=rep(as.numeric(plotEnd),2),y=c(min(ratioEpiDif),max(ratioEpiDif)),col="red")
      
      } else if (cutEndMode == "lastMax"){
      cat("Taking threshold as first maximum value...\n")
      endCut <- as.numeric(gsub("end","",names(ratioEpiDif)[which(ratioEpiDif==max(ratioEpiDif))])[length(which(ratioEpiDif==max(ratioEpiDif)))])
      
      plotEnd <- which(ratioEpiDif==max(ratioEpiDif))[length(which(ratioEpiDif==max(ratioEpiDif)))]
      lines(x=rep(as.numeric(plotEnd),2),y=c(min(ratioEpiDif),max(ratioEpiDif)),col="red")
      } else { 
        cat("Value in cutEndMode does not exist. Using fit...\n")
        
        cat("Taking threshold as maximum of fit...\n")
        endCut <- as.numeric(gsub("end","",names(ratioEpiDif)[idx]))
        
        plotEnd <- idx
        lines(x=rep(as.numeric(plotEnd),2),y=c(min(ratioEpiDif),max(ratioEpiDif)),col="red")
      }
    
  } else {endCut<-endIdx}
  
  cellsToTake <- names(d)[1:endCut]
  cellsStart <- names(sort(ValuesAll[cellsToTake],decreasing=FALSE)[1:50])
  cellsEnd <- names(sort(ValuesAll[cellsToTake],decreasing=TRUE)[1:50])
  
  
  cat("Plotting start (blue) and end (red) cells after filtering...\n")
  plot(coorPlot[,1],coorPlot[,2],pch=20,col="gray",
       xlab="Coordinate 1",ylab="Coordinate 2",main="After Filter - Start (blue); end (red)")
  points(x=coorPlot[as.character(cellsStart),1],coorPlot[as.character(cellsStart),2],col="blue",pch=20)
  points(x=coorPlot[cellsEnd,1],coorPlot[cellsEnd,2],col="red",pch=20)

  
  plotsOutput[["plot2DstartEnd"]] <- function(){
    plot(coorPlot[,1],coorPlot[,2],pch=20,col="gray",
         xlab="Coordinate 1",ylab="Coordinate 2",main="After Filter - Start (blue); end (red)")
    points(x=coorPlot[as.character(cellsStart),1],coorPlot[as.character(cellsStart),2],col="blue",pch=20)
    points(x=coorPlot[cellsEnd,1],coorPlot[cellsEnd,2],col="red",pch=20)
  }
  
  if(strengthen==TRUE){

    
    cellsStart <- names(sort(ValuesAll[cellsToTake],decreasing=FALSE)[1:strenCells])
    cellsEnd <- names(sort(ValuesAll[cellsToTake],decreasing=TRUE)[1:strenCells])
    
    weight1 <- do.call("rbind", replicate(strenPower, manifold[cellsStart,], simplify = FALSE))
    weight2 <- do.call("rbind", replicate(strenPower, manifold[cellsEnd,], simplify = FALSE))
    cellsPrinCurve <- c(cellsToTake,rownames(weight1),rownames(weight2))
    manifold1.0 <- manifold[cellsToTake,]
    rownames(weight1) <- rownames(weight2)<- rownames(manifold1.0) <- NULL
    manifoldWeight <- rbind(manifold1.0,weight1,weight2)
    
  } else{
    cellsPrinCurve <- c(cellsToTake)
    manifoldWeight <- manifold[cellsToTake,]
    
  }
  cat("Calculating principal curve...\n")
  
  prinCurve <- princurve::principal.curve(manifoldWeight,plot=FALSE)
  
  orderCells <- unique(cellsPrinCurve[prinCurve$tag])
  if (sum(which(orderCells%in%cellsEnd))<sum(which(orderCells%in%cellsStart))){
    orderCells <- rev(orderCells)
  }
  cat("Plotting cells participating in the calculation...\n")
  
  plot(coorPlot[,1],coorPlot[,2],pch=20,col="gray",
       xlab="Coordinate 1",ylab="Coordinate 2",main="After filter: Cells in trajectory")
  points(x=coorPlot[orderCells,1],coorPlot[orderCells,2],col="blue",pch=20)

  
  plotsOutput[["plot2DCellsPrinCurve"]] <- function(){
    plot(coorPlot[,1],coorPlot[,2],pch=20,col="gray",
         xlab="Coordinate 1",ylab="Coordinate 2",main="After filter: Cells in trajectory")
    points(x=coorPlot[orderCells,1],coorPlot[orderCells,2],col="blue",pch=20)
    
  }
  cellsLin <- orderCells
  values <- sort(seq(1,length(cellsLin),1),decreasing=FALSE)
  names(values)<-cellsLin
  ValuesAll <- rep(0,nrow(coorPlot))
  names(ValuesAll)<- rownames(coorPlot)
  ValuesAll[names(values)] <- values
  valuesMock <- ValuesAll[rownames(coorPlot)]

  
  
  cols <- colIntensity
  redRamp <- colorRampPalette(cols)
  df <- data.frame(x = coorPlot[,1], y = coorPlot[,2], exp = valuesMock)
  df <- df[order(df$exp,decreasing=F),]
  dfsub <- df[df$exp>0,]
  
  interval <- findInterval(dfsub$exp, seq(min(dfsub$exp), 
                                          max(dfsub$exp), 
                                          (max(dfsub$exp)-min(dfsub$exp))/10))
  
  interval[interval==0]<-1
  colorsPlot <- redRamp(11)[interval]
  
  plot(df$x, df$y, col=cols[1], pch=16, cex=0.75,
       xlab="", ylab="", main="Order Differentiation",
       axes=F, cex.main=1.5)
  box(bty="l")
  points(dfsub$x, dfsub$y, pch=20, cex=0.5, 
         col=colorsPlot)
  
  
  plotsOutput[["plot2DTrajectory"]] <- function(){
    
    plot(df$x, df$y, col=cols[1], pch=16, cex=0.75,
         xlab="", ylab="", main="Order Differentiation",
         axes=F, cex.main=1.5)
    box(bty="l")
    points(dfsub$x, dfsub$y, pch=20, cex=0.75, 
           col=colorsPlot)
    
  }
  return(list(cellsOrdered = orderCells,
              prinCurve=prinCurve,
              ratioOriDif = ratioEpiDif,
              fitLoess = lo,
              endCutIdx=endCut,
              CellsTrajtrall=cellsLinSort,
              plotsOutput=plotsOutput)
         )
  }


