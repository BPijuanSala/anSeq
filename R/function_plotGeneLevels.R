########################################################################################
## Title: function_plotGeneLevels.R
## Author: Blanca Pijuan Sala
## Description: Plot gene expression levels by color.
## Date: 21 February 2017
## **anSeq package**
########################################################################################


#' @title Plot gene expression
#' @description It plots the level of gene expression of a particular gene
#'  within the plot you want.
#' @param data Dataframe or matrix containing the counts (genes x cells). The cells
#' should correspond to the ones used for the plot. This matrix is only used to take 
#' the values of the row of the gene to plot.
#' @param x x coordinate. It needs to have entries per each cell in the counts matrix. 
#' e.g. tsne$Y[,1] when performing tSNE on those cells.The plot does not need to be based
#' all genes.
#' @param y y coordinate. It needs to have entries per each cell in the counts matrix. 
#' e.g. tsne$Y[,2] when performing tSNE on those cells.The plot does not need to be based
#' all genes.
#' @param gene Gene that you want the expression of. This should be the row ID in 
#' your data table.
#' @param cols Colour range for the expression values. Default: blue,gold.
#' @param xlab x-label. Default: "x".
#' @param ylab y-label. Default: "y".
#' @param titlePlot plot title. Default: the gene name introduced in "gene".
#' @param boxplot Logical. TRUE: It will plot a boxplot next to it with the gene 
#' expression levels per cluster.
#' @param clusters If boxplot=TRUE, this should be specified. It should be a vector of
#' cells as names and colors defining the cluster the cell belongs to. The colours will
#' be used for colouring the boxplots.
#' @param pointSize cex size for points. Default: 1.
#' @note The code where boxplot=TRUE is adapted from Ximena Ibarra-Soria.
#' @author Blanca Pijuan-Sala.
#' @export
#' @rdname plotGeneLevels
#' @import ggplot2
#' @importFrom graphics par box points rect segments text plot layout layout.show
#' @importFrom utils head tail
#' @importFrom stats median
plotGeneLevels <- function(data, x, y, gene, cols=c("blue","gold"),
                           pointSize=1,
                             xlab="x",ylab="y",titlePlot=gene, 
                           clusters=NULL, boxplot=FALSE)
{
  if (boxplot==FALSE){
    

    redRamp <- colorRampPalette(cols)
    df <- data.frame(x = x, y = y, exp = data[gene,])
    df <- df[order(df$exp,decreasing=F),]
    dfsub <- df[df$exp>0,]
    
    interval <- findInterval(dfsub$exp, seq(min(dfsub$exp), 
                                            max(dfsub$exp), 
                                            (max(dfsub$exp)-min(dfsub$exp))/10))
    
    interval[interval==0]<-1
    colorsPlot <- redRamp(11)[interval]
    
    par(mar=c(8,4,8,4),xpd=NA)
    plot(df$x, df$y, col=cols[1], pch=20, cex=pointSize,
         xlab="", ylab="", main=titlePlot, axes=F, cex.main=1.5)
    box(bty="l")
    points(dfsub$x, dfsub$y, pch=20, cex=pointSize, 
           col=colorsPlot)
    
    xl <- min(df$x)- (min(df$x)*(-5.57*10^(-4))); yb <- (min(df$y))-(-0.34*min(df$y)); xr <- max(df$x); yt <- (min(df$y))-(min(df$y)*(-0.17))
    rect(xleft=head(seq(xl,xr,(xr-xl)/10),-1), ybottom = yb, 
         xright = tail(seq(xl,xr,(xr-xl)/10),-1), ytop = yt,
         col=redRamp(10), border=redRamp(10), lwd=0.5, xpd=NA)
    rect(xl,yb,xr,yt, xpd=NA)
    segments(seq(xl,xr,((xr-xl)/5)),yb,seq(xl,xr,((xr-xl)/5)),yb-(-yb*0.0381), xpd=NA)
    text(seq(xl,xr,((xr-xl)/5)), yb-(-yb*0.089), labels = round(seq(min(dfsub$exp), max(dfsub$exp), (max(dfsub$exp)-min(dfsub$exp))/5),1), cex=1, xpd=NA)
    text(stats::median(seq(xl,xr,((xr-xl)/5))), yb-(-yb*0.216), labels = expression('log'[10]*' normalized counts + 1'),cex=1.2, xpd=NA)
    
    #data.order <- order(data[gene,],decreasing=FALSE)
    #plot = qplot(x=x[data.order], y=y[data.order], colour=data[gene,][data.order])
    #print(plot + theme_bw() + labs(x = xlab, y=ylab) + scale_color_gradient("Expression",low=col.low, high=col.high) + 
#            ggtitle(titlePlot) + theme(legend.text = element_text(size = 15),
 #                                      legend.title = element_text(size = 20, face="bold", vjust=1),
  #                                     axis.title = element_text(size = 15),axis.text=element_text(size=15), 
   #                                    plot.title = element_text(vjust= 1,hjust=0.5, size = 20 ,face = "bold"), 
    #                                   panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
     #                                  panel.background = element_blank())) 
    
  } else {
    if (is.null(clusters)){
      stop("Please specify the 'clusters' section.")
    }
    #bluefunc <- colorRampPalette(c("lightblue", "darkblue"))
    redRamp <- colorRampPalette(cols)
    clust = clusters[rownames(x)]
    df <- data.frame(x = x, y = y, exp = data[gene,],clust=clusters)
    df <- df[order(df$exp,decreasing=F),]
    #dfsub <- df
    dfsub <- df[df$exp>0,]
    interval <- findInterval(dfsub$exp, seq(min(dfsub$exp), 
                                                       max(dfsub$exp), 
                                                       (max(dfsub$exp)-min(dfsub$exp))/10))
    
    interval[interval==0]<-1
    colorsPlot <- redRamp(11)[interval]
    #dfsub <- df[df$exp>0,]
    nf <- layout(matrix(c(1,2),1,2,byrow = TRUE), c(3,3),c(3), TRUE)
    layout.show(nf)
    par(mar=c(4,2,4,2),xpd=NA)
    #bottom, left, top and right
    
    plot(df$x, df$y, col=cols[1], pch=20, cex=pointSize,
         xlab="", ylab="", main=titlePlot, axes=F, cex.main=1.5)
    box(bty="l")
    points(dfsub$x, dfsub$y, pch=20, cex=pointSize, 
           col=colorsPlot )
    
    xl <- min(df$x)- (min(df$x)*(-5.57*10^(-4))); 
    yb <- (min(df$y))-(-0.64*min(df$y)); 
    xr <- max(df$x); 
    yt <- (min(df$y))-(min(df$y)*(-0.47))
    rect(xleft=head(seq(xl,xr,(xr-xl)/10),-1), ybottom = yb, 
         xright = tail(seq(xl,xr,(xr-xl)/10),-1), ytop = yt,
         col=redRamp(10), border=redRamp(10), lwd=0.5, xpd=NA)
    rect(xl,yb,xr,yt, xpd=NA)
    segments(seq(xl,xr,((xr-xl)/5)),yb,seq(xl,xr,((xr-xl)/5)),yb-(-yb*0.1381), xpd=NA)
    text(seq(xl,xr,((xr-xl)/5)), yb-(-yb*0.289), labels = round(seq(min(dfsub$exp), max(dfsub$exp), (max(dfsub$exp)-min(dfsub$exp))/5),1), cex=1, xpd=NA)
    text(stats::median(seq(xl,xr,((xr-xl)/5))), yb-(-yb*0.616), labels = expression('log'[10]*' normalized counts + 1'),cex=1.2, xpd=NA)

    par(mar=c(4,3,4,1))
    boxplot(df$exp~df$clust, col=sort(as.character(unique(df$clust))), 
            main=titlePlot, ylab=expression('Expression'), las=2)
    
    
    
  }
  
  
}  

#plotGeneLevels(data = counts.endo[,cell.order], tsne[cell.order,1], tsne[cell.order,2],hbb_genes[1])
