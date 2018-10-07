########################################################################################
## Title: external_data_save.R
## Author: Blanca Pijuan-Sala
## Description: Save external data.
## Date: 30 September 2017
## **anSeq package**
########################################################################################

###SAVE EXTERNAL DATASET.
#http://r-pkgs.had.co.nz/data.html
wd <- "/Users/blancap/Documents/PhD_CAMBRIDGE/PhD/My_R_package/00.PACKAGES/anSeq_development/anSeq/"


geneTable <- read.csv(paste0(wd,"inst/extdata/mm10_mart_export.txt"))
save(geneTable, file=paste0(wd,"data/geneTable.rda"))

########################################################################################
#Surface markers
########################################################################################
##====MOUSE


tableRef="/Users/blancap/Documents/PhD_CAMBRIDGE/PhD/Experiments/bin_data/surfaceomeproteins_wlabethz_File.xlsx"
species="Mouse"
tab = "Table B"

tableRef <- xlsx::read.xlsx(tableRef,sheetName = tab)

MouseSurface <- tableRef$CSPA.category
names(MouseSurface) <- tableRef$ENTREZ.gene.symbol

save(MouseSurface,file="/Users/blancap/Documents/repos_bitbucket/r-packages/00.PACKAGES/anSeq_development/anSeq/data/MouseSurface.rda")
MouseSurface <- tableRef$CSPA.category

##====HUMAN

tableRef="/Users/blancap/Documents/PhD_CAMBRIDGE/PhD/Experiments/bin_data/surfaceomeproteins_wlabethz_File.xlsx"
species="Human"
tab = "Table A"

tableRef <- xlsx::read.xlsx(tableRef,sheetName = tab)

HumanSurface <- tableRef$CSPA.category
names(HumanSurface) <- tableRef$ENTREZ.gene.symbol

save(HumanSurface,file="/Users/blancap/Documents/repos_bitbucket/r-packages/00.PACKAGES/anSeq_development/anSeq/data/HumanSurface.rda")
