


tableRef="/Users/blancap/Documents/PhD_CAMBRIDGE/PhD/Experiments/bin_data/surfaceomeproteins_wlabethz_File.xlsx"
species="Mouse"
tab = "Table B"

tableRef <- xlsx::read.xlsx(tableRef,sheetName = tab)

MouseSurface <- tableRef$CSPA.category
names(MouseSurface) <- tableRef$ENTREZ.gene.symbol

save(MouseSurface,file="/Users/blancap/Documents/repos_bitbucket/r-packages/00.PACKAGES/anSeq_development/anSeq/data/MouseSurface.rda")
MouseSurface <- tableRef$CSPA.category
#========================

tableRef="/Users/blancap/Documents/PhD_CAMBRIDGE/PhD/Experiments/bin_data/surfaceomeproteins_wlabethz_File.xlsx"
species="Human"
tab = "Table A"

tableRef <- xlsx::read.xlsx(tableRef,sheetName = tab)

HumanSurface <- tableRef$CSPA.category
names(HumanSurface) <- tableRef$ENTREZ.gene.symbol

save(HumanSurface,file="/Users/blancap/Documents/repos_bitbucket/r-packages/00.PACKAGES/anSeq_development/anSeq/data/HumanSurface.rda")
