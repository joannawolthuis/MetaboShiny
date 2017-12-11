biocLite("topGO")
install.packages("visNetwork")
biocLite("ALL")

library(topGO)
library(visNetwork)
library(ALL)

data(ALL)
data(geneList)

ALL@assayData$exprs

geneList

# Jelle Goeman, gene ontology
# fisNetwork
# 'Application note' subgedeelte van bioinformatics / BMC
