help(jorissen.colon.ras)
library(data.table)
library(predictionet)


## load gene expression data for colon cancer data, list of genes related to RAS signaling pathway and the corresponding priors
data(expO.colon.ras)
## create matrix of perturbations (no perturbations in this dataset)
pert <- matrix(0, nrow=nrow(data.ras), ncol=ncol(data.ras), dimnames=dimnames(data.ras))

## number of genes to select for the analysis
genen <- 10
## select only the top genes
goi <- dimnames(annot.ras)[[1]][order(abs(log2(annot.ras[ ,"fold.change"])), decreasing=TRUE)[1:genen]]
mydata <- data.ras[ , goi, drop=FALSE]
myannot <- annot.ras[goi, , drop=FALSE]
mypriors <- priors.ras[goi, goi, drop=FALSE]
mydemo <- demo.ras
mypert <- pert[ , goi, drop=FALSE]

mydata
myannot
priors.ras
mypriors
## regression-based network inference
########################
## infer global network from data and priors
mynet <- netinf(data=mydata, 
                perturbations=mypert, 
                priors=mypriors, 
                priors.count=TRUE, 
                priors.weight=0.5, 
                maxparents=3, 
                method="regrnet", 
                seed=54321)

mynet
## plot network topology
mytopo <- mynet$topology
library(network)
xnet <- network(x=mytopo, matrix.type="adjacency", directed=TRUE, loops=FALSE, vertex.attrnames=dimnames(mytopo)[[1]])
plot.network(x=xnet, displayisolates=TRUE, displaylabels=TRUE, boxed.labels=FALSE, label.pos=0, arrowhead.cex=2, vertex.cex=4, vertex.col="royalblue", jitter=FALSE, pad=0.5)

## export network as a 'gml' file that you can import into Cytoscape
## Not run: rr <- netinf2gml(object=mynet, file="/predictionet_regrnet")


load("testMset.RData")

# --- OTHER OPTION...---

source("http://bioconductor.org/biocLite.R")
biocLite("pathview")
biocLite("rsbml")
biocLite("visNetwork")

nodes = data.frame(id=colnames(mytopo))
topomelt <- as.data.table(melt(mytopo))
colnames(topomelt) <- c("from", "to", "value")

nodes
edges = data.frame(topomelt[value == 1, -"value"])
edges

require(visNetwork, quietly = TRUE)
# minimal example
nodes <- data.frame(id = 1:3)
edges <- data.frame(from = c(1,2), to = c(1,3))
edges
visNetwork(nodes, edges)

library(MetaboAnalystR)

# === METABOANALHYST ENRICHMENT ===
mSet<-Setup.MapData(mSet);
mSet<-CrossReferencing(mSet, "name");
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "predicted", 0);
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_0_", "bar", "png", 72, width=NA)

# !! needs DB

# ===== GAGE ENRICHMENT ===

library(piano)
library(RSQLite)
library(DBI)
db <- "backend/db/smpdb.base.db"
conn <- dbConnect(RSQLite::SQLite(), db) # change this to proper var later

dbGetQuery(conn, "SELECT distinct * FROM pathways limit 10;")
gset <- dbGetQuery(conn,"SELECT DISTINCT  c.baseformula AS cpd,
                                          p.name as name
                              FROM pathways p
                              JOIN base c
                              ON c.pathway = p.identifier 
                               ")

gset_proc <- loadGSC(gset)
gsaRes <- runGSA(pvals, 
                 gsc=gset_proc)
summary <- GSAsummaryTable(gsaRes) 

