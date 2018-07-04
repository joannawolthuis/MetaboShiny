# BiocInstaller::biocLite("NetPathMiner")
# require(NetPathMiner)
# SBML2igraph(filename, parse.as = c("metabolic", "signaling"),
#             miriam.attr = "all", gene.attr, expand.complexes, verbose = TRUE)
# shiny::runGitHub('MetaMapR','dgrapov') 
# install_github("dgrapov/MetaMapR")
# #https://www.bioconductor.org/packages/release/bioc/vignettes/graphite/inst/doc/graphite.pdf

BiocInstaller::biocLite("graphite")

library(graphite)

pigPWs <- pathways("sscrofa", "kegg")
names(pigPWs)[1]

#graphs <- lapply(1:5, function(i) igraph::igraph.from.graphNEL(pathwayGraph(pigReactome[[i]])))

edges <- lapply(1:length(pigPWs), function(i) edges(pigPWs[[i]], which="metabolites"))

all_edges <- unique(data.table::rbindlist(edges)[,1:2])

edges_reformat <- unlist(lapply(1:nrow(all_edges), function(j){
  unlist(all_edges[j,])
}))

g <- igraph::make_graph(edges = edges_reformat, directed = FALSE)

install.packages("MCL")
library(MCL)

giant.component <- function(graph) {
  cl <- igraph::clusters(graph)
  igraph::induced_subgraph(graph, which(cl$membership == which.max(cl$csize)))
}

largest_cmp <- giant.component(g)

#g <- largest_cmp

adjacency <- igraph::as_adjacency_matrix(largest_cmp, 
                                         type = "both",
                                         sparse = igraph::igraph_opt("sparsematrices"))

mcl_clust <- mcl(x = adjacency, addLoops = TRUE)

#coords <- igraph::layout_nicely(g)
#plot(g, layout = coords)

# =========== CLUSTERING (DAY 4) ==============

nmols = 100
expr <- mSet$dataSet$norm
m = apply(expr,1,mad)
ind = order(m,decreasing = TRUE,na.last = NA)
expr = expr[,ind[1:nmols]]

Dexpr <- 0.5*(1-cor(expr))
heatmap(Dexpr)

# =========== PETRI NETS (DAY 5) =============

install.packages("petrinetR")

library(petrinetR)

# create_PN(places, transitions, flows, marking)

places = c("STAT1", 
           "T-bet", 
           "GATA3", 
           "STAT6", 
           "IL-4R", 
           "IL-12R", 
           "IL-12",
           "SOCS1",
           "IFN-yR",
           "IFN-y",
           "IL-4",
           "STAT4")

transitions = c("STAT1|T-bet",
                "STAT1|IL-4",
                "STAT1|SOCS1",
                "IFN-yR|STAT1",
                "SOCS1|IFN-yR",
                "SOCS1|IL-4R",
                "T-bet|T-bet",
                "T-bet|GATA3",
                "T-bet|IFN-y",
                "T-bet|SOCS1",
                "STAT4|IFN-y",
                "STAT6|GATA3",
                "STAT6|IL-12R",
                "GATA3|IL-4",
                "GATA3|T-bet",
                "GATA3|STAT4",
                "IL-12|IL-12R",
                "IL-12R|STAT4",
                "IL-4|IL-4R")

flows <- data.table::rbindlist(lapply(transitions, function(pair){
  print(pair)
  split <- strsplit(x = pair, split = "\\|")[[1]]
  if(grepl(x = split[[1]], pattern = "_"))
    {
    split[[1]] <- gsub(x = split[[1]],
                       pattern = "inv_|tau_",
                       replacement = "")
    }
  data.frame(from = c(split[[1]], pair),
             to = c(pair, split[[2]]
                    ))
}))

PN <- create_PN(places, transitions, flows, marking=places)
visNetwork_from_PN(PN)

steps = 10

start_trans = "GATA3|T-bet"
PN_first_act <- execute(PN, start_trans)
visNetwork_from_PN(PN_first_act)
parsel(PN_first_act, "T-bet|SOCS1")

for(i in steps){
  ...
}

post_set(PN, "IL-4")

# =====================================================================

# - - - infer from mz - - -

source <- mSet$dataSet$norm[,1:20]

#install.packages("bnlearn")

library(bnlearn)
library(Rgraphviz)

hist(as.matrix(source))
plot(source[, 1], source[,2])

#heatmap(t(as.matrix(source)), col=colorRampPalette(c("blue", "white", "red"))(100), scale="row")

dsachs = discretize(source, method = "hartemink", breaks = 3, ibreaks = 60, idisc = "quantile")
## model average
boot = boot.strength(data = dsachs, R = 500, algorithm = "tabu", algorithm.args=list(tabu = 50))
boot[(boot$strength >= 0.85) & (boot$direction >= 0.5), ]
plot(boot)
avg.boot = averaged.network(boot, threshold = 0.85)
graphviz.plot(avg.boot)
