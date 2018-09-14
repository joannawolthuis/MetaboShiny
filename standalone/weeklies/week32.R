# - - - network analysis - - -

paper.files <- "/Users/jwolthuis/Downloads/1-s2.0-S1570023217305688-mmc12/"
setwd(paper.files)

# need a matrix to begin with
# use HMDB as source
# conn <- RSQLite::dbConnect(RSQLite::SQLite(), patdb)
# chosen.db <- "/Users/jwolthuis/Google Drive/MetaboShiny/backend/db/hmdb.full.db"
# adducts <- c(pos_adducts$Name, neg_adducts$Name)
# 
# res <- get_all_matches(conn, db, adducts)
# 

# - - ggm finding - - 

source <- mSet$dataSet$norm

require(GeneNet)

inferred.pcor <- ggm.estimate.pcor(source)
xtest.results <- network.test.edges(inferred.pcor)
xtest.results[1:20,]

# extract network containing edges with prob > 0.9 (i.e. local fdr < 0.1)

net <- extract.network(xtest.results, cutoff.ggm=0.9)

head(net)

# try to visualize??

named.net <- net
orig_order <- 1:ncol(source)
named.net$node1 = as.numeric(colnames(source)[match(named.net$node1, orig_order)])
named.net$node2 = as.numeric(colnames(source)[match(named.net$node2, orig_order)])

# show/plot subset ...

all_edges <- named.net[1:100,]

get.nodes <- function(netobj){
  unlist(lapply(1:nrow(netobj), function(j){
    unlist(netobj[j,c("node1", "node2")])
  }))
}

g <- igraph::make_graph(edges = as.character(get.nodes(all_edges)), directed = FALSE)

plot(g,
     edge.arrow.size=.5,
     vertex.size=5)

# - - 

require(data.table)

score.netw.add <- function(cpd, patdb, ppm){

  conn <- RSQLite::dbConnect(RSQLite::SQLite(), patdb)
  
  results <- unique(DBI::dbGetQuery(conn, "SELECT * FROM adducts"))
  
  keep.cpds <- pbapply::pbsapply(unique(results$fullmz), function(cpd){
    res <- RSQLite::dbGetQuery(conn, gsubfn::fn$paste(strwrap(
      "SELECT DISTINCT count(*)
      FROM mzvals mz
      JOIN mzranges rng ON rng.ID = mz.ID
      WHERE $cpd BETWEEN rng.mzmin AND rng.mzmax",
      width=10000, 
      simplify=TRUE)))[,1]  
    if(res > 0) cpd else NA
  })
  
  keepos <- keep.cpds[!is.na(keep.cpds)]
  
  results <- data.table::as.data.table(results)[fullmz %in% keepos,]
  
  split.by.formula <- split(results, list(results$baseformula))
  split.by.formula <- split.by.formula[which(sapply(split.by.formula, function(x) nrow(x) > 0))]
  
  involved.rows <- which(named.net$node1 == cpd | named.net$node2 == cpd)
  network.of.mz <- as.data.table(named.net[involved.rows,])
  
  # make all node1 the orig mz
  need.swap <- which(network.of.mz$node2 == cpd)
  
  for(i in need.swap){
    store.1 = network.of.mz[i, "node1"][[1]]
    store.2 = network.of.mz[i, "node2"][[1]]
    # - - - -
    network.of.mz$node1[i] <- store.2
    network.of.mz$node2[i] <- store.1
  }
  
  g <- igraph::make_graph(edges = as.character(get.nodes(network.of.mz)), directed = FALSE)

  # - - - - - - - - - 

  scores <- pbapply::pblapply(split.by.formula, function(current){
    current = unique(current[,-9])
    print(paste("Scoring", unique(current$baseformula)))
    matches <- pbapply::pblapply(current$fullmz, function(mz){
      #print(paste("---", mz, "---"))
      bounds <- c(mz - mz*(ppm*1e-6), mz + mz*(ppm*1e-6))
      node.matches <- network.of.mz[as.numeric(node2) %between% bounds,]
      # - - -
      node.matches
    })
    score = sum(sapply(matches, function(x) if(nrow(x)>0) 1 else 0))
    print(paste("Total rows:", nrow(current)))
    print(paste("Rows with match:", score))
    percscore = round(score/nrow(current)*100.00, digits=2)
    print(paste("Fractional score (intra-adduct):", percscore,"%"))
    
    data.table(baseformula = unique(current$baseformula), score = percscore)
  })
  # - - - - - - -
  rbindlist(scores)
}

cpd <- curr_cpd
patdb <- global$paths$patdb
ppm <- global$constants$ppm

score <- score.netw.add(cpd, patdb, ppm)

View(global$tables$last_matches[unique(score), on = c("baseformula")])

g <- igraph::make_graph(edges = as.character(get.nodes(network.of.mz)), directed = FALSE)
plot(g,edge.arrow.size=.5,vertex.size=5)


# === METHOD 2 == CHECK REACTIONS ===

require(enviPat)

# loop through all pairs

# molecule modification list

adduct = c()
deduct = c()

# - - paper files - - 

source("configurationFile.R")
source("R/constructNetworkModel.R")
source("R/getMetabolitePairs.R")
source("R/assignReactionsByMass.R")
source("R/filterReactions.R")
source("R/summarizeUnknownAnnotation.R")
source("R/filterColumnByNA.R")
source("R/learnPathways.R")
source("R/predictPathways.R")
source("R/exportGraph.R")