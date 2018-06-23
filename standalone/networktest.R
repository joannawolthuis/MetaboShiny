BiocInstaller::biocLite("NetPathMiner")
require(NetPathMiner)

SBML2igraph(filename, parse.as = c("metabolic", "signaling"),
            miriam.attr = "all", gene.attr, expand.complexes, verbose = TRUE)