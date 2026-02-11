# install_metaboshiny_local.R
# Mirrors the posted Dockerfile (rocker/rstudio:4.1.0) as closely as possible.

req_major <- 4
req_minor <- 1

if (!(getRversion()[1] == req_major && getRversion()[2] == req_minor)) {
  warning(sprintf(
    "This install script mirrors rocker/rstudio:4.1.0 (R 4.1.x). You are running R %s. Expect breakage.",
    getRversion()
  ), call. = FALSE)
}

cran <- "https://cloud.r-project.org"

install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) install.packages(missing, repos = cran)
}

cat("== Core tools ==\n")
install_if_missing(c("devtools", "remotes", "BiocManager", "pacman"))

cat("== Bioc org.Hs.eg.db ==\n")
BiocManager::install("GenomeInfoDbData", ask = FALSE, update = FALSE)
BiocManager::install("org.Hs.eg.db", ask = FALSE, update = FALSE)

cat("== mvtnorm (source) ==\n")
if (!requireNamespace("mvtnorm", quietly = TRUE)) {
  install.packages("mvtnorm", repos = cran, type = "source")
}

cat("== METADBPARSE REQUIREMENTS (pacman::p_load) ==\n")
pacman::p_load(
  pacman, rcdk, rJava, parallel, pbapply, enviPat, data.table,
  RSQLite, DBI, gsubfn, utils, RCurl, XML, base,
  stringr, WikidataQueryServiceR, webchem, openxlsx, jsonlite,
  R.utils, KEGGREST, zip, ChemmineR, rvest, xml2, stringi, reshape2,
  Hmisc, httr, RJSONIO, readxl, cmmr, progress, Rdisop, rlist
)

install.packages(c("sparsevctrs", "S7"), type = "source")

cat("== METABOSHINY REQUIREMENTS (batch 1) ==\n")
pacman::p_load(
  lava,
  ggplot2, data.table, plotly, shinyBS, shinyjs, caret, grDevices,
  RColorBrewer, colorRamps, tidytext, qdapDictionaries, tm, shiny, htmltools,
  BiocManager, pacman, devtools, classyfireR, httr, jsonlite, RCurl, shinyFiles,
  DT, RSQLite, pbapply, stringr, gsubfn, shinyWidgets, parallel,
  mice, sva, limma, tools, plyr, heatmaply, wordcloud2, shinyjqui, rmarkdown,
  enviPat, ROCR, tsne, e1071
)

cat("== METABOSHINY REQUIREMENTS (batch 2) ==\n")
pacman::p_load(
  pls, rhandsontable, testthat, shinytest, showtext, sysfonts, colourpicker,
  reshape, ggdark, ECharts2Shiny, shinyalert, shinybusy, rcdk, RISmed, dplyr,
  InterpretMSSpectrum, DBI, qdap, reshape2, Hmisc, ggbeeswarm, Rmisc, rgl,
  stats, pROC, car, doParallel, missForest, ggfortify, fdrtool, plsdepot,
  vroom, umap, ica, svglite, beepr, showtextdb
)

cat("== extra loads ==\n")
pacman::p_load(ctc, gdata, glasso, huge, ppcor, crmn)

pacman::p_load(
  BiocParallel, IRanges, plyr, preprocessCore, vsn,
  grid, stats4, affy, impute, pcaMethods, MALDIquant, mzID,
  digest, lattice, ggplot2, XML, scales, MASS, Rcpp
)

pacman::p_load(multtest, siggenes, KEGGgraph, SSPA, preprocessCore,
               Rgraphviz, GlobalAncova, globaltest, pcaMethods, impute)

pacman::p_load(pdftools, magick)

pacman::p_load(MSnbase)
pacman::p_load(fgsea)

BiocManager::install("RBGL", ask = FALSE, update = FALSE)
pacman::p_load(crmn)

# Dockerfile repeats devtools + mvtnorm reinstall; we mirror it.
install_if_missing("devtools")

if (requireNamespace("mvtnorm", quietly = TRUE)) remove.packages("mvtnorm")
install.packages("mvtnorm", repos = cran)

install.packages("mutoss", repos = cran, type = "source")

pacman::p_load(char = c(
  "genefilter", "sva", "limma", "siggenes", "BiocParallel", "MSnbase",
  "multtest", "RBGL", "edgeR", "fgsea", "httr", "gdata", "glasso",
  "huge", "robustbase", "qqconf", "ppcor", "crmn", "plotly"
))

pacman::p_load(char = c("anytime", "rJava", "ggraph", "tidygraph", "WikidataR", "xmlparsedata"))
pacman::p_load(char = c("rJava", "rcdk", "webchem", "KEGGREST", "ChemmineR", "Rdisop"))

install_if_missing("igraph")
install.packages("rJava", repos = cran, type = "source")
install.packages("latticeExtra", repos = cran, type = "source")

cat("== GitHub installs (non-MetaboAnalystR) ==\n")
# NOTE: keeping the same repos/refs you used
devtools::install_github("joannawolthuis/ggVennDiagram")
devtools::install_github("dengkuistat/WaveICA")

# PINS exactly like your Dockerfile
cat("== Pin rlang to 1.1.0 (as in Dockerfile) ==\n")
if (requireNamespace("rlang", quietly = TRUE)) remove.packages("rlang")
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/rlang/rlang_1.1.7.tar.gz",clean = T,
  repos = NULL, type = "source"
)

cat("== Pin httr2 to 0.2.3 (as in Dockerfile) ==\n")
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/httr2/httr2_0.2.3.tar.gz",
  repos = NULL, type = "source"
)

devtools::install_github("lvaudor/glitter", ref = "674418b")
devtools::install_github("joannawolthuis/MetaDBparse")

pacman::p_load(ggpp, pathview, ggplot2)

remotes::install_github("deepanshu88/shinyDarkmode")

pacman::p_load(kohonen, ucminf, mclust, modeltools, scoringRules, ordinal,
               kernlab, diptest, prabclus, flexmix, crch, fpc)

devtools::install_github("joannawolthuis/MetaboShiny", ref = "dev", upgrade = "always")

  # metadbparse fixes
devtools::install_github("Bioconductor/KEGGREST", ref = "devel")
remotes::install_github("WMBEdmands/CompMS2miner")

# misc pins (same as Dockerfile)
install.packages(
  "https://www.bioconductor.org/packages/3.12/bioc/src/contrib/limma_3.46.0.tar.gz",
  repos = NULL, type = "source"
)
install.packages(
  "https://www.bioconductor.org/packages/3.12/bioc/src/contrib/qvalue_2.22.0.tar.gz",
  repos = NULL, type = "source"
)
install.packages(
  "https://bioconductor.org/packages/3.12/bioc/src/contrib/SSPA_2.30.0.tar.gz",
  repos = NULL, type = "source"
)
install.packages("mvtnorm", INSTALL_opts="--no-multiarch")
install.packages("mutoss",  INSTALL_opts="--no-multiarch")
install.packages("metap",   INSTALL_opts="--no-multiarch")

options(pkgType = "win.binary")
options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages("Rserve")

options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2022-06-11"))
options(pkgType = "source")
# --------------------------
# MetaboAnalystR MUST BE LAST
# --------------------------
cat("== MetaboAnalystR LAST ==\n")
# Your Dockerfile installs MetaboAnalystR earlier, but you asked to put it late.
# Also: it needs rlang >= 1.1.1 now, so we temporarily bump rlang right before it.

cat("Bumping rlang for MetaboAnalystR...\n")
if (requireNamespace("rlang", quietly = TRUE)) remove.packages("rlang")
install.packages("rlang", repos = cran)

devtools::install_github("xia-lab/MetaboAnalystR", ref = "0d61192")

cat("\n??? Done.\n")
