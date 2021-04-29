install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
install.packages("pacman")
library(pacman)

# METADBPARSE REQUIREMENTS
pacman::p_load(pacman, rcdk, rJava, parallel, pbapply, enviPat, data.table,
               RSQLite, DBI, gsubfn, utils, RCurl, XML, base, 
               stringr, WikidataQueryServiceR, webchem, openxlsx, jsonlite,
               R.utils, KEGGREST, zip, ChemmineR, rvest, xml2, stringi, reshape2,
               Hmisc, httr, RJSONIO, readxl, cmmr, progress, Rdisop, rlist,
               # metaboshiny
               ggplot2, data.table, plotly, shinyBS, shinyjs, caret, grDevices,
               RColorBrewer, colorRamps, tidytext, qdapDictionaries, tm, shiny, htmltools,
               BiocManager, pacman, devtools, classyfireR, httr, jsonlite, RCurl, shinyFiles,
               DT, RSQLite, pbapply, stringr, gsubfn, shinyWidgets, parallel,
               mice, sva, limma, tools, plyr, heatmaply, wordcloud2, shinyjqui, rmarkdown,
               enviPat, ROCR, tsne, e1071, genefilter,
               pls, rhandsontable, testthat, shinytest, showtext, sysfonts, colourpicker,
               reshape, ggdark, ECharts2Shiny, shinyalert, shinybusy, rcdk, RISmed, dplyr,
               InterpretMSSpectrum, DBI, qdap, reshape2, Hmisc, ggbeeswarm, Rmisc, rgl,
               stats, pROC, car, doParallel, missForest, ggfortify, fdrtool, plsdepot,
               vroom, umap, ica, svglite, beepr, showtextdb,
               ctc, gdata, glasso, huge, ppcor,crmn,
               BiocParallel, IRanges, plyr, preprocessCore, vsn,
               grid, stats4, affy, impute, pcaMethods, MALDIquant, mzID,
               digest, lattice, ggplot2, XML, scales, MASS, Rcpp,
               multtest, siggenes, KEGGgraph, SSPA, preprocessCore,
               Rgraphviz, GlobalAncova, globaltest, pcaMethods, impute,
               multiROC, ggpubr,
               MSnbase,
               fgsea,
               crmn,
               smotefamily,
               ROSE,
               DGCA,
               heatmap3
               )

BiocManager::install("RBGL")
devtools::install_github("xia-lab/MetaboAnalystR", "0d61192")
devtools::install_github("yixuan/showtext")
devtools::install_github("gaospecial/ggVennDiagram")
devtools::install_github("joannawolthuis/WaveICA")

install.packages("rJava", type="source")

devtools::install_github("UMCUGenetics/MetaDBparse")
devtools::install_github("UMCUGenetics/MetaboShiny", "dev")

