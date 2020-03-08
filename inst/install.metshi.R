install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
install.packages("pacman")
library(pacman)

pacman::p_load(rJava, ggplot2, data.table, plotly, caret, shinybusy, svglite,
               grDevices, RColorBrewer, colorRamps, tidytext, qdapDictionaries, 
               tm, shinyBS, shiny, htmltools, BiocManager, pacman, 
               devtools, classyfireR, httr, jsonlite, RCurl, shinyjqui, 
               shinyWidgets, DT, heatmaply, wordcloud2, shinyFiles, 
               RSQLite, pbapply, stringr, gsubfn, shinyjs, 
               parallel, mice, sva, limma, tools, rmarkdown, enviPat, tsne, e1071, 
               pls, ROCR, RISmed, rhandsontable, sysfonts, 
               colourpicker, reshape, ggdark, ECharts2Shiny, shinyalert, 
               shinycssloaders, rcdk, dplyr, InterpretMSSpectrum, DBI, 
               qdap, reshape2, Hmisc, ggbeeswarm, Rmisc, rgl, VennDiagram, 
               shadowtext, stats, pROC, car, doParallel, missForest,
               impute, pcaMethods, globaltest, GlobalAncova, Rgraphviz, 
               preprocessCore, genefilter, SSPA, sva, limma, 
               KEGGgraph, siggenes, BiocParallel, multtest, CAMERA, Rdisop,
               qdap, fgsea, ChemmineR, rpubchem)

devtools::install_github("xia-lab/MetaboAnalystR")
devtools::install_github("rwehrens/BatchCorrMetabolomics")
devtools::install_github("yixuan/showtext")
devtools::install_github("UMCUGenetics/MetaDBparse")
devtools::install_github("UMCUGenetics/MetaboShiny")
devtools::install_github("joannawolthuis/ggVennDiagram")

library(MetaboShiny)
start.metshi(inBrowser=T)

update.packages(ask = F)
