install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
install.packages("pacman")

pacman::p_load(pacman, rcdk, rJava, parallel, pbapply, enviPat, data.table, 
               #MetaDBparse, 
               RSQLite, DBI, gsubfn, utils, RCurl, XML, base, 
               stringr, WikidataQueryServiceR, webchem, openxlsx, jsonlite,
               R.utils, KEGGREST, zip, ChemmineR, rvest, xml2, stringi, reshape2, 
               Hmisc, httr, RJSONIO, readxl, cmmr, progress, Rdisop, rlist)

# METABOSHINY REQUIREMENTSa
pacman::p_load(ggplot2, data.table, plotly, shinyBS, shinyjs, caret, grDevices, 
               RColorBrewer, colorRamps, tidytext, qdapDictionaries, tm, shiny, htmltools, 
               BiocManager, pacman, devtools, classyfireR, httr, jsonlite, RCurl, shinyFiles, 
               DT, RSQLite, pbapply, stringr, gsubfn, shinyWidgets, parallel, 
               mice, sva, limma, tools, plyr, heatmaply, wordcloud2, shinyjqui, rmarkdown, 
               enviPat, ROCR, tsne, e1071)

pacman::p_load(pls, rhandsontable, testthat, shinytest, showtext, sysfonts, colourpicker, 
               reshape, ggdark, ECharts2Shiny, shinyalert, shinybusy, rcdk, RISmed, dplyr, 
               InterpretMSSpectrum, DBI, qdap, reshape2, Hmisc, ggbeeswarm, Rmisc, rgl,
               stats, pROC, car, doParallel, missForest, ggfortify, fdrtool, plsdepot,
               vroom, umap, ica, svglite, beepr, showtextdb)

pacman::p_load(BiocParallel, IRanges, plyr, preprocessCore, vsn, 
               grid, stats4, affy, impute, pcaMethods, MALDIquant, mzID, 
               digest, lattice, ggplot2, XML, scales, MASS, Rcpp)

pacman::p_load(ctc, gdata, glasso, huge, ppcor,crmn, MSnbase)

devtools::install_github("joannawolthuis/ggVennDiagram")
devtools::install_github("dengkuistat/WaveICA")
devtools::install_github("joannawolthuis/MetaDBparse")
devtools::install_github("joannawolthuis/MetaboShiny", "dev")
MetaboShiny::start_metshi(inBrowser=T)
