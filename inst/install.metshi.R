# METADBPARSE REQUIREMENTS
pacman::p_load(pacman, rcdk, rJava, parallel, pbapply, enviPat, data.table,
               RSQLite, DBI, gsubfn, utils, RCurl, XML, base, 
               stringr, WikidataQueryServiceR, webchem, openxlsx, jsonlite,
               R.utils, KEGGREST, zip, ChemmineR, rvest, xml2, stringi, reshape2,
               Hmisc, httr, RJSONIO, readxl, cmmr, progress, Rdisop, rlist)

# METABOSHINY REQUIREMENTS
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

pacman::p_load(ctc, gdata, glasso, huge, ppcor,crmn)

pacman::p_load(BiocParallel, IRanges, plyr, preprocessCore, vsn,
                         grid, stats4, affy, impute, pcaMethods, MALDIquant, mzID,
                         digest, lattice, ggplot2, XML, scales, MASS, Rcpp)

pacman::p_load(multtest, siggenes, KEGGgraph, SSPA, preprocessCore,
                         Rgraphviz, GlobalAncova, globaltest, pcaMethods, impute, pmp)

pacman::p_load(MSnbase)
pacman::p_load(fgsea)
BiocManager::install("RBGL")
pacman::p_load(crmn)

# 0d61192
devtools::install_github("xia-lab/MetaboAnalystR","0d61192")
#devtools::install_github("yixuan/showtext")
devtools::install_github("joannawolthuis/ggVennDiagram")
devtools::install_github("dengkuistat/WaveICA")

# covbat
devtools::install_github("andy1764/CovBat_Harmonization/R")

install.packages("rJava", type="source", configure.args=c("--disable-jri"))

devtools::install_github("joannawolthuis/MetaDBparse")
devtools::install_github("joannawolthuis/MetaboShiny", "dev")


MetaboShiny::start_metshi(inBrowser=T)
