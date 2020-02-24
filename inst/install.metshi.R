install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
install.packages("pacman")
library(pacman)

pacman::p_load(Rserve, ellipse, scatterplot3d, Cairo, randomForest, caTools, e1071, som, impute, 
               caMethods, RJSONIO, ROCR, globaltest, GlobalAncova, Rgraphviz, preprocessCore, 
               genefilter, pheatmap, SSPA, sva, Rcpp, pROC, data.table, limma, car, fitdistrplus, 
               lars, Hmisc, magrittr, methods, xtable, pls, caret, lattice, igraph, gplots, 
               KEGGgraph, reshape, RColorBrewer, tibble, siggenes, plotly, xcms, 
               CAMERA, fgsea, MSnbase, BiocParallel, metap, reshape2, scales,
               mice, shinydashboard, shinyFiles, shinyBS, 
               rhandsontable, colorRamps, shinyalert, shinyWidgets, colourpicker, 
               ECharts2Shiny, shinyjqui, shinycssloaders, qdapDictionaries, wordcloud2,
               rmarkdown, minval, sqldf, heatmaply, xlsx, manhattanly,
               rgl, glmnet, TSPred, VennDiagram, missForest, InterpretMSSpectrum,
               tm, RISmed, qdap, extrafont, gmp, shadowtext, Rmisc, BioMedR,
               ggbeeswarm, plotROC, ROSE, DMwR, tidytext, TSP, qap, gclus, 
               registry, miniUI, webshot, fracdiff, lmtest, tseries, urca, 
               quadprog, quantmod, xts, TTR, fields, locfit, spam, maps,
               dotCall64, strucchange, chron, NLP, openNLPdata, slam, xlsxjars,
               SnowballC, ISOcodes, broom, shinyjs, tinytex, seriation, egg, 
               manipulateWidget, shape, forecast, KFAS, MuMIn, EMD, wavelets, 
               vars, qdapRegex, qdapTools, gender, openNLP, reports, stringdist, 
               venneuler, wordcloud, extrafontdb, Rttf2pt1, beeswarm,
               vipor, gridSVG, hunspell, tokenizers, janeaustenr, stopwords,
               R.utils, WikidataQueryServiceR, 
               gsubfn, pbapply, rJava, rcdk, enviPat, 
               webchem, rpubchem, rlist, rvest, rapport, SPARQL,
               R.methodsS3, Matrix, selectr, pander, R.oo, cluster, proto,
               rcdklibs, fingerprint, itertools, rapportools,
               sysfonts, showtextdb,ChemometricsWithR, fpc, crch, modeltools,
               mclust, ucminf, RcppArmadillo, kohonen, flexmix, prabclus, diptest, 
               kernlab, ordinal, scoringRules,globaltest, pcaMethods,impute,
               GlobalAncova, Rgraphviz, preprocessCore, genefilter, SSPA, sva,
               limma,KEGGgraph,siggenes,BiocParallel,multtest,spls,ChemometricsWithR, 
               fpc, crch, modeltools)

devtools::install_github("xia-lab/MetaboAnalystR")
devtools::install_github("rwehrens/BatchCorrMetabolomics")
devtools::install_github("yixuan/showtext")
BiocManager::install("ChemmineR")
BiocManager::install("KEGGREST")
devtools::install_github("UMCUGenetics/MetaDBparse")
devtools::install_github("UMCUGenetics/MetaboShiny")
library(MetaboShiny)
start.metshi(inBrowser=T)