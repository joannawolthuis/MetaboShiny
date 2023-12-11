# Install R version 4.0.3
#from ubuntu:22.04
#ENV R_BASE_VERSION=4.0.3
FROM rocker/tidyverse:4.0.2

# Install Ubuntu packages
RUN apt-get update
#RUN apt-get update && apt-get install -y \
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
    sudo \
    gdebi-core \
    xclip \
    pandoc \
    openjdk-11-jdk \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libxml2-dev \
    libgmp-dev \
    libnetcdf-dev \
    mesa-common-dev \
    apt-utils \
    libgit2-dev \
    libglu1-mesa \
    libglu1-mesa-dev \
    libpng-dev \
    libudunits2-dev \
    librsvg2-dev \
    ca-certificates \
    curl \
    sqlite3 \
    libsqlite3-dev \
    libgtk2.0-0 \
    dconf-gsettings-backend \
    xvfb \
    fuse \
    desktop-file-utils \
    libgdal-dev \
    gnupg
    
RUN apt update && apt install -y libcurl4-openssl-dev libssl-dev

RUN sudo R CMD javareconf

RUN R -e 'install.packages("devtools")'
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install("org.Hs.eg.db")'
RUN R -e 'install.packages("pacman")'
RUN R -e 'library(pacman)'
RUN R -e 'install.packages("mvtnorm", type="source")'

# METADBPARSE REQUIREMENTS
RUN R -e 'pacman::p_load(pacman, rcdk, rJava, parallel, pbapply, enviPat, data.table,\
               RSQLite, DBI, gsubfn, utils, RCurl, XML, base, \
               stringr, WikidataQueryServiceR, webchem, openxlsx, jsonlite,\
               R.utils, KEGGREST, zip, ChemmineR, rvest, xml2, stringi, reshape2,\
               Hmisc, httr, RJSONIO, readxl, cmmr, progress, Rdisop, rlist)'

# METABOSHINY REQUIREMENTS
RUN R -e 'pacman::p_load(ggplot2, data.table, plotly, shinyBS, shinyjs, caret, grDevices,\
               RColorBrewer, colorRamps, tidytext, qdapDictionaries, tm, shiny, htmltools,\
               BiocManager, pacman, devtools, classyfireR, httr, jsonlite, RCurl, shinyFiles,\
               DT, RSQLite, pbapply, stringr, gsubfn, shinyWidgets, parallel,\
               mice, sva, limma, tools, plyr, heatmaply, wordcloud2, shinyjqui, rmarkdown,\
               enviPat, ROCR, tsne, e1071)'
               
RUN R -e 'pacman::p_load(pls, rhandsontable, testthat, shinytest, showtext, sysfonts, colourpicker,\
               reshape, ggdark, ECharts2Shiny, shinyalert, shinybusy, rcdk, RISmed, dplyr,\
               InterpretMSSpectrum, DBI, qdap, reshape2, Hmisc, ggbeeswarm, Rmisc, rgl,\
               stats, pROC, car, doParallel, missForest, ggfortify, fdrtool, plsdepot,\
               vroom, umap, ica, svglite, beepr, showtextdb)'

RUN R -e 'pacman::p_load(ctc, gdata, glasso, huge, ppcor,crmn)'

RUN R -e 'pacman::p_load(BiocParallel, IRanges, plyr, preprocessCore, vsn,\
                         grid, stats4, affy, impute, pcaMethods, MALDIquant, mzID,\
                         digest, lattice, ggplot2, XML, scales, MASS, Rcpp)'
                         
RUN R -e 'pacman::p_load(multtest, siggenes, KEGGgraph, SSPA, preprocessCore,\
                         Rgraphviz, GlobalAncova, globaltest, pcaMethods, impute)'
                        
RUN R -e 'pacman::p_load(pdftools, magick)'

RUN R -e 'pacman::p_load(MSnbase)'
RUN R -e 'pacman::p_load(fgsea)'
RUN R -e 'BiocManager::install("RBGL")'
RUN R -e 'pacman::p_load(crmn)'

RUN R -e 'install.packages("devtools")'
RUN R -e 'remove.packages("mvtnorm")' 
RUN R -e 'install.packages("mvtnorm")'
RUN R -e 'install.packages("mutoss", type="source")'

RUN R -e 'pacman::p_load(char=c("genefilter", "sva", "limma", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea","httr","gdata", "glasso", "huge","robustbase","qqconf","ppcor","crmn","plotly"))'

RUN R -e 'pacman::p_load(char = c("anytime", "rJava", "ggraph", "tidygraph", "WikidataR", "xmlparsedata"))'
RUN R -e 'pacman::p_load(char=c("rJava","rcdk", "webchem","KEGGREST", "ChemmineR", "Rdisop"))'

RUN R -e 'install.packages("igraph")' 

RUN R -e 'install.packages("rJava", type="source")'
RUN R -e 'install.packages("latticeExtra", type="source")'

RUN DEBIAN_FRONTEND=noninteractive apt-get install -y libglpk-dev

RUN R -e 'devtools::install_github("xia-lab/MetaboAnalystR", "0d61192")'
RUN R -e 'devtools::install_github("yixuan/showtext")'
RUN R -e 'devtools::install_github("joannawolthuis/ggVennDiagram")'
RUN R -e 'devtools::install_github("dengkuistat/WaveICA")'

RUN R -e 'remove.packages("rlang")'
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/rlang/rlang_1.1.0.tar.gz", repo=NULL, type="source")'

RUN R -e 'devtools::install_github("r-lib/httr2")'
RUN R -e 'devtools::install_github("lvaudor/glitter", "674418b")'

RUN R -e 'devtools::install_github("joannawolthuis/MetaDBparse")'

RUN R -e 'pacman::p_load(ggpp, pathview, ggplot2)'
RUN R -e 'devtools::install_github("aphalo/ggpp")'

RUN R -e 'devtools::install_github("joannawolthuis/MetaboShiny", "dev", upgrade="always")'

# Make the ShinyApp available at port 8080
EXPOSE 8080

