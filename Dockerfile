# Install R version 3.6
FROM r-base:4.0.2

RUN pwd

# Install Ubuntu packages
RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    xclip \
    pandoc \
    openjdk-11-jdk \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev/unstable \
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
    gconf-gsettings-backend \
    xvfb \
    fuse \
    desktop-file-utils \
    libgdal-dev \
    gnupg

RUN R -e 'install.packages("devtools")'
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install("org.Hs.eg.db")'
RUN R -e 'install.packages("pacman")'
RUN R -e 'library(pacman)'
RUN R -e 'pacman::p_load(Rserve, RSclient, ellipse, scatterplot3d, Cairo, randomForest, caTools, e1071, som, impute, pcaMethods, RJSONIO, ROCR, globaltest, GlobalAncova, Rgraphviz, preprocessCore, genefilter, pheatmap, SSPA, sva, Rcpp, pROC, data.table, limma, car, fitdistrplus, lars, Hmisc, magrittr, methods, xtable, pls, caret, lattice, igraph, gplots, KEGGgraph, reshape, RColorBrewer, tibble, siggenes, plotly, fgsea, metap, reshape2, scales, CAMERA)'
RUN R -e 'pacman::p_load(ctc, gdata, glasso, huge, ppcor,crmn)'
RUN R -e 'devtools::install_github("xia-lab/MetaboAnalystR")'
RUN R -e 'devtools::install_github("yixuan/showtext")'
RUN R -e 'devtools::install_github("joannawolthuis/ggVennDiagram")'
RUN R -e 'devtools::install_github("dengkuistat/WaveICA")'
RUN R -e 'install.packages("rJava", type="source")'
RUN R -e 'pacman::p_load(pacman, rcdk, rJava, parallel, pbapply, enviPat, data.table, \
               RSQLite, DBI, gsubfn, utils, RCurl, XML, base, \
               stringr, WikidataQueryServiceR, webchem, openxlsx, jsonlite,\
               R.utils, KEGGREST, zip, ChemmineR, rvest, xml2, stringi, reshape2, \
               Hmisc, httr, RJSONIO, readxl, cmmr, progress, Rdisop, rlist, SPARQL)'
 
RUN sudo R CMD javareconf

RUN R -e 'devtools::install_github("UMCUGenetics/MetaDBparse")'
RUN R -e 'pacman::p_load(ggplot2, data.table, plotly, shinyBS, shinyjs, caret, grDevices, \
               RColorBrewer, colorRamps, tidytext, qdapDictionaries, tm, shiny, htmltools, \
               BiocManager, pacman, devtools, classyfireR, httr, jsonlite, RCurl, shinyFiles, \
               DT, RSQLite, pbapply, stringr, gsubfn, shinyWidgets, parallel, \
               mice, sva, limma, tools, plyr, heatmaply, wordcloud2, shinyjqui, rmarkdown, enviPat, ROCR, tsne, e1071,\
               pls, rhandsontable, testthat, shinytest, showtext, sysfonts, colourpicker, \
               reshape, ggdark, ECharts2Shiny, shinyalert, shinybusy, rcdk, RISmed, dplyr, \
               InterpretMSSpectrum, DBI, qdap, reshape2, Hmisc, ggbeeswarm, Rmisc, rgl,\
               stats, pROC, car, doParallel, missForest, ggfortify, fdrtool, plsdepot, vroom, qs, WaveICA)'
RUN R -e 'devtools::install_github("UMCUGenetics/MetaboShiny", "dev")'

# Make the ShinyApp available at port 8080
CMD ['R -e "MetaboShiny::start_metshi(inBrowser=F, port=8080, runmode='docker')"']
EXPOSE 8080

