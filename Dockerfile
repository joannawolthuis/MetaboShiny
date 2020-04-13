# Install R version 3.6
FROM r-base:3.6.0

RUN pwd

# Install Ubuntu packages
RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
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

# Download orca AppImage, extract it, and make it executable under xvfb
RUN wget https://github.com/plotly/orca/releases/download/v1.1.1/orca-1.1.1-x86_64.AppImage -P /home
RUN chmod 777 /home/orca-1.1.1-x86_64.AppImage 

# To avoid the need for FUSE, extract the AppImage into a directory (name squashfs-root by default)
RUN cd /home && /home/orca-1.1.1-x86_64.AppImage --appimage-extract
RUN printf '#!/bin/bash \nxvfb-run --auto-servernum --server-args "-screen 0 640x480x24" /home/squashfs-root/app/orca "$@"' > /usr/bin/orca
RUN chmod 777 /usr/bin/orca
RUN chmod -R 777 /home/squashfs-root/

RUN R -e 'install.packages("devtools")'
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install("org.Hs.eg.db")'
RUN R -e 'install.packages("pacman")'
RUN R -e 'library(pacman)'
RUN R -e 'pacman::p_load(Rserve, RSclient, ellipse, scatterplot3d, Cairo, randomForest, caTools, e1071, som, impute, pcaMethods, RJSONIO, ROCR, globaltest, GlobalAncova, Rgraphviz, preprocessCore, genefilter, pheatmap, SSPA, sva, Rcpp, pROC, data.table, limma, car, fitdistrplus, lars, Hmisc, magrittr, methods, xtable, pls, caret, lattice, igraph, gplots, KEGGgraph, reshape, RColorBrewer, tibble, siggenes, plotly, fgsea, metap, reshape2, scales, CAMERA)'
RUN R -e 'devtools::install_github("xia-lab/MetaboAnalystR")'
RUN R -e 'devtools::install_github("rwehrens/BatchCorrMetabolomics")'
RUN R -e 'devtools::install_github("yixuan/showtext")'
RUN R -e 'devtools::install_github("joannawolthuis/ggVennDiagram")'
RUN R -e 'install.packages("rJava", type="source")'
RUN R -e 'pacman::p_load(pacman, rcdk, rJava, parallel, pbapply, enviPat, data.table, \
               RSQLite, DBI, gsubfn, utils, RCurl, XML, base, \
               stringr, WikidataQueryServiceR, webchem, openxlsx, jsonlite,\
               R.utils, KEGGREST, zip, ChemmineR, rvest, xml2, stringi, reshape2, \
               Hmisc, httr, RJSONIO, readxl, cmmr, progress, Rdisop, rlist, SPARQL)'
RUN R -e 'devtools::install_github("UMCUGenetics/MetaDBparse")'
RUN R -e 'pacman::p_load(ggplot2, data.table, plotly, shinyBS, shinyjs, caret, grDevices, \
               RColorBrewer, colorRamps, tidytext, qdapDictionaries, tm, shiny, htmltools, \
               BiocManager, pacman, devtools, classyfireR, httr, jsonlite, RCurl, shinyFiles, \
               DT, RSQLite, pbapply, stringr, gsubfn, shinyWidgets, parallel, \
               mice, sva, limma, tools, plyr, heatmaply, wordcloud2, shinyjqui, rmarkdown, enviPat, ROCR, tsne, e1071,\
               pls, rhandsontable, testthat, shinytest, showtext, sysfonts, colourpicker, \
               reshape, ggdark, ECharts2Shiny, shinyalert, shinybusy, rcdk, RISmed, dplyr, \
               InterpretMSSpectrum, DBI, qdap, reshape2, Hmisc, ggbeeswarm, Rmisc, rgl,\
               stats, pROC, car, doParallel, missForest)'
RUN R -e 'devtools::install_github("UMCUGenetics/MetaboShiny", "dev")'
RUN sudo R CMD javareconf

# Make the ShinyApp available at port 8080
CMD ['R -e "MetaboShiny::start.metshi(inBrowser=F, port=8080, runmode='docker')"']
EXPOSE 8080

