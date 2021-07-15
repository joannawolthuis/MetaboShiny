# Install R version 4.0.2
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

RUN sudo R CMD javareconf

RUN R -e 'install.packages("devtools")'
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install("org.Hs.eg.db")'
RUN R -e 'install.packages("pacman")'
RUN R -e 'library(pacman)'

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

                         
RUN R -e 'devtools::install_github("xia-lab/MetaboAnalystR", "0d61192")'
RUN R -e 'devtools::install_github("yixuan/showtext")'
RUN R -e 'devtools::install_github("joannawolthuis/ggVennDiagram")'
RUN R -e 'devtools::install_github("dengkuistat/WaveICA")'

RUN R -e 'install.packages("rJava", type="source")'

RUN R -e 'devtools::install_github("UMCUGenetics/MetaDBparse")'
RUN R -e 'devtools::install_github("UMCUGenetics/MetaboShiny", "dev")'

# Make the ShinyApp available at port 8080
EXPOSE 8080

