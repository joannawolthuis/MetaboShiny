# Install R version 3.5
FROM r-base:3.5.2

RUN pwd

# Install Ubuntu packages
RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
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
    librsvg2-dev \
    default-jre \
    default-jdk \
    ca-certificates \
    curl \
    sqlite3 \
    libsqlite3-dev

RUN sudo R CMD javareconf

# Copy files into the Docker image
COPY  . /

RUN R -e "install.packages('BiocManager')"

RUN R -e "BiocManager::install(c('rJava', 'shiny', 'shinydashboard', 'httr', 'curl', 'git2r', 'devtools', 'pacman', 'gsubfn', 'DT', 'R.utils'))"
RUN R -e "BiocManager::install(c('data.table', 'shinyFiles', 'shinyBS', 'rhandsontable', 'XML', 'colorRamps', 'enviPat', 'shinyalert', 'shinyWidgets', 'colourpicker'))"
RUN R -e "BiocManager::install(c('ECharts2Shiny', 'shinyjqui', 'later', 'shinycssloaders', 'qdapDictionaries', 'sysfonts', 'showtext', 'wordcloud2', 'Rserve'))"
RUN R -e "BiocManager::install(c('RColorBrewer', 'xtable', 'som', 'ROCR', 'RJSONIO', 'gplots', 'e1071', 'caTools', 'igraph'))"
RUN R -e "BiocManager::install(c('Cairo', 'pls', 'lattice', 'rmarkdown', 'knitr', 'pROC', 'Rcpp', 'caret', 'ellipse', 'scatterplot3d'))"
RUN R -e "BiocManager::install(c('impute', 'pcaMethods', 'siggenes', 'globaltest', 'GlobalAncova', 'Rgraphviz', 'KEGGgraph', 'preprocessCore', 'genefilter', 'SSPA'))"
RUN R -e "BiocManager::install(c('sva', 'DBI', 'RSQLite', 'ggplot2', 'minval', 'plotly', 'pbapply', 'sqldf', 'plyr', 'ChemmineR'))"
RUN R -e "BiocManager::install(c('stringr', 'heatmaply', 'reshape2', 'xlsx', 'pheatmap', 'rJava', 'KEGGREST', 'manhattanly', 'rgl', 'glmnet'))"
RUN R -e "BiocManager::install(c('TSPred', 'VennDiagram', 'rcdk', 'SPARQL', 'webchem', 'WikidataQueryServiceR', 'openxlsx', 'doParallel', 'missForest', 'InterpretMSSpectrum'))"
RUN R -e "BiocManager::install(c('tm', 'RISmed', 'qdap', 'extrafont', 'gmp', 'shadowtext', 'xcms', 'CAMERA', 'fgsea', 'MSnbase', 'Rmisc', 'ggbeeswarm'))"

# Make the ShinyApp available at port 8080
EXPOSE 8080

#CMD ./start.sh
