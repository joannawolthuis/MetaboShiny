# Install R version 3.6
FROM r-base:3.6.0

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

RUN R -e 'install.packages("devtools")'
RUN R -e 'devtools::install_github("UMCUGenetics/MetaboShiny")'

# Make the ShinyApp available at port 8080
CMD ['R -e "MetaboShiny::start.metshi(inBrowser=F, port=8080)"']
EXPOSE 8080

