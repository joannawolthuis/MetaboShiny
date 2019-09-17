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
    libsqlite3-dev \
    libgtk2.0-0 \
    gconf-gsettings-backend \
    xvfb \
    fuse \
    desktop-file-utils \
    gnupg

RUN sudo R CMD javareconf

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
RUN R -e 'devtools::install_github("UMCUGenetics/MetaDBparse")'
RUN R -e 'devtools::install_github("UMCUGenetics/MetaboShiny")'

# Make the ShinyApp available at port 8080
CMD ['R -e "MetaboShiny::start.metshi(inBrowser=F, port=8080, runmode='docker')"']
EXPOSE 8080

