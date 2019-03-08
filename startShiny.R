# - - helloo - - 
"
This script is currently used to start MetaboShiny.
It takes care of installing packages necessary for 
MetaboShiny to even start.
"


#' Function to install packages, either through regular method or through downloading from git directly
#' @param package package name to install, either CRAN or bioconductor
install.if.not <- function(package){
  if(package %in% rownames(installed.packages())){
    print(paste("Already installed base package", package))
  }else{
    if(package %in% c("MetaboAnalystR", "BatchCorrMetabolomics", "showtext")){
      metanr_packages() # Installs MetaboAnalyst-specific packages
      # Step 1: Install devtools
      install.if.not("devtools")
      # Step 2: Install MetaboAnalystR with documentation
      gitfolder <- switch(package, MetaboAnalystR = "xia-lab/MetaboAnalystR",
                          BatchCorrMetabolomics = "rwehrens/BatchCorrMetabolomics",
                          showtext = "yixuan/showtext")
      devtools::install_github(gitfolder)#, build_vignettes=TRUE)
    }else{
      install.packages(package)
    }
  }
}

# installs packages that metaboanalyst needs to run
metanr_packages <- function(){
  # - - - - - - - - - - - - - -
  metr_pkgs <- c("Rserve", "RColorBrewer",
                 "xtable", "som", "ROCR", "RJSONIO", "gplots",
                 "e1071", "caTools", "igraph", "randomForest", "Cairo",
                 "pls", "pheatmap", "lattice", "rmarkdown", "knitr",
                 "data.table", "pROC", "Rcpp", "caret", "ellipse",
                 "scatterplot3d", "impute", "rhandsontable", "pcaMethods", 
                 "siggenes", "globaltest", "GlobalAncova", "Rgraphviz", 
                 "KEGGgraph", "preprocessCore", "genefilter", "SSPA",
                 "sva", "showtext", "wordcloud2")
  list_installed <- installed.packages()
  
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  
  if(length(new_pkgs) != 0){
    source("https://bioconductor.org/biocLite.R")
    BiocInstaller::biocLite(new_pkgs, dependencies = TRUE, ask = FALSE)
    print(paste0(new_pkgs, " added..."))
  }
  if((length(new_pkgs) < 1)){
    print("No new packages added...")
  }
}

# -- IF YOU HAVE PROBLEMS STARTING UP --
# requires openssl, libcurl, libxt, 
# libxml2, libnetcdf, cairo, java, mesa-common, 
# libgit2, (LIKELY -devel versions)
# libgl1, libglu1, libpng installed
# IF RJAVA DOESN'T WANT TO INSTALL
# sudo R CMD javareconf # FIXES EVERYTHING!!! JUST NEED TO USE ADMIN MODE...
# install.packages("rJava", type="source")

# packages needed to start up
base.packs <<- c('httr', 'curl', 'git2r', 'devtools', 
                'pacman', 'gsubfn', 'shiny', 'DT', 
                'R.utils', 'data.table', 'shinyFiles', 
                'shinyBS', 'rhandsontable', 'XML', 
                'MetaboAnalystR', 'BatchCorrMetabolomics', 
                'colorRamps', 'enviPat', 'shinyalert',
                'shinyWidgets', 'colourpicker', 'here',
                'ECharts2Shiny', 'shinyjqui', 'later',
                'shinycssloaders', 'qdapDictionaries',
                'sysfonts', 'showtext', 'wordcloud2')

# install the base packages needed to start up
for(package in base.packs){
  install.if.not(package)
}

# set working directory to where the rstudio file currently is
wdir <<- dirname(rstudioapi::getSourceEditorContext()$path) # TODO: make this not break when not running from rstudio
setwd(wdir)

# create options file if it doesnt exist yet
if(!file.exists(file.path(wdir, "user_options.txt"))){
  default_options = 
gsubfn::fn$paste(
  "db_dir = $wdir/backend/db
  work_dir = $wdir/analysis
  proj_name = MY_PROJECT
  ppm = 2
  packages_installed = Y
  font1 = Supermercado One
  font2 = Supermercado One
  font3 = Open Sans
  font4 = Open Sans
  col1 = #000000
  col2 = #DBDBDB
  col3 = #FFFFFF
  col4 = #FFFFFF
  size1 = 50
  size2 = 20
  size3 = 15
  size4 = 12
  taskbar_image = gemmy_rainbow.png
  gtheme = classic
  gcols = #000000&#FFA1C3&#FFC914&#2E282A&#8A00ED&#00E0C2&#95C200&#FF6BE4
  gspec = RdGy")
con = file("./user_options.txt", "w")
writeLines(text = default_options, con = con)
close.connection(con)
}

# go run it! :-)
shiny::runApp(".")#, launch.browser = T)
