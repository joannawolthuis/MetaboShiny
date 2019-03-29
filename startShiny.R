# - - helloo - -
"
This script is currently used to start MetaboShiny.
It takes care of installing packages necessary for
MetaboShiny to even start.
"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function to install packages, either through regular method or through downloading from git directly
#' @param package package name to install, either CRAN or bioconductor
install.if.not <- function(package){
  if(package %in% rownames(installed.packages())){
    print(paste("Already installed base package", package))
  }else{
    if(package %in% c("MetaboAnalystR", "BatchCorrMetabolomics")){
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

# Install R packages that are required
# TODO: add further package if you need!

needed.packages <- c("BiocManager", "shiny", "shinydashboard", "httr", "curl", "git2r", "devtools",
                     "pacman", "gsubfn", "DT", "R.utils", "data.table", "shinyFiles",
                     "shinyBS", "rhandsontable", "XML", "colorRamps", "enviPat", "shinyalert",
                     "shinyWidgets", "colourpicker", "here", "ECharts2Shiny", "shinyjqui",
                     "later", "shinycssloaders", "qdapDictionaries", "sysfonts", "showtext",
                     "wordcloud2", "Rserve", "RColorBrewer", "xtable", "som", "ROCR",
                     "RJSONIO", "gplots", "e1071", "caTools", "igraph", "randomForest",
                     "Cairo", "pls", "lattice", "rmarkdown", "knitr", "pROC", "Rcpp",
                     "caret", "ellipse", "scatterplot3d", "impute", "pcaMethods",
                     "siggenes", "globaltest", "GlobalAncova", "Rgraphviz", "KEGGgraph",
                     "preprocessCore", "genefilter", "SSPA", "sva", "DBI", "RSQLite",
                     "ggplot2", "minval", "plotly", "pbapply", "sqldf", "plyr", "ChemmineR",
                     "stringr", "heatmaply", "reshape2", "xlsx", "pheatmap", "rJava",
                     "KEGGREST", "manhattanly", "rgl", "glmnet", "TSPred", "VennDiagram",
                     "rcdk", "SPARQL", "webchem", "WikidataQueryServiceR", "openxlsx",
                     "doParallel", "missForest", "InterpretMSSpectrum", "tm", "RISmed",
                     "qdap", "extrafont", "gmp", "shadowtext", "fgsea")

missing.packages <- setdiff(needed.packages,rownames(installed.packages()))

if("BiocManager" %in% missing.packages){
  install.packages('BiocManager')
}

print(missing.packages)

if(length(missing.packages)>0){
  BiocManager::install(missing.packages)
}


## used to create the docker run statements
# groupies <- split(unique(to.install), ceiling(seq_along(unique(to.install))/10))
# for(group in groupies){
#   packages <- paste0("'",paste0(group, collapse="', '"), "'")
#   cat(gsubfn::fn$paste('RUN R -e "BiocManager::install(c($packages))"\n'))
# }
## docker internet
# docker run --net=host -it metaboshiny/test1
## docker with volume mounted
# docker run -p 8080:8080 --net=host -it --mount src='~/MetaboShiny',target='/databases/,type=bind metaboshiny/test1

# attempt 1: WORKS FOR NOW
# had to saet memory to max lol

# make metaboshiny_storage dir in home first..
# docker run -p 8080:8080 --mount src=~/MetaboShiny_storage,target=/userfiles/,type=bind --rm -it metaboshiny/test1 sh

# packages needed to start up
git.packages <<- c("MetaboAnalystR", 
                   "BatchCorrMetabolomics")

# install the base packages needed to start up
for(package in git.packages){
  install.if.not(package)
}

library(httr)

# does this work?
# https://stackoverflow.com/questions/3160909/how-do-i-deal-with-certificates-using-curl-while-trying-to-access-an-https-url
# nope..
#options('download.file.method' = 'curl')

options('unzip.unzip' = getOption("unzip"), 'download.file.extra' = "--insecure", 'download.file.method' = 'curl')

runmode <- if(file.exists(".dockerenv")) 'docker' else 'local'

opt.loc <<- if(runmode == 'local') '~/Documents/MetaboShiny/user_options_local.txt' else '/userfiles/user_options_docker.txt'

optfolder <- dirname(opt.loc)

if(!dir.exists(optfolder)) dir.create(optfolder,recursive = T)

if(!file.exists(opt.loc)){
  # write options file if it doesn't exist yet
  contents <- switch(runmode,
         local = {
'db_dir = ~/Documents/MetaboShiny/databases
work_dir = ~/Documents/MetaboShiny/userfiles
proj_name = MY_METSHI
ppm = 2
packages_installed = Y
font1 = Pacifico
font2 = Pacifico
font3 = Open Sans
font4 = Open Sans
col1 = #000000
col2 = #DBDBDB
col3 = #FFFFFF
col4 = #FFFFFF
size1 = 50
size2 = 20
size3 = 15
size4 = 11
taskbar_image = gemmy_rainbow.png
gtheme = classic
gcols = #FF0004&#38A9FF&#FFC914&#2E282A&#8A00ED&#00E0C2&#95C200&#FF6BE4
gspec = RdBu'}, docker = {
'db_dir = /userfiles/databases
work_dir = /userfiles/userfiles
proj_name = MY_METSHI
ppm = 2
packages_installed = Y
font1 = Pacifico
font2 = Pacifico
font3 = Open Sans
font4 = Open Sans
col1 = #000000
col2 = #DBDBDB
col3 = #FFFFFF
col4 = #FFFFFF
size1 = 40
size2 = 20
size3 = 15
size4 = 11
taskbar_image = gemmy_rainbow.png
gtheme = classic
gcols = #FF0004&#38A9FF&#FFC914&#2E282A&#8A00ED&#00E0C2&#95C200&#FF6BE4
gspec = RdBu'})
  writeLines(contents, opt.loc)
}


switch(runmode,
       local = {
         wdir <<- dirname(rstudioapi::getSourceEditorContext()$path) # TODO: make this not break when not running from rstudio
         setwd(wdir)
         shiny::runApp(".")
       }, 
       docker = {
         shiny::runApp(".",
                       port = 8080,
                       host = "0.0.0.0",
                       launch.browser = FALSE)
       })
# go run it! :-)

