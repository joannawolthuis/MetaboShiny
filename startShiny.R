
# --- install base base packages ---

install.if.not <- function(package){
  if(package %in% rownames(installed.packages())){
    print(paste("Already installed base package", package))
  }else{
    install.packages(package)
  }
}

# ----------------------------------

install.if.not("pacman")

base.packs <- c("pacman", "shiny", "DT", "data.table", "shinyFiles", "plotly")

p_load(char=base.packs)

# bioconductor 

source("https://bioconductor.org/biocLite.R")
biocLite(suppressUpdates = T)

#wdir <<- "C:/Users/joby/Software/MetaboShiny" #windows
#wdir <<- "/Users/jwolthuis/Google Drive/MetaboShiny/" #mac

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #rstudio only...

# ---------------------------------

runApp(".")
