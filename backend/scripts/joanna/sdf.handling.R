# require(parallel)
# require(pbapply)
# require(ChemmineR)
# require(stringr)

# --- install stuff ---

# pkgs <- c(
#   "S4Vectors", "IRanges", "GenomicRanges", "DelayedArray",
#   "XVector", "GenomicAlignments", "ShortRead",
#   "VariantAnnotation", "AnnotationHub", "GGtools",
#   "ggbio", "ChemmineR", "InteractionSet", "flowCore",
#   "GenomicTuples", "CNEr", "MultiAssayExperiment",
#   "genomeIntervals", "TFBSTools", "IWTomics", "spliceSites",
#   "podkat", "kebabs", "matter", "dada2",
#   "ClusterSignificance", "gespeR", "HiTC", "tigre", "soGGi"
#          )
#
# update <- intersect(rownames(installed.packages()), pkgs)
# BiocInstaller::biocLite(update, type="source", ask=FALSE)

# ------------------------------------------
# 
# hpcDir <- "/Users/jwolthuis/MountPoint/Joanna/Databases/PubChem"
# sdfDir <- "/Users/jwolthuis/MountPoint/Joanna/Databases/PubChem/"
# 
# 
# #hpcDir <- "/hpc/shared/dbg_mz/Joanna/Databases/PubChem"
# sdfDir <- "/hpc/shared/dbg_mz/Joanna/Databases/PubChem/xls"
# 
# # --- spider function puts stuff in the xls file ---
# 
# # Calculate the number of cores
# no_cores <- Sys.getenv( "NSLOTS" )
# 
# print(paste(no_cores, "cores found!", sep=" "))
# 
# cl <- makeCluster(3, type="FORK")
# 
# # --- LEGGO! ---
# 
# setwd(hpcDir)
# 
# in.files <- list.files(pattern = "\\.sdf\\.gz$")
# print(in.files)
# print("Building base database... ")
# 
# pbsapply(cl=cl, in.files, FUN=function(input){
#   output <- file.path(sdfDir, gsub("\\.sdf.gz$", ".xls", x = input))
#   if(file.exists(output)) return(NA)
#   # --- check for file permissions ---
#   if(file.access(input, mode=4) == -1){stopCluster(cl);
#     stop(paste("no read permission on input file"));}
#   if(file.access(output, mode=4) == -1)
#     stop(paste("no read permission on output file"))
#   # -----
#   sdfStream(input=input,
#             output=output,
#             fct=function(sdfset, test){
#               blockmatrix <- datablock2ma(datablocklist=datablock(sdfset)) # Converts data block to matrix
#               name.col <- if("PUBCHEM_IUPAC_TRADITIONAL_NAME" %not in% colnames(blockmatrix)) "PUBCHEM_MOLECULAR_FORMULA" else("PUBCHEM_IUPAC_TRADITIONAL_NAME")
#               # --------------
#               barebones.table <- data.frame(
#                 CompoundName = blockmatrix[, name.col],
#                 BaseFormula = blockmatrix[, "PUBCHEM_MOLECULAR_FORMULA"],
#                 Identifier = as.numeric(blockmatrix[, "PUBCHEM_COMPOUND_CID"]),
#                 Source = "PubChem"
#               )
#               # --- return ---
#               barebones.table
#             },
#             append = FALSE,
#             silent = TRUE,
#             Nlines = 100000 )
# })

# ----------------------------------
# no_cores <- Sys.getenv( "NSLOTS" )
# 
# print(paste(no_cores, "cores found!", sep=" "))
# 
# cl <- makeCluster(no_cores, type="FORK")
# 
# sdfDir <- "/Users/jwolthuis/MountPoint/Joanna/Databases/PubChem/"
# dbDir <-
# setwd(sdfDir)
# 
# in.files <- list.files(pattern = "\\.xls$")
# print(in.files)
# 
# pubchem <- pblapply(in.files, FUN=function(file){
#   print(paste("Current file:", file, sep=" "))
#   tab <- fread(file)
#   fread(file)[,c(4:7)]
# })
# 
# pubchem.base <- unique(rbindlist(pubchem))
# 
# save(file=file.path(sdfDir, "PubChem_base.RData")

# # ---------------------------------
# 
# dbDir <- "/Users/jwolthuis/Downloads/PubChem"
# setwd(dbDir)
# in.files <- list.files(pattern = "\\.csv$")
# in.files
# 
# require(feather)
# require(parallel)
# require(data.table)
# 
# #pubchem_base <- unique(read_feather(in.files, columns = c("CompoundName", "BaseFormula", "Identifier", "Source")))
# pubchem_base <- fread("PubChem_base.csv",sep = "\t",colClasses = "character",showProgress = TRUE)
# pubchem_uncharged <- subset(pubchem_base, !grepl("[+-]", BaseFormula))
# pubchem_clean_formulae <- subset(pubchem_uncharged, !is.na(BaseFormula))
# 
# names(pubchem_clean_formulae)
# pubchem.pos <- db.add.adduct.isotope(pubchem_clean_formulae[c(1:100),], wkz.adduct.confirmed, "positive", "pubchem_pos")
# pubchem.neg <- db.add.adduct.isotope(pubchem_uncharged$BaseFormula, wkz.adduct.confirmed, "negative", "pubchem_neg")
# 
# save(file="PubChem.RData", pubchem.pos, pubchem.neg)
# save(wkz.adduct.confirmed, file="WKZ_adduct_table.RData")
# # ----------------------------------
# 
