
replicates = 3
sheet1 = "~/Documents/umc/data/Data/Project 2017-021 DSM feed-5 (Brazil-1) - Saskia v Mil/Brazil - injection list no 1.xlsx"
sheet2 = "~/Documents/umc/data/Data/Project 2017-022 DSM feed-6 (Brazil-2 & UVESA-1) - Saskia v Mil/Brazil(p2) & UVESA(p1) - injection list.xlsx"
sheet3 = "~/Documents/umc/data/Data/Project 2017-023 DSM feed-7 (Spain Uvesa-2 & Bonarea-1) - Saskia v Mil/Injection list UVESA (p2) Bonarea (p1).xlsx"
#file4 = "~/Documents/umc/data/Data/Project 2017-021 DSM feed-5 (Brazil-1) - Saskia v Mil/Brazil - injection list no 1.xlsx"

excel_paths = c(sheet1, sheet2, sheet3)

excel = paths[1]

samples <- lapply(excel_paths, FUN=function(excel){
  
  # get paths
  inj_path <- excel
  
  parent_folder <- dirname(excel)
  
  parent_dirs <- list.dirs(parent_folder) 
  
  mzpath <- grep("RES",
                 parent_dirs,
                 value = T)
  
  # find mz files
  
  mzfiles <- list.files(mzpath,
                        pattern="\\.mzXML$")
  
  # find necessary sheet
  wb <- xlsx::loadWorkbook(inj_path)
  sheets <- xlsx::getSheets(wb)
  inj_sheet_name <- grep(pattern = "inj",
                                 x = names(sheets),
                                 ignore.case = T,
                                 value = T)
  
  # load sheet
  injfile <- openxlsx::read.xlsx(inj_path,
                                 inj_sheet_name,
                                 skipEmptyRows=T)
  
  # injfile <- xlsx::read.xlsx(file = inj_path,
  #                            sheetName = inj_sheet_name)
  df_1 <- as.data.frame(injfile[1:(length(mzfiles)/replicates), 
                                1:3])
  
  # get correct rows and columns (remove empty ones)
  colnames(df_1) <- c("Injection","Sample_Name", 
                      "Identifier")
  
  filepfx <- gsub(mzfiles[[1]], 
                  pattern = "_\\d*\\.mzXML", 
                  replacement = "")
  
  df_2 <- df_1[rep(1:nrow(df_1),
                   each = replicates),] # replicados

  first_inj <- min(df_1$Injection)
  
  filenames <- paste(filepfx, 
                     1:nrow(df_2),  # THIS NEEDS FIXED - NOT CORRECT NOW
                     sep = "_")

  filenames <- gsub(filenames,
                    pattern="(_)([0-9])$",
                    replacement="\\100\\2")
  filenames <- gsub(filenames,
                    pattern="(_)([0-9][0-9])$",
                    replacement="\\10\\2")
  
  df_3 <- df_2
  df_3$File_Name <- filenames
  dt <- data.table::as.data.table(df_3)
  
  print(dt)
  
  sampleNames <- dt[,c("Injection","Sample_Name", "File_Name")]
  # === RETURN ===
  sampleNames
})

  # rename...
library(data.table)
sampleNames_all <- rbindlist(samples)

View(sampleNames_all)
out <- split( sampleNames_all , f = sampleNames_all$Injection )

qcs = 1
library(pbapply)
cl=parallel::makeCluster(3, "FORK")
stopCluster(cl)
qc_fixed <- pblapply(out, cl=NULL,FUN=function(triple){
  if(unique(triple$Sample_Name) == "QC"){
    #
    print(qcs)
    triple$Sample_Name <- paste0("QC", qcs)
    qcs <<- qcs + 1
    print(triple)
  }
  # === RETURN ===
  res = triple[,2:3]
  res
})

final_sample_list <- rbindlist(qc_fixed)

# RENAME SPANISH SAMPLES LATER...

# MOVE TO LOC
dir.create("~/Documents/umc/data/Data/BrazilAndSpain")

groupdir = "~/Documents/umc/data/Data/BrazilAndSpain/MZXML"

dir.create(groupdir)

# SYMLINK STUFF
for(excel in excel_paths){
  parent_folder <- dirname(excel)
  parent_dirs <- list.dirs(parent_folder)
  mzpath <- grep("RES",
                 parent_dirs,
                 value = T)
  # find mz files
  mzfiles <- list.files(mzpath,
                        pattern="\\.mzXML$",
                        full.names = T)
  sapply(mzfiles, FUN=function(file){
    fn = basename(file)
    print(fn)
    file.symlink(file, file.path(groupdir, fn))
  })
}

# reorder sample names based on symlinks
file_order <- sapply(list.files(groupdir), FUN=function(file){
  gsub(basename(file),pattern = "\\.mzXML", replacement = "")
})

reordered_sample_list <- final_sample_list[match(as.character(file_order), final_sample_list$File_Name),]

View(reordered_sample_list)
# --------------------------------------
fwrite(reordered_sample_list,file = file.path(groupdir,"sampleNames.txt"),sep = "\t")

