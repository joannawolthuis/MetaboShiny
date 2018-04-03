
replicates = 3

excel_paths = list(
  sheet1 = "~/Documents/umc/data/Data/Project 2017-021 DSM feed-5 (Brazil-1) - Saskia v Mil/Brazil - injection list no 1.xlsx",
  sheet2 = "~/Documents/umc/data/Data/Project 2017-022 DSM feed-6 (Brazil-2 & UVESA-1) - Saskia v Mil/Brazil(p2) & UVESA(p1) - injection list.xlsx",
  sheet3 = "~/Documents/umc/data/Data/Project 2017-023 DSM feed-7 (Spain Uvesa-2 & Bonarea-1) - Saskia v Mil/Injection list UVESA (p2) Bonarea (p1).xlsx",
  sheet4 = "~/Documents/umc/data/Data/Project 2018-009 DSM feed-11 (Veronesi) - Saskia v Mil_part1/DBS Veronesi.xlsm",
  sheet5 = "~/Documents/umc/data/Data/Project 2018-009 DSM feed-11 (Veronesi) - Saskia v Mil_part2/DBS Veronesi part2.xlsm"
)

excel_paths = list(
  sheet1 = "~/Documents/umc/data/Data/Project 2017-029 analysis feces (pig)_DI-MS untargeted - Myrthe Gilbert.WUR/Injection list of piglets feces (Gilbert WUR).xlsx"
)


#file4 = "~/Documents/umc/data/Data/Project 2017-021 DSM feed-5 (Brazil-1) - Saskia v Mil/Brazil - injection list no 1.xlsx"

paths = sapply(excel_paths, function(x) x)

excel = paths[[1]]

samples <- lapply(paths, FUN=function(excel){
  
  # get paths
  inj_path <- excel
  
  parent_folder <- dirname(excel)
  
  parent_dirs <- list.dirs(parent_folder) 
  
  mzpath <- grep("RES",
                 parent_dirs,
                 value = T)
  
  # find mz files
  
  mzfiles <- list.files(mzpath,
                        pattern="\\.raw$")
  
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
  df_0 <- as.data.frame(injfile[1:(length(mzfiles)/replicates), 
                                1:3])
  
  injcol <- 1:nrow(df_0)
  print(head(df_0))
  sncol_idx <- readline(prompt = "Which column has your sample names? (enter integer): ")
  sncol = df_0[,as.numeric(sncol_idx)]

  df_1 <- as.data.table(cbind(injcol, sncol))
  
  # get correct rows and columns (remove empty ones)
  colnames(df_1) <- c("Injection","Sample_Name")
  
  filepfx <- gsub(mzfiles[[1]], 
                  pattern = "_\\d*(\\.mzXML|\\.raw)", 
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

sampleNames_all$Batch <- as.numeric(as.factor(gsub(sampleNames_all$File_Name,
                                      pattern = "_\\d\\d\\d$",
                                      replacement="")))
#inj = rep(1:(nrow(sampleNames_all)/replicates), each = replicates)
#sampleNames_all$Injection <- inj
out <- split( sampleNames_all[Sample_Name == "QC"] , f = sampleNames_all[Sample_Name == "QC", c("Injection", "Batch")],drop = T)
out 

qcs = 1

library(pbapply)
cl=parallel::makeCluster(3, "FORK")
parallel::stopCluster(cl)

qc_fixed <- pblapply(out, cl=NULL,FUN=function(triple){
  triple$Sample_Name <- paste0("QC", qcs)
    qcs <<- qcs + 1
    print(triple)
  # === RETURN ===
  res = triple[,1:4]
  res
})

final_sample_list <- rbind(sampleNames_all[Sample_Name != "QC"], rbindlist(qc_fixed))

# RENAME SPANISH SAMPLES LATER...

# MOVE TO LOC
dir.create("~/Documents/umc/data/Data/BrSpIt")

groupdir = "~/Documents/umc/data/Data/BrSpIt/MZXML"

groupdir = "~/Documents/umc/data/Data/Project 2017-029 analysis feces (pig)_DI-MS untargeted - Myrthe Gilbert.WUR/RES_2017-12-07 piglets feces (Gilbert WUR)"

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
                        pattern="\\.raw$",
                        full.names = T)
  sapply(mzfiles, FUN=function(file){
    fn = basename(file)
    print(fn)
    file.symlink(file, file.path(groupdir, fn))
  })
}

# reorder sample names based on symlinks
file_order <- sapply(list.files(groupdir,pattern = "(\\.mzXML|\\.raw)"), FUN=function(file){
  gsub(basename(file),pattern = "(\\.mzXML|\\.raw)", replacement = "")
})

reordered_sample_list <- final_sample_list[match(as.character(file_order), final_sample_list$File_Name),]

# --------------------------------------

fwrite(reordered_sample_list[complete.cases(reordered_sample_list),],file = file.path(groupdir,"sampleNames.txt"),sep = "\t")

# DO PIPELINE
NULL

# BATCH EFFECT CORRECTION
library(data.table)
library(sva)
library(devtools)
library(Biobase)
library(limma)

groupdir = "~/Documents/umc/data/Data/BrazilAndSpain/MZXML"

for(mode in c("positive", "negative")){
  csv <- fread(file.path(groupdir,"results",paste0("outlist_",mode,".csv")),
               data.table = F)
  sn <- fread(file.path(groupdir,"sampleNames.txt"))
  
  sn$batch <- as.numeric(as.factor(gsub(sn$File_Name, 
                                        pattern = "_\\d\\d\\d$", 
                                        replacement="")))
  sn
  # remove QCs and 0.. for now
  keep_rows <- grepl(pattern = "BR\\d|S\\d",
                     x = sn$Sample_Name)
  keep_cols <- grepl(pattern = "BR\\d|S\\d",
                     x = colnames(csv))

  sn_adj <- sn[keep_rows,]
  csv_adj = csv[,keep_cols]
  
  # --------------------------
  csv_edata <- as.matrix(csv_adj)
  rownames(csv_edata) <- csv$mzmed
  
  # sample  outcome batch cancer
  # GSM71019.CEL      1   Normal     3 Normal
  # GSM71020.CEL      2   Normal     2 Normal
  # GSM71021.CEL      3   Normal     2 Normal
  # GSM71022.CEL      4   Normal     3 Normal
  # GSM71023.CEL      5   Normal     3 Normal
  
  sn_pheno <- unique(sn_adj[,c(1,3)])
  sn_pheno$farm <- as.numeric(as.factor(gsub(sn_pheno$Sample_Name, 
                                             pattern = "-\\d*$", 
                                             replacement = "")))
  # remove missing ones
  missing <- sn_pheno$Sample_Name[!sn_pheno$Sample_Name %in% colnames(csv_adj)]
  if(!length(missing) == 0) sn_pheno <- sn_pheno[sn_pheno$Sample_Name != missing,]
  
  csv_pheno <- data.frame(sample = sn_pheno$Sample_Name,
                          outcome = sn_pheno$farm,
                          batch = sn_pheno$batch)
  # ================
  
  mod = model.matrix(~as.factor(outcome), data=csv_pheno)
  # parametric adjustment
  
  combat_corrected = ComBat(dat=csv_edata, 
                            batch=as.factor(csv_pheno$batch))
  print("here")
  final_csv <- cbind(csv[,1:6], combat_corrected)
  fwrite(file = file.path(groupdir, paste0("outlist_", mode, "_combat.RData")),
         x = final_csv)  
}


# ==== JOIN EXCEL SHEETS ===
sheets = list.files("~/Documents/umc/metaboshiny_all/ExcelSheetsDSM/BR",
                    pattern = "BR[1-9]\\.xlsx|BR10\\.xlsx",
                    full.names = T)

all = lapply(sheets, FUN=function(sheet){
  tab <- openxlsx::read.xlsx(sheet,
                             "Individual Data",
                             skipEmptyRows=T)
  tab
})

tab_total = data.table::rbindlist(all)
tab_total

wb <- openxlsx::loadWorkbook(sheets[1])
wb
openxlsx::writeData(wb, sheet = "Individual Data", tab_total)
