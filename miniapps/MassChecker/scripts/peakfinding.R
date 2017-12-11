library(data.table)
library(xlsx)

mzfiles <- list.files("/Users/jwolthuis/MountPoint/Data/Metabolomics/DSM/RES-2017-10-31_DSM DBS Brazil part 1/MZXML/")
injfile <- xlsx::read.xlsx(file = "/Users/jwolthuis/Downloads/Data/Project 2017-021 DSM feed-5 (Brazil-1) - Saskia v Mil/Brazil - Sample names VS ID.xlsx",sheetName = "randomized injection list")
df <- as.data.frame(injfile[c(1,2)])
df <- df[-c(87:89),]
filepfx <- "RES_20171031"

df.expanded <- df[rep(1:nrow(df),each=3),] # replicados
samplecount = nrow(df)
filenames <- paste(filepfx, c(1:nrow(df.expanded)), sep = "_")
filenames <- gsub(filenames,
                  pattern="(_)([0-9])$",
                  replacement="\\100\\2")
filenames <- gsub(filenames,
                  pattern="(_)([0-9][0-9])$",
                  replacement="\\10\\2")
df.w.filenames <- df.expanded
df.w.filenames$filename <- filenames
dt <- as.data.table(df.w.filenames)

dt$farm <- gsub(dt$farm.sample, pattern = "\\-\\d*$", replacement = "")

dt
wantedFarms <- c("BR1", "BR2")
wantedFiles <- dt[farm %in% wantedFarms, "filename"]

expandedFileList <- sapply(c(wantedFiles), FUN=function(f){
  paste("/Users/jwolthuis/MountPoint/Data/Metabolomics/DSM/RES-2017-10-31_DSM DBS Brazil part 1/MZXML/", f, ".mzXML", sep="")
})

# copy somewhere else...

# file.copy(from=expandedFileList, to='/Users/jwolthuis/MountPoint/Data/Metabolomics/DSM/BrazilFarm1and2/', copy.mode = TRUE)

to.xlsx <- df.w.filenames[,-1]

to.xlsx <- data.table(to.xlsx,keep.rownames = F)
to.xlsx

colnames(to.xlsx) <- c("Card ID", "File Name")

folder <- "/Users/jwolthuis/Downloads/Data/Project 2017-021 DSM feed-5 (Brazil-1) - Saskia v Mil"

file <- paste(file.path(folder, "File names vs Card IDs"), "xlsx", sep=".")
write.xlsx(to.xlsx, file=file,showNA = F)

# --------------------------------------
folder <- "/Users/jwolthuis/Downloads/Data/Project 2017-021 DSM feed-5 (Brazil-1) - Saskia v Mil"

sheetpath = file.path(folder, "File names vs Card IDs.xlsx")

sample_tab <- as.data.table(read.xlsx(file = sheetpath, 
                                      sheetIndex = 1))[, 2:3]

# === FOR HPC... ===

sampleNames <- sample_tab
colnames(sampleNames) <- c("Sample_Name", "File_Name")
write.table(sampleNames, file = "sampleNames.txt",sep = "\t",quote = F,row.names = F)

# ==================

# colnames(sample_tab) <- c("Card.ID")
# colnames(sample_tab) <- c("File.Name", "Card.ID")

qcloc <- grep("QC|Quality",
              x = levels(sample_tab$Card.ID)
              )

levels(sample_tab$Card.ID)[qcloc] <- "QC"

rmv.cols = c("fq.best", "fq.worst", "nrsamples", "avg.int")

load("~/Downloads/Data/Project 2017-021 DSM feed-5 (Brazil-1) - Saskia v Mil/outlist2identify_negative.RData")
outlist_neg <- outlist
load("~/Downloads/Data/Project 2017-021 DSM feed-5 (Brazil-1) - Saskia v Mil/outlist2identify_positive.RData")
outlist_pos <- outlist

## --------------------------------------------

colnames(outlist_pos) <- gsub(colnames(outlist_pos), pattern = "\\.pgrp", replacement = "")
to_replace_name_idx <- grep(colnames(outlist_pos), pattern = "RES_")
to_replace_names <- gsub(colnames(outlist_pos)[to_replace_name_idx], pattern = "\\.mzXML", replacement = "")
new_names <- as.character(factor(to_replace_names, levels = sample_tab$File.Name, labels = sample_tab$Card.ID))
colnames(outlist_pos)[to_replace_name_idx] <- new_names
unique_new_names <- ave(new_names, new_names, FUN=function(x) if (length(x)>1) paste0(x[1], "_", seq_along(x)) else x[1])
colnames(outlist_pos)[to_replace_name_idx] <- unique_new_names
rmv.ints <- which(colnames(outlist_pos) %in% rmv.cols)
outlist_pos <- outlist_pos[,-rmv.ints]

colnames(outlist_pos)
## --------------------------------------------

colnames(outlist_neg) <- gsub(colnames(outlist_neg), pattern = "\\.pgrp", replacement = "")
to_replace_name_idx <- grep(colnames(outlist_neg), pattern = "RES_")
to_replace_names <- gsub(colnames(outlist_neg)[to_replace_name_idx], pattern = "\\.mzXML", replacement = "")
new_names <- as.character(factor(to_replace_names, levels = sample_tab$File.Name, labels = sample_tab$Card.ID))
colnames(outlist_neg)[to_replace_name_idx] <- new_names
unique_new_names <- ave(new_names, new_names, FUN=function(x) if (length(x)>1) paste0(x[1], "_", seq_along(x)) else x[1])
colnames(outlist_neg)[to_replace_name_idx] <- unique_new_names
rmv.ints <- which(colnames(outlist_neg) %in% rmv.cols)
outlist_neg <- outlist_neg[,-rmv.ints]

setdiff(colnames(outlist_neg), colnames(outlist_pos))
outlist_neg <- outlist_neg[,-which(colnames(outlist_neg) == "BR8-8_3")]

# --------------------------------------------

save(file = file.path(folder, "poslist_summed_wkz.RData"), outlist_pos)
save(file = file.path(folder, "neglist_summed_wkz.RData"), outlist_neg)

# ============================================
# filter just farm 1...

farm1 <- grep(colnames(outlist_pos), pattern = "BR1")
outlist_pos <- cbind(outlist_pos[,1:4], outlist_pos[,..farm1])

farm1 <- grep(colnames(outlist_neg), pattern = "BR1")
outlist_neg <- cbind(outlist_neg[,1:4], outlist_neg[,..farm1])

ncol(outlist_pos)
ncol(outlist_neg)


outlist_neg <- outlist_neg[,-"BR1-9_3"]
