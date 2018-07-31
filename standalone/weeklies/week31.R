# - - miranda toasted status - - -
f <- "/Users/jwolthuis/Desktop/xls/miranda_toasted.csv"
tab <- fread(f)
n <- tab[toasted == "n"]$litter
t <- tab[toasted == "t"]$litter

t.test(n, t)
# p-value = 0.1474

m <- ggplot(tab, aes(x = litter, group = toasted))
m + geom_histogram(binwidth = 0.1)

# - - - COMPARISONS FOR JEROEN - - - 

# === COMPARE SPLIT TABLES W/ WHOLE BATCH CORRECTED TABLE ===

# split csv
csv_loc <- "/Users/jwolthuis/Analysis/SP/CHICKENS.csv"

csv <- fread(csv_loc, header = T)

csv[1:10,1:10]

require(data.table)

split.by.country <- split(csv, by = "Country")
current.exp <- split.by.country

results <- pbapply::pblapply(1:length(current.exp), function(i){
  name <- names(current.exp)[i]
  table <- current.exp[[i]]
  file.name <- paste0(name, "_", basename(csv_loc))
  new.loc <- file.path(dirname(csv_loc), file.name)
  print(new.loc)
  # - get QCs -
  batches <- unique(table$Batch)
  print(batches)
  QCs <- csv[Sample %like% "^QC" & Batch %in% batches]
  result <- rbind(table, QCs)
  print(dim(result))
  # - - -
  fwrite(x = result, file = new.loc)
})
# === CHICKENS ===

week31 <- list(chicken=list(), pig=list())

save(week31, file="week31.RData")

#week31$chicken$all <- mSet$analSet
#week31$chicken$brasil <- mSet$analSet
#week31$chicken$spain <- mSet$analSet
#week31$chicken$italy <- mSet$analSet #italy is bugged RN

#week31$pig$all <- mSet$analSet
#week31$pig$italy <- mSet$analSet
#week31$pig$brasil <- mSet$analSet
#week31$pig$netherlands <- mSet$analSet



# country 1 tt fc ls
# country 2 tt fc ls
# country 3 tt fc ls
# merge tt fc ls

# === PIGS ===
# country 1 tt fc ls
# country 2 tt fc ls
# country 3 tt fc ls
# merge tt fc ls

# === COMPARE SCORING METHODS FOR INTERNAL STANDARDS ===
# mean squared error
# mscore
# sirius?
# chisquare for proportions

# === COMPARE MIRANDA TOAST VS NON TOAST W/ HITS FROM THE ABOVE (if time left) ===