trim <- function (x) gsub("^\\s+|\\s+$", "", x)

trim.trailing <- function (x) sub("\\s+$", "", x)

#' @export
list.to.df <- function(lst,
                       type){
  dfStorage <- pblapply(seq_along(lst), FUN=function(y, n, i){
    matches <- as.character(y[[i]])
    if((length(matches) == 0) && (typeof(matches) == "character")){ matches = NA }
    temp.df <- data.frame(stringsAsFactors = FALSE, "mz"=rep(n[[i]], length(matches)), "CompoundName"=matches[[1]], "Adduct"=matches[[2]], "Isotope"=matches[[3]], "Source"=rep(type, length(matches)))
    temp.df
  },y=lst, n=names(lst))
  df <- rbindlist(dfStorage)
  # --- return ---
  df
}

#' @export
factorize <- function(vector){
  items <- as.factor(vector)
  new.levels <- c(0:(length(levels(items)) - 1 ))
  levels(items) <- as.numeric(new.levels)
  # --- return ---
  items
}

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

ppm_range <- function(mz, ppm) c((mz - (ppm/1000000 * mz)), (mz + (ppm/1000000 * mz)))

# --- isotopey testing stuff ---

# lactate <- "C3H5O3"
# glucose <- "C6H12O6"
# 
# glucose.2H <- mergeform.joanna(glucose, "H2")
# 
# lact.isos <- isopattern(isotopes, lactate, charge = -1)
# gluc.isos <- isopattern(isotopes, glucose.2H, charge = 2)
# 
# plot.isos <- function(isopat.table){
#   plot(isopat.table[, 1], isopat.table[, 2], type = "h", xlab = "m/z", 
#        ylab = "Relative abundance")
# }
# 
# plot.isos(lact.isos[[1]])
# plot.isos(gluc.isos[[1]])
# 
# # --- testy stuff ---
# 
# plot(lact.isos[[1]][, 1], lact.isos[[1]][, 2], type = "h", xlab = "m/z", 
#      ylab = "Relative abundance", col='blue')
# lines(gluc.isos[[1]][, 1], gluc.isos[[1]][, 2], type = "h", xlab = "m/z", 
#       ylab = "Relative abundance", col='red')
# 
# lact.isos[[1]]

# --- UNMATCHED FRACTION ---
"select(
((select distinct count([mzmed.pgrp]) 
from matches_hmdb 
where baseformula ISNULL) - (select distinct count([mzmed.pgrp]) 
from matches_hmdb))) as Unmatched_percentage"

"select (u.cnt / t.cnt) * 100 as UnmatchedFraction
from (select distinct count([mzmed.pgrp]) * 1.0 as cnt
from matches_chebi
where baseformula ISNULL
) u cross join
(select distinct count([mzmed.pgrp]) * 1.0 as cnt
from matches_chebi) t;"

# --- ADD ADDUCTS ---
# add new adduct : - )

# wkz.adduct.confirmed <- wkz.adduct.confirmed[1:(nrow(wkz.adduct.confirmed)-2)]
# colnames(wkz.adduct.confirmed)
# wkz.adduct.confirmed <- rbind(wkz.adduct.confirmed,
#                               data.table(Name=c("M"),
#                                          Ion_mode=c("positive", "negative"),
#                                          Charge=c(0),
#                                          Multi=c(1),
#                                          Formula_add=FALSE,
#                                          Formula_ded=FALSE,
#                                          Source="naturally charged",
#                                          Note="helps detect molecules that are naturally charged"))

# write adduct table
#save(wkz.adduct.confirmed, file="AdductTableWKZ.RData")

# -- some testing for internal db ---
# internal.base.db <- read.csv(file.path(dbDir, "NeededFiles", "InternalDB", "internal.db.txt"), sep="\t")
# # --- ROUTE 1: ASSUME ALL CHARGES 0 ---
# 
# internal.base.db$charge <- c(0)
# internal.base.db <- internal.base.db[,1:3]
# 
# #--------------------------------------
# conn <- dbConnect(RSQLite::SQLite(), file.path(dbDir, "hmdb.full.db"))
# # ------------------------------------------------
# temp.int.db <- data.frame(compoundname = 1:nrow(internal.base.db),
#                           formula = 1:nrow(internal.base.db),
#                           charge = 1:nrow(internal.base.db))
# # ------------------------------------------------
# for(x in 1:nrow(internal.base.db)){
#   print(x)
#   row <- internal.base.db[x,]
#   if("Dimeric" %in% row[[1]]) NA # skip dimers, should be in adduct table instead
#   name <- row[1]
#   formula <- row[2]
#   # --- check chemical formula ---
#   name <- gsub(x=name, pattern='( \\(.*\\))|( $)', replacement="")
#   if(name == "") next
#   print(paste("Current target:", name))
#   checked.formula <- check.chemform.joanna(isotopes, gsub(x=formula, pattern="( )|(\xa0)|(\xe1)", replacement=""))
#   row[2] <- checked.formula$new_formula
#   # --- find charge (some rest API?) ---
#   query <- fn$paste('SELECT compoundname, charge FROM base WHERE compoundname LIKE "%$name%" LIMIT 10')
#   matches <- dbGetQuery(conn, query)
#   print(matches)
#   match.row <- readline("Please select a row number: ")
#   row[3] <- if(match.row != "") matches[match.row, 'charge'] else "unknown"
#   print(row[3])
#   if(row[3] == "stop" | row[3] == "STOP") break
#   # use hmdb to get some answers
#   temp.int.db[x,] <- row[1:3]
# }
# 
# 
# internal.dt <- as.data.table(temp.int.db)
# ### SAVE MATCHED ONES O_O"" then do this again for chebi
# 
# int.unmatched <- internal.dt[charge == "unknown"]
# int.matched <- internal.dt[charge != "unknown"]  
# # -------------------
# # repeat for chebi?
# 
# conn2 <- dbConnect(RSQLite::SQLite(), file.path(dbDir, "chebi.full.db"))
# # ------------------------------------------------
# temp.int.db2 <- data.frame(compoundname = 1:nrow(int.unmatched),
#                           formula = 1:nrow(int.unmatched),
#                           charge = 1:nrow(int.unmatched))
# # ------------------------------------------------
# for(x in 1:nrow(int.unmatched)){
#   print(x)
#   row <- int.unmatched[x,]
#   print(row)
#   name <- row[,1]
#   formula <- row[,2]
#   # --- check chemical formula ---
#   name <- gsub(x=name, pattern='( \\(.*\\))|( $)|(ic acid)', replacement="")
#   if(name == "") next
#   print(paste("Current target:", name))
#   # --- find charge (some rest API?) ---
#   query <- fn$paste('SELECT compoundname, charge FROM base WHERE compoundname LIKE "%$name%" LIMIT 10')
#   matches <- dbGetQuery(conn2, query)
#   print(matches)
#   match.row <- readline("Please select a row number: ")
#   row[,3] <- if(match.row != "") matches[match.row, 'charge'] else "unknown"
#   print(row[,3])
#   if(row[,3] == "stop" | row[,3] == "STOP") break
#   # use hmdb to get some answers
#   temp.int.db2[x,] <- row[,1:3]
# }
# 
