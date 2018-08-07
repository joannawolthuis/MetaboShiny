# - - - network analysis - - -

paper.files <- "/Users/jwolthuis/Downloads/1-s2.0-S1570023217305688-mmc12/"
setwd(paper.files)

## NEED RECONMAP METABOLITES (about 3000 of em)

api_url <- "https://vmh.uni.lu/_api/metabolites/"

pagerange = 110
# get the first page

stop = FALSE

table_list <- pbapply::pblapply(1:pagerange, function(i){
  tbl = NA
  try({
    url = gsubfn::fn$paste("http://vmh.uni.lu/_api/metabolites/?page=$i")
    print(url)
    r <- httr::GET(url, httr::accept(".json"))
    lst <- jsonlite::fromJSON(httr::content(r, "text"))
    tbl <- lst[[4]]
    Sys.sleep(.1)
  })
  # - - return - - 
  tbl
})

table_main <- data.table::rbindlist(table_list[!is.na(table_list)])

db.formatted <- data.table::data.table(compoundname = table_main$fullName,
                                       description = table_main$description,
                                       baseformula = table_main$chargedFormula, 
                                       identifier= table_main$reconMap,
                                       charge= table_main$charge,
                                       structure= table_main$smile)

missing.desc <- which(db.formatted$description == "<NA>" | db.formatted$description == "")
replacements <- table_main$synonyms # use synonum instead
db.formatted$description[missing.desc] <- replacements[missing.desc]
missing.desc <- which(db.formatted$description == "<NA>" | db.formatted$description == "")
db.formatted$description[missing.desc] <- c("No further intel")


# need a matrix to begin with
# use HMDB as source
conn <- RSQLite::dbConnect(RSQLite::SQLite(), patdb)
db <- "/Users/jwolthuis/Google Drive/MetaboShiny/backend/db/hmdb.full.db"
adducts <- c(pos_adducts$Name, neg_adducts$Name)

res <- get_all_matches(conn, db, adducts)



# - - paper files - - 

source("configurationFile.R")
source("R/constructNetworkModel.R")
source("R/getMetabolitePairs.R")
source("R/assignReactionsByMass.R")
source("R/filterReactions.R")
source("R/summarizeUnknownAnnotation.R")
source("R/filterColumnByNA.R")
source("R/learnPathways.R")
source("R/predictPathways.R")
source("R/exportGraph.R")