#' Clean HMDB DB from previous pipeline
#'
#' @param db
#'
#' @return A cleaned DB with adducts and compoundnames and isotopes in seperate columns.
#' @export
#'
#' @examples clean.db.extended(HMDB.pos)
clean.db.extended <- function(db){
  require(Hmisc)
  db.clean <- db
  # create 'Adduct' column
  db.clean$Adduct <- c(NA)
  db.clean$Isotope <- c(0)
  # isomers
  adj.loc <- grep('(\\[M[+-].*])|( iso .)', db.clean$CompoundName)
  # get part of db with adducts
  db.clean.adduct.list <- pblapply(adj.loc, FUN=function(x){
    row <- db.clean[x,]
    compound.full <- row$CompoundName
    # -- split elements ---
    isotope <- ifelse(grepl(' iso ', compound.full), as.numeric(strsplit(compound.full, ' iso ')[[1]][[2]]), 0)
    compound <- gsub(' (\\[M[+-].*\\][+-])|( iso .*)', '', compound.full)
    adduct <- ifelse(grepl('\\[M[+-].*\\]', compound.full),
                     gsub('^(.*\\[M[+-])|(\\][+-])|( iso .*)', '', compound.full),
                     NA)
    # --- clean compound name ---
    compound <- capitalize(gsub('[\\\"]', '', compound))
    # -- assign ---
    row$CompoundName <- compound
    row$Adduct <- adduct
    row$Isotope<- isotope
    row})
  db.clean.adducts <- rbindlist(db.clean.adduct.list)
  cleaned.db <- rbind(db.clean[-adj.loc,], db.clean.adducts)
  # --- adjust composition for later isotope calculation ---
  cleaned.db
}

#' @export
clean.db.base <- function(db){
  # ---
  require(enviPat)
  data(isotopes)
  # ---
  db.clean <- db
  # create 'Adduct' column
  db.clean$Adduct <- c(NA)
  db.clean$Isotope <- c(0)
  adj.loc <- grep('(\\[M[+-].*])|( iso .)', db.clean$CompoundName)
  cleaned.db <- db.clean[-adj.loc,]
  # --- adjust composition for later isotope calculation ---
  checked <- check_chemform(isotopes,cleaned.db$Composition);
  cleaned.db$Composition <- checked$new_formula
  cleaned.db$Identifier <- row.names(cleaned.db)
  cleaned.db
}
