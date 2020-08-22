#' @title Get experimental variables from a MetShi csv
#' @description Sees which columns are numeric and which aren't, the ones that are not are the metadata
#' @param patcsv Source CSV
#' @return Head of table metadata
#' @seealso 
#'  \code{\link[data.table]{fread}}
#' @rdname get_exp_vars
#' @export 
#' @importFrom data.table fread
get_exp_vars <- function(patcsv){
  header = data.table::fread(patcsv, nrows = 5, header=T)
  return(header[getColDistribution(header)$meta])
}

#' @title Browse a base database
#' @description Deprecated, should be moved to MetaDBparse. Returns full database (can be memory heavy)
#' @param chosen.db Which base db to browse (full path).
#' @return Data table of whole database
#' @seealso 
#'  \code{\link[RSQLite]{SQLite}}
#' @rdname browse_db
#' @export 
#' @importFrom RSQLite dbConnect SQLite dbGetQuery
browse_db <- function(chosen.db){
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), chosen.db) # change this to proper var later
  # --- browse ---
  result <- RSQLite::dbGetQuery(conn, "SELECT DISTINCT compoundname, baseformula as formula, description as description, charge as charge FROM base")
  # --- result ---
  result
}

#' @title Get results of a m/z search
#' @description Search results are stored in a SQLITE table, which can be returned through this function with filters.
#' @param who Which m/z to return, Default: NA
#' @param what Which column are we matching the 'who' to?, Default: 'query_mz'
#' @param patdb SQLITE database file of current project
#' @param showdb Which databases to include, Default: c()
#' @param showadd Which adducts to include, Default: c()
#' @param showiso Which isotope category to include, Default: c()
#' @return Data table with results
#' @seealso 
#'  \code{\link[RSQLite]{SQLite}}
#'  \code{\link[gsubfn]{fn}}
#' @rdname get_prematches
#' @export 
#' @importFrom RSQLite dbConnect SQLite dbGetQuery dbDisconnect
#' @importFrom gsubfn fn
get_prematches <- function(who = NA,
                           what = "query_mz",
                           patdb,
                           showdb=c(),
                           showadd=c(),
                           showiso=c()){
  
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), patdb)
  
  firstpart = strwrap("SELECT DISTINCT map.query_mz, compoundname,
                               map.baseformula as baseformula,
                               map.adduct as adduct,
                               con.identifier,
                               `%iso`,
                               fullformula,
                               finalcharge,
                               dppm,
                               description, con.structure as structure,
                               source
                               FROM match_mapper map
                               JOIN match_content con
                               ON map.baseformula = con.baseformula
                               AND map.adduct = con.adduct
                               AND map.query_mz = con.query_mz", simplify=T, width=1000)
  
  showadd <- if(is.null(showadd)) c() else if(length(showadd) > 0) paste0(showadd, collapse=" OR map.adduct = '", "'") else c()
  showdb <- if(is.null(showdb)) c() else if(length(showdb) > 0) paste0(showdb, collapse=" OR source = '", "'") else c()
  showiso <- if(is.null(showiso)) c() else{
    if(length(showiso) > 0){
      if(length(showiso) == 2) c() else showiso
    } else c()
  }
  
  dbfrag = if(length(showdb)>0) gsubfn::fn$paste("AND (source = '$showdb)") else ""
  addfrag = if(length(showadd)>0) gsubfn::fn$paste("AND (map.adduct = '$showadd)") else ""
  isofrag = if(length(showiso)>0) switch(showiso, 
                                         main = "AND `%iso` > 99.9999", 
                                         minor = "AND `%iso` < 99.9999") else ""
  
  who = gsub("0+$", "", who)
  query = gsubfn::fn$paste("$firstpart WHERE $what = '$who' $dbfrag $addfrag $isofrag")

  res = RSQLite::dbGetQuery(conn, query)
 
  if(any(grepl(pattern = "iso", colnames(res)))){
    res$isocat <- sapply(res$`%iso`, function(perc) if(perc == 100) "main" else "minor")
  }
  
  if(nrow(res) > 0){
    has.no.struct = which(trimws(res$structure) == "")
    if(length(has.no.struct) > 0){
      res[has.no.struct,]$structure <- paste0("[", 
                                             res[has.no.struct,]$baseformula, "]", 
                                             res[has.no.struct,]$finalcharge, "_", 
                                             res[has.no.struct,]$identifier)
    }  
  }
  
  
  RSQLite::dbDisconnect(conn)
  return(res)
}

#' @title Perform isotope scoring
#' @description For a given compound, check patient data to see if the theoretical isotope pattern is close to the one seen.
#' @param table Result table from match finding
#' @param mSet mSet item
#' @param method Method to perform scoring, currently only mScore from 'InterpretMSSpectrum', Default: 'mscore'
#' @param inshiny Are we running in shiny or outside?, Default: TRUE
#' @param session Shiny session, Default: 0
#' @param intprec Percentage the intensity is expected to be off
#' @param ppm Parts per million m/z error margin allowed
#' @param dbdir Where is your database stored?
#' @return Input 'table' with score column added.
#' @seealso 
#'  \code{\link[shiny]{showNotification}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[InterpretMSSpectrum]{mScore}}
#'  \code{\link[data.table]{rbindlist}}
#' @rdname score.isos
#' @export 
#' @importFrom shiny showNotification
#' @importFrom pbapply pblapply
#' @importFrom InterpretMSSpectrum mScore
#' @importFrom data.table data.table rbindlist
score.isos <- function(qmz, table, mSet, method="mscore", inshiny=TRUE, 
                       session=0, intprec, ppm, dbdir,
                       rtmode=F, rtperc=0.1, useint=T){
  
  if(is.null(rtmode)) rtmode = F
  
  if(inshiny) shiny::showNotification("Scoring isotopes...")

  formulas = unique(table$fullformula)
  
  repr.smiles <- sapply(formulas, function(form){
    table[fullformula == form][1,]$structure
  })
  
  mini.table <- table[structure %in% repr.smiles]

  isotopies = lapply(1:nrow(mini.table), function(i){
    smi=mini.table$structure[i]
    add=mini.table$adduct[i]
    form=mini.table$baseformula[i]
    revres = data.table::as.data.table(MetaDBparse::searchRev(smi, "extended", dbdir))
    isotopes_to_find = revres[adduct == add]
    isotopes_to_find$form = c(form)
    isotopes_to_find
  })
  
  isotopies <- unique(isotopies)
  ionMode = if(grepl("\\-", qmz)) "neg" else "pos"
  
  if(rtmode){
    rt = as.numeric(gsub("(.*RT)", "", qmz))
    rtRange = c(rt - rtperc/100 * rt, rt + rtperc/100 * rt)  
    splitMzRt=stringr::str_split(colnames(mSet$dataSet$proc), "RT")
  }else{
    if(grepl("RT", qmz)){
      splitMzRt=stringr::str_split(colnames(mSet$dataSet$proc), "RT")
    }else{
      splitMzRt=lapply(colnames(mSet$dataSet$proc), function(mz) list(mz, NA))
    }
  }
  
  mzRtTable = data.table::as.data.table(do.call("rbind", splitMzRt))
  
  colnames(mzRtTable) = c("mz", "rt")
  mzRtTable$mzfull = colnames(mSet$dataSet$proc)
  mzRtTable$mzmode = sapply(mzRtTable$mz, function(mz) if(grepl("\\-", mz)) "neg" else "pos")
  mzRtTable <- mzRtTable[mzmode == ionMode]
  mzRtTable$mz = gsub("\\+|\\-", "", mzRtTable$mz)
  mzRtTable$mz <- as.numeric(mzRtTable$mz) 
  mzRtTable$rt <- as.numeric(mzRtTable$rt) 
  
  if(rtmode){
    mzRtTable <- mzRtTable[mzRtTable$rt %between% rtRange,]
  }
  
  if(nrow(mzRtTable) == 0){
    if(inshiny) shiny::showNotification("No isotopes within RT range...")
    return(data.table::data.table(fullformula = mini.table$fullformula, score = c(0)))
  }
  
  score_rows = pbapply::pblapply(isotopies, function(l){
    
    mzs = l$fullmz
    formula = unique(l$fullformula)
    
    per_mz_cols = lapply(mzs, function(mz){
      
      matches = mzRtTable$mzfull[which(mzRtTable$mz %between% MetaboShiny::ppm_range(mz, ppm))]
      
      if(length(matches) > 0){
        int = data.table::as.data.table(mSet$dataSet$proc)[, ..matches]
        int[is.na(int)] <- 0
        int = rowMeans(int)
      }else{
        int = rep(0, nrow(mSet$dataSet$proc))
      }
      l = list(values = int)
      names(l) = mz
      l[[1]]
    })
    
    bound = do.call("cbind", per_mz_cols)
    colnames(bound) = mzs
    
    
    theor = matrix(c(l$fullmz, l$isoprevalence), nrow=2, byrow = T)

    scores_persamp = apply(bound, MARGIN=1, FUN = function(row){
      
      foundiso = which(row != 0)

      if(length(foundiso) <= 1){
        return(0)
      }

      if(useint){
        row <- sapply(row, function(x) if(x == max(row)) 100 else x/max(row))
        obs = matrix(c(as.numeric(names(row)), row), nrow=2, byrow = T)
        switch(method,
               mape={
                 actual = obs[2,]
                 theor = theor[2,]
                 deltaSignal = abs(theor - actual)
                 percentageDifference = deltaSignal / actual * 100 # Percent by element.
                 # - - -
                 mean(percentageDifference) #Average percentage over all elements.
               },
               mscore={
                 score = InterpretMSSpectrum::mScore(obs = obs,
                                                     the = theor,
                                                     dppm = ppm,
                                                     int_prec = intprec/100)
                 score
               },
               sirius={NULL},
               chisq={
                 test <- chisq.test(obs[2,],
                                    p = theor[2,],
                                    rescale.p = T)
                 # - - -
                 as.numeric(test$p.value)
               }
        ) 
      }else{
        return(length(foundiso) / ncol(theor) * 100) # percentage of possible isotopes found
      }
    })
    meanScore = mean(scores_persamp, na.rm = T)
    data.table::data.table(fullformula = formula, score = meanScore)
  })
  data.table::rbindlist(score_rows)
}


#' @title Filter patient database
#' @description Remove samples from DB file that do not have metadata
#' @param patdb Patient database file (full path)
#' @seealso 
#'  \code{\link[RSQLite]{SQLite}}
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[gsubfn]{fn}}
#' @rdname filterPatDB
#' @export 
#' @importFrom RSQLite dbConnect SQLite dbGetQuery dbExecute dbDisconnect
#' @importFrom pbapply pblapply
#' @importFrom gsubfn fn
filterPatDB <- function(patdb){
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), normalizePath(patdb))
  # which samples to remove?
  cat("Removing samples without metadata from new DB file...\n")
  to_remove <- RSQLite::dbGetQuery(conn, "SELECT DISTINCT filename FROM mzintensities WHERE filename
                                          NOT IN (SELECT DISTINCT sample FROM individual_data)")[,1]
  
  pbapply::pblapply(to_remove, function(sample){
    RSQLite::dbExecute(conn, gsubfn::fn$paste("DELETE FROM mzintensities WHERE filename='$sample'"))
  })
  
  # drop mz values that are not in mzintensities anymore
  cat("Removing mz values without samples from new DB file...\n")
  RSQLite::dbExecute(conn, "DELETE FROM mzvals WHERE mzmed
                     NOT IN (SELECT DISTINCT mzmed FROM mzintensities)")
  RSQLite::dbExecute(conn, "VACUUM")
  RSQLite::dbDisconnect(conn)
}

#' @title Checks for integrating metadata and peak table
#' @param conn Database connection
#' @rdname prepDatabase
#' @export 
#' @importFrom RSQLite dbExecute
prepDatabase <- function(conn){
  cat("Checking for mismatches between peak tables and metadata... \n")
  
  fn_meta <- allSampInMeta(conn)
  fn_int <- allSampInPeaktable(conn)
  
  cat(paste0("-- in peaklist, not in metadata: --- \n", 
             paste0(setdiff(fn_int,
                            fn_meta), 
                    collapse=", "), 
             "\n"))
  cat(paste0("-- in metadata, not in peaklist: --- \n", 
             paste0(setdiff(fn_meta,
                            fn_int), 
                    collapse=", "), 
             "\n\n"))
  
  RSQLite::dbExecute(conn, "PRAGMA journal_mode=WAL;")
  RSQLite::dbExecute(conn, "CREATE INDEX IF NOT EXISTS filenames ON mzintensities(filename)")
}

#' @title Generate query to get long format MetaboShiny table from database file
#' @description Checks if 'setup' table exists (original metadata format) and adjusts query based on that.
#' @param conn Database connection
#' @return SQLITE query to get long format MetaboShiny table.
#' @seealso 
#'  \code{\link[DBI]{dbExistsTable}}
#'  \code{\link[gsubfn]{fn}}
#' @rdname getCSVquery
#' @export 
#' @importFrom DBI dbExistsTable
#' @importFrom gsubfn fn
getCSVquery <- function(conn){
  if(DBI::dbExistsTable(conn, "setup")){
    query <- strwrap(gsubfn::fn$paste("select distinct d.sample as sample, d.*, s.*
                                        from mzintensities i
                                        join individual_data d
                                        on i.filename = d.sample
                                        join setup s on d.[Group] = s.[Group]"),
                     width=10000,
                     simplify=TRUE)   
  }else{
    query <- strwrap(gsubfn::fn$paste("select distinct d.sample as sample, d.*
                                        from mzintensities i
                                        join individual_data d
                                        on i.filename = d.sample"),
                     width=10000,
                     simplify=TRUE)
  }
}

#' @title Get all m/z values from database
#' @description Collects all m/z values from given database.
#' @param conn Database connection.
#' @return Vector of m/z values.
#' @rdname allMZ
#' @export 
#' @importFrom RSQLite dbGetQuery
allMZ <- function(conn){
  RSQLite::dbGetQuery(conn, "select distinct i.mzmed
                             from mzintensities i
                             join individual_data d
                             on i.filename = d.sample")[,1]  
}

#' @title Get all samples in metadata table
#' @param conn Database connection
#' @return Vector of sample names.
#' @rdname allSampInMeta
#' @export 
#' @importFrom RSQLite dbGetQuery
allSampInMeta <- function(conn){
  RSQLite::dbGetQuery(conn, "SELECT DISTINCT sample FROM individual_data")[,1]
}

#' @title Get all samples in peaktable
#' @param conn Database connection
#' @return Vector of sample names.
#' @rdname allSampInPeaktable
#' @export 
#' @importFrom RSQLite dbGetQuery
allSampInPeaktable <- function(conn){
  RSQLite::dbGetQuery(conn, "SELECT DISTINCT filename FROM mzintensities")[,1]
}

#' @title Get metadata for given sample
#' @description Given a file name, grabs metadata for that specific sample.
#' @param conn Database connection
#' @param filename File name as seen in peak table
#' @param query Base query to fetch relevant sample data
#' @return Metadata vector
#' @seealso 
#'  \code{\link[gsubfn]{fn}}
#'  \code{\link[data.table]{as.data.table}}
#' @rdname getSampMeta
#' @export 
#' @importFrom gsubfn fn
#' @importFrom data.table as.data.table
#' @importFrom RSQLite dbGetQuery
getSampMeta <- function(conn, filename, query){
  # adjust query
  query_add = gsubfn::fn$paste(" WHERE i.filename = '$filename'")
  
  # get results for sample
  z.meta = data.table::as.data.table(RSQLite::dbGetQuery(conn, paste0(query, query_add)))
  colnames(z.meta) <- tolower(colnames(z.meta))
  z.meta$sample <- gsub(z.meta$sample, pattern=" |\\(|\\)|\\+", replacement="")
  
  if(nrow(z.meta)==0) return(NA) else return(z.meta)
}

#' @title Get m/z intensities for given sample
#' @description Given a filename, fetches all m/z value intensities
#' @param conn Database connection
#' @param filename File name as seen in peak table.
#' @param all_mz All m/z values (needed for proper reordering)
#' @return Vector of intensities
#' @seealso 
#'  \code{\link[gsubfn]{fn}}
#'  \code{\link[data.table]{as.data.table}},\code{\link[data.table]{dcast.data.table}}
#' @rdname getSampInt
#' @export 
#' @importFrom gsubfn fn
#' @importFrom data.table as.data.table dcast.data.table
#' @importFrom RSQLite dbGetQuery
getSampInt <- function(conn, filename, all_mz){
  query_add = gsubfn::fn$paste(" WHERE i.filename = '$filename'")
  z.int = data.table::as.data.table(RSQLite::dbGetQuery(conn, 
                                                        paste0("SELECT DISTINCT
                                                                i.mzmed,
                                                                i.intensity
                                                                FROM mzintensities i", query_add)))
  if(nrow(z.int)==0) return(NA)
  
  missing_mz <- setdiff(all_mz, z.int$identifier)
  
  # cast to wide
  cast.dt <- data.table::dcast.data.table(z.int,
                                          formula = ... ~ identifier,
                                          fun.aggregate = sum,
                                          value.var = "intensity")
  suppressWarnings({
    complete = cast.dt[1,]
  })
  names(complete) = colnames(cast.dt)
  missing = rep(NA, length(missing_mz))
  names(missing) <- missing_mz
  complete.row = c(complete[-1], missing)
  reordered <- order(names(complete.row))
  complete.row <- complete.row[reordered]    
  complete.row.dt <- data.table::as.data.table(t(data.table::as.data.table(complete.row)))
  colnames(complete.row.dt) <- names(complete.row)  
  complete.row.dt
}

#' @export
score.add <- function(qmz, table, mSet, inshiny=TRUE, 
                     session=0, mzppm, rtperc=0.1, rtmode=F, dbdir, 
                     adducts_considered = adducts$Name){
  
  if(is.null(rtmode)) rtmode = F
  
  maptable = pbapply::pblapply(1:nrow(table), function(i){
    smi=table$structure[i]
    form=table$baseformula[i]
    revres = data.table::as.data.table(MetaDBparse::searchRev(smi, 
                                                              "extended", 
                                                              dbdir))
    to_find = revres[adduct %in% adducts_considered & isoprevalence == 100]
    to_find$form = c(form)
    to_find$smi = smi
    to_find
  })
  
  maptable <- unique(maptable)
  maptable = maptable[sapply(maptable, function(x) nrow(x) > 0)]
  
  ionMode = if(grepl("\\-", qmz)) "neg" else "pos"
  
  if(rtmode){
    rt = as.numeric(gsub("(.*RT)", "", qmz))
    rtRange = c(rt - rtperc/100 * rt, rt + rtperc/100 * rt)
    splitMzRt=stringr::str_split(colnames(mSet$dataSet$proc), "RT")
    mzRtTable = data.table::as.data.table(do.call("rbind", splitMzRt))
    colnames(mzRtTable) = c("mz", "rt")
    mzRtTable$mzfull = colnames(mSet$dataSet$proc)
    mzRtTable <- mzRtTable[mzmode == ionMode]
    mzRtTable$mz = gsub("\\+|\\-", "", mzRtTable$mz)
    mzRtTable$mz <- as.numeric(mzRtTable$mz) 
    mzRtTable$rt <- as.numeric(mzRtTable$rt) 
    mzRtTable <- mzRtTable[mzRtTable$rt %between% rtRange,] 
  }else{
    mzRtTable <- data.table::data.table(mz = gsub("(RT.*)-?$", "",
                                                  colnames(mSet$dataSet$proc)))
  }
  
  mzRtTable$mzmode = sapply(mzRtTable$mz, function(mz) if(grepl("\\-", mz)) "neg" else "pos")
  mzRtTable$mz = as.numeric(gsub("\\+|\\-", "", mzRtTable$mz))
  
  score_rows = pbapply::pblapply(maptable, function(l){
    mzs = l$fullmz
    formula = unique(l$fullformula)
    per_mz_cols = sapply(mzs, function(mz){
      matches_mz = which(mzRtTable$mz %between% MetaboShiny::ppm_range(mz, mzppm) &
                         mzRtTable$mzmode == ionMode)
      matchTable = mzRtTable[matches_mz,]
      nrow(matchTable) > 0
    })
    list(structure = unique(l$smi), 
         score = sum(per_mz_cols)/length(adducts_considered) * 100)
  })
  data.table::rbindlist(score_rows)
}
