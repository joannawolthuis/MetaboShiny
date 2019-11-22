get_exp_vars <- function(from, patdb){
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), patdb) # change this to proper var later
  RSQLite::dbGetQuery(conn, gsubfn::fn$paste("PRAGMA table_info($from)"))$name
}

#' @export
browse_db <- function(chosen.db){
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), chosen.db) # change this to proper var later
  # --- browse ---
  result <- RSQLite::dbGetQuery(conn, "SELECT DISTINCT compoundname as name, baseformula as formula, description as description, charge as charge FROM base")
  # --- result ---
  result
}

#' @export
get_prematches <- function(who = NA,
                           what = "query_mz",
                           patdb,
                           showdb=c(),
                           showadd=c(),
                           showiso=c()){
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), patdb)
  
  firstpart = "SELECT DISTINCT
               map.query_mz as query_mz, lower(name) as name,baseformula,adduct,
                              fullformula,finalcharge,`%iso`,
                              dppm,
               description, map.structure as structure,
               GROUP_CONCAT(source) as source
               FROM match_mapper map
               JOIN match_content con
               ON map.structure = con.structure
               AND map.query_mz = con.query_mz"
  
  dbfrag = if(length(showdb)>0) gsubfn::fn$paste("AND source = '$showdb'") else ""
  addfrag = if(length(showadd)>0) gsubfn::fn$paste("AND adduct = '$showadd'") else ""
  isofrag = if(length(showiso)>0) switch(showiso, 
                                         main = "AND `%iso` > 99.9999", 
                                         minor = "AND `%iso` < 99.9999") else ""
  
  query = gsubfn::fn$paste("$firstpart WHERE $what = '$who' $dbfrag $addfrag $isofrag")

  query = paste0(query, " GROUP BY map.query_mz, name, baseformula, fullformula, finalcharge, adduct, `%iso`, dppm, map.structure, description")
  
  res = RSQLite::dbGetQuery(conn, query)
 
  if(any(grepl(pattern = "iso", colnames(res)))){
    res$isocat <- sapply(res$`%iso`, function(perc) if(perc == 100) "main" else "minor")
  }
  
  RSQLite::dbDisconnect(conn)
  return(res)
}

score.isos <- function(table, mSet, patdb, method="mscore", inshiny=TRUE, session=0, intprec, ppm){
  
  shiny::showNotification("Scoring isotopes...")
  require(InterpretMSSpectrum)
  
  formulas = unique(table$fullformula)
  
  repr.smiles <- sapply(formulas, function(form){
    table[fullformula == form][1,]$structure
  })
  
  mini.table <- table[structure %in% repr.smiles]

  isotopies = lapply(1:nrow(mini.table), function(i){
    smi=mini.table$structure[i]
    add=mini.table$adduct[i]
    form=mini.table$baseformula[i]
    revres = as.data.table(MetaDBparse::searchRev(smi, "extended", lcl$paths$db_dir))
    isotopes_to_find = revres[adduct == add]
    isotopes_to_find$form = c(form)
    isotopes_to_find
  })
  
  isotopies <- unique(isotopies)
  
  score_rows = pbapply::pblapply(isotopies, function(l){
    
    mzs = l$fullmz
    formula = unique(l$fullformula)
    
    per_mz_cols = lapply(mzs, function(mz){
      matches = which(as.numeric(colnames(mSet$dataSet$orig)) %between% MetaboShiny::ppm_range(mz, ppm))
      if(length(matches) > 0){
        int = as.data.table(mSet$dataSet$orig)[,..matches]
        int[is.na(int)] <- 0
        int = rowMeans(int)
      }else{
        int = rep(0, nrow(mSet$dataSet$orig))
      }
      l = list(values = int)
      names(l) = mz
      l[[1]]
    })
    
    bound = do.call("cbind", per_mz_cols)
    colnames(bound) = mzs
    
    theor = matrix(c(l$fullmz, l$isoprevalence), nrow=2, byrow = T)

    scores_persamp = apply(bound, MARGIN=1, FUN = function(row){
      foundiso = which(row > 0)
      if(length(foundiso) <= 1){
        return(0)
      }

      row <- sapply(row, function(x) if(x == max(row)) 100 else x/max(row))
      obs = matrix(c(as.numeric(names(row)), row), nrow=2, byrow = T)
      
      switch(method,
             mape={
               actual = obs[2,]
               theor = theor[2,]
               deltaSignal = abs(theor - actual)
               percentageDifference = deltaSignal / actual * 100# Percent by element.
               # - - -
               mean(percentageDifference) #Average percentage over all elements.
             },
             mscore={
               score = InterpretMSSpectrum::mScore(obs = obs,
                                                   the = theor,
                                                   dppm = ppm,
                                                   int_prec = input$int_prec/100)
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
    })
    meanScore = mean(scores_persamp, na.rm = T)
    data.table::data.table(fullformula = formula, score = meanScore)
  })
  data.table::rbindlist(score_rows)
}

get_user_role <- function(username, password){
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), "users.db") # change this to proper var later
  role = RSQLite::dbGetQuery(conn, gsubfn::fn$paste(
    "SELECT role FROM users WHERE username = '$username' AND password = '$password'"))
  if(nrow(role) == 0){
    return(NULL)
  }else{
    return(role[1,1])
  }
  RSQLite::dbDisconnect(conn)
}

getIonMode = function(mzs, patdb){
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), normalizePath(patdb))
  if(length(mzs) == 1){
    mode = RSQLite::dbGetQuery(conn, gsubfn::fn$paste("SELECT foundinmode FROM mzvals WHERE mzmed LIKE $mzs"))[,1]
  }else{
    temp.tbl = data.table::data.table(mzmed = mzs)
    RSQLite::dbExecute(conn, "CREATE TEMP TABLE search_query(mzmed INT)")
    RSQLite::dbWriteTable(conn, "search_query", temp.tbl, append=T)
    mode_tbl = RSQLite::dbGetQuery(conn, gsubfn::fn$paste("SELECT sq.mzmed, foundinmode FROM mzvals
                                                            JOIN search_query sq
                                                            ON mzvals.mzmed LIKE sq.mzmed"))
    mode = mode_tbl$foundinmode[match(mzs, mode_tbl$mzmed)]
    RSQLite::dbDisconnect(conn)
  }
  mode
}

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