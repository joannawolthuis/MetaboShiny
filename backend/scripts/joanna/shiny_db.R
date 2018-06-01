get_exp_vars <- function(from){
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), patdb) # change this to proper var later
  RSQLite::dbGetQuery(conn, gsubfn::fn$paste("PRAGMA table_info($from)"))$name
}

get_times <- function(chosen.db){
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), chosen.db) # change this to proper var later
  # --- browse ---
  result <- RSQLite::dbGetQuery(conn, "SELECT DISTINCT sampling_date as Date FROM individual_data WHERE sampling_date != ''")
  print(result$Date)
  times <- as.numeric(as.factor(as.Date(result$Date)))
  # --- result ---
  times
}

#' @export
browse_db <- function(chosen.db){
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), chosen.db) # change this to proper var later
  # --- browse ---
  result <- RSQLite::dbGetQuery(conn, "SELECT DISTINCT compoundname as Compound, baseformula as Formula, description as Description, charge as Charge FROM base")
  # --- result ---
  result
}

#' @export
get_matches <- function(cpd = NA, 
                        chosen.db, 
                        search_formula=F,
                        searchid)
{
  # --- connect to db -=
  req("patdb")
  # change this to proper var later
  #print(searchid)
  # 0. Attach db
  if(!is.null(searchid) && searchid != "mz"){
    conn <- RSQLite::dbConnect(RSQLite::SQLite(), chosen.db)
    cpd <- gsub(cpd, pattern = "http\\/\\/", replacement = "http:\\/\\/")
    query <- if(searchid == "pathway"){
      gsubfn::fn$paste(strwrap("SELECT DISTINCT name as Name
                               FROM pathways
                               WHERE identifier = '$cpd'"
                               , width=10000, simplify=TRUE))
    } else{
      gsubfn::fn$paste(strwrap(
        "SELECT DISTINCT compoundname as Name, baseformula as 'Mol. Formula', identifier as Identifier, description as Description, structure as Structure
        FROM base indexed by b_idx1
        WHERE $searchid = '$cpd'"
        , width=10000, simplify=TRUE))
    }
    res <- RSQLite::dbGetQuery(conn, query)
    return(res)
  }else{
    
    #patdb = "/Users/jwolthuis/Analysis/SP/BrazilAndSpain_W.db"
    #chosen.db <- "/Users/jwolthuis/Google Drive/MetaboShiny/backend/db/wikidata.full.db"
    
    conn <- RSQLite::dbConnect(RSQLite::SQLite(), patdb)
    query.zero <- gsubfn::fn$paste("ATTACH '$chosen.db' AS db")

    print("here??")
    #cpd = 178.01787
    RSQLite::dbExecute(conn, query.zero)
    print("here??")
    
    # --- OLD (don't touch!!!) ---
    
    #RSQLite::dbExecute(conn, "drop table unfiltered")
    
    query.one <- gsubfn::fn$paste(strwrap(
      "CREATE TEMP TABLE unfiltered AS
      SELECT cpd.baseformula, cpd.adduct
      FROM mzvals mz
      JOIN mzranges rng ON rng.ID = mz.ID
      JOIN db.extended cpd indexed by e_idx2
      ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
      AND mz.foundinmode = cpd.foundinmode
      WHERE ABS(mz.mzmed - $cpd) < 0.000000000001",width=10000, simplify=TRUE))
    # 1. Find matches in range (reasonably fast <3)
    RSQLite::dbExecute(conn, query.one)
    #  2. get isotopes for these matchies (reverse search)
    
    #RSQLite::dbExecute(conn,"drop table isotopes")
    
    query.two <- gsubfn::fn$paste(strwrap(
      "SELECT cpd.baseformula, cpd.adduct, cpd.isoprevalence, cpd.basecharge, int.* 
      FROM db.extended cpd indexed by e_idx1
      JOIN unfiltered u
      ON u.baseformula = cpd.baseformula
      AND u.adduct = cpd.adduct
      JOIN mzranges rng
      ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
      JOIN mzintensities int
      ON int.mzmed BETWEEN rng.mzmin AND rng.mzmax"
      , width=10000, simplify=TRUE))
    
    table <- RSQLite::dbGetQuery(conn,query.two)
    
    table <- table[complete.cases(table),]
    
    p.cpd <- split(x = table, 
                   f = table$baseformula)
    
    if(nrow(table) > 0){
      res_rows <- pbapply::pblapply(p.cpd, cl=session_cl, function(cpd_tab){
        if(any(cpd_tab$isoprevalence > 99.999999)){
          formula = unique(cpd_tab$baseformula)
          adduct = unique(cpd_tab$adduct)
          # - - - - - - - - -
          sorted <- unique(cpd_tab[order(cpd_tab$isoprevalence, 
                                         decreasing = FALSE),])
          split.by.samp <- split(sorted, 
                                 sorted$filename)
          # - - - - - - - - - 
          score <- sapply(split.by.samp, function(samp_tab){
            samp_tab <- data.table::as.data.table(samp_tab)
            if(samp_tab[isoprevalence == max(samp_tab$isoprevalence),intensity] == max(samp_tab$intensity)) 1 else 0
          })
          confidence = sum(unlist(score))/length(split.by.samp) * 100.00
          data.table::data.table(baseformula = formula,
                                 adduct = adduct,
                                 confidence = confidence)  
        }else{
          data.table::data.table()
        }
      })
      results <- data.table::rbindlist(res_rows[!is.na(res_rows)])
      
      RSQLite::dbGetQuery(conn, "SELECT * FROM base limit 20;")
      
      formula.matches <- data.table::rbindlist(lapply(unique(results$baseformula), function(formula){
        query <- gsubfn::fn$paste(strwrap("SELECT DISTINCT compoundname as name, baseformula, identifier, description, structure
                                          FROM base indexed by b_idx1
                                          WHERE baseformula = '$formula'",width=10000, simplify=TRUE))
        RSQLite::dbGetQuery(conn, query)
      }))
      res=data.table::data.table()
      try({
        res <- merge(formula.matches, results, by = "baseformula") 
        })
    res      
    }else{
      data.table::data.table() 
    }
    
    # --- TESTING ----
    
    # RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS unfiltered")
    # RSQLite::dbExecute(conn, gsubfn::fn$paste(strwrap(
    #   "CREATE TEMP TABLE unfiltered AS
    #   SELECT mz.mzmed, cpd.fullmz, cpd.baseformula, cpd.adduct, cpd.isoprevalence
    #   FROM mzvals mz
    #   JOIN mzranges rng ON rng.ID = mz.ID
    #   JOIN db.extended cpd indexed by e_idx2
    #   ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
    #   AND mz.foundinmode = cpd.foundinmode
    #   WHERE ABS(mz.mzmed - $cpd) < 0.000001",width=10000, simplify=TRUE)))
    # 
    # table <- RSQLite::dbGetQuery(conn, 
    # "SELECT cpd.baseformula, cpd.adduct, cpd.isoprevalence, int.mzmed, int.filename, int.intensity
    # FROM unfiltered u
    # JOIN db.extended cpd
    # ON u.baseformula = cpd.baseformula
    # AND u.adduct = cpd.adduct
    # JOIN mzranges rng
    # ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
    # JOIN mzvals mz
    # on rng.id = mz.id
    # JOIN mzintensities int
    # ON mz.mzmed = int.mzmed")
    # 
    # table <- table[complete.cases(table),]
    # p.cpd <- split(x = table, 
    #                f = table$baseformula)
    
    # keep <- pbapply::pblapply(p.cpd, function(cpd_tab){
    #   print(cpd_tab)
    #   if(any(cpd_tab$isoprevalence > 99.999999)){
    #     cpd_tab
    #   }else(
    #     NA
    #     )})
    #continue <- keep[!is.na(keep)]
    
    # res_rows <- lapply(p.cpd, function(cpd_tab){
    #   formula = unique(cpd_tab$baseformula)
    #   adduct = unique(cpd_tab$adduct)
    #   # - - - - - - - - -
    #   sorted <- unique(cpd_tab[order(cpd_tab$isoprevalence, 
    #                                  decreasing = FALSE),])
    #   split.by.samp <- split(sorted, 
    #                          sorted$filename)
    #   # - - - - - - - - - 
    #   score <- sapply(split.by.samp, function(samp_tab){
    #     samp_tab <- data.table::as.data.table(samp_tab)
    #     if(samp_tab[isoprevalence == max(samp_tab$isoprevalence),intensity] == max(samp_tab$intensity)) 1 else 0
    #   })
    #   confidence = sum(unlist(score))/length(split.by.samp) * 100.00
    #   data.table::data.table(baseformula = formula,
    #                          adduct = adduct,
    #                          confidence = confidence)
    # })
    # data.table::rbindlist(res_rows)
  }}

#' @export
get_mzs <- function(baseformula, charge, chosen.db){
  # --- connect to db ---
  req("patdb")
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), patdb) # change this to proper var later
  query.zero <- gsubfn::fn$paste("ATTACH '$chosen.db' AS db")
  RSQLite::dbExecute(conn, query.zero)
  # search combo of baseformula and charge matching your choice and find all possible mzvals and adducts
  query.one <-  gsubfn::fn$paste(strwrap(
    "CREATE TEMP TABLE possible_options AS
    SELECT DISTINCT e.fullmz, e.adduct, e.isoprevalence
    FROM db.extended e
    WHERE e.baseformula = '$baseformula' 
    AND e.basecharge = $charge"
    , width=10000, simplify=TRUE))
  RSQLite::dbExecute(conn, query.one)
  
  # join with patdb
  query.two <- gsubfn::fn$paste(strwrap(
    "CREATE TEMP TABLE isotopes AS
    SELECT DISTINCT mz.mzmed, o.*
    FROM possible_options o
    JOIN mzranges rng
    ON o.fullmz BETWEEN rng.mzmin AND rng.mzmax
    JOIN mzvals mz
    ON rng.ID = mz.ID", width=10000, simplify=TRUE))
  print(query.two)
  
  RSQLite::dbExecute(conn, query.two)
  # isofilter and only in
  query.three <-  strwrap(
    "SELECT DISTINCT iso.mzmed, adduct
    FROM isotopes iso
    GROUP BY iso.adduct
    HAVING COUNT(iso.isoprevalence > 99.99999999999) > 0"
    , width=10000, simplify=TRUE)
  results <- RSQLite::dbGetQuery(conn,query.three)
  # --------
  results
}

#' @export
get_all_matches <- function(#exp.condition=NA, 
  pat.conn=NA, 
  which_dbs=NA, 
  which_adducts=c("M+H", "M-H", "M"),
  group_by="baseformula"
  #,var_table="setup"
  #,batches = NULL
){
  # --- connect to db ---
  # req("patdb")
  # available.RSQLite::dbs <- list.files(options$db_dir,pattern = ".full.db$",full.names = T)
  # which_dbs <- available.dbs
  # group_by="mz"
  # exp.condition="stool_condition"
  # pat.conn =RSQLite::dbConnect(RSQLite::SQLite(),patdb)
  # RSQLite::dbDisconnect(pat.conn)
  join.query <- c()
  # --- GET POOL OF ALL DBS ---
  for(i in seq_along(which_dbs)){
    chosen.db <- which_dbs[i]
    if(is.na(chosen.db)) next
    # --------------------------
    dbshort <- paste0("db", i)
    print(chosen.db)
    RSQLite::dbExecute(pat.conn, gsubfn::fn$paste("ATTACH '$chosen.db' AS $dbshort"))
    # --- extended ---
    dbext <- paste0(dbshort, ".extended")
    extpfx <- paste0("dbext", i)
    mzquery <- paste0(extpfx, ".fullmz")
    modequery <- paste0(extpfx, ".foundinmode")
    extchargequery <- paste0(extpfx, ".basecharge")
    extformquery <- paste0(extpfx, ".baseformula")
    # --- base ---
    dbbase <- paste0(dbshort, ".base")
    basepfx <- paste0("dbbase", i)
    basechargequery <- paste0(basepfx, ".charge")
    baseformquery <- paste0(basepfx, ".baseformula")
    # --- VERY OPTIONAL ---
    pathway = if(group_by == "pathway") ",pathway" else{""}
    # ------------------------
    join.query <- c(join.query, 
                    gsubfn::fn$paste(strwrap("SELECT DISTINCT mz.mzmed as mz, 
                                             compoundname,
                                             identifier,
                                             $baseformquery, 
                                             adduct, 
                                             isoprevalence, 
                                             basecharge
                                             $pathway
                                             FROM mzvals mz
                                             JOIN mzranges rng ON rng.ID = mz.ID
                                             JOIN $dbext $extpfx indexed by e_idx2
                                             ON $mzquery BETWEEN rng.mzmin AND rng.mzmax
                                             AND mz.foundinmode = $modequery
                                             JOIN $dbbase $basepfx ON $extchargequery = $basechargequery
                                             AND $extformquery = $baseformquery",
                                             width=10000, simplify=TRUE)
                    )
                    )
  }
  union <- paste(join.query, collapse = " UNION ")
  # ------------
  query.one <- gsubfn::fn$paste(strwrap(
    "CREATE TEMP TABLE isotopes AS
    $union",width=10000, simplify=TRUE))
  print("Exec query one")
  RSQLite::dbExecute(pat.conn, query.one)
  adductfilter <- paste0("WHERE adduct = '", paste(which_adducts, collapse= "' OR adduct = '"), "'")
  idquery <- paste0("iso.", group_by)
  # adductfilter <- paste0("WHERE adduct NOT REGEXP '[", paste(which_adducts, collapse="|"), "]'")
  query.two <- gsubfn::fn$paste(strwrap(
    "CREATE temp TABLE results AS
    SELECT DISTINCT iso.mz as mz, $idquery as identifier, iso.adduct as adduct
    FROM isotopes iso
    $adductfilter
    GROUP BY mz
    HAVING COUNT(iso.isoprevalence > 99.99999999999) > 0"
    , width=10000, simplify=TRUE))
  print("Exec query two")
  RSQLite::dbExecute(pat.conn, query.two)
  # a1 = RSQLite::dbGetQuery(pat.conn, gsubfn::fn$paste("SELECT distinct filename FROM avg_intensities"))
  # a2 = RSQLite::dbGetQuery(pat.conn, gsubfn::fn$paste("SELECT distinct card_id FROM individual_data"))
  summary = if(group_by == "pathway") "sum(abs(i.intensity)) as intensity" else{"sum(i.intensity) as intensity"}
  
  # --- batch ---
  
  query.collect <-  strwrap(gsubfn::fn$paste("select distinct i.filename, 
                                             d.*,
                                             s.*,
                                             b.*,
                                             r.identifier as identifier, 
                                             $summary
                                             from mzintensities i
                                             join individual_data d
                                             on i.filename = d.card_id
                                             join setup s on d.[Group] = s.[Group]
                                             join results r
                                             on r.mz = i.mzmed
                                             join batchinfo b
                                             on b.sample = d.card_id
                                             group by d.card_id, 
                                             d.sampling_date, 
                                             r.identifier"),
                            width=10000,
                            simplify=TRUE)
  print(query.collect)
  print("Exec query three")
  res <- RSQLite::dbGetQuery(pat.conn, query.collect)
  # ---------------------------
  res
}

multimatch <- function(cpd, dbs, searchid){
  match_list <- lapply(dbs, FUN=function(match.table){
    dbname <- gsub(basename(match.table), pattern = "\\.full\\.db", replacement = "")
    print(dbname)
    res <- get_matches(cpd, 
                       match.table, 
                       searchid=searchid)
    if(nrow(res) > 0){
      res = cbind(res, source = c(dbname))
    }
  })
  
  if(is.null(unlist(match_list))) return(data.table(name = "None",
                                                    description = "Unknown compound"))
  
  match_table <- (as.data.table(rbindlist(match_list))[name != ""])
  # --- sort ---
  match_table <- match_table[order(match(match_table$adduct, sort_order))]
  # --- merge description and identifier ---
  merged_cols <- paste(match_table$source, paste0("(",match_table$identifier,")."), match_table$description)
  match_table$description <- merged_cols
  # ----------------------------------------
  if("isopercentage" %in% colnames(match_table)){
    match_table$isopercentage <- round(match_table$IsoPerc, digits = 2)
  }
  # --- return ---
  match_table[,-c("source", "identifier")]
}
