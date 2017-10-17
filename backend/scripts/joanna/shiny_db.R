get_exp_vars <- function(){
  conn <- dbConnect(RSQLite::SQLite(), patdb) # change this to proper var later
  dbGetQuery(conn, "PRAGMA table_info(setup)")$name
}

get_times <- function(chosen.db){
  conn <- dbConnect(RSQLite::SQLite(), chosen.db) # change this to proper var later
  # --- browse ---
  result <- dbGetQuery(conn, "SELECT DISTINCT sampling_date as Date FROM individual_data WHERE sampling_date != ''")
  print(result$Date)
  times <- as.numeric(as.factor(as.Date(result$Date)))
  # --- result ---
  times
}

#' @export
browse_db <- function(chosen.db){
  conn <- dbConnect(RSQLite::SQLite(), chosen.db) # change this to proper var later
  # --- browse ---
  result <- dbGetQuery(conn, "SELECT DISTINCT compoundname as Compound, baseformula as Formula, description as Description, charge as Charge FROM base")
  # --- result ---
  result
}

#' @export
get_matches <- function(cpd = NA, chosen.db, search_formula=F){
  # --- connect to db -=
  req("patdb")
   # change this to proper var later
  # 0. Attach db
  print(search_formula)
  if(search_formula){
    conn <- dbConnect(RSQLite::SQLite(), chosen.db)
    query <- fn$paste(strwrap(
        "SELECT DISTINCT compoundname as Compound, identifier as Identifier, description as Description 
        FROM base indexed by b_idx1
        WHERE baseformula = '$cpd'"
        , width=10000, simplify=TRUE))
    print(query)
    dbGetQuery(conn, query)
  }else{
    conn <- dbConnect(RSQLite::SQLite(), patdb)
    query.zero <- fn$paste("ATTACH '$chosen.db' AS db")
    dbExecute(conn, query.zero)
    query.one <- fn$paste(strwrap(
      "CREATE TEMP TABLE unfiltered AS
      SELECT cpd.baseformula, cpd.adduct
      FROM mzvals mz
      JOIN mzranges rng ON rng.ID = mz.ID
      JOIN db.extended cpd indexed by e_idx2
      ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
      AND mz.foundinmode = cpd.foundinmode
      WHERE ABS(mz.[mzmed.pgrp] - $cpd) < 0.000000000001",width=10000, simplify=TRUE))
    # 1. Find matches in range (reasonably fast <3)
    dbExecute(conn, query.one)
    #  2. get isotopes for these matchies (reverse search)
    query.two <- fn$paste(strwrap(
      "CREATE TEMP TABLE isotopes AS
      SELECT cpd.baseformula, cpd.adduct, cpd.isoprevalence, cpd.basecharge 
      FROM db.extended cpd indexed by e_idx1
      JOIN unfiltered u
      ON u.baseformula = cpd.baseformula
      AND u.adduct = cpd.adduct
      JOIN mzranges rng
      ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax"
      , width=10000, simplify=TRUE))
    dbExecute(conn, query.two)
    query.three <-strwrap (
        "SELECT DISTINCT base.compoundname as Compound, base.identifier as Identifier, iso.adduct as Adduct, base.description as Description 
        FROM isotopes iso
        JOIN db.base base indexed by b_idx1
        ON base.baseformula = iso.baseformula AND
        base.charge = iso.basecharge
        GROUP BY iso.baseformula, iso.adduct
        HAVING COUNT(iso.isoprevalence > 99.99999999999) > 0"
        , width=10000, simplify=TRUE)
    # 3. get the info you want
    results <- dbGetQuery(conn,query.three)
    # --- results ---
    results
  }
  }

#' @export
get_mzs <- function(baseformula, charge, chosen.db){
  # --- connect to db ---
  req("patdb")
  conn <- dbConnect(RSQLite::SQLite(), patdb) # change this to proper var later
  print(baseformula)
  print(charge)
  query.zero <- fn$paste("ATTACH '$chosen.db' AS db")
  print(query.zero)
  dbExecute(conn, query.zero)
  # search combo of baseformula and charge matching your choice and find all possible mzvals and adducts
  query.one <-  fn$paste(strwrap(
    "CREATE TEMP TABLE possible_options AS
    SELECT DISTINCT e.fullmz, e.adduct, e.isoprevalence
    FROM db.extended e
    WHERE e.baseformula = '$baseformula' 
    AND e.basecharge = $charge"
    , width=10000, simplify=TRUE))
  print(query.one)
  
  dbExecute(conn, query.one)
  
  # join with patdb
  query.two <- fn$paste(strwrap(
    "CREATE TEMP TABLE isotopes AS
    SELECT DISTINCT mz.[mzmed.pgrp], o.*
    FROM possible_options o
    JOIN mzranges rng
    ON o.fullmz BETWEEN rng.mzmin AND rng.mzmax
    JOIN mzvals mz
    ON rng.ID = mz.ID", width=10000, simplify=TRUE))
  print(query.two)
  
  dbExecute(conn, query.two)
  # isofilter and only in
  query.three <-  strwrap(
    "SELECT DISTINCT iso.[mzmed.pgrp], adduct
    FROM isotopes iso
    GROUP BY iso.adduct
    HAVING COUNT(iso.isoprevalence > 99.99999999999) > 0"
    , width=10000, simplify=TRUE)
  print(query.three)
  
  results <- dbGetQuery(conn,query.three)
  print(results)
}

#' @export
get_all_matches <- function(exp.condition=NA, pat.conn=NA){
  # --- connect to db ---
  req("patdb")
  available.dbs <- list.files(options$db_dir,pattern = ".full.db$",full.names = T)
  join.query <- c()
  # --- GET POOL OF ALL DBS ---
  for(i in seq_along(available.dbs)){
     chosen.db <- available.dbs[i]
     dbshort <- paste0("db", i)
     dbExecute(pat.conn, fn$paste("ATTACH '$chosen.db' AS $dbshort"))
     dbpfx <- paste0(dbshort, ".extended")
     tabpfx <- paste0("dbext", i)
     mzquery <- paste0(tabpfx, ".fullmz")
     modequery <- paste0(tabpfx, ".foundinmode")
     # ------------------------
     join.query <- c(join.query, 
      fn$paste(strwrap("SELECT DISTINCT
                        mz.[mzmed.pgrp] as mz, 
                        baseformula, 
                        adduct, 
                        isoprevalence, 
                        basecharge
                        FROM mzvals mz
                        JOIN mzranges rng ON rng.ID = mz.ID
                        JOIN $dbpfx $tabpfx indexed by e_idx2
                        ON $mzquery BETWEEN rng.mzmin AND rng.mzmax
                        AND mz.foundinmode = $modequery",
                       width=10000, simplify=TRUE)
               )
      )
   }
  union <- paste(join.query, collapse = " UNION ")
  # ------------
  query.one <- fn$paste(strwrap(
    "CREATE TEMP TABLE isotopes AS
    $union",width=10000, simplify=TRUE))
    print("Exec query one")
  dbExecute(pat.conn, query.one)
  query.two <- strwrap(
    "CREATE temp TABLE results AS
    SELECT DISTINCT iso.mz as mz, iso.baseformula as Baseformula, iso.adduct as Adduct
    FROM isotopes iso
    GROUP BY mz
    HAVING COUNT(iso.isoprevalence > 99.99999999999) > 0"
    , width=10000, simplify=TRUE)
  print("Exec query two")
  dbExecute(pat.conn, query.two)

  query.collect <-  strwrap(fn$paste("select distinct i.filename, d.animal_internal_id, s.[$exp.condition] as label, d.sampling_date, r.baseformula as identifier, sum(i.intensity) as intensity
                                     from avg_intensities i
                                     join individual_data d
                                     on i.filename = d.card_id
                                     join setup s
                                     on d.[group] = s.[group]
                                     join results r
                                     on r.mz = i.mz
                                     group by d.animal_internal_id, 
                                     d.sampling_date, 
                                     r.baseformula"),
                            width=10000,
                            simplify=TRUE)
  print("Exec query three")
  # ---------------------------
  dbGetQuery(pat.conn, query.collect)
}