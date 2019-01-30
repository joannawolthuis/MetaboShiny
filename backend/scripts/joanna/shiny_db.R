get_exp_vars <- function(from){
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), global$paths$patdb) # change this to proper var later
  RSQLite::dbGetQuery(conn, gsubfn::fn$paste("PRAGMA table_info($from)"))$name
}

get_times <- function(chosen.db){
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), chosen.db) # change this to proper var later
  # --- browse ---
  result <- RSQLite::dbGetQuery(conn, "SELECT DISTINCT sampling_date as Date FROM individual_data WHERE sampling_date != ''")
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
                        search_formula = F,
                        searchid = "mz",
                        inshiny=TRUE,
                        append=FALSE){
  
  # --- connect to db ---
  
  req("global$paths$patdb")
  
  # 0. Attach db
  
  if(!exists(searchid) && searchid != "mz"){
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
    
    conn <- RSQLite::dbConnect(RSQLite::SQLite(), global$paths$patdb)
    
    #cpd = curr_cpd
    
    RSQLite::dbExecute(conn, gsubfn::fn$paste("ATTACH '$chosen.db' AS db"))
    
    func <- function(){
      if(!DBI::dbExistsTable(conn, "unfiltered") | !append)
      {
        RSQLite::dbExecute(conn, 'DROP TABLE IF EXISTS unfiltered')
        # create
        RSQLite::dbExecute(conn, gsubfn::fn$paste(strwrap(
          "CREATE TABLE unfiltered AS
          SELECT DISTINCT cpd.baseformula as baseformula, 
          cpd.fullformula, 
          cpd.adduct as adduct, 
          cpd.isoprevalence as isoprevalence, 
          cpd.basecharge as charge, 
          (ABS($cpd - cpd.fullmz) / $cpd)/1e6 AS dppm
          FROM mzvals mz
          JOIN mzranges rng ON rng.ID = mz.ID
          JOIN db.extended cpd indexed by e_idx2
          ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
          AND mz.foundinmode = cpd.foundinmode
          WHERE $cpd BETWEEN rng.mzmin AND rng.mzmax",width=10000, simplify=TRUE)))
      } else{
        # append
        RSQLite::dbExecute(conn, gsubfn::fn$paste(strwrap(
          "INSERT INTO unfiltered
          SELECT DISTINCT cpd.baseformula as baseformula, 
          cpd.fullformula, 
          cpd.adduct as adduct, 
          cpd.isoprevalence as isoprevalence, 
          cpd.basecharge as charge, 
          (ABS($cpd - cpd.fullmz) / $cpd)/1e6 AS dppm
          FROM mzvals mz
          JOIN mzranges rng ON rng.ID = mz.ID
          JOIN db.extended cpd indexed by e_idx2
          ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
          AND mz.foundinmode = cpd.foundinmode
          WHERE $cpd BETWEEN rng.mzmin AND rng.mzmax",width=10000, simplify=TRUE)))
      }     
      
      # 1. Find matches in range (reasonably fast <3)
      if(inshiny) shiny::setProgress(value = 0.4)
      
      if(!DBI::dbExistsTable(conn, "isotopes") | !append){
        RSQLite::dbExecute(conn, 'DROP TABLE IF EXISTS isotopes')
        # create
        RSQLite::dbExecute(conn,
                           "CREATE TABLE isotopes AS
                           SELECT DISTINCT cpd.*
                           FROM db.extended cpd
                           JOIN unfiltered u
                           ON u.baseformula = cpd.baseformula
                           AND u.adduct = cpd.adduct")
      }else{
        # append
        RSQLite::dbExecute(conn,
                           "INSERT INTO isotopes
                           SELECT cpd.*
                           FROM db.extended cpd
                           JOIN unfiltered u
                           ON u.baseformula = cpd.baseformula
                           AND u.adduct = cpd.adduct")
      }
      
      if(!DBI::dbExistsTable(conn, "adducts") | !append){
        RSQLite::dbExecute(conn, 'DROP TABLE IF EXISTS adducts')
        # - - - - - - - - -
        RSQLite::dbExecute(conn,
                           "CREATE TABLE adducts AS
                           SELECT DISTINCT cpd.*
                           FROM db.extended cpd
                           JOIN unfiltered u
                           ON u.baseformula = cpd.baseformula")
      }else{
        # append
        RSQLite::dbExecute(conn,
                           "INSERT INTO adducts
                           SELECT cpd.*
                           FROM db.extended cpd
                           JOIN unfiltered u
                           ON u.baseformula = cpd.baseformula")
      }
      }
    
    if(inshiny) withProgress({func()}) else func()
    
    # - - - save adducts too - - - 
    
    
    results <- DBI::dbGetQuery(conn, "SELECT 
                               b.compoundname as name, 
                               b.baseformula,
                               u.adduct,
                               u.isoprevalence as perciso,
                               u.dppm,
                               b.identifier,
                               b.description, 
                               b.structure 
                               FROM unfiltered u
                               JOIN db.base b
                               ON u.baseformula = b.baseformula
                               AND u.charge = b.charge")
    
    results$perciso <- round(results$perciso, 2)
    results$dppm <- signif(results$dppm, 2)
    
    #colnames(results)[which(colnames(results) == "dppm")] <- "Δppm"
    colnames(results)[which(colnames(results) == "perciso")] <- "%iso"
    
    round_df <- function(df, digits) {
      nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
      
      df[,nums] <- round(df[,nums], digits = digits)
      
      (df)
    }
    
    
    DBI::dbDisconnect(conn)
    
    # - - return - - 
    
    results
    
      }}

score.isos <- function(patdb, method="mscore", inshiny=TRUE){
  
  func <- function(){
    conn <- RSQLite::dbConnect(RSQLite::SQLite(), global$paths$patdb)
    
    if(inshiny) shiny::setProgress(value = 0.6)
    
    table <- RSQLite::dbGetQuery(conn,gsubfn::fn$paste(strwrap(
      "SELECT int.mzmed, iso.baseformula, iso.adduct, iso.fullmz, iso.fullformula, iso.isoprevalence, int.filename, int.intensity
      FROM isotopes iso
      JOIN mzranges rng
      ON iso.fullmz BETWEEN rng.mzmin AND rng.mzmax
      JOIN mzintensities int
      ON int.mzmed BETWEEN rng.mzmin AND rng.mzmax
      WHERE int.filename NOT LIKE 'QC%'"
      , width=10000, simplify=TRUE)))
    
    if(inshiny) shiny::setProgress(value = 0.8)
    
    table <- as.data.table(table[complete.cases(table),])
    
    table <- data.table::setDT(table)[, .(mzmed = mean(mzmed), intensity = sum(intensity)),
                                      by=.(baseformula, fullmz, fullformula, adduct, isoprevalence, filename)]
    p.cpd <- split(x = table, 
                   f = list(table$fullformula))
    
    i <<- 0
    
    res_rows <- pbapply::pblapply(p.cpd, cl = session_cl, function(cpd_tab, method = method, i = i, inshiny = inshiny){
      
      formula = unique(cpd_tab$baseformula)
      adduct = unique(cpd_tab$adduct)
      
      # https://assets.thermofisher.com/TFS-Assets/CMD/Reference-Materials/pp-absoluteidq-qexactive-ms-targeted-metabolic-lipid-metabolomics2017-en.pdf
      
      if(any(cpd_tab$isoprevalence > 99.999999)){
        
        # - - - - - - - - -
        
        sorted <- data.table::as.data.table(unique(cpd_tab[order(cpd_tab$isoprevalence, 
                                                                 decreasing = TRUE),]))
        
        split.by.samp <- split(sorted, 
                               sorted[,"filename"])
        # - - - - - - - - - 
        
        score <- sapply(split.by.samp, function(samp_tab){
          
          samp_tab <- data.table::as.data.table(samp_tab)
          
          if(nrow(samp_tab) == 1){
            res = 0
          }else{
            
            theor_mat <- samp_tab[,c("fullmz", "isoprevalence")]
            theor <- matrix(ncol = nrow(theor_mat), nrow = 2, data = c(theor_mat$fullmz, theor_mat$isoprevalence),byrow = T)
            
            obs_mat <- samp_tab[,c("mzmed", "intensity")]
            obs <- matrix(ncol = nrow(obs_mat), nrow = 2, data = c(obs_mat$mzmed, obs_mat$intensity),byrow = T)
            
            theor[2,] <- theor[2,]/sum(theor[2,])
            obs[2,] <- obs[2,]/sum(obs[2,])
            
            res <- switch(method,
                          mape={
                            actual = obs[2,]
                            theor = theor[2,]
                            deltaSignal = abs(theor - actual)
                            percentageDifference = deltaSignal / actual * 100# Percent by element.
                            # - - -
                            mean(percentageDifference) #Average percentage over all elements.
                          },
                          mscore={
                            InterpretMSSpectrum::mScore(obs=obs, the=theor, dppm = 1, int_prec = 0.225)
                          },
                          sirius={NULL},
                          chisq={
                            test <- chisq.test( obs[2,], p = theor[2,], rescale.p = T)
                            # - - -
                            as.numeric(test$p.value)
                          }
            )
            
          }
          res
        })
        
        #i <<- i + 1
        
        #if(inshiny) shiny::setProgress(value = 0.8 + i * (.2/length(cpd_tab)))
        
        mean_error <- round(mean(score, na.rm=TRUE), digits = 1)
        
        data.table::data.table(baseformula = formula,
                               adduct = adduct,
                               score = as.numeric(mean_error))
      }else{
        data.table::data.table()
      } 
    }, method = method, i = i, inshiny = inshiny)
    
    # - - - - - - - - -
    
    rbindlist(res_rows)
  }
  
  func()
  
  #if(inshiny) func() else func()
  
}

#' @export
get_mzs <- function(baseformula, charge, chosen.db){
  # --- connect to db ---
  req("patdb")
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), global$paths$patdb) # change this to proper var later
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
  
  RSQLite::dbExecute(conn, query.two)
  # isofilter and only in
  query.three <-  strwrap(
    "SELECT DISTINCT iso.mzmed, adduct
    FROM isotopes iso
    GROUP BY iso.adduct"
    #HAVING COUNT(iso.isoprevalence > 99.99999999999) > 0"
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
  join.query <- c()
  
  # --- GET POOL OF ALL DBS ---
  for(i in seq_along(which_dbs)){
    chosen.db <- which_dbs[i]
    
    print(paste("Looking for matches in", chosen.db))
    
    if(is.na(chosen.db)) next
    # --------------------------
    dbshort <- paste0("db", i)
    
    try({
      RSQLite::dbExecute(pat.conn, gsubfn::fn$paste("DETACH $dbshort"))
    })
    
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
                    gsubfn::fn$paste(strwrap("SELECT mz.mzmed as mz, 
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
  RSQLite::dbExecute(pat.conn, gsubfn::fn$paste("DROP TABLE IF EXISTS isotopes"))
  RSQLite::dbExecute(pat.conn, gsubfn::fn$paste("DROP TABLE IF EXISTS results"))
  
  query.one <- gsubfn::fn$paste(strwrap(
    "CREATE TEMP TABLE isotopes AS
    $union",width=10000, simplify=TRUE))
  
  RSQLite::dbExecute(pat.conn, query.one)
  adductfilter <- paste0("WHERE adduct = '", paste(which_adducts, collapse= "' OR adduct = '"), "'")
  idquery <- paste0("iso.", group_by)
  
  query.two <- gsubfn::fn$paste(strwrap(
    "CREATE temp TABLE results AS
    SELECT DISTINCT iso.mz as mz, $idquery as identifier, iso.adduct as adduct
    FROM isotopes iso
    $adductfilter
    GROUP BY mz"
    #HAVING COUNT(iso.isoprevalence > 99.99999999999) > 0"
    , width=10000, simplify=TRUE))
  
  RSQLite::dbExecute(pat.conn, query.two)
  
  summary = if(group_by == "pathway") "sum(abs(i.intensity)) as intensity" else{"sum(i.intensity) as intensity"}
  
  # --- batch ---
  
  query.collect <-  strwrap(gsubfn::fn$paste("select distinct i.filename, 
                                             d.*,
                                             s.*,
                                             b.batch, b.injection,
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
  
  res <- RSQLite::dbGetQuery(pat.conn, query.collect)
  # ---------------------------
  res
}

multimatch <- function(cpd, dbs, searchid="mz", inshiny=T, search_pubchem=F){
  
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), global$paths$patdb) # change this to proper var later
  DBI::dbExecute(conn, "DROP TABLE IF EXISTS unfiltered")
  DBI::dbExecute(conn, "DROP TABLE IF EXISTS isotopes")
  DBI::dbExecute(conn, "DROP TABLE IF EXISTS adducts")
  
  DBI::dbDisconnect(conn)
  
  if(exists("match_table")){
    match_table <<- match_table[,-"score"]
  }
  
  i <<- 1
  
  count = length(dbs)
  
  shiny::withProgress({
    
    match_list <- lapply(dbs, FUN=function(match.table){
      
      shiny::setProgress(i/count)
      
      dbname <- gsub(basename(match.table), pattern = "\\.full\\.db", replacement = "")

      if(dbname == "magicball"){
        # get ppm from db...
        conn <- RSQLite::dbConnect(RSQLite::SQLite(), global$paths$patdb) # change this to proper var later
        A = DBI::dbGetQuery(conn, "SELECT * FROM mzvals WHERE ID = 1")
        B = DBI::dbGetQuery(conn, "SELECT * FROM mzranges WHERE ID = 1")
        DBI::dbDisconnect(conn)
        ppm = round((abs(A$mzmed - B$mzmin) / A$mzmed) * 1e6, digits = 0)
        res <- get_predicted(cpd, ppm = ppm, search_pubchem = search_pubchem)
        
      }else{
        res <- get_matches(cpd, 
                           match.table, 
                           searchid=searchid,
                           inshiny=inshiny,
                           append = if(i == 1) F else T)
      }
      
      #print(res)
      
      i <<- i + 1
      
      if(nrow(res) > 0){
        res = cbind(res, source = c(dbname))
      }
      
      res
    })
    
  })
  
  if(is.null(unlist(match_list))) return(data.table(name = "None",
                                                    description = "Unknown compound",
                                                    source = "None"))
  
  match_table <- (as.data.table(rbindlist(match_list, fill=T))[name != ""])
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
  match_table[,-c("identifier")]
}

get_predicted <- function(mz, 
                          charge = NULL, 
                          ppm = 2, 
                          scanmode = "positive", 
                          checkdb = T, 
                          elements = "CHNOPSNaClKILi",
                          search_pubchem = T){
  
cat(" 
      _...._
    .`      `.
   / ***      \\            MagicBall
  : **         :         is searching...
  :            :        
  \\           /       
    `-.,,,,.-'              
     _(    )_
    )        (
   (          )
    `-......-`
")
  # find which mode we are in
  if(checkdb){
    conn <- RSQLite::dbConnect(RSQLite::SQLite(), global$paths$patdb)
    scanmode <- DBI::dbGetQuery(conn, paste0("SELECT DISTINCT foundinmode FROM mzvals WHERE mzmed = ", mz))[,1]
  }
  
  # get which formulas are possible
  predicted = Rdisop::decomposeMass(as.numeric(mz), 
                                    ppm = ppm,
                                    elements = elements)
  # charged
  charged = which( (predicted$DBE %% 1)  == 0.5 )
  posdbe = which( predicted$DBE > 0 )
  
  candidates = predicted$formula[intersect(charged, posdbe)]
  candidates <- candidates[!is.na(candidates)]
  
  res = lapply(candidates, function(formula){
    
    checked <- check_chemform(isotopes, formula)
    new_formula <- checked[1,]$new_formula
    # switch between positive and negative mode
    check_adducts <- adducts[Ion_mode == scanmode]
    # check which adducts are possible
    adductvars = lapply(1:nrow(check_adducts), function(i){
      row = check_adducts[i,]
      theor_orig_formula = new_formula
      # if there's an adduct, remove it from the original formula
      if(row$Formula_add != FALSE){
        add_possible <- !as.logical(enviPat::check_ded(theor_orig_formula, row$Formula_add))
        if(add_possible){
          theor_orig_formula <- Rdisop::subMolecules(theor_orig_formula, row$Formula_add)$formula
        }else{
          NULL # if not possible, skip this adduct
        }
      }
      if(row$Formula_ded != FALSE){
        theor_orig_formula <- Rdisop::addMolecules(theor_orig_formula, row$Formula_ded)$formula
      }
      if(!is.null(theor_orig_formula)){
        #if(theor_orig_formula == "C9H17N3O6") print("here!!")
        # check params of the original formula (DBE must be integer)
        calc <- Rdisop::getMolecule(theor_orig_formula)
        
        # if(new_formula == "C9H18N3O6"){
        #   print(row$Name)
        #   print(theor_orig_formula)
        #   print(calc$DBE)
        # }
        #print(calc$DBE)
        if(calc$DBE %% 1 != 0){
          NULL
        }else{
          data.table(name = theor_orig_formula, 
                     baseformula = theor_orig_formula, 
                     adduct = row$Name, 
                     `%iso` = 100,
                     structure = NA, 
                     identifier = "???",
                     #                     description = "Predicted possible formula for this m/z value.",
                     source = "magicball")
        }
      }
    })
  })
  res_proc = flattenlist(res)
  tbl = rbindlist(res_proc[!sapply(res_proc, is.null)])
  
  # do a pubchem search
  
  uniques <- unique(tbl$baseformula)
  
  if(search_pubchem){
    i = 0
    count = length(uniques)
    
    shiny::withProgress({
      pc_rows <<- pbapply::pblapply(uniques, function(formula){
        i <<- i + 1
        url = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/", formula, "/cids/JSON")
        description = "No PubChem hits for this predicted formula."
        try({
          pc_res <- jsonlite::read_json(url,simplifyVector = T)
          cids <- pc_res$IdentifierList$CID
          description <- paste0("PubChem found these Compound IDs (check ChemSpider or PubChem): ", paste0(cids, collapse = ", "))
        }) 
        shiny::setProgress(value = i/count)
        row = data.table(baseformula = formula, description = description)
      })
    })
    pc_tbl <- rbindlist(pc_rows)
    tbl <- merge(tbl, pc_tbl, by = "baseformula")
  }else{
    tbl$description = c("Predicted possible formula for this m/z value.")
  }
  tbl
}