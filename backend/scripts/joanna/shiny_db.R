get_exp_vars <- function(from, patdb){
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), patdb) # change this to proper var later
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
  result <- RSQLite::dbGetQuery(conn, "SELECT DISTINCT compoundname as name, baseformula as formula, description as description, charge as charge FROM base")
  # --- result ---
  result
}

#' @export
get_matches <- function(cpd = NA,
                        base.db,
                        search_formula = F,
                        searchid = "mz",
                        inshiny = TRUE,
                        append = FALSE,
                        patdb,db_dir){
  
  # 0. Attach db
  
  
  if(!exists("searchid") && searchid != "mz"){
    
    conn <- RSQLite::dbConnect(RSQLite::SQLite(), base.db)
    
    query <- if(searchid == "pathway"){
      cpd <- gsub(cpd, pattern = "http\\/\\/", replacement = "http:\\/\\/")
      gsubfn::fn$paste(strwrap("SELECT DISTINCT name as Name
                               FROM pathways
                               WHERE identifier = '$cpd'"
                               , width=10000, simplify=TRUE))
    } else{
      gsubfn::fn$paste(strwrap(
        "SELECT DISTINCT compoundname as name, baseformula as formula, identifier, description, structure
        FROM base
        WHERE $searchid = '$cpd'"
        , width=10000, simplify=TRUE))
    }
    res <- RSQLite::dbGetQuery(conn, query)
    return(res)
    
  }else{
    
    conn <- RSQLite::dbConnect(RSQLite::SQLite(), patdb)
    
    query = gsubfn::fn$paste("ATTACH '$base.db' AS base")
    RSQLite::dbExecute(conn, query)
    
    print(query)
    
    ext.db <- file.path(db_dir, 'extended.db')
    
    query = gsubfn::fn$paste("ATTACH '$ext.db' AS full")
    
    RSQLite::dbExecute(conn, query)
    
    print(query)
    
    func <- function(){
      if(!DBI::dbExistsTable(conn, "unfiltered") | !append)
      {
        RSQLite::dbExecute(conn, 'DROP TABLE IF EXISTS unfiltered')
        query <- gsubfn::fn$paste(strwrap(
          "CREATE TABLE unfiltered AS
          SELECT DISTINCT cpd.fullformula,
          cpd.adduct as adduct,
          cpd.isoprevalence as isoprevalence,
          struc.smiles as structure,
          (ABS($cpd - cpd.fullmz) / $cpd)/1e6 AS dppm
          FROM mzvals mz
          JOIN mzranges rng ON rng.ID = mz.ID
          JOIN full.extended cpd indexed by e_idx2
          ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
          AND mz.foundinmode = cpd.foundinmode
          JOIN full.structures struc
          ON cpd.struct_id = struc.id
          WHERE $cpd BETWEEN rng.mzmin AND rng.mzmax",width=10000, simplify=TRUE)) 
        # create
        print(query)
        
        RSQLite::dbExecute(conn, query)
      } else{
        # append
        query <- gsubfn::fn$paste(strwrap(
          "INSERT INTO unfiltered
          SELECT DISTINCT cpd.fullformula,
          cpd.adduct as adduct,
          cpd.isoprevalence as isoprevalence,
          struc.smiles as structure,
          (ABS($cpd - cpd.fullmz) / $cpd)/1e6 AS dppm
          FROM mzvals mz
          JOIN mzranges rng ON rng.ID = mz.ID
          JOIN full.extended cpd indexed by e_idx2
          ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
          AND mz.foundinmode = cpd.foundinmode
          JOIN full.structures struc
          ON cpd.struct_id = struc.id
          WHERE $cpd BETWEEN rng.mzmin AND rng.mzmax",width=10000, simplify=TRUE))
        print(query)
        
        RSQLite::dbExecute(conn, query)
      }
      # if(!DBI::dbExistsTable(conn, "isotopes") | !append){
      #   
      #   RSQLite::dbExecute(conn, 'DROP TABLE IF EXISTS isotopes')
      #   # create
      #   query = strwrap("CREATE TABLE isotopes AS
      #                      SELECT DISTINCT cpd.*
      #                      FROM full.extended cpd indexed by e_idx3
      #                      JOIN unfiltered u
      #                      ON u.baseformula = cpd.baseformula
      #                      AND u.adduct = cpd.adduct",width=10000, simplify=TRUE)
      #   print(query)
      #   RSQLite::dbExecute(conn,
      #                      query)
      # }else{
      #   # append
      #   query = strwrap("INSERT INTO isotopes
      #                    SELECT cpd.*
      #                    FROM full.extended cpd indexed by e_idx3
      #                    JOIN unfiltered u
      #                    ON u.baseformula = cpd.baseformula
      #                    AND u.adduct = cpd.adduct",width=10000, simplify=TRUE)
      #   print(query)
      #   RSQLite::dbExecute(conn,
      #                      query)
      # }
    }
    
    func()
    
    # - - - save adducts too - - -
    
    results <- DBI::dbGetQuery(conn, strwrap("SELECT
                                     b.compoundname as name,
                                     b.baseformula,
                                     u.adduct,
                                     u.isoprevalence as perciso,
                                     u.dppm,
                                     b.identifier,
                                     b.description,
                                     u.structure
                                     FROM unfiltered u
                                     JOIN base.base b
                                     ON u.structure = b.structure",width=10000, simplify=TRUE))
    
    results$perciso <- round(results$perciso, 2)
    results$dppm <- signif(results$dppm, 2)
    
    colnames(results)[which(colnames(results) == "perciso")] <- "%iso"
    
    DBI::dbDisconnect(conn)
    
    # - - return - -
    
    results
    
  }}

get.ppm <- function(patdb){
  # get ppm error retrospectively
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), patdb) # change this to proper var later
  A = DBI::dbGetQuery(conn, "SELECT * FROM mzvals WHERE ID = 1")
  B = DBI::dbGetQuery(conn, "SELECT * FROM mzranges WHERE ID = 1")
  DBI::dbDisconnect(conn)
  ppm = round((abs(A$mzmed - B$mzmin) / A$mzmed) * 1e6, digits = 0)
  # - - -
  ppm
}

score.isos <- function(table, mSet, patdb, method="mscore", inshiny=TRUE, intprec){
  
  func <- function(){
    
    ppm <- get.ppm(patdb=patdb)
    
    conn <- RSQLite::dbConnect(RSQLite::SQLite(), patdb)
    
    if(inshiny) shiny::setProgress(value = 0.2)
    
    RSQLite::dbSendQuery(conn, "DROP TABLE IF EXISTS selection")
    RSQLite::dbWriteTable(conn, "selection", table)
    
    RSQLite::dbSendQuery(conn, "DROP TABLE IF EXISTS selected_isos")
    RSQLite::dbSendQuery(conn,"CREATE TABLE selected_isos AS
      SELECT DISTINCT iso.baseformula, iso.adduct, iso.fullformula, iso.isoprevalence, iso.fullmz
      FROM isotopes iso
      JOIN selection sel
      ON iso.baseformula = sel.baseformula
      AND iso.adduct = sel.adduct")
    
    mzmatches <- RSQLite::dbGetQuery(conn,gsubfn::fn$paste(strwrap(
      "SELECT mz.mzmed, iso.baseformula, iso.adduct, iso.fullformula, iso.isoprevalence, iso.fullmz
      FROM selected_isos iso
      JOIN mzranges rng
      ON iso.fullmz BETWEEN rng.mzmin AND rng.mzmax
      JOIN mzvals mz
      ON rng.ID = mz.ID"
      , width=10000, simplify=TRUE)))
    
    if(inshiny) shiny::setProgress(value = 0.4)
    
    mapper = as.data.table(unique(mzmatches[,2:4]))
    
    mapper = mapper[baseformula %in% table$baseformula & adduct %in% table$adduct]
    
    mzmatches <- mzmatches[,-c(2:3)]
    mzmatches <- as.data.table(unique(mzmatches[complete.cases(mzmatches),]))
    mzmatches$mzmed <- as.factor(mzmatches$mzmed)
    
    sourceTable <- as.data.table(mSet$dataSet$orig, keep.rownames = T)
    longints <- melt(sourceTable, id.var="rn")
    
    colnames(longints) <- c("filename", "mzmed", "intensity")
    
    table <- merge(mzmatches, longints)
    table <- unique(table[,c("fullformula", "isoprevalence", "filename", "mzmed", "fullmz", "intensity")])
    
    # table <- data.table::setDT(table)[, .(mzmed = mean(mzmed), intensity = sum(intensity)),
    #                                   by=.(fullmz, fullformula, isoprevalence, filename)]
    
    p.cpd <- split(x = table,
                   f = list(table$fullformula))
    
    if(inshiny) shiny::setProgress(value = 0.6)
    
    i <<- 0
    
    #cpd_tab <- p.cpd$C9H18N3O6
    
    res_rows <- pbapply::pblapply(p.cpd, cl = session_cl, function(cpd_tab, method = method, i = i, inshiny = inshiny, ppm = ppm, intprec = intprec){
      
      formula = unique(cpd_tab$fullformula)
      #adduct = unique(cpd_tab$adduct)
      
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
                            InterpretMSSpectrum::mScore(obs=obs,
                                                        the=theor,
                                                        dppm = ppm,
                                                        int_prec = intprec)#, int_prec = 0.225)
                          },
                          sirius={NULL},
                          chisq={
                            test <- chisq.test( obs[2,],
                                                p = theor[2,],
                                                rescale.p = T)
                            # - - -
                            as.numeric(test$p.value)
                          }
            )
            
          }
          res
        })
        
        #i <<- i + 1
        #if(inshiny) shiny::setProgress(value = 0.6 + i * (.4/length(cpd_tab)))
        
        mean_error <- round(mean(score, na.rm=TRUE), digits = 1)
        
        data.table::data.table(fullformula = formula,
                               score = as.numeric(mean_error))
      }else{
        data.table::data.table()
      }
    }, method = method, i = i, inshiny = inshiny, ppm = ppm, intprec=intprec)
    
    # - - - - - - - - -
    
    score_tbl <- rbindlist(res_rows)
    merged_tbl <- merge(score_tbl, mapper)
    unique(merged_tbl[,-"fullformula"])
    
  }
  
  tbl <<- func()
  
  tbl
  
  #if(inshiny) func() else func()
  
}

#' @export
get_mzs <- function(baseformula, charge, chosen.db, patdb){
  # --- connect to db ---
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
  
  RSQLite::dbExecute(conn, query.two)
  # isofilter and only in
  query.three <-  strwrap(
    "SELECT DISTINCT iso.mzmed, adduct, iso.isoprevalence
    FROM isotopes iso"
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

multimatch <- function(cpd, dbs, searchid="mz",
                       inshiny=T, calc_adducts = c("M+H", "M-H"),
                       search_pubchem=F, pubchem_detailed=F, patdb, db_dir){
  
  
  patdb <- normalizePath(patdb)
  db_dir <- normalizePath(db_dir)
  dbs <- lapply(dbs, normalizePath)
  
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), patdb) # change this to proper var later
  DBI::dbExecute(conn, "DROP TABLE IF EXISTS unfiltered")
  DBI::dbExecute(conn, "DROP TABLE IF EXISTS isotopes")
  DBI::dbExecute(conn, "DROP TABLE IF EXISTS adducts")
  
  DBI::dbDisconnect(conn)
  
  if(exists("match_table")){
    match_table <<- match_table[,-"score"]
  }
  
  i <<- 1
  
  # check which dbs are even available
  
  avail.dbs <- list.files(db_dir, pattern = "\\.base\\.db",full.names = T)
  
  # fix magicball name
  
  which.magic <- which(grepl(dbs, pattern = "magicball"))
  if(length(which.magic) == 1){
    dbs[[which.magic]] <- "magicball"
  }
  
  # - - - - - - -
  
  keep.dbs <- c(intersect(unlist(dbs), c(avail.dbs, "magicball")))
  
  # - - - - - - -
  
  count = length(keep.dbs)
  
  f <- function(){
    lapply(keep.dbs, FUN=function(match.table){
      
      if(inshiny) shiny::setProgress(i/count)
      
      dbname <- gsub(basename(match.table), pattern = "\\.base\\.db", replacement = "")
      
      if(dbname == "magicball"){
        # get ppm from db...
        
        res <- get_predicted(cpd,
                             ppm = get.ppm(patdb=patdb),
                             search_pubchem = search_pubchem,
                             pubchem_detailed = pubchem_detailed,
                             calc_adducts = calc_adducts, inshiny=inshiny)
        
      }else{
        res <- get_matches(cpd,
                           normalizePath(match.table),
                           searchid=searchid,
                           inshiny=inshiny,
                           append = if(i == 1) F else T,
                           patdb = patdb,
                           db_dir=db_dir)
      }
      
      i <<- i + 1
      
      if(nrow(res) > 0){
        res = cbind(res, source = c(dbname))
      }
      
      res
    })
  }
  
  match_list <- if(inshiny){
    shiny::withProgress({
      f()
    })
  }else{
    f()
  }
  
  if(is.null(unlist(match_list))) return(data.table())
  
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
                          #scanmode = "positive",
                          checkdb = T,
                          elements = "CHNOPSNaClKILi",
                          search_pubchem = T,
                          pubchem_detailed = T,
                          calc_adducts ,#= c("M+H", "M-H"),
                          inshiny=F){
  
  cat("
      _...._
    .`      `.    *             *
   / ***      \\            MagicBall   *
  : **         :    *     is searching...
  :            :
   \\          /
*   `-.,,,,.-'        *
     _(    )_             *
  * )        (                  *
   (          ) *
    `-......-`
")
  # find which mode we are in
  
  per_adduct_results <- pbapply::pblapply(calc_adducts, function(add_name){
    
    print(add_name)
    
    row <- adducts[Name == add_name]
    
    if(!is.na(row$AddEx)){
      add.ele = unlist(strsplit(row$AddEx,
                                split = "\\d*"))
      add.ele <<- add.ele[add.ele != ""]
    }else{
      add.ele <<- c()
    }
    
    settings = list(
      CHNOPS = c("C","H","N","O","P","S"),
      CHNOP = c("C","H","N","O","P"),
      CHNO = c("C","H","N","O"),
      CHO = c("C","H","O")
    )
    
    filter="."
    
    temp_res <- pbapply::pblapply(settings, function(def.ele){
      
      print(def.ele)
      
      add.only.ele <- setdiff(add.ele,
                              def.ele)
      
      total.ele <- unique(c(def.ele,
                            add.ele))
      
      total.ele <- total.ele[total.ele != ""]
      
      # get which formulas are possible
      predicted = Rdisop::decomposeMass(as.numeric(mz),
                                        ppm = ppm,
                                        elements = Rdisop::initializeElements(names = total.ele),
                                        z = row$Charge
      )
      
      # 
      corrected = enviPat::check_chemform(isotopes = isotopes, chemforms = predicted$formula)
      
      deconstructed = data.table::data.table(
        nrC = as.numeric(str_match(corrected$new_formula, pattern = "C(\\d*)")[,2]),
        nrH = as.numeric(str_match(corrected$new_formula, pattern = "H(\\d*)")[,2]),
        nrBr = as.numeric(str_match(corrected$new_formula, pattern = "Br(\\d*)")[,2]),
        nrCl = as.numeric(str_match(corrected$new_formula, pattern = "Cl(\\d*)")[,2]),
        nrF = as.numeric(str_match(corrected$new_formula, pattern = "F(\\d*)")[,2]),
        nrN = as.numeric(str_match(corrected$new_formula, pattern = "N(\\d*)")[,2]),
        nrO = as.numeric(str_match(corrected$new_formula, pattern = "O(\\d*)")[,2]),
        nrP = as.numeric(str_match(corrected$new_formula, pattern = "P(\\d*)")[,2]),
        nrS = as.numeric(str_match(corrected$new_formula, pattern = "S(\\d*)")[,2]),
        nrSi = as.numeric(str_match(corrected$new_formula, pattern = "Si(\\d*)")[,2])  
      )
      
      electron.per.atom <- data.table(
        nrC = 12,
        nrH = 1,
        nrBr = 79,
        nrCl = 35,
        nrF = 19, 
        nrN = 14,
        nrO = 16,
        nrP =31,
        nrS = 32,
        nrSi = 28
      )
      
      deconstructed[is.na(deconstructed)] <- 0
      deconstructed$nrAtoms <- rowSums(deconstructed)
    
      #deconstructed$nrElectrons <- pbapply::pbsapply(1:nrow(deconstructed), function(i){
      #   row = deconstructed[i,]
      #   total = 0
      #   for(atom.type in colnames(row)){
      #     if(atom.type == "nrAtoms") next
      #     nAtom = as.numeric(row[,..atom.type])
      #     nElec = as.numeric(electron.per.atom[1,..atom.type]) * nAtom
      #     total = total + nElec
      #   }
      #   total
      # })
      
      deconstructed$senior <- pbapply::pbsapply(1:nrow(deconstructed), function(i){
        row = deconstructed[i,]
        with(row,{
          if((4*nrC+1*nrH+1*nrBr+1*nrCl+1*nrF+5*nrN+2*nrO+5*nrP+6*nrS+4*nrSi)>=(2*(nrAtoms-1))) TRUE else FALSE
        })
      })
      
      deconstructed$eminus <- pbapply::pbsapply(1:nrow(deconstructed), function(i){
        row = deconstructed[i,]
        with(row,{
          # e- : =4*C4+D4+7*E4+7*F4+7*G4+5*H4+6*I4+5*J4+6*K4+4*L4
          4*nrC+nrH+7*nrBr+7*nrCl+7*nrF+5*nrN+6*nrO+5*nrP+6*nrS+4*nrSi
        })
      })
      
      deconstructed$lewis <- pbapply::pbsapply(1:nrow(deconstructed), function(i){
        row = deconstructed[i,]
        with(row,{
          # 4*C4+ 1*D4 +1*E4 +1*F4 +1*G4 +3*H4 +2*I4 +3*J4 +2*K4 +4*L4
          lewis.sum = 4* nrC + 1* nrH +1* nrBr +1* nrCl +1* nrF +3* nrN +2* nrO +3* nrP +2* nrS +4* nrSi
          if(lewis.sum %% 2 == 0 & eminus >7) TRUE else FALSE
        })
      })
      
      ch_nops_chnops_rows <- pbapply::pblapply(1:nrow(deconstructed), function(i){
        row = deconstructed[i,]
        res= data.table(hc = F, chnops = F, nops = F)
        # T4, V4, X4, X4, Y4 -> H/C N/C	O/C	P/C	S/C
        # chnops
        #IF(AND(T4>=0.2,T4<=3,V4>=0,V4<=2,W4>=0,W4<=1.2,X4>=0,X4<=0.32,Y4>=0,Y4<=0.65),"YES","NO")
        # nops
        # IF(AND(V4>=0,V4<=4,W4>=0,W4<=3,X4>=0,X4<=2,Y4>=0,Y4<=3),"YES","NO"
        try({
        res = with(row,{
            HC = nrH/nrC
            NC = nrN/nrC
            OC = nrO/nrC
            PC = nrP/nrC
            SC = nrS/nrC
            hc = HC %between% c(0, 6)
            chnops = HC %between% c(0.2, 3) & NC %between% c(0, 2) & OC %between% c(0, 1.2) & PC %between% c(0, 0.32) & SC %between% c(0,65)  
            nops = NC %between% c(0, 4) & OC %between% c(0,3) & PC %between% c(0,2) & SC %between% c(0,3)
            data.table(hc = hc, chnops = chnops, nops = nops)
            })  
        })
        res
      })
      
      checks = rbindlist(ch_nops_chnops_rows)
      deconstructed <- cbind(deconstructed, checks)
      
      passes.checks <- with(deconstructed, {
        which(senior & lewis & hc & chnops & nops)
      })
      
      keep.candidates <- predicted$formula[passes.checks]
      
      res = lapply(keep.candidates, function(formula, row){
        
        checked <- check_chemform(isotopes, formula)
        new_formula <- checked[1,]$new_formula
        # check which adducts are possible
        theor_orig_formula = new_formula
        
        # remove last adduct
        if(!is.na(row$AddEx)){
          theor_orig_formula <- Rdisop::subMolecules(theor_orig_formula, row$AddEx)$formula
        }
        if(!is.na(row$RemEx)){
          theor_orig_formula <- Rdisop::addMolecules(theor_orig_formula, row$RemEx)$formula
        }
        
        # remove multiplic
        if(row$xM > 1){
          theor_orig_formula <- multiform.joanna(theor_orig_formula, 1/row$xM)
        }
        
        # remove initial adduct
        if(!is.na(row$AddAt)){
          theor_orig_formula <- Rdisop::subMolecules(theor_orig_formula, row$AddAt)$formula
        }
        if(!is.na(row$RemAt)){
          theor_orig_formula <- Rdisop::addMolecules(theor_orig_formula, row$RemAt)$formula
        }
        data.table(name = theor_orig_formula,
                   baseformula = theor_orig_formula,
                   adduct = row$Name,
                   `%iso` = 100,
                   structure = NA,
                   identifier = "???",
                   # description = "Predicted possible formula for this m/z value.",
                   source = "magicball")
          
      }, row = row)
      
      if(length(res) > 0){
        
        res_proc = flattenlist(res)
        
        tbl <<- rbindlist(res_proc[!sapply(res_proc, is.null)])
        
        if(nrow(tbl) == 0) return(NULL)
        
      }
      
      tbl
      
    })
    
    tbl <- unique(rbindlist(temp_res[!sapply(temp_res, is.null)]))
    uniques <- unique(tbl$baseformula)
    
    if(is.null(uniques)) return(NULL)
    
    if(search_pubchem){
      i = 0
      count = length(uniques)
      if(count == 0) return(NULL)
      
      f <- function(){
        pbapply::pblapply(uniques, function(formula){
          i <<- i + 1
          url = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/", formula, "/cids/JSON")
          description = "No PubChem hits for this predicted formula."
          rows = data.table(name = formula,
                            baseformula = formula,
                            structure = NA,
                            description = description)
          
          try({
            pc_res <- jsonlite::read_json(url,simplifyVector = T)
            cids <- pc_res$IdentifierList$CID
            
            if(pubchem_detailed){ # SLOW!!
              rows <- info_from_cids(cids, charge = if(row$Name == "[M1+.]1+") row$Charge else 0)
            }else{
              rows$description <- paste0("PubChem found these Compound IDs (check ChemSpider or PubChem): ", 
                                         paste0(cids, collapse = ", "))
            }
          })
          if(inshiny) shiny::setProgress(value = i/count)
          rows
        })
      }
      
      pc_rows <- if(inshiny){
        shiny::withProgress({
          f()
        })
      }else{
        f()
      }
      
      if(is.null(pc_rows)) return(NULL)
      
      if(length(pc_rows) > 0){
        pc_tbl <- rbindlist(flattenlist(pc_rows), fill=T)
        
        tbl.merge <- merge(pc_tbl, tbl, by = "baseformula")
        
        checked <- check.chemform.joanna(chemforms = tbl.merge$baseformula, isotopes = isotopes)
        tbl.merge$baseformula <- checked$new_formula
        
        tbl <- tbl.merge[, list(name = name.x, baseformula, adduct, `%iso`, structure = structure.x, description = description)]
      }        
    }else{
      tbl$description = c("Predicted possible formula for this m/z value.")  
    }
    
    tbl_uniq = unique(tbl)
    
    # check SMARTS
    valid.struct = which(!is.na(tbl_uniq$structure))
    
    mols <- lapply(1:nrow(tbl_uniq), function(i) return(NA))
    
    if(length(valid.struct) > 0){
      saveme <<- tbl_uniq
      mols[valid.struct] <- rcdk::parse.smiles(tbl_uniq$structure[valid.struct])#sapply(smiles, rcdk::parse.smiles)
      backtrack <- data.table::data.table(structure = tbl_uniq$structure[valid.struct])

      print(valid.struct)
      
      backtrack_molinfo <- lapply(valid.struct, function(i,row){
          mol = mols[[i]]
          # get molecular formula
          mf = rcdk::get.mol2formula(mol)@string
          # get charge
          ch = rcdk::get.total.formal.charge(mol)
          data.table(baseformula = mf, charge = ch)
      }, row = row)
      backtrack <- cbind(backtrack, rbindlist(backtrack_molinfo))
    }else{
      backtrack <- rbindlist(lapply(1:nrow(tbl_uniq), function(i,row){
        r = tbl_uniq[i,]
        data.table(baseformula = r$baseformula, 
                   structure = c(NA),
                   charge = 0)
      }, row = row))
    }
    
    tbl_fin <- merge(tbl_uniq, backtrack, by = c("baseformula", "structure"))
    print(tbl_fin)
    # - - - 
    tbl_fin
   })
  
  total_tbl <- rbindlist(per_adduct_results[sapply(per_adduct_results, function(x)!is.null(x))], fill=T)
  
  # get more info
  has.struct <- which(!is.na(total_tbl$structure))
  if(length(has.struct) > 0){
    iatoms = rcdk::parse.smiles(total_tbl$structure[has.struct])
    tpsas <- sapply(iatoms, rcdk::get.tpsa)
    total_tbl$tpsa <- c(NA)
    total_tbl$tpsa[has.struct] <- tpsas
  }else{
    total_tbl$tpsa <- c(NA)
  }
  
  # - - - - - - -
  
  total_tbl
  
}

info_from_cids <- function(cids,
                           charge = 0,
                           maxn = 30,
                           write2db=F){
  # structural info
  split.cids = split(cids,
                     ceiling(seq_along(cids) / maxn))
  
  chunk.row.list <- lapply(split.cids, function(cidgroup){
    
    dput(cidgroup)
    url_struct = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                        paste0(cidgroup, collapse=","),
                        "/property/MolecularFormula,CanonicalSMILES,Charge/JSON")
    
    struct_res <- jsonlite::fromJSON(url_struct,
                                     simplifyVector = T)
    
    # CHECK IF ORIGINAL CHARGE IS ZERO
    
    keep.cids <- which(struct_res$PropertyTable$Properties$Charge == charge)
    
    if(length(keep.cids)==0) return(NULL)
    cidgroup = cidgroup[keep.cids]
    
    rows <- struct_res$PropertyTable$Properties[keep.cids,]
    
    # descriptions
    url_desc = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                      paste0(cidgroup, collapse=","),
                      "/description/JSON")
    
    try({
      
      desc_res <- jsonlite::fromJSON(url_desc,
                                     simplifyVector = T)
      
      descs <- as.data.table(desc_res$InformationList$Information)
      
      if("Description" %in% colnames(descs)){
        descs.adj <- descs[, list(name = Title[!is.na(Title)], Description = paste(Description[!is.na(Description)], collapse=" ")), by = CID]
      }else{
        descs.adj <- descs[, list(name = Title[!is.na(Title)], Description = c("No further description available")), by = CID]
      }
      
      descs.adj$Description <- gsub(descs.adj$Description,
                                    pattern = "</?a(|\\s+[^>]+)>",
                                    replacement = "",
                                    perl = T)
      
      if(any(descs.adj$Description == "")){
        descs.adj[Description == ""]$Description <- c("No further description available")
      }
      
      rows <- unique(merge(rows, descs.adj, by.x="CID", by.y="CID"))
      
    })
    
    if(is.null(rows)) return(NULL)
    if(nrow(rows) == 0) return(NULL)
    
    print(rows)
    
    colnames(rows) <- c("identifier", "baseformula", "structure", "charge","name","description")
    
    url_syn = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                     paste0(cidgroup, collapse=","),
                     "/synonyms/JSON")
    
    syn.adj = data.table()
    
    print("..")
    
    try({
      synonyms <- jsonlite::fromJSON(url_syn,
                                     simplifyVector = T)
      syn.adj = synonyms$InformationList$Information
    })
    
    if(nrow(syn.adj) > 0){
      rows.adj <- merge(rows, syn.adj, by.x="identifier", by.y="CID")
      rows.renamed <- lapply(1:nrow(rows.adj), function(i){
        row = rows.adj[i,]
        synonyms = row$Synonym[[1]]
        old.name <- row$name
        new.name <- synonyms[1]
        
        if(is.null(new.name)) new.name <- old.name
        
        desc.names <- synonyms[-1]
        
        row$name <- new.name
        
        row$description <- paste0(paste0("PubChem(", row$identifier, "). ",
                                         "Other names: ",
                                         paste0(if(length(desc.names) > 0) c(old.name, desc.names) else old.name, collapse="; "),
                                         ". "),
                                  row$description)
        row <- as.data.table(row)
        row[,-"Synonym"]
        
      })
      tbl.renamed <- rbindlist(rows.renamed, fill=T)
    }else{
      tbl.renamed <- rows
    }
    
    print("got here")
    tbl.fin <- as.data.table(tbl.renamed)[,-"charge"]
    
    tbl.fin$source <- "PubChem"
    
    # - - - return rows - - -
    
    result <- tbl.fin[,c("name", "baseformula", "structure", "description", "source")]
    result
    
  })
  
  chunk.row.list <<- chunk.row.list
  res <- chunk.row.list[sapply(chunk.row.list, function(x){
    if(is.null(nrow(x))) TRUE else FALSE})]
  
  res
  
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
