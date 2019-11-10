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
    
    row <- adducts[Name == add_name]
    
    add.ele <- if(!is.na(row$AddEx)){
      ele = unlist(strsplit(row$AddEx,
                                split = "\\d*"))
      ele[ele != ""]
    }else{
      c()
    }
    
    settings = list(
      CHNOPS = c("C","H","N","O","P","S"),
      CHNOP = c("C","H","N","O","P"),
      CHNO = c("C","H","N","O"),
      CHO = c("C","H","O")
    )
    
    filter="."
    
    temp_res <- pbapply::pblapply(settings, function(def.ele){
      
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
        nrC = as.numeric(stringr::str_match(corrected$new_formula, pattern = "C(\\d*)")[,2]),
        nrH = as.numeric(stringr::str_match(corrected$new_formula, pattern = "H(\\d*)")[,2]),
        nrBr = as.numeric(stringr::str_match(corrected$new_formula, pattern = "Br(\\d*)")[,2]),
        nrCl = as.numeric(stringr::str_match(corrected$new_formula, pattern = "Cl(\\d*)")[,2]),
        nrF = as.numeric(stringr::str_match(corrected$new_formula, pattern = "F(\\d*)")[,2]),
        nrN = as.numeric(stringr::str_match(corrected$new_formula, pattern = "N(\\d*)")[,2]),
        nrO = as.numeric(stringr::str_match(corrected$new_formula, pattern = "O(\\d*)")[,2]),
        nrP = as.numeric(stringr::str_match(corrected$new_formula, pattern = "P(\\d*)")[,2]),
        nrS = as.numeric(stringr::str_match(corrected$new_formula, pattern = "S(\\d*)")[,2]),
        nrSi = as.numeric(stringr::str_match(corrected$new_formula, pattern = "Si(\\d*)")[,2])  
      )
      
      electron.per.atom <- data.table::data.table(
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
        res= data.table::data.table(hc = F, chnops = F, nops = F)
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
            data.table::data.table(hc = hc, chnops = chnops, nops = nops)
            })  
        })
        res
      })
      
      checks = data.table::rbindlist(ch_nops_chnops_rows)
      deconstructed <- cbind(deconstructed, checks)
      
      passes.checks <- with(deconstructed, {
        which(senior & lewis & hc & chnops & nops)
      })
      
      keep.candidates <- predicted$formula[passes.checks]
      
      res = lapply(keep.candidates, function(formula, row){
        
        checked <- enviPat::check_chemform(isotopes, formula)
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
          theor_orig_formula <- multiform(theor_orig_formula, 1/row$xM)
        }
        
        # remove initial adduct
        if(!is.na(row$AddAt)){
          theor_orig_formula <- Rdisop::subMolecules(theor_orig_formula, row$AddAt)$formula
        }
        if(!is.na(row$RemAt)){
          theor_orig_formula <- Rdisop::addMolecules(theor_orig_formula, row$RemAt)$formula
        }
        data.table::data.table(name = theor_orig_formula,
                   baseformula = theor_orig_formula,
                   adduct = row$Name,
                   `%iso` = 100,
                   structure = NA,
                   identifier = "???",
                   # description = "Predicted possible formula for this m/z value.",
                   source = "magicball")
          
      }, row = row)
      
      if(length(res) > 0){
        
        res_proc = MetaboShiny::flattenlist(res)
        
        tbl <- data.table::rbindlist(res_proc[!sapply(res_proc, is.null)])
        
        if(nrow(tbl) == 0) return(NULL)
        
      }
      
      tbl
      
    })
    
    tbl <- unique(data.table::rbindlist(temp_res[!sapply(temp_res, is.null)]))
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
          rows = data.table::data.table(name = formula,
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
        pc_tbl <- data.table::rbindlist(MetaboShiny::flattenlist(pc_rows), fill=T)
        
        tbl.merge <- merge(pc_tbl, tbl, by = "baseformula")
        
        checked <- check.chemform(chemforms = tbl.merge$baseformula, isotopes = isotopes)
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
      mols[valid.struct] <- rcdk::parse.smiles(tbl_uniq$structure[valid.struct])#sapply(smiles, rcdk::parse.smiles)
      backtrack <- data.table::data.table(structure = tbl_uniq$structure[valid.struct])
      backtrack_molinfo <- lapply(valid.struct, function(i,row){
          mol = mols[[i]]
          # get molecular formula
          mf = rcdk::get.mol2formula(mol)@string
          # get charge
          ch = rcdk::get.total.formal.charge(mol)
          data.table(baseformula = mf, charge = ch)
      }, row = row)
      backtrack <- cbind(backtrack, data.table::rbindlist(backtrack_molinfo))
    }else{
      backtrack <- data.table::rbindlist(lapply(1:nrow(tbl_uniq), function(i,row){
        r = tbl_uniq[i,]
        data.table(baseformula = r$baseformula, 
                   structure = c(NA),
                   charge = 0)
      }, row = row))
    }
    
    tbl_fin <- merge(tbl_uniq, backtrack, by = c("baseformula", "structure"))
    # - - - 
    tbl_fin
   })
  
  total_tbl <- data.table::rbindlist(per_adduct_results[sapply(per_adduct_results, function(x)!is.null(x))], fill=T)
  
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
    
    colnames(rows) <- c("identifier", "baseformula", "structure", "charge","name","description")
    
    url_syn = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                     paste0(cidgroup, collapse=","),
                     "/synonyms/JSON")
    
    syn.adj = data.table()
    
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
      tbl.renamed <- data.table::rbindlist(rows.renamed, fill=T)
    }else{
      tbl.renamed <- rows
    }
    
    tbl.fin <- as.data.table(tbl.renamed)[,-"charge"]
    
    tbl.fin$source <- "PubChem"
    
    # - - - return rows - - -
    
    result <- tbl.fin[,c("name", "baseformula", "structure", "description", "source")]
    result
    
  })
  
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