#' @export
build.base.db <- function(dbname=NA, 
                          outfolder=getOptions("user_options.txt")$db_dir, 
                          cl=FALSE){
  # --- check if user chose something ---
  if(is.na(dbname)) return("~ Please choose one of the following options: HMDB, ChEBI, PubChem, MetaCyc, internal, noise, KEGG! d(>..w.<)b ~")
  # --- make dir if not exists ---
  if(!dir.exists(outfolder)) dir.create(outfolder)
  # --- create .db file ---
  db <- file.path(outfolder, paste0(dbname, ".base.db"))
  if(file.exists(db)) file.remove(db)
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), db)
  RSQLite::dbExecute(conn, statement = "create table base(compoundname text, description text, baseformula text, identifier text, charge text, structure text)")
  # --- create base ---
  function.of.choice <- switch(tolower(dbname),
                               internal = function(dbname){ # BOTH NOISE AND NORMAL
                                 # --- uses csv package ---
                                 int.loc <- file.path(".", "backend","umcfiles", "internal")
                                 # --- non-noise ---
                                 internal.base.db <- read.csv(file.path(int.loc, 
                                                                        "TheoreticalMZ_NegPos_noNoise.txt"),
                                                              sep="\t")[,1:2]
                                 db.formatted <- data.table::data.table(compoundname = internal.base.db$CompoundName,
                                                                        description = c("Internal"),
                                                                        baseformula = internal.base.db$Composition, 
                                                                        identifier=c("Internal"),
                                                                        charge=c(0),
                                                                        structure=c(NA))
                                 
                                 db.rows <- pbapply::pblapply(1:nrow(db.formatted), cl = cl, function(i, db.formatted){
                                   try({
                                     row = db.formatted[i,]
                                     name = row$compoundname
                                     
                                     if(grepl(name, pattern="\\(")){
                                       cpds.1 = gsub(" $", "", gsub("\\(.*$", "", name))
                                       cpds.2 = gsub("[\\(\\)]", "", regmatches(name, gregexpr("\\(.*?\\)", name))[[1]])
                                       if(grepl("_", cpds.2)){
                                         cpds.3 = unlist(strsplit(cpds.2, "_", fixed = FALSE, perl = FALSE, useBytes = FALSE))
                                         cpds = c(cpds.1, cpds.3)
                                       }else{
                                         cpds = c(cpds.1, cpds.2)
                                       }
                                     }else{
                                       cpds = name
                                     }
                                     
                                     structs = NA
                                     charges = 0
                                     
                                     try({
                                       structs = sapply(cpds, function(cpd){
                                         struc = webchem::cir_query(utils::URLencode(cpd))
                                         res = struc[[1]][[1]]
                                         res
                                       })
                                       charges <- sapply(structs, function(smi){
                                         iatom <- rcdk::parse.smiles(smi)[[1]]
                                         charge = rcdk::get.total.charge(iatom)
                                         charge[[1]]
                                       })
                                     })
                                     
                                     if(grepl(name, pattern="\\(IS\\)")) cpds <- name
                                     
                                     print("done!")
                                     result = data.table::data.table(compoundname = c(cpds),
                                                                     description = c("Internal"),
                                                                     baseformula = c(row$baseformula), 
                                                                     identifier=c("Internal"),
                                                                     charge=c(charges),
                                                                     structure=c(structs)
                                     )
                                     # ----------
                                     result
                                   })
                                 }, db.formatted = db.formatted)
                                 
                                 db.formatted = rbindlist(db.rows[!is.na(db.rows)])
                                 
                                 db.formatted$baseformula <-  gsub(x=db.formatted$baseformula, 
                                                                   pattern="( )|(\xa0)|(\xe1)|( .*$)", 
                                                                   replacement="")
                                 db.formatted$baseformula <-  gsub(db.formatted$baseformula, 
                                                                   pattern = "\\((\\d*)([A-Z]*)\\)", 
                                                                   replacement = "\\[\\1\\]\\2")
                                 
                                 checked <- data.table::as.data.table(check.chemform.joanna(isotopes,
                                                                                            db.formatted$baseformula))
                                 db.formatted$baseformula <- checked$new_formula
                                 keep <- checked[warning == FALSE, which = TRUE]
                                 db.formatted <- db.formatted[keep]
                                 # --- write ---
                                 RSQLite::dbWriteTable(conn, "base", db.formatted, append=TRUE)},
                               noise = function(dbname){
                                 int.loc <- file.path(".", "backend","umcfiles", "internal")
                                 # --- noise ---
                                 noise.base.db <- read.csv(file.path(int.loc, 
                                                                     "TheoreticalMZ_NegPos_yesNoise.txt"), 
                                                           sep="\t")
                                 # - - structures? - - -
                                 
                                 names <- noise.base.db$CompoundName
                                 names_stripped <- unique(gsub(pattern = " \\(.*$", replacement = "", x = names))
                                 
                                 # - - - - - - - - - - -
                                 
                                 db.formatted <- data.table::data.table(compoundname = as.character(noise.base.db[,1]),
                                                                        description = as.character(str_match(noise.base.db[,1], "(?<=\\().+?(?=\\))")),
                                                                        baseformula = as.character(noise.base.db[,2]), 
                                                                        identifier=c("Noise"),
                                                                        charge=c(0),
                                                                        structure=c(NA))
                                 RSQLite::dbWriteTable(conn, "base", db.formatted, append=TRUE)
                                 ### do extended table and write to additional db (exception...)
                                 db.full <- file.path(outfolder, paste0(dbname, ".full.db"))
                                 if(file.exists(db.full)) file.remove(db.full)
                                 full.conn <- RSQLite::dbConnect(RSQLite::SQLite(), db.full)
                                 # --- create base table here too ---
                                 RSQLite::dbExecute(full.conn, statement = "create table base(compoundname text, description text, baseformula text, identifier text, charge text, structure text)")
                                 # --- create extended table ---
                                 sql.make.meta <- strwrap("CREATE TABLE IF NOT EXISTS extended(
                                                          baseformula text,
                                                          fullformula text, 
                                                          basemz decimal(30,13), 
                                                          fullmz decimal(30,13), 
                                                          adduct text,
                                                          basecharge int,
                                                          totalcharge int,
                                                          isoprevalence float,
                                                          foundinmode text)", width=10000, simplify=TRUE)
                                 RSQLite::dbExecute(full.conn, sql.make.meta)
                                 # -- reformat noise table ---
                                 db.formatted.list <- pbapply::pblapply(1:nrow(noise.base.db), FUN=function(x){
                                   row <- c(noise.base.db[x,])
                                   print(row[1])
                                   # --- find adduct ---
                                   split.col <- strsplit(unlist(row[1]), split = "\\[|\\]")
                                   compoundname <- split.col[[1]][1]
                                   adduct <- split.col[[1]][2]
                                   # --------------
                                   mz.neg <- as.numeric(row[3])
                                   mz.pos <- as.numeric(row[4])
                                   # --------------
                                   result <- data.table::data.table(
                                     baseformula = as.character(row[2]),
                                     fullformula = as.character(row[2]), 
                                     basemz = c(mz.pos, mz.neg), 
                                     fullmz = c(mz.pos, mz.neg), 
                                     adduct = if(is.na(adduct)) c("M+H", "M-H") else c(as.character(adduct)),
                                     basecharge = c(0),
                                     totalcharge = c(1, -1),
                                     isoprevalence = c(100),
                                     foundinmode = c("positive", "negative")
                                   )
                                   result
                                 })
                                 db.formatted.full <- rbindlist(db.formatted.list[!is.na(db.formatted.list)])
                                 db.formatted.full <- db.formatted.full[basemz > 0] # remove not-found ions
                                 # --- write to db ---
                                 db.formatted$compoundname <- as.character(strsplit(db.formatted$compoundname, 
                                                                                    split=", \\[.*\\][\\+-]"))
                                 RSQLite::dbWriteTable(full.conn, 
                                                       "base", 
                                                       db.formatted, 
                                                       append=TRUE)
                                 
                                 RSQLite::dbWriteTable(full.conn, 
                                                       "extended", 
                                                       db.formatted.full, 
                                                       append=TRUE)
                                 # --- index ---
                                 RSQLite::dbExecute(full.conn, "create index b_idx1 on base(baseformula, charge)")
                                 RSQLite::dbExecute(full.conn, "create index e_idx1 on extended(baseformula, basecharge)")
                                 RSQLite::dbExecute(full.conn, "create index e_idx2 on extended(fullmz, foundinmode)")# -------------
                                 RSQLite::dbDisconnect(full.conn)
                               },
                               hmdb = function(dbname){
                                 print("Downloading XML database...")
                                 file.url <- "http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip"
                                 # ----
                                 #base.loc <- getOptions("user_options.txt")$db_dir
                                 base.loc <- file.path(getOptions("user_options.txt")$db_dir, "hmdb_source")
                                 if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
                                 zip.file <- file.path(base.loc, "HMDB.zip")
                                 utils::download.file(file.url, zip.file,mode = "w")
                                 utils::unzip(zip.file, exdir = base.loc)
                                 # --- go through xml ---
                                 print("Parsing XML...")
                                 data <- XML::xmlParse(file.path(base.loc,"hmdb_metabolites.xml"), useInternalNodes = T)
                                 # --- xpath magic! :-) --- [NOTE: leaves out anything that doens't have formal charge available, another option would be to default to zero]
                                 compoundnames <- XML::getNodeSet(data, 
                                                                  "/*/pf:metabolite[pf:predicted_properties/pf:property[pf:kind='formal_charge']]/pf:name", 
                                                                  c(pf = "http://www.hmdb.ca"))
                                 formulae <- XML::getNodeSet(data, 
                                                             "/*/pf:metabolite[pf:predicted_properties/pf:property[pf:kind='formal_charge']]/pf:chemical_formula", 
                                                             c(pf = "http://www.hmdb.ca"))
                                 identifiers <- XML::getNodeSet(data, 
                                                                "/*/pf:metabolite[pf:predicted_properties/pf:property[pf:kind='formal_charge']]/pf:accession", 
                                                                c(pf = "http://www.hmdb.ca"))
                                 charges <- XML::getNodeSet(data, 
                                                            "/*/pf:metabolite/pf:predicted_properties/pf:property[pf:kind='formal_charge']/pf:value", 
                                                            c(pf = "http://www.hmdb.ca"))
                                 description <- XML::getNodeSet(data, 
                                                                "/*/pf:metabolite[pf:predicted_properties/pf:property[pf:kind='formal_charge']]/pf:description", 
                                                                c(pf = "http://www.hmdb.ca"))
                                 structures <- XML::getNodeSet(data, 
                                                               "/*/pf:metabolite[pf:predicted_properties/pf:property[pf:kind='formal_charge']]/pf:smiles", 
                                                               c(pf = "http://www.hmdb.ca"))
                                 # --- make nice big table ---
                                 db.formatted <- data.table::data.table(
                                   compoundname = sapply(compoundnames, XML::xmlValue),
                                   description = sapply(description, XML::xmlValue),
                                   baseformula = sapply(formulae, XML::xmlValue),
                                   identifier =  gsub(sapply(identifiers, XML::xmlValue), pattern = "(HMDB0*)", replacement = ""),
                                   charge = sapply(charges, XML::xmlValue),
                                   structure = sapply(structures, XML::xmlValue)
                                 )
                                 # --- check formulae ---
                                 checked <- data.table::as.data.table(check.chemform.joanna(isotopes,
                                                                                            db.formatted$baseformula))
                                 db.formatted$baseformula <- checked$new_formula
                                 keep <- checked[warning == FALSE, which = TRUE]
                                 db.formatted <- db.formatted[keep]
                                 # ----------------------
                                 RSQLite::dbWriteTable(conn, "base", db.formatted, append=TRUE)
                               },
                               chebi = function(dbname, ...){
                                 
                                 # ----------------------
                                 
                                 db.full <- data.table::as.data.table(download.chebi.joanna(release = "latest", 
                                                                                            woAssociations = FALSE))
                                 db.formatted <- unique(db.full[, list(compoundname = ChEBI, 
                                                                       description = DEFINITION,
                                                                       baseformula = FORMULA,
                                                                       identifier = ID, 
                                                                       charge = gsub(CHARGE,pattern = "$\\+\\d", replacement = ""),
                                                                       structure = toupper(STRUCTURE)
                                 )])
                                 
                                 # --- check formulae ---
                                 checked <- data.table::as.data.table(check.chemform.joanna(isotopes,
                                                                                            db.formatted$baseformula))
                                 db.formatted$baseformula <- checked$new_formula
                                 keep <- checked[warning == FALSE, which = TRUE]
                                 db.formatted <- db.formatted[keep]
                                 # ----------------------
                                 RSQLite::dbWriteTable(conn, "base", db.formatted, append=TRUE)
                               }, 
                               wikipathways = function(dbname){
                                 chebi.loc <- file.path(getOptions("user_options.txt")$db_dir, "chebi.full.db")
                                 # ---------------------------------------------------
                                 chebi <- SPARQL::SPARQL(url="http://sparql.wikipathways.org/",
                                                         query='prefix wp:      <http://vocabularies.wikipathways.org/wp#>
                                                         prefix rdfs:    <http://www.w3.org/2000/01/rdf-schema#>
                                                         prefix dcterms: <http://purl.org/dc/terms/>
                                                         prefix xsd:     <http://www.w3.org/2001/XMLSchema#>
                                                         PREFIX wdt: <http://www.wikidata.org/prop/direct/>
                                                         
                                                         select  ?mb 
                                                         (group_concat(distinct str(?labelLit);separator=", ") as ?label )
                                                         ?idurl as ?csid 
                                                         (group_concat(distinct ?pwTitle;separator=", ") as ?description)
                                                         ?pathway
                                                         where {
                                                         ?mb a wp:Metabolite ;
                                                         rdfs:label ?labelLit ;
                                                         wp:bdbChEBI ?idurl ;
                                                         dcterms:isPartOf ?pathway .
                                                         ?pathway a wp:Pathway ;
                                                         dc:title ?pwTitle .
                                                         FILTER (BOUND(?idurl))
                                                         }
                                                         GROUP BY ?mb ?wp ?idurl ?pathway')
                                 chebi.ids <- gsub(chebi$results$csid, pattern = ".*:|>", replacement = "")
                                 conn.chebi <- RSQLite::dbConnect(RSQLite::SQLite(), chebi.loc)
                                 chebi.join.table <- data.table::data.table(identifier = chebi.ids,
                                                                            description = chebi$results$description,
                                                                            widentifier = chebi$results$mb,
                                                                            pathway = chebi$results$pathway)
                                 RSQLite::dbWriteTable(conn.chebi, "wikipathways", chebi.join.table, overwrite=TRUE)
                                 db.formatted <- RSQLite::dbGetQuery(conn.chebi, "SELECT DISTINCT  b.compoundname, 
                                                                     w.description,
                                                                     b.baseformula, 
                                                                     w.widentifier as identifier, 
                                                                     b.charge,
                                                                     w.pathway,
                                                                     b.structure
                                                                     FROM base b
                                                                     JOIN wikipathways w
                                                                     ON b.identifier = w.identifier")
                                 # --- get pathway info ---
                                 # -- R PACKAGE EXISTS: source("https://bioconductor.org/biocLite.R"); biocLite("rWikiPathways"); --
                                 sparql.pathways <- SPARQL::SPARQL(url="http://sparql.wikipathways.org/",
                                                                   query='PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
                                                                   PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
                                                                   PREFIX dc:  <http://purl.org/dc/elements/1.1/>
                                                                   PREFIX foaf: <http://xmlns.com/foaf/0.1/>
                                                                   PREFIX schema: <http://schema.org/>
                                                                   PREFIX wp:      <http://vocabularies.wikipathways.org/wp#>
                                                                   PREFIX dcterms:  <http://purl.org/dc/terms/>
                                                                   
                                                                   SELECT DISTINCT str(?titleLit) as ?name ?identifier ?ontology
                                                                   WHERE {
                                                                   ?pathway dc:title ?titleLit .
                                                                   ?pathway dc:identifier ?identifier .
                                                                   OPTIONAL {?pathway wp:ontologyTag ?ontology .}
                                                                   } ')
                                 db.pathways <- sparql.pathways$results
                                 db.formatted$identifier <- gsub(db.formatted$identifier, pattern = "<|>", replacement = "")
                                 db.formatted$pathway <- gsub(db.formatted$pathway, pattern = "<|_.*", replacement = "")
                                 db.pathways$identifier <- gsub(db.pathways$identifier, pattern = "<|>", replacement = "")
                                 # ---------------------------------
                                 RSQLite::dbWriteTable(conn, "base", db.formatted, overwrite=TRUE)
                                 RSQLite::dbWriteTable(conn, "pathways", db.pathways, overwrite=TRUE)
                                 RSQLite::dbRemoveTable(conn.chebi, "wikipathways")
                                 # ---------------------------------
                                 # SOMETHING FOR THE UNIDENTIFIED ONES HERE...
                                 # no.chebi <- SPARQL(url="http://sparql.wikipathways.org/",
                                 #                    query='prefix wp:      <http://vocabularies.wikipathways.org/wp#>
                                 #                          prefix rdfs:    <http://www.w3.org/2000/01/rdf-schema#>
                                 #                          prefix dcterms: <http://purl.org/dc/terms/>
                                 #                          prefix xsd:     <http://www.w3.org/2001/XMLSchema#>
                                 #                          PREFIX wdt: <http://www.wikidata.org/prop/direct/>
                                 #                          
                                 #                          select distinct ?mb str(?labelLit) as ?label ?idurl as ?csid ?pwTitle
                                 #                          where {
                                 #                          ?mb a wp:Metabolite ;
                                 #                          rdfs:label ?labelLit ;
                                 #                          wp:bdbChEBI ?idurl ;
                                 #                          dcterms:isPartOf ?pathway .
                                 #                          ?pathway a wp:Pathway ;
                                 #                          dc:title ?pwTitle .
                                 #                          FILTER (!BOUND(?idurl))
                                 #                          }')
                                 # ---------------------------------
                               },
                               smpdb = function(dbname){
                                 # ---------------
                                 file.url <- "http://smpdb.ca/downloads/smpdb_metabolites.csv.zip"
                                 # ----
                                 base.loc <- file.path(getOptions("user_options.txt")$db_dir, "smpdb_source")
                                 if(!dir.exists(base.loc)) dir.create(base.loc)
                                 zip.file <- file.path(base.loc, "SMPDB.zip")
                                 utils::download.file(file.url, zip.file)
                                 utils::unzip(zip.file, exdir = base.loc)
                                 # -------------------------------
                                 
                                 smpdb.paths <- list.files(path = base.loc, pattern = "\\.csv$", full.names = T)
                                 smpdb.tabs <- pbapply::pblapply(smpdb.paths, fread)
                                 smpdb.tab <- unique(rbindlist(smpdb.tabs, fill=TRUE))
                                 
                                 # --- get charges ---
                                 
                                 smi_mapper <- unique(smpdb.tab[,c("Metabolite ID", "SMILES")])
                                 
                                 smi_mapper$charge <- pbapply::pbsapply(smi_mapper$SMILES, cl = 0, function(smi){
                                   charge = 0 # set default charge to zero
                                   try({
                                     iatom <- rcdk::parse.smiles(smi)[[1]]
                                     charge = rcdk::get.total.charge(iatom)
                                   })
                                   #print(charge)
                                   charge
                                 })
                                 
                                 db.formatted <- unique(data.table::data.table(compoundname = smpdb.tab$`Metabolite Name`,
                                                                               description = c("See HMDB for more info :-)"),
                                                                               baseformula = smpdb.tab$Formula,
                                                                               identifier = smpdb.tab$`Metabolite ID`,
                                                                               structure = smpdb.tab$SMILES))
                                 db.formatted <- merge(smi_mapper, db.formatted, by.x = "Metabolite ID", by.y = "identifier")
                                 colnames(db.formatted)[which(colnames(db.formatted) == 'Metabolite ID')] <- "identifier"
                                 db.mol2pathway <- unique(data.table::data.table(identifier = smpdb.tab$`Metabolite ID`,
                                                                          pathway = smpdb.tab$`SMPDB ID`))
                                 db.pathways <- unique(data.table::data.table(name = smpdb.tab$`Pathway Name`,
                                                                              identifier = smpdb.tab$`SMPDB ID`,
                                                                              subject = smpdb.tab$`Pathway Subject`))
                                 checked <- data.table::as.data.table(check.chemform.joanna(isotopes,
                                                                                            db.formatted$baseformula))
                                 db.formatted$baseformula <- checked$new_formula
                                 keep <- checked[warning == FALSE, which = TRUE]
                                 db.formatted <- db.formatted[keep]
                                 # --- create ---
                                 RSQLite::dbWriteTable(conn, "base", db.formatted, overwrite=TRUE)
                                 RSQLite::dbWriteTable(conn, "mol2pathway", db.mol2pathway, overwrite=TRUE)
                                 RSQLite::dbWriteTable(conn, "pathways", db.pathways, overwrite=TRUE)
                               },
                               kegg = function(dbname){
                                 # ---------------
                                 batches <- split(0:2000, ceiling(seq_along(0:2000)/100))
                                 cpds <- pbapply::pblapply(batches, cl=session_cl, FUN=function(batch){
                                   names(KEGGREST::keggFind("compound", batch, "mol_weight"))
                                 })
                                 cpd.ids <- Reduce(c, cpds)
                                 id.batches <- split(cpd.ids, ceiling(seq_along(cpd.ids)/10))
                                 # ---------------
                                 # KEGGREST:keggGet("C00032") # cpd
                                 # KEGGREST:keggGet("map00860") # pathway
                                 # KEGGREST:keggGet("M00121") # module
                                 # KEGGREST:keggGet("H00201") # disease
                                 # --- GET COMPOUNDS ---
                                 
                                 #cl = parallel::makeCluster(3)
                                 parallel::clusterExport(session_cl, c("kegg.charge", "rbindlist"))
                                 
                                 kegg.cpd.list <- pbapply::pblapply(id.batches, cl=session_cl, FUN=function(batch){
                                   rest.result <- KEGGREST::keggGet(batch)
                                   # ---------------------------
                                   base.list <- lapply(rest.result, FUN=function(cpd){
                                     cpd$NAME_FILT <- gsub(cpd$NAME, pattern = ";", replacement = "")
                                     data.table::data.table(compoundname = c(paste(cpd$NAME_FILT, collapse=", ")),
                                                            description = c(if("BRITE" %in% names(cpd)) paste(cpd$BRITE, collapse=", ") else{NA}),
                                                            baseformula = c(cpd$FORMULA),
                                                            identifier = c(cpd$ENTRY),
                                                            widentifier = c(paste(cpd$DBLINKS, collapse=";")),
                                                            charge = c(kegg.charge(cpd$ATOM)),
                                                            structure = c(NA),
                                                            pathway = if("PATHWAY" %in% names(cpd)) names(cpd$PATHWAY) else{NA}
                                     )
                                   })
                                   rbindlist(base.list)
                                 })
                                 
                                 # webchem::cir_query("glycine")
                                 
                                 db.cpds <- rbindlist(kegg.cpd.list)
                                 
                                 # --- get most from chebi ---
                                 #
                                 # w.chebi.rows <- which(grepl(db.cpds$widentifier, pattern = "ChEBI")) 
                                 #
                                 # --- and rest through pubchem ---
                                 #
                                 # no.chebi <- db.cpds[-w.chebi.rows,]
                                 # length(which(grepl(no.chebi$widentifier, pattern = "PubChem")))
                                 # --------------------------------
                                 
                                 db.cpds$baseformula <- gsub(pattern = " |\\.", replacement = "", db.cpds$baseformula) # fix salt and spacing issues                                 
                                 
                                 # -------------------------------
                                 checked <- data.table::as.data.table(check.chemform.joanna(isotopes,
                                                                                            db.cpds$baseformula))
                                 db.cpds$baseformula <- checked$new_formula
                                 keep <- checked[warning == FALSE, which = TRUE]
                                 db.cpds <- db.cpds[keep]
                                 # --- GET PATHWAYS ---
                                 pathways <- unique(db.cpds$pathway)
                                 pw.batches <- split(pathways, ceiling(seq_along(pathways)/10))
                                 kegg.pw.list <- pbapply::pblapply(pw.batches, cl=session_cl, FUN=function(batch){
                                   rest.result <- KEGGREST::keggGet(batch)
                                   # ---------------------------
                                   base.list <- lapply(rest.result, FUN=function(pw){
                                     data.table::data.table(identifier = as.character(pw$ENTRY),
                                                            name = as.character(pw$NAME),
                                                            class = if("CLASS" %in% names(pw)) pw$CLASS else{NA},
                                                            module = if("MODULE" %in% names(pw)) pw$MODULE else{NA},
                                                            disease = if("DISEASE" %in% names(pw)) as.character(pw$DISEASE) else{NA})
                                   })
                                   rbindlist(base.list)
                                 })
                                 db.pathways <- rbindlist(kegg.pw.list)
                                 # --- GET MODULES & DISEASES? (necessary for future?) ---
                                 # modules <- unique(db.cpds$module)
                                 # md.batches <- split(modules, ceiling(seq_along(modules)/10))
                                 # RSQLite::dbWriteTable(conn, "modules", db.modules, overwrite=TRUE)
                                 # -------------------------------------------------------
                                 RSQLite::dbWriteTable(conn, "base", db.cpds, overwrite=TRUE)
                                 RSQLite::dbWriteTable(conn, "pathways", db.pathways, overwrite=TRUE)
                               },
                               metacyc = function(dbname, ...){
                                 # NOTE: Requires downloading this SmartTable as delimited file: https://metacyc.org/group?id=biocyc17-31223-3729417004
                                 # May need to remake smartTable if anything on the website changes unfortunately
                                 # TODO: download file directly from link, will need a javascript. Maybe Rselenium??
                                 source.file = "./backend/db/metacyc_source/All_compounds_of_MetaCyc.txt"
                                 if(!file.exists(source.file)){
                                   message("Please download SmartTable from 'https://metacyc.org/group?id=biocyc17-31223-3729417004' as .txt and save in the backend/db/metacyc_source folder.")
                                   return(NULL)
                                 }
                                 metacyc.raw = read.table(source.file, fill = T, header = T)
                                 #metacyc.raw = fread(source.file, sep = "\t", quote = '\"', header = T, fill=T)
                                 charges = pbapply::pbsapply(metacyc.raw$SMILES, cl=session_cl, FUN=function(smile){
                                   m <- rcdk::parse.smiles(smile)
                                   #print(m)
                                   ## perform operations on this molecule
                                   try({
                                     charge = rcdk::get.total.formal.charge(m[[1]])
                                     # --- return ---
                                     return(charge)
                                   })
                                 }) 
                                 charges[grep("Error", charges)] <- 0
                                 charges = as.numeric(charges)
                                 
                                 compounds <- pbapply::pbsapply(metacyc.raw$Compound, cl=session_cl, FUN=function(pw){
                                   pw <- pw[pw != " // "]
                                   pw <- gsub(pw, pattern = "&", replacement="")
                                   pw <- gsub(pw, pattern = ";", replacement="")
                                   res <- gsub(pw, pattern = "<((i|\\/i)|sub)>|\\/|\\|", replacement = "",perl = T)
                                   paste0(res, collapse=" --- ")
                                 })
                                 
                                 pathways <- pbapply::pbsapply(metacyc.raw$Pathways.of.compound, cl=session_cl, FUN=function(pw){
                                   pw <- unlist(strsplit(pw, split = '\\"'))
                                   pw <- pw[pw != " // "]
                                   pw <- gsub(pw, pattern = "&", replacement="")
                                   pw <- gsub(pw, pattern = ";", replacement="")
                                   res <- gsub(pw, pattern = "<((i|\\/i)|sub)>|\\/|\\|", replacement = "",perl = T)
                                   paste0(res, collapse=" --- ")
                                 })
                                 
                                 # give the pathwys some identifiers
                                 uniq.pws <- unique(unlist(pbapply::pbsapply(pathways, cl=session_cl, function(x){unlist(strsplit(x, split = " --- "))})))
                                 
                                 db.pathways <- data.table::data.table(name = uniq.pws,
                                                                       identifier = paste0("METACYC_PW_", 1:length(uniq.pws)))
                                 
                                 db.formatted <- data.table(compoundname = compounds, 
                                                            description = metacyc.raw$Summary,metacyc.raw,
                                                            baseformula = metacyc.raw$Chemical.Formula,
                                                            identifier = paste0("METACYC_CP_", 1:length(charges)), 
                                                            charge = charges,
                                                            structure = metacyc.raw$SMILES,
                                                            pathway = pathways
                                 )
                                 # --- check formulae ---
                                 checked <- data.table::as.data.table(check.chemform.joanna(isotopes,
                                                                                            db.formatted$baseformula))
                                 db.formatted$baseformula <- checked$new_formula
                                 keep <- checked[warning == FALSE, which = TRUE]
                                 db.formatted <- db.formatted[keep]
                                 # --- fix the multiple pathway thingy ---
                                 db.fixed.rows <- pbapply::pblapply(1:nrow(db.formatted), cl = session_cl, FUN=function(i, db.cpd = db.cpd, db.pw = db.pw){
                                   row <- db.cpd[i,]
                                   pathways <- unlist(strsplit(row$pathway, split = " --- "))
                                   if(length(pathways) == 0) pathways <- c(1)
                                   pw.ids <- sapply(pathways, function(pw) db.pw[name == pw]$identifier)
                                   res <- data.table(compoundname = rep(row$compoundname, length(pathways)),
                                                     description = rep(row$description, length(pathways)),
                                                     baseformula = rep(row$baseformula, length(pathways)),
                                                     identifier = rep(row$identifier, length(pathways)),
                                                     charge = rep(row$charge, length(pathways)),
                                                     structure = rep(row$structure, length(pathways)),
                                                     pathway = pw.ids)
                                   # - - - - - - - - - -
                                   res
                                 }, db.cpd=db.formatted, db.pw = db.pathways)
                                 
                                 db.formatted <- rbindlist(db.fixed.rows)
                                 
                                 #db <- file.path(outfolder, paste0(dbname, ".base.db"))
                                 #if(file.exists(db)) file.remove(db)
                                 #conn <- RSQLite::dbConnect(RSQLite::SQLite(), db)
                                 #DBI::dbDisconnect(conn)
                                 db.formatted$pathway <- as.character(db.formatted$pathway)
                                 
                                 # ---------------------------------------
                                 
                                 RSQLite::dbWriteTable(conn, "base", db.formatted, overwrite=TRUE)
                                 RSQLite::dbWriteTable(conn, "pathways", db.pathways, overwrite=TRUE)
                               },
                               pubchem = function(dbname, ...){
                                 # --- create working space ---
                                 baseLoc <- file.path(getOptions("user_options.txt")$db_dir, "pubchem_source")
                                 sdf.loc <- file.path(baseLoc, "sdf")
                                 csv.loc <- file.path(baseLoc, "csv")
                                 # --- check the user's determination ---
                                 # continue <- readline("You chose PubChem, this will take a while and requires at least 30gb of space.Continue? (yes/no): ")
                                 # if(continue == "no" | continue == "n") return(NA)
                                 # --- download the 1 million sdf files ---
                                 if(!dir.exists(sdf.loc)) dir.create(sdf.loc)
                                 if(!dir.exists(csv.loc)) dir.create(csv.loc)
                                 print("Downloading SDF files...")
                                 folder.url = "ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/"
                                 ftp.handle = new_handle(dirlistonly=TRUE)
                                 ftp.conn = curl(folder.url, "r", ftp.handle)
                                 file.table = read.table(ftp.conn, stringsAsFactors=TRUE, fill=TRUE)
                                 close(ftp.conn)
                                 file.urls <- paste0(folder.url, file.table[,1])
                                 # ---------------------------------------
                                 counter  = 0
                                 max_counter = length(file.urls)
                                 # --- start downloady ---
                                 pool <- new_pool()
                                 pbapply::pbsapply(file.urls, cl=session_cl, FUN=function(url){
                                   print(url)
                                   counter <<- counter + 1
                                   print(paste(counter, length(file.urls), sep=" of "))
                                   fn <- file.path(sdf.loc, basename(url))
                                   # -------------------------------
                                   if(file.exists(fn)) return(NA)
                                   # -------------------------------
                                   print(fn)
                                   download.file(url = url, destfile = fn, method = "auto")
                                 })
                                 # lotsa files, need tiem.
                                 # ------------------------------
                                 print("Converting SDF files to tab delimited matrices...")
                                 sdf.files <- list.files(path = sdf.loc, pattern = "\\.sdf\\.gz$")
                                 # -----------------------
                                 pbapply::pbsapply(cl=session_cl, sdf.files, FUN=function(sdf.file){
                                   input <- file.path(sdf.loc, sdf.file)
                                   output <- file.path(csv.loc, gsub("\\.sdf.gz$", "\\.csv", sdf.file))
                                   # -------------------------------
                                   if(file.exists(output)) return(NA)
                                   # -------------------------------
                                   sdfStream.joanna(input=input,
                                                    output=output,
                                                    fct = function(sdfset, test){
                                                      valid <- validSDF(sdfset)
                                                      sdfset <- sdfset[valid]
                                                      print(head(datablock(sdfset)))
                                                      blockmatrix <- datablock2ma(datablock(sdfset)) # Converts data block to matrix
                                                      # --------------
                                                      db.formatted <- data.table::data.table(
                                                        compoundname = blockmatrix[, if("PUBCHEM_IUPAC_TRADITIONAL_NAME" %not in% colnames(blockmatrix)) "PUBCHEM_MOLECULAR_FORMULA" else("PUBCHEM_IUPAC_TRADITIONAL_NAME")],
                                                        description = c("PubChem"),
                                                        baseformula = gsub(x = blockmatrix[, "PUBCHEM_MOLECULAR_FORMULA"], 
                                                                           pattern="[\\+\\-]\\d*$", 
                                                                           replacement=""),
                                                        identifier = as.numeric(blockmatrix[, "PUBCHEM_COMPOUND_CID"]),
                                                        charge = blockmatrix[,"PUBCHEM_TOTAL_CHARGE"]
                                                      )
                                                      # --- check formulae ---
                                                      checked <- data.table::as.data.table(check.chemform.joanna(isotopes,
                                                                                                                 db.formatted$baseformula))
                                                      db.formatted$baseformula <- checked$new_formula
                                                      keep <- checked[warning == FALSE, which = TRUE]
                                                      db.keep <- db.formatted[keep]
                                                      db.keep},
                                                    append = FALSE,
                                                    silent = FALSE,
                                                    Nlines = 100000 )
                                 })
                                 # --- assemble ---
                                 csv.files <- list.files(path = csv.loc, pattern = "\\.csv$", full.names = TRUE)
                                 print("Assembling and putting in db file...")
                                 sapply(csv.files, FUN=function(file){
                                   first.row <- read.csv(file, nrows=3, header=TRUE, sep="\t")
                                   if("description" %not in% colnames(first.row)){print("NOPE!"); file.remove(file); return(NULL)}
                                   read.csv.sql(file, sep="\t", sql = c(gsubfn::fn$paste("insert into base select compoundname, description, baseformula, identifier, charge from file")), dbname = db)
                                 })
                                 print("Done!")
                                 # --- cleanup --- downloading takes forever, would not recommend unless compressing... ---
                                 # remove.residuals <- readline("Done! Remove downloaded SDF / XLS files? (yes/no): ")
                                 # if(remove.residuals == "yes" | remove.residuals == "y") unlink(sdf.loc, recursive=TRUE) else return("Alrighty, files kept.")
                               },
                               metabolights = function(dbname, ...){
                                 require(webchem)
                                 require(rcdk)
                                 require(curl)
                                 
                                 all_ids <- read.table("ftp://ftp.ebi.ac.uk/pub/databases/metabolights/eb-eye/unichem.tsv", sep="\t")
                                 #all_ids <- all_ids
                                 colnames(all_ids) <- c("identifier", "inchi", "inchikey")
                                 #query <- all_ids$inchi[[1]]
                                 
                                 #token = "b7e19fad-46bd-48d3-80c1-765140acace1"
                                 
                                 # w.csid.rows <- pbapply::pblapply(1:nrow(all_ids), cl=0, function(i, all_ids){
                                 #   row = all_ids[i,]
                                 #   inchi = row$inchikey
                                 #   print(i/nrow(all_ids)*100)
                                 #   try({
                                 #     csid = webchem::cs_convert(inchi, from = c("inchikey"),
                                 #                                to = c("csid"))[[1]]
                                 #     info = webchem::cs_extcompinfo(csid,token=token)
                                 #     smi = info$smiles
                                 #     iatom <- rcdk::parse.smiles(smi)[[1]]
                                 #     charge = rcdk::get.total.charge(iatom)
                                 #     formula = rcdk::get.mol2formula(iatom, charge = charge)@string
                                 #     # - - - - -
                                 #     data.table::data.table(compoundname = info$common_name,
                                 #                            baseformula = formula,
                                 #                            description = "MetaboLights compound",
                                 #                            identifier = row$identifier,
                                 #                            structure = info$smiles)
                                 #   })
                                 # }, all_ids = all_ids)
                                 # 
                                 # # metabs = jsonlite::fromJSON(txt = "http://ftp.ebi.ac.uk/pub/databases/metabolights/compounds/MetabolitesReport.json",
                                 # #                             simplifyDataFrame = TRUE)
                                 # # pathways = jsonlite::fromJSON(txt = "http://ftp.ebi.ac.uk/pub/databases/metabolights/compounds/reactome.json", 
                                 # #                               simplifyDataFrame = TRUE)
                                 # url = "ftp://ftp.ebi.ac.uk/pub/databases/metabolights/eb-eye/eb-eye_metabolights_complete.xml"
                                 # result <- XML::xmlParse(url)
                                 # tbl <- XML::xmlToList(url)
                                 # rootnode <- XML::xmlRoot(result)
                                 # 
                                 # # Find number of nodes in the root.
                                 # rootsize <- XML::xmlSize(rootnode)
                                 # 
                                 # db.rows.from.study <- pbapply::pblapply(tbl$entries, cl=session_cl, function(study){
                                 #   ids = sapply(study$cross_references, function(x) x[[1]])
                                 #   ids = ids[grepl(x=ids, pattern = "MTBLC")]
                                 #   names = sapply(study$additional_fields, function(x){
                                 #     #print(x)
                                 #     try(
                                 #       {if(x[[2]] == "metabolite_name") x[[1]]
                                 #       })
                                 #     names <- unlist(names)
                                 #     names <- grep(x = names, pattern = "Error", value = TRUE, invert = TRUE)
                                 #   })
                                 #   #print(names)
                                 #   #print(unique(names)[,1])
                                 #   if(length(ids)>0){
                                 #     res <- data.table::data.table(identifier = unique(ids),
                                 #                                   compoundname = unique(names)[,1])
                                 #     print(res)
                                 #     res
                                 #   }
                                 # })
                                 # 
                                 # joined <- data.table::rbindlist(db.rows.from.study)
                                 # 
                                 # # Print the result.
                                 # print(rootsize)
                                 # 
                                 # rootnode
                                 # compoundnames <- XML::getNodeSet(metabs_all, 
                                 #                                  "/cross_references/*")
                                 
                                 # file.url <- "ftp://ftp.ebi.ac.uk/pub/databases/metabolights/xml_feeds/unichem.tsv"
                                 # # ----
                                 # base.loc <- file.path(getOptions("user_options.txt")$db_dir, "metabolights_source")
                                 # if(!dir.exists(base.loc)) dir.create(base.loc)
                                 # csv.file <- file.path(base.loc, basename(file.url))
                                 # utils::download.file(file.url, csv.file)
                                 # # -------------------------------
                                 # metabol.tab <- fread(csv.file)
                                 # cpds <- read.table(file =  sep="\t")
                                 # met_info <- jsonlite::fromJSON(txt = "https://www.ebi.ac.uk/metabolights/webservice/beta/compound/MTBLC11449", simplifyDataFrame = TRUE)
                                 # - - - - - -
                                 
                                 metabs <- all_ids$identifier
                                 #cl = parallel::makeCluster(3)
  
                                 metabs <- pbapply::pblapply(metabs, cl = session_cl, FUN=function(id){
                                   met_info = NA
                                   try({
                                     url <- paste0("https://www.ebi.ac.uk/metabolights/webservice/beta/compound/", id)
                                     tries = 4
                                     while(is.na(met_info) & tries > 0){
                                       #if(tries < 4) print(paste0("retrying ",id))
                                       met_info <- jsonlite::fromJSON(txt = url)
                                       tries = tries - 1
                                     }
                                   })
                                   met_info
                                 })
                                 
                                 print(length(metabs))
                                 
                                 db_rows <- pbapply::pblapply(metabs[!is.na(metabs)], cl = 0, FUN=function(met_info){
                                   formula <- NA
                                   charge <- NA
                                   try({
                                     formula <- met_info$formula
                                     charge = met_info$charge
                                   })
                                   smi <- toupper(gsub(" ", "", met_info$smiles))
                                   # - - - - - - - - - - - - - - --
                                   if(is.na(formula) | is.na(charge)){
                                     # inchi or smiles
                                     worked = FALSE
                                     try({
                                       smi <- met_info$smiles
                                       iatom <- rcdk::parse.smiles(smi)[[1]]
                                       charge = rcdk::get.total.charge(iatom)
                                       formula = rcdk::get.mol2formula(iatom, charge = charge)@string
                                       worked = TRUE
                                     })
                                     if(!worked){
                                       print("no smiles..")
                                       smi <- toupper(webchem::cs_inchi_smiles(inchi = met_info$inchi,verbose = TRUE))
                                       iatom <- rcdk::parse.smiles(smi)[[1]]
                                       charge <- rcdk::get.total.charge(iatom)
                                       formula <- rcdk::get.mol2formula(iatom, charge = charge)@string
                                     }
                                   }
                                   met_dt <- data.table(compoundname = met_info$name,
                                                        description = if(!is.null(met_info$definition)) met_info$definition else "Unknown",
                                                        baseformula = formula,
                                                        identifier = met_info$id,
                                                        structure = smi,
                                                        charge = charge)
                                   # -------
                                   list(metabs=met_dt, pathway=met_info$pathways)
                                 })
                                 
                                 db_rows_a <- lapply(db_rows, function(x) x$metabs)
                                 db.formatted <- rbindlist(db_rows_a[!is.na(db_rows_a)])
                                 
                                 #db.formatted$pathway = c(NA)
                                 #has_pathway = sapply(metabs, function(x){ as.logical(x$flags$hasPathways)})
                                 #metab_has_pathway = names(has_pathway[has_pathway])
                                 #db.formatted.has.pw <- db.formatted[identifier %in% metab_has_pathway]
                                 #db.pathway.info <- lapply(db_rows, function(x) if(x$metabs$identifier %in% metab_has_pathway){ list(id = x$metabs$identifier, pw = x$pathway);  })
                                 #db.pathway.info <- Filter(Negate(is.null), db.pathway.info)
                                 
                                 # pathway.tb.list <- pbapply::pblapply(db.pathway.info, cl=0, FUN=function(item){
                                 #   rows <- lapply(item$pw, function(group){
                                 #     rows <- lapply(group, function(subgroup){
                                 #       data.table(subgroup)
                                 #     })
                                 #     rbindlist(rows)
                                 #   })
                                 #   rbindlist(rows, fill=TRUE)
                                 # })
                                 
                                 #names(pathway.tb.list) <- sapply(db.pathway.info, function(x) x$id)
                                 
                                 # with.pw.info <- pbapply::pblapply(1:length(pathway.tb.list), cl=session_cl, function(i){
                                 #   tb = pathway.tb.list[[i]]
                                 #   id = names(pathway.tb.list)[i]
                                 #   namecol = which(colnames(tb) == "name")
                                 #   reordered <- lapply(grep(colnames(tb),pattern = "id|ID|Id"), function(idcol){
                                 #     subsection = data.table::as.data.table(unique(tb[,c(idcol, namecol),with=FALSE]))
                                 #     colnames(subsection) <- c("id", "name")
                                 #     # - - - - - 
                                 #     subsection
                                 #   })
                                 #   repasted <- data.table::rbindlist(reordered)
                                 #   keep <- unique(repasted[complete.cases(repasted),"name"])
                                 #   keep$cpdid <- id
                                 #   # - - - - -
                                 #   keep
                                 # })
                                 # 
                                 # pathways <- rbindlist(with.pw.info)
                                 # 
                                 # cpds.w.pws <- merge(db.formatted, pathways, by.x = "identifier", by.y = "cpdid",all.x = TRUE)
                                 # 
                                 BACKUP <<- db.formatted
                                 RSQLite::dbWriteTable(conn, "base", db.formatted, overwrite=TRUE)
                                 
                               }, dimedb = function(dbname, ...){
                                 files = c(#"structures.zip", 
                                   "dimedb_pathways.zip", 
                                   "dimedb_sources.zip",
                                   "dimedb_pc_info.zip",
                                   "dimedb_id_info.zip")
                                 file.url <- "https://dimedb.ibers.aber.ac.uk/help/downloads/"
                                 file.urls <- paste0(file.url, files)
                                 # ----
                                 print("Downloading files...")
                                 base.loc <- file.path(getOptions("user_options.txt")$db_dir, "dimedb_source")
                                 if(!dir.exists(base.loc)) dir.create(base.loc)
                                 pbapply::pbsapply(file.urls, function(url){
                                   zip.file <- file.path(base.loc, basename(url))
                                   utils::download.file(url, zip.file)
                                   utils::unzip(zip.file, exdir = base.loc)  
                                 })
                                 atom <- fread(file.path(base.loc,"dimedb_pc_info.tsv"))
                                 ids <- fread(file.path(base.loc,"dimedb_id_info.tsv"))
                                 source <- fread(file.path(base.loc,"dimedb_sources.tsv"))
                                 pathway <- fread(file.path(base.loc, "dimedb_pathways.tsv"))
                                 
                                 unique.inchi <- unique(ids$InChIKey)
                                 
                                 joined <- rbind(ids, atom)
                                 casted <- reshape2::dcast(joined, InChIKey ~ Property, function(vec) paste0(vec, collapse=","))
                                 
                                 db.formatted <- data.table(compoundname = Hmisc::capitalize(tolower(casted$Name)), 
                                                            description = do.call(paste0, c(casted[,c("IUPAC Name", "Synonym")], col="-")),
                                                            baseformula = casted$`Molecular Formula`,
                                                            identifier =  casted$InChIKey, 
                                                            charge = casted$`Formal Charge`,
                                                            structure = casted$SMILES,
                                                            pathway = c(NA)
                                 )
                                 
                                 # - - -
                                 
                                 db.formatted[which(db.formatted$description == "-")] <- c("Unknown")
                                 db.formatted$description <- gsub(db.formatted$description, pattern = "^-|-$", replacement = "")
                                 
                                 # - - -
                                 
                                 RSQLite::dbWriteTable(conn, "base", db.formatted, overwrite=TRUE)
                               }, wikidata = function(dbname, ...){
                                 
                                 sparql_query <- 'PREFIX wd: <http://www.wikidata.org/entity/>
                                 PREFIX wds: <http://www.wikidata.org/entity/statement/>
                                 PREFIX wdv: <http://www.wikidata.org/value/>
                                 PREFIX wdt: <http://www.wikidata.org/prop/direct/>
                                 PREFIX wikibase: <http://wikiba.se/ontology#>
                                 PREFIX p: <http://www.wikidata.org/prop/>
                                 PREFIX ps: <http://www.wikidata.org/prop/statement/>
                                 PREFIX pq: <http://www.wikidata.org/prop/qualifier/>
                                 PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
                                 PREFIX bd: <http://www.bigdata.com/rdf#>
                                 
                                 SELECT ?chemical_compound ?chemical_compoundLabel ?chemical_formula ?chemical_compoundDescription ?canonical_SMILES ?roleLabel WHERE {
                                 SERVICE wikibase:label { bd:serviceParam wikibase:language "en, de". }
                                 ?chemical_compound wdt:P31 wd:Q11173.
                                 ?chemical_compound wdt:P274 ?chemical_formula.
                                 ?chemical_compound wdt:P233 ?canonical_SMILES.
                                 OPTIONAL {?chemical_compound wdt:P2868 ?role.}}'
                                 
                                 db.1 <- WikidataQueryServiceR::query_wikidata(sparql_query, 
                                                                               format = "simple")
                                 
                                 db.1$chemical_compound <- basename(db.1$chemical_compound)
                                 db.1$chemical_compoundDescription[db.1$chemical_compoundDescription == "chemical compound"] <- NA
                                 db.1$roleLabel[db.1$roleLabel == ""] <- NA
                                 
                                 # https://spark.apache.org/ for speed increases? is it useful locally? more an HPC thing?
                                 
                                 #cl = parallel::makeCluster(3, "FORK")
                                 
                                 charge.rows <- pbapply::pblapply(unique(db.1$canonical_SMILES), 
                                                                  cl=session_cl, 
                                                                  function(smi){
                                                                    row = data.table(canonical_SMILES = smi,
                                                                                     charge = 0)
                                                                    # - - - - -
                                                                    try({
                                                                      iatom <- rcdk::parse.smiles(smi)[[1]]
                                                                      charge = rcdk::get.total.charge(iatom)
                                                                      row = data.table(canonical_SMILES = smi,
                                                                                       charge =charge[[1]])
                                                                    },silent = TRUE)
                                                                    # - return -
                                                                    row
                                                                  })
                                 
                                 charge.dt <- rbindlist(charge.rows)
                                 
                                 db.2 <- merge(db.1, 
                                               charge.dt, 
                                               by = "canonical_SMILES", 
                                               all.y = TRUE)
                                 
                                 # NOTE: tool - myFAIR, SOAR, EUDAT, JUNIPER (ML, can use R, jupyter notebooks?), GALAXY'S GUI/API CAN BE CHANGED??
                                 
                                 db.3 <- as.data.table(aggregate(db.2, by=list(db.2$chemical_compoundLabel), function(x) c(unique((x)))))
                                 
                                 db.3$description = apply(db.3[,c("roleLabel", "chemical_compoundDescription")], 1, FUN=function(x){
                                   x <- unlist(x)
                                   print(x)
                                   paste(x[!is.na(x)], collapse=", ")
                                 })
                                 
                                 db.formatted <<- data.table::data.table(compoundname = db.3$chemical_compoundLabel,
                                                                         description = db.3$description,
                                                                         baseformula = db.3$chemical_formula, 
                                                                         identifier= db.3$chemical_compound,
                                                                         charge= db.3$charge,
                                                                         structure=db.3$canonical_SMILES)
                                 
                                 from = "\u2080\u2081\u2082\u2083\u2084\u2085\u2086\u2087\u2088\u2089"
                                 to = "0123456789"
                                 db.formatted$baseformula <- chartr(from, to, db.formatted$baseformula)
                                 
                                 db.formatted$description[db.formatted$description == ""] <<- "Unknown"
                                 
                                 # - - write - -
                                 
                                 RSQLite::dbWriteTable(conn, "base", db.formatted[, lapply(.SD, as.character),], overwrite=TRUE)
                                 
                               },vmh = function(dbname, ...){
                                 
                                 api_url <- "https://vmh.uni.lu/_api/metabolites/"
                                 
                                 pagerange = 150
                                 # get the first page
                                 
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
                                                                        identifier= table_main$reconMap3,
                                                                        charge= table_main$charge,
                                                                        structure= table_main$smile,
                                                                        isHuman = table_main$isHuman,
                                                                        isMicrobe = table_main$isMicrobe)
                                 
                                 #filter for weights that are way too much w/ weird formulas
                                 
                                 mzlimit=6000 # for now...
                                 
                                 checked <- check.chemform.joanna(isotopes, db.formatted$baseformula)
                                 keep <- which(checked$monoisotopic_mass %between% c(0, mzlimit))
                                 
                                 db.formatted <- db.formatted[keep,]
                                 
                                 # - - - - - - - - - - - - - - - - - - - 
                                 
                                 missing.id <- which(is.na(db.formatted$identifier))
                                 replacements <- table_main$reconMap # use synonum instead
                                 db.formatted$identifier[missing.id] <- replacements[missing.id]
                                 
                                 missing.desc <- which(db.formatted$description == "<NA>" | db.formatted$description == "" | is.na(db.formatted$description))
                                 replacements <- table_main$synonyms # use synonum instead
                                 db.formatted$description[missing.desc] <- replacements[missing.desc]
                                 missing.desc <- which(db.formatted$description == "<NA>" | db.formatted$description == "" | is.na(db.formatted$description))
                                 db.formatted$description[missing.desc] <- c("Unknown")
                                 
                                 descriptions <- sapply(1:nrow(db.formatted), function(i){
                                   
                                   row = db.formatted[i,]
                                   
                                   if(row$isHuman & row$isMicrobe){
                                     suffix = "Found in humans and microbes."
                                   }else if(row$isHuman & !row$isMicrobe){
                                     suffix = "Found in humans."
                                   }else{
                                     suffix = "Found in microbes"
                                   }
                                   
                                   if(length(row$description) > 1){
                                     if(substring(row$description, nchar(row$description)) == "."){
                                       paste0(row$description," -- ", suffix, " -- ")
                                     }else{
                                       paste0(row$description, ". -- ", suffix, " -- ")
                                     }
                                   }else{
                                     paste0(row$description," -- ", suffix, " -- ")
                                   }
                                 })
                                 
                                 db.formatted$description <- descriptions
                                 
                                 db.formatted <- db.formatted[,-c("isHuman", "isMicrobe")]
                                 
                                 # find formula and charge of some molecules through SMILES
                                 
                                 structs <- db.formatted$structure
                                 
                                 charges_formulae <- pbapply::pblapply(structs, function(smi){
                                   res = list(charge = NA, formula = NA)
                                   try({
                                     iatom <- rcdk::parse.smiles(smi)[[1]]
                                     charge = rcdk::get.total.charge(iatom)
                                     formula = rcdk::get.mol2formula(iatom)
                                     res = list(charge = charge[[1]], formula = formula@string)
                                   })
                                   # - - return - -
                                   res
                                 })
                                 
                                 formulae <- sapply(charges_formulae, function(x) x$formula)
                                 charges <- sapply(charges_formulae, function(x) x$charge)
                                 
                                 found.formula <- which(!is.na(formulae))
                                 found.charges <- which(!is.na(charges))
                                 
                                 db.formatted$baseformula[found.formula] <- formulae[found.formula]
                                 db.formatted$charge[found.charges] <- charges[found.charges]
                                 
                                 # check duplicated formulae and give same ID (no better solution for now..)
                                 
                                 duplicates <- split(db.formatted, by=c("baseformula", "charge"))
                                 
                                 no.duplicates <- which(sapply(duplicates, nrow) == 1)
                                 has.duplicates <- which(sapply(duplicates, nrow) > 1)
                                 
                                 single.groups <- duplicates[no.duplicates]
                                 duplicate.groups <- duplicates[has.duplicates]
                                 
                                 filled.duplicate.groups <- pbapply::pblapply(duplicate.groups, function(group){
                                   if(length(unique(group$identifier)) == 1){print(group)}
                                 })
                                 # check integrity of formulae
                                 
                                 checked <- data.table::as.data.table(check.chemform.joanna(isotopes,
                                                                                            db.formatted$baseformula))
                                 db.formatted$baseformula <- checked$new_formula
                                 
                                 # wrong <- checked[warning == TRUE, which = TRUE]
                                 # fine.compounds <- db.formatted[-wrong,]
                                 # wack.compounds <- db.formatted[wrong,]
                                 # 
                                 # # fix halogens
                                 # 
                                 # halogens <- c("F1"="fluoride", "Cl1"="chloride", "Br1"="bromide", "I1"="iodide", "At1"="astatide")
                                 # w.halogens <- grep(x=wack.compounds$baseformula, "X")
                                 # 
                                 # halogen.fixed.list <- pbapply::pblapply(w.halogens, function(i){
                                 #   
                                 #   row <- last.wack.compounds[i,]
                                 #   
                                 #   if(row$baseformula == "X"){
                                 #     print("yeah skip this one")
                                 #     return(NA)
                                 #   }else{
                                 #     rows <- lapply(names(halogens), function(X){
                                 #       new.row <- row
                                 #       new.row$baseformula = gsub(x = new.row$baseformula, pattern="X", replacement=X)
                                 #       new.row$compoundname = paste0(new.row$compoundname, " + ", halogens[[X]])
                                 #       new.row$structure = NA
                                 #       # - - return - - 
                                 #       new.row
                                 #     })
                                 #     # - - -
                                 #     rows
                                 #   }
                                 #   data.table::rbindlist(rows)
                                 # })
                                 # 
                                 # halogen.compounds <- data.table::rbindlist(halogen.fixed.list[!is.na(halogen.fixed.list)])
                                 # 
                                 # # - - residual compounds - - 
                                 # 
                                 # w.residual <- grep(x=last.wack.compounds$baseformula, "FULL")
                                 # 
                                 # residuals <- last.wack.compounds[w.residual,]
                                 # 
                                 # # how to fix these?? they're always attached to something... skip for now
                                 # 
                                 # # - - totals - - 
                                 # 
                                 # subtables <- list(
                                 #   fine.compounds,
                                 #   fixed.compounds,
                                 #   halogen.compounds
                                 # )
                                 # 
                                 # db.formatted <- unique(rbindlist(subtables))
                                 # 
                                 # - - write - - 
                                 
                                 RSQLite::dbWriteTable(conn, "base", db.formatted[, lapply(.SD, as.character),], overwrite=TRUE)
                                 
                               }, knapsack = function(dbname, ...){
                                 # TOO SLOW 
                                 # url <- "http://kanaya.naist.jp/KNApSAcK/"
                                 # response <- XML::htmlParse(url)
                                 # tab <- as.data.table(XML::readHTMLTable(response)[[1]])
                                 # metab_count <- as.numeric(gsub(tab[V1 == "metabolite", 2], pattern = " entries", replacement = ""))
                                 # 
                                 # url2 = "http://kanaya.naist.jp/knapsack_jsp/information.jsp?sname=C_ID&word=C00001001"
                                 # response <- XML::htmlParse(url2)
                                 # 
                                 # 
                                 # url <- "http://kanaya.naist.jp/knapsack_jsp/information.jsp?sname=C_ID&word="
                                 # hrefs <- list()
                                 # 
                                 # max_digits = 8 #C00000001
                                 # ids <- sapply(1:metab_count, function(i){
                                 #   paste0("C", str_pad(i, max_digits, pad = "0"))
                                 #   })
                                 # cl = makeCluster(3, "FORK")
                                 # responses = pbapply::pblapply(ids, cl=session_cl, function(id){
                                 #   #print(id)
                                 #   response <- XML::htmlParse(paste0(url,id))
                                 #   # - - - return - - -
                                 #   response
                                 # })
                               }, respect = function(dbname, ...){
                                 # - - download reSpect database, phytochemicals - - 
                                 
                                 file.url <- "http://spectra.psc.riken.jp/menta.cgi/static/respect/respect.zip"
                                 
                                 # ----
                                 #base.loc <- getOptions("user_options.txt")$db_dir
                                 base.loc <- file.path(getOptions("user_options.txt")$db_dir, "respect_source")
                                 if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
                                 zip.file <- file.path(base.loc, "respect.zip")
                                 utils::download.file(file.url, zip.file,mode = "w")
                                 utils::unzip(zip.file, exdir = base.loc)
                                 
                                 
                                 cpd_files <- list.files(base.loc, 
                                                         full.names = T)
                                 
                                 db_rows <- pbapply::pblapply(cpd_files, function(fn){
                                   row = NA
                                   try({
                                     lines <- readLines(fn,skipNul = T, n = 100)   
                                     split.lines <- sapply(lines, strsplit, ": ")
                                     names(split.lines) <- sapply(split.lines, function(x) x[1])
                                     split.lines <- lapply(split.lines, function(x) x[2:length(x)])
                                     row <- data.table(
                                       compoundname = split.lines$`CH$NAME`, 
                                       description = split.lines$RECORD_TITLE,
                                       baseformula = split.lines$`CH$FORMULA`,
                                       identifier = split.lines$ACCESSION, 
                                       charge = {
                                         smi = split.lines$`CH$SMILES`
                                         charge=0
                                         try({
                                           iatom <- rcdk::parse.smiles(smi)[[1]]
                                           charge = rcdk::get.total.charge(iatom)
                                         })
                                         charge
                                       },
                                       structure = split.lines$`CH$SMILES`
                                     )
                                     if(row$structure == "N/A") row$structure <- split.lines$`CH$INCHI`
                                   })
                                   row
                                 })
                                 
                                 db.formatted <- rbindlist(db_rows[!is.na(db_rows)])
                                 db.formatted <- unique(db.formatted[!is.na(baseformula),])
                                 # --- check formulae ---
                                 checked <- data.table::as.data.table(check.chemform.joanna(isotopes,
                                                                                            db.formatted$baseformula))
                                 db.formatted$baseformula <- checked$new_formula
                                 
                                 RSQLite::dbWriteTable(conn, "base", db.formatted, overwrite=TRUE)
                               },
                               foodb = function(dbname, ...){
                                 
                                 file.url <- "http://www.foodb.ca/system/foodb_2017_06_29_csv.tar.gz"
                                 
                                 # ----
                                 #base.loc <- getOptions("user_options.txt")$db_dir
                                 base.loc <- file.path(getOptions("user_options.txt")$db_dir, "foodb_source")
                                 if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
                                 zip.file <- file.path(base.loc, "foodb.zip")
                                 utils::download.file(file.url, zip.file,mode = "w")
                                 utils::untar(zip.file, exdir = base.loc,list =  T)
                                 
                                 base.table <- data.table::fread(file = file.path(base.loc, "foodb_2017_06_29_csv", "compounds.csv"))
                                 
                                 db.formatted <- data.table::data.table(compoundname = base.table$name,
                                                                        description = base.table$description,
                                                                        baseformula = base.table$moldb_formula, 
                                                                        identifier= base.table$id,
                                                                        charge= base.table$charge,
                                                                        structure= base.table$moldb_smiles)
                                 
                                 db.formatted <- unique(db.formatted)
                                 missing.charges <- which(is.na(db.formatted$charge))
                                 
                                 charges <- pbapply::pbsapply(missing.charges, cl=0, function(i, db){
                                   smi = db[i, "structure"][[1]]
                                   charge = 0 # set default placeholder to zero, seems like a decent assumption
                                   try({
                                     iatom <- rcdk::parse.smiles(smi)[[1]]
                                     charge = rcdk::get.total.charge(iatom)
                                   })
                                   charge
                                 }, db = db.formatted)
                                 
                                 db.formatted$charge[missing.charges] <- charges
                                 # --- check formulae ---
                                 checked <- data.table::as.data.table(check.chemform.joanna(isotopes,
                                                                                            db.formatted$baseformula))
                                 
                                 db.formatted$baseformula <- checked$new_formula
                                 
                                 missing.formula <- which(checked$warning)
                                 
                                 db.formatted <- db.formatted[-missing.formula,]
                                 
                                 RSQLite::dbWriteTable(conn, "base", db.formatted, overwrite=TRUE)
                               },
                               massbank = function(dbname, ...){
                                 file.url <- "https://github.com/MassBank/MassBank-data/archive/master.zip"
                                 base.loc <- file.path(getOptions("user_options.txt")$db_dir, "massbank_source")
                                 if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
                                 zip.file <- file.path(base.loc, "massbank.zip")
                                 utils::download.file(file.url, zip.file,mode = "w")
                                 utils::unzip(zip.file, exdir = base.loc)
                                 cpd_files <- list.files(base.loc, 
                                                         pattern = ".txt$",
                                                         full.names = T,
                                                         recursive = T)
                                 
                                 # - - - - - - - - - - - - - - -
                                 
                                 db_rows <- pbapply::pblapply(cpd_files, cl=session_cl, function(fn){
                                   row = NA
                                   try({
                                     lines <- readLines(fn)   
                                     split.lines <- sapply(lines, strsplit, ": ")
                                     names(split.lines) <- sapply(split.lines, function(x) x[1])
                                     split.lines <- lapply(split.lines, function(x) x[2:length(x)])
                                     row <- data.table(
                                       compoundname = split.lines$`CH$NAME`,
                                       description = split.lines$RECORD_TITLE,
                                       baseformula = split.lines$`CH$FORMULA`,
                                       identifier = split.lines$ACCESSION,
                                       charge = {
                                         smi = split.lines$`CH$SMILES`
                                         charge=0
                                         try({
                                           iatom <- rcdk::parse.smiles(smi)[[1]]
                                           charge = rcdk::get.total.charge(iatom)
                                         })
                                         charge
                                       },
                                       structure = {
                                         struct = "N/A"
                                         try({
                                           struct = split.lines$`CH$SMILES`
                                         })
                                         struct
                                       }
                                     )
                                     if(row$structure[[1]] == "N/A"){
                                       row$structure <- split.lines$`CH$INCHI`
                                       }
                                   })
                                   row
                                 })
                                 
                                 db.formatted <- rbindlist(db_rows[!is.na(db_rows)], fill=TRUE)
                                 db.formatted <- db.formatted[!is.na(baseformula),]
                                 # --- check formulae ---
                                 checked <- data.table::as.data.table(check.chemform.joanna(isotopes,
                                                                                            db.formatted$baseformula))
                                 db.formatted$baseformula <- checked$new_formula
                                 
                                 RSQLite::dbWriteTable(conn, "base", db.formatted, overwrite=TRUE)
                               })
  # --- execute ;) ---
  function.of.choice(dbname)
  RSQLite::dbDisconnect(conn)
}

#' @export
build.extended.db <- function(dbname,
                              outfolder,
                              adduct.table,
                              continue = F,
                              cl = 0,
                              fetch.limit = -1,
                              cpd.limit = -1){
  # ------------------------
  
  # build.base.db("internal")
  # dbname = "noise"
  # outfolder = getOptions("user_options.txt")$db_dir
  # adduct.table = adducts
  # continue = F
  # cpd.limit = 100
  # fetch.limit = 100
  
  # ------------------------
  data(isotopes, package = "enviPat")
  base.db <- file.path(outfolder, paste0(dbname, ".base.db"))
  full.db <- file.path(outfolder, paste0(dbname, ".full.db"))
  # ------------------------
  
  full.conn <- RSQLite::dbConnect(RSQLite::SQLite(), full.db)
  base.conn <- RSQLite::dbConnect(RSQLite::SQLite(), base.db)
  
  if(!RSQLite::dbExistsTable(full.conn, "done") | continue == FALSE){
    continue <- FALSE
    RSQLite::dbDisconnect(full.conn)
    file.remove(full.db)
    full.conn <- RSQLite::dbConnect(RSQLite::SQLite(), full.db)
    # boot.null(boot.null(boot.null(boot.null(boot.null(boot.null(boot.null(boot)))))))  
  } 
  
  # ------------------------
  RSQLite::dbExecute(full.conn, "PRAGMA journal_mode=wal") # THIS SHOULD HELP DB LOCKING if continue = ON
  # ------------------------
  limit.query <- if(cpd.limit == -1) "" else gsubfn::fn$paste("LIMIT $cpd.limit")
  if(continue){
    RSQLite::dbExecute(base.conn, gsubfn::fn$paste("ATTACH '$full.db' as db"))
    continue.query <- strwrap("SELECT DISTINCT baseformula, charge FROM base b
                              WHERE NOT EXISTS(SELECT DISTINCT baseformula, charge 
                              FROM db.done d
                              WHERE b.baseformula = d.baseformula
                              AND b.charge = d.basecharge)", width=10000, simplify=TRUE)
    total.formulae <- RSQLite::dbGetQuery(base.conn, gsubfn::fn$paste("SELECT Count(*)
                                                                      FROM ($continue.query)"))
    formula.count <- total.formulae[1,]
    results <- RSQLite::dbSendQuery(base.conn, continue.query)
  } else{
    # --- add base db to the new one ---
    print("Attaching base...")
    RSQLite::dbExecute(full.conn, gsubfn::fn$paste("ATTACH '$base.db' AS tmp"))
    RSQLite::dbExecute(full.conn, gsubfn::fn$paste("CREATE TABLE IF NOT EXISTS done(baseformula text, basecharge text)"))
    RSQLite::dbExecute(full.conn, gsubfn::fn$paste("CREATE TABLE IF NOT EXISTS base AS SELECT * FROM tmp.base"))
    if(RSQLite::dbExistsTable(base.conn, "pathways")){
      RSQLite::dbExecute(full.conn, gsubfn::fn$paste("CREATE TABLE IF NOT EXISTS pathways AS SELECT * FROM tmp.pathways"))
    }
    print("Indexing base...")
    RSQLite::dbExecute(full.conn, "CREATE INDEX IF NOT EXISTS b_idx1 on base(baseformula, charge)")
    # ----------------------
    sql.make.meta <- strwrap("CREATE TABLE IF NOT EXISTS extended(
                             baseformula text,
                             fullformula text, 
                             basemz decimal(30,13), 
                             fullmz decimal(30,13), 
                             adduct text,
                             basecharge int,
                             totalcharge int,
                             isoprevalence float,
                             foundinmode text)", width=10000, simplify=TRUE)
    RSQLite::dbExecute(full.conn, sql.make.meta)
    # --------------------
    get.query <- gsubfn::fn$paste("SELECT DISTINCT b.baseformula, b.charge FROM base b $limit.query")
    total.formulae <- RSQLite::dbGetQuery(base.conn, gsubfn::fn$paste("SELECT Count(*)
                                                                      FROM ($get.query)"))
    formula.count <- total.formulae[1,]
    results <- RSQLite::dbSendQuery(base.conn, get.query)
  }
  # --- start pb ---
  pb <- pbapply::startpb(0, formula.count)
  print("Starting DB generation.")
  print(paste("Approximate batches:", formula.count / fetch.limit ))
  # --- waow, my first while in R ---
  while(!RSQLite::dbHasCompleted(results)){
    # --- fetch part of results ---
    partial.results <- data.table::as.data.table(RSQLite::dbFetch(results, fetch.limit))
    if(length(partial.results$baseformula) == 0) next
    # -----------------------
    
    #print(partial.results$baseformula)
    
    print(paste(RSQLite::dbGetRowCount(results), formula.count, sep=" / "))
    # -----------------------
    
    checked.formulae <- data.table::as.data.table(check.chemform.joanna(isotopes, 
                                                                        partial.results$baseformula))
    # -----------------------------
    # keep.rows <- checked.formulae[warning == FALSE & monoisotopic_mass %between% c(60, 600), which=TRUE]
    # -----------------------------
    keep.rows <- checked.formulae[warning == FALSE, which=TRUE]
    # -----------------------------
    if(length(keep.rows) == 0) next
    backtrack <- data.table::data.table(baseformula = checked.formulae[keep.rows, new_formula],
                                        basemz = checked.formulae[keep.rows, monoisotopic_mass],
                                        charge = c(partial.results$charge)[keep.rows])
    # -- go through each adduct --
    do.calc <- function(x){
      row <- adduct.table[x,]
      name <- row$Name
      
      #print(paste("--- CURRENT ADDUCT:", name, "---"))
      
      adduct <- row$Formula_add
      deduct <- row$Formula_ded
      # --- fix booleans ---
      adduct <- if(adduct == "FALSE") FALSE else adduct
      deduct <- if(deduct == "FALSE") FALSE else deduct
      # --------------------
      mode <- row$Ion_mode
      multiplier <- as.numeric(row$Multi)
      # --- multiplication ---
      formulae.mult <- multiform.joanna(backtrack$baseformula, multiplier)
      backtrack$multiform <- formulae.mult
      backtrack$adducted <- backtrack$multiform
      # --- is adduction necessary? ---
      if(adduct != FALSE){
        formulae.add <- mergeform.joanna(formula1 = backtrack$baseformula, 
                                         formula2 = adduct)
        backtrack$adducted <- formulae.add
      }
      # --- placeholder ---
      backtrack$final <- backtrack$adducted
      # --- is deduction necessary? ---
      if(deduct != FALSE){
        formulae <- backtrack$final
        can.deduct <- which(!check.ded.joanna(formulas = backtrack$adducted,
                                              deduct = deduct))
        if(length(can.deduct) == 0) return(NA)
        deductibles <- formulae[can.deduct]
        formulae.ded <- subform.joanna(deductibles, deduct)
        backtrack$final[can.deduct] <- formulae.ded
        backtrack <- backtrack[can.deduct]
      }
      if(nrow(backtrack) == 0) return(NA)
      backtrack$final.charge <- c(as.numeric(backtrack$charge)) + c(as.numeric(row$Charge))
      backtrack <- backtrack[final.charge != 0]
      if(nrow(backtrack) == 0) return(NA)
      
      # --- get isotopes ---
      
      # base_mzs <- enviPat::isopattern(
      #   isotopes,
      #   backtrack$baseformula, #backtrack$final,
      #   threshold = 0.1,
      #   plotit = FALSE,
      #   charge = FALSE, #backtrack$final.charge,
      #   algo = 2,
      #   verbose = FALSE
      # )
      
      isotables <- enviPat::isopattern(
        isotopes,
        backtrack$final,
        threshold = 0.1,
        plotit = FALSE,
        charge = backtrack$final.charge,
        algo = 2,
        verbose = FALSE
      )
      
      # for(i in 1:length(full_mzs)){
      #   match = all(full_mzs[[i]][,"abundance"] == base_mzs[[i]][,"abundance"])
      #   # - - -
      #   if(!match){
      #     print(nrow(base_mzs[[i]]))
      #     # print(gsubfn::fn$paste("--- M ---"))
      #     # print(base_mzs[[i]])
      #     # print(gsubfn::fn$paste("--- $name ---"))
      #     # print(full_mzs[[i]])
      #     # replace
      #     full_mzs[[i]] <- full_mzs[[i]][1:nrow(base_mzs[[i]]),]
      #     full_mzs[[i]][,"abundance"] <- full_mzs[[i]][,"abundance"]
      #   }
      #   #print(full_mzs[[i]])
      # }
      # print("here???")
      # -------------------------
      
      isolist <- lapply(isotables, function(isotable){
        #print(isotable)
        if(isotable[[1]] == "error") return(NA)
        iso.dt <- data.table::data.table(isotable, fill=TRUE)
        result <- iso.dt[,1:2]
        names(result) <- c("fullmz", "isoprevalence")
        # --- return ---
        result
      })
      
      isolist.nonas <- isolist[!is.na(isolist)]
      isotable <- rbindlist(isolist.nonas)
      
      
      keep.isos <- names(isolist.nonas)
      
      # --- remove 'backtrack' rows that couldn't be calculated ---
      backtrack.final <- backtrack[final %in% keep.isos]
      
      # --------------------
      repeat.times <- c(unlist(lapply(isolist.nonas, FUN=function(list) nrow(list)), use.names = FALSE))
      meta.table <- data.table::data.table(baseformula = rep(backtrack.final$baseformula, repeat.times),
                                           fullformula = rep(backtrack.final$final, repeat.times),
                                           basemz = rep(backtrack.final$basemz, repeat.times),
                                           fullmz =isotable$fullmz,
                                           adduct = c(name),
                                           basecharge = rep(backtrack.final$charge, repeat.times),
                                           totalcharge = rep(backtrack.final$final.charge, repeat.times),
                                           isoprevalence = isotable$isoprevalence,
                                           foundinmode = c(mode)
      )
      # --- return ---
      meta.table
    }
    tab.list <- pbapply::pblapply(cl=session_cl, 1:nrow(adduct.table), FUN=function(x) do.calc(x))
    #tab.list <- parallel::parLapply(cl=session_cl, 1:nrow(adduct.table), fun=function(x) do.calc(x))
    # --- progress bar... ---
    pbapply::setpb(pb, RSQLite::dbGetRowCount(results))
    total.table <- rbindlist(tab.list[!is.na(tab.list)])
    done.table <- unique(total.table[,c("baseformula", "basecharge")])
    # --- filter on m/z ===
    RSQLite::dbWriteTable(full.conn, "extended", total.table, append=T)
    RSQLite::dbWriteTable(full.conn, "done", done.table, append=T)
  }
  RSQLite::dbClearResult(results)
  # --- indexy ---
  print("Indexing extended table...")
  RSQLite::dbExecute(full.conn, "create index e_idx1 on extended(baseformula, basecharge)")
  RSQLite::dbExecute(full.conn, "create index e_idx2 on extended(fullmz, foundinmode)")
  # ====== range table =====
  print("Adding range tables for raw peak grouping...")
  # ppm = 2 ### PPM HARD CODED HERE FOR NOW
  # # ------------------------
  # mzvals <- data.table::data.table(mzmed = total.table$fullmz,
  #                                  foundinmode = total.table$foundinmode)
  # mzranges <- data.table::data.table(mzmin = pbapply::pbsapply(total.table$fullmz,cl=session_cl,
  #                                                              FUN=function(mz, ppm){
  #                                                                mz - mz * (ppm / 1E6)}, ppm=ppm),
  #                                    mzmax = pbapply::pbsapply(total.table$fullmz,cl=session_cl,
  #                                                              FUN=function(mz, ppm){
  #                                                                mz + mz * (ppm / 1E6)}, ppm=ppm))
  # sql.make.meta <- strwrap("CREATE TABLE mzvals(
  #                          ID INTEGER PRIMARY KEY AUTOINCREMENT,
  #                          mzmed decimal(30,13),
  #                          foundinmode text)", width=10000, simplify=TRUE)
  # RSQLite::dbExecute(full.conn, sql.make.meta)
  # RSQLite::dbExecute(full.conn, "create index mzfind on mzvals(mzmed, foundinmode);")
  # # --- write vals to table ---
  # RSQLite::dbWriteTable(full.conn, "mzvals", mzvals, append=TRUE) # insert into
  # # --- make range table (choose if R*tree or not) ---
  # sql.make.rtree <- strwrap("CREATE VIRTUAL TABLE mzranges USING rtree(
  #                           ID INTEGER PRIMARY KEY AUTOINCREMENT,
  #                           mzmin decimal(30,13),
  #                           mzmax decimal(30,13));"
  #                           , width=10000, simplify=TRUE)
  # sql.make.normal <- strwrap("CREATE TABLE mzranges(
  #                            ID INTEGER PRIMARY KEY AUTOINCREMENT,
  #                            mzmin decimal(30,13),
  #                            mzmax decimal(30,13));", width=10000, simplify=TRUE)
  # RSQLite::dbExecute(full.conn, sql.make.rtree)
  # # --- write ranges to table ---
  # RSQLite::dbWriteTable(full.conn, "mzranges", mzranges, append=TRUE) # insert into
  # --- cleanup ---
  RSQLite::dbExecute(full.conn, "VACUUM")
  # ==========================
  print("Disconnecting...")
  RSQLite::dbDisconnect(base.conn)
  RSQLite::dbDisconnect(full.conn)
  on.exit(pbapply::closepb(pb))
  print("Done! :-)")
}

#' @export
build.pat.db <- function(db.name, 
                         pospath, 
                         negpath, 
                         overwrite=FALSE,
                         rtree=TRUE,
                         make.full = TRUE,
                         ppm=2){
  
  
  # ------------------------
  ppm = as.numeric(ppm)
  # ------------------------
  
  poslist <- fread(pospath,header = T)
  neglist <- fread(negpath,header = T)
  
  # ------------------------------------
  
  keepcols <- intersect(colnames(poslist), colnames(neglist))
  
  poslist <- poslist[,..keepcols]
  neglist <- neglist[,..keepcols]
  
  # - - - fix QCs - - -
  
  which.qc <- grep(colnames(poslist), pattern = "^QC")
  qc.i = 1
  
  for(qc in which.qc){
    print(qc.i)
    new.qc.name <- paste0("QC", qc.i)
    new.qc.name <- gsub(colnames(poslist)[qc], pattern = "(^QC[\\d|\\d\\d])", replacement = new.qc.name,perl = T)
    print(new.qc.name)
    colnames(poslist)[qc] <- new.qc.name
    colnames(neglist)[qc] <- new.qc.name
    # - - -
    qc.i = qc.i + 1
  }
  # - - - - - - - - - -
  setProgress(.20)
  
  gc()
  
  mzvals <- data.table::data.table(mzmed = c(as.numeric(poslist$mzmed), as.numeric(neglist$mzmed)),
                                   foundinmode = c(rep("positive", nrow(poslist)), rep("negative", nrow(neglist))))
  
  mzranges <- data.table::data.table(mzmin = sapply(c(as.numeric(poslist$mzmed), as.numeric(neglist$mzmed)), 
                                                    FUN=function(mz, ppm){
                                                      mz - mz * (ppm / 1E6)}, ppm=ppm),
                                     mzmax = sapply(c(as.numeric(poslist$mzmed), as.numeric(neglist$mzmed)), 
                                                    FUN=function(mz, ppm){
                                                      mz + mz * (ppm / 1E6)}, ppm=ppm))
  
  mzvals$foundinmode <- trimws(mzvals$foundinmode)
  
  # --- SAVE BATCH INFO (kinda ugly...  ; _;") --- 
  
  if(any(grepl("\\*", x = colnames(poslist)))){
    samp_split = strsplit(colnames(poslist)[2:ncol(poslist)], "\\*")
    batch_split = strsplit(unlist(lapply(samp_split, function(x) x[2])), "\\_")
    batch_info = data.table::data.table(sample = sapply(samp_split, function(x) x[1]),
                                        batch = sapply(samp_split, function(x) x[2]),
                                        injection = sapply(samp_split, function(x) x[3]))
    print(batch_info)
    colnames(poslist) = gsub(colnames(poslist), pattern = "(\\*.*$)", replacement = "")
    colnames(neglist) = gsub(colnames(poslist), pattern = "(\\*.*$)", replacement = "")
  } else{
    batch_info = data.table::data.table(sample = colnames(poslist)[2:ncol(poslist)], batch = c(1))
  }
  
  gc()
  
  setProgress(.30)
  
  poslist <- data.table::melt(poslist,#[,(rmv.cols) := NULL], 
                              id.vars="mzmed",
                              variable.name="filename",
                              value.name="intensity",
                              variable.factor=TRUE
  )
  neglist <- data.table::melt(neglist,#[,(rmv.cols) := NULL], 
                              id.vars="mzmed",
                              variable.name="filename",
                              value.name="intensity",
                              variable.factor=TRUE
  )
  setProgress(.40)
  
  #poslist$filename <- trimws(poslist$filename)
  #neglist$filename <- trimws(neglist$filename)
  
  # COMPUTER CANT HANDLE THE RBIND, WRITE THEM SEPERATELY
  
  #mzintensities = rbind(poslist, neglist)
  
  # ------------------------
  
  if(overwrite==TRUE & file.exists(db.name)) file.remove(db.name)
  
  # --- reconnect / remake ---
  
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), db.name)
  
  # ------------------------
  
  #RSQLite::dbCreateTable()
  RSQLite::dbWriteTable(conn, "batchinfo", batch_info, overwrite=T) # insert into
  
  # ------------------------
  
  sql.make.int <- strwrap("CREATE TABLE mzintensities(
                          ID INTEGER PRIMARY KEY AUTOINCREMENT,
                          mzmed decimal(30,13),
                          filename text,
                          intensity float)", width=10000, simplify=TRUE)
  
  RSQLite::dbExecute(conn, sql.make.int)
  
  setProgress(.60)
  
  # --- write intensities to table and index ---
  RSQLite::dbWriteTable(conn, "mzintensities", poslist, append=TRUE) # insert into
  RSQLite::dbWriteTable(conn, "mzintensities", neglist, append=TRUE) # insert into
  
  RSQLite::dbExecute(conn, "CREATE INDEX intindex ON mzintensities(filename,'mzmed',intensity)")
  RSQLite::dbExecute(conn, "CREATE INDEX intindex2 ON mzintensities('mzmed')")
  
  # ------------------------
  
  sql.make.meta <- strwrap("CREATE TABLE mzvals(
                           ID INTEGER PRIMARY KEY AUTOINCREMENT,
                           mzmed decimal(30,13),
                           foundinmode text)", width=10000, simplify=TRUE)
  RSQLite::dbExecute(conn, sql.make.meta)
  RSQLite::dbExecute(conn, "create index mzfind on mzvals(mzmed, foundinmode);")
  
  setProgress(.70)
  
  # --- write vals to table ---
  RSQLite::dbWriteTable(conn, "mzvals", mzvals, append=TRUE) # insert into
  
  # --- make range table (choose if R*tree or not) ---
  sql.make.rtree <- strwrap("CREATE VIRTUAL TABLE mzranges USING rtree(
                            ID INTEGER PRIMARY KEY AUTOINCREMENT,
                            mzmin decimal(30,13),
                            mzmax decimal(30,13));"
                            , width=10000, simplify=TRUE)
  sql.make.normal <- strwrap("CREATE TABLE mzranges(
                             ID INTEGER PRIMARY KEY AUTOINCREMENT,
                             mzmin decimal(30,13),
                             mzmax decimal(30,13));", width=10000, simplify=TRUE)
  RSQLite::dbExecute(conn, if(rtree) sql.make.rtree else sql.make.normal)
  
  setProgress(.80)
  
  # --- write ranges to table ---
  RSQLite::dbWriteTable(conn, "mzranges", mzranges, append=TRUE) # insert into
  # --- cleanup ---
  RSQLite::dbExecute(conn, "VACUUM")
  
  setProgress(.90)
  
  # ----------------
  RSQLite::dbDisconnect(conn)
  print("Made!")}

#' @export
load.excel <- function(path.to.xlsx, 
                       path.to.patdb,
                       tabs.to.read = c(
                         #"General",
                         "Setup",
                         "Individual Data"
                         #,"Pen Data",
                         #"Admin"
                       )){
  
  #path.to.xlsx <- "~/Desktop/xls/DSM_NL_BR_IT.xlsx"
  #path.to.patdb <- global$paths$patdb
  
  print(path.to.patdb)
  print(path.to.xlsx)
  
  # --- connect to sqlite db ---
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), path.to.patdb)
  # -------------------------------
  data.store <- pbapply::pblapply(tabs.to.read, FUN=function(tab.name){
    #tab <- data.table::as.data.table(xlsx::read.xlsx(path.to.xlsx, sheetName = tab.name))
    tab <- data.table::as.data.table(openxlsx::read.xlsx(path.to.xlsx, sheet = tab.name))
    # --- reformat colnames ---
    colnames(tab) <- tolower(gsub(x=colnames(tab), pattern = "\\.$|\\.\\.$", replacement = ""))
    colnames(tab) <- gsub(x=colnames(tab), pattern = "\\.|\\.\\.", replacement = "_")
    colnames(tab)[grep(x=colnames(tab), pattern= "*date*")] <- "sampling_date"
    print(colnames(tab))
    # -----------------------------------------------------------------------------------------
    data.table::as.data.table(tab, keep.rownames=F)
  })
  
  # --- convert to data table --- ## make this nicer loooking in the future
  #general <- data.store[[1]]
  setup <- data.store[[1]]
  individual.data <- data.store[[2]]
  
  # --- fill empty cells w/ na ---
  
  indx <- which(sapply(setup, is.character)) 
  for (j in indx) set(setup, i = grep("^$|^ $", setup[[j]]), j = j, value = NA_character_)
  
  # indx <- which(sapply(general, is.character)) 
  # for (j in indx) set(general, i = grep("^$|^ $", general[[j]]), j = j, value = NA_character_)
  
  indx <- which(sapply(individual.data, is.character)) 
  for (j in indx) set(individual.data, i = grep("^$|^ $", individual.data[[j]]), j = j, value = NA_character_)
  
  # --- remove empty lines ---
  
  #general <- general[rowSums(is.na(general)) != ncol(general),]
  setup <- setup[rowSums(is.na(setup)) != ncol(setup),]
  individual.data <- individual.data[rowSums(is.na(individual.data)) != ncol(individual.data),]
  
  # --------------------------
  
  if(any(is.na(as.numeric(individual.data$sampling_date)))){
    individual.data$sampling_date <- as.factor(as.Date(as.character(individual.data$sampling_date), 
                                                       format = "%d-%m-%y"))
  }else{
    individual.data$sampling_date <- as.factor(as.Date(as.numeric(individual.data$sampling_date), 
                                                       origin = "1899-12-30"))
    
  }
  
  individual.data$card_id <- as.character(individual.data$card_id)
  individual.data$animal_internal_id <- as.character(individual.data$animal_internal_id)
  
  if(is.na(individual.data$sampling_date[1])) levels(individual.data$sampling_date) <- factor(1) 

  setup <- data.table::as.data.table(apply(setup, MARGIN=2, trimws))
  individual.data <- data.table::as.data.table(apply(individual.data, MARGIN=2, trimws))
  #general <- data.table::as.data.table(apply(general, MARGIN=2, trimws))
  
  # --- add the QC samples ---
  
  qc_samps = RSQLite::dbGetQuery(conn, "SELECT * FROM batchinfo WHERE sample LIKE '%QC%'")
  
  placeholder_date <- individual.data$sampling_date[[1]]
    
  qc_ind_data <- lapply(qc_samps$sample, function(qc) {
    data.table(label = c(1),
               card_id = qc,
               animal_internal_id = qc,
               sampling_date = placeholder_date,
               sex = "qc",
               group = "qc",
               farm = "QcLand")
  })
  
  qc_tab_setup = data.table(group = "qc",
                            stool_condition = "qc")
  qc_tab_ind = unique(rbindlist(qc_ind_data))
  
  # --- join to existing ---
  
  setup <- rbind(setup, qc_tab_setup, fill=TRUE)
  individual.data <- rbindlist(list(individual.data, qc_tab_ind), fill=TRUE)
  individual.data$label <- 1:nrow(individual.data)
  
  #pen.data <- data.table::as.data.table(apply(pen.data, MARGIN=2, trimws))
  #admin <- data.table::as.data.table(apply(admin, MARGIN=2, trimws))
  
  # --- import to patient sql file ---
  #RSQLite::dbWriteTable(conn, "general", general, overwrite=TRUE) # insert into BUGGED FIX LATER
  RSQLite::dbWriteTable(conn, "setup", setup, overwrite=TRUE) # insert into
  RSQLite::dbWriteTable(conn, "individual_data", individual.data, overwrite=TRUE) # insert into
  #RSQLite::dbWriteTable(conn, "pen_data", pen.data, overwrite=TRUE) # insert into
  #RSQLite::dbWriteTable(conn, "admin", admin, overwrite=TRUE) # insert into
  # --- disconnect ---
  RSQLite::dbDisconnect(conn)
}
