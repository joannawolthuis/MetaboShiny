#' @export
build.base.db <- function(dbname=NA, 
                          outfolder="/Users/jwolthuis/Google Drive/Metabolomics/Databases", 
                          cl=FALSE){
  # --- check if user chose something ---
  if(is.na(dbname)) return("~ Please choose one of the following options: HMDB, ChEBI, PubChem! d(>..w.<)b ~")
  # --- make dir if not exists ---
  if(!dir.exists(outfolder)) dir.create(outfolder)
  # --- create .db file ---
  db <- file.path(outfolder, paste0(dbname, ".base.db"))
  if(file.exists(db)) file.remove(db)
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), db)
  RSQLite::dbExecute(conn, statement = "create table base(compoundname text, description text, baseformula text, identifier text, charge text)")
  # --- create base ---
  function.of.choice <- switch(tolower(dbname),
                               internal = function(dbname){ # BOTH NOISE AND NORMAL
                                 # --- uses csv package ---
                                 int.loc <- file.path(".", "backend","umcfiles", "internal")
                                 # --- non-noise ---
                                 internal.base.db <- read.csv(file.path(int.loc, 
                                                                        "TheoreticalMZ_NegPos_noNoise.txt"),
                                                              sep="\t")
                                 db.formatted <- data.table(compoundname = internal.base.db[,1],
                                                            description = c("Internal"),
                                                            baseformula = internal.base.db[,2], 
                                                            identifier=c("Internal"),
                                                            charge=c(0))
                                 # fix some unicode stuff
                                 db.formatted <- db.formatted[baseformula != ""]
                                 db.formatted$baseformula <-  gsub(x=db.formatted$baseformula, 
                                                                   pattern="( )|(\xa0)|(\xe1)|( .*$)", 
                                                                   replacement="")
                                 checked <- as.data.table(check.chemform.joanna(isotopes,
                                                                                db.formatted$baseformula))
                                 db.formatted$baseformula <- checked$new_formula
                                 keep <- checked[warning == FALSE, which = TRUE]
                                 db.formatted <- db.formatted[keep]
                                 # --- write ---
                                 RSQLite::dbWriteTable(conn, "base", db.formatted, append=TRUE)},
                               noise = function(dbname){
                                 print("noisee")
                                 int.loc <- file.path(".", "backend","umcfiles", "internal")
                                 # --- noise ---
                                 noise.base.db <- read.csv(file.path(int.loc, 
                                                                     "TheoreticalMZ_NegPos_yesNoise.txt"), 
                                                           sep="\t")
                                 db.formatted <- data.table(compoundname = as.character(noise.base.db[,1]),
                                                            description = as.character(str_match(noise.base.db[,1], "(?<=\\().+?(?=\\))")),
                                                            baseformula = as.character(noise.base.db[,2]), 
                                                            identifier=c("Noise"),
                                                            charge=c(0))
                                 RSQLite::dbWriteTable(conn, "base", db.formatted, append=TRUE)
                                 ### do extended table and write to additional db (exception...)
                                 db.full <- file.path(outfolder, paste0(dbname, ".full.db"))
                                 if(file.exists(db.full)) file.remove(db.full)
                                 full.conn <- RSQLite::dbConnect(RSQLite::SQLite(), db.full)
                                 # --- create base table here too ---
                                 RSQLite::dbExecute(full.conn, statement = "create table base(compoundname text, description text, baseformula text, identifier text, charge text)")
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
                                   result <- data.table(
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
                                 print("heree")
                                 
                                 RSQLite::dbWriteTable(full.conn, 
                                              "extended", 
                                              db.formatted.full, 
                                              append=TRUE)
                                 print("heree")
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
                                 base.loc <- file.path(options$db_dir, "hmdb_source")
                                 if(!dir.exists(base.loc)) dir.create(base.loc)
                                 zip.file <- file.path(base.loc, "HMDB.zip")
                                 utils::download.file(file.url, zip.file)
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
                                 # --- make nice big table ---
                                 db.formatted <- data.table(
                                   compoundname = sapply(compoundnames, xmlValue),
                                   description = sapply(description, xmlValue),
                                   baseformula = sapply(formulae, xmlValue),
                                   identifier =  gsub(sapply(identifiers, xmlValue), pattern = "(HMDB0*)", replacement = ""),
                                   charge = sapply(charges, xmlValue)
                                 )
                                 # --- check formulae ---
                                 checked <- as.data.table(check.chemform.joanna(isotopes,
                                                                                db.formatted$baseformula))
                                 db.formatted$baseformula <- checked$new_formula
                                 keep <- checked[warning == FALSE, which = TRUE]
                                 db.formatted <- db.formatted[keep]
                                 # ----------------------
                                 RSQLite::dbWriteTable(conn, "base", db.formatted, append=TRUE)
                               },
                               chebi = function(dbname, ...){
                                 db.full <- as.data.table(download.chebi.joanna(release = "latest", 
                                                                                woAssociations = FALSE))
                                 db.formatted <- unique(db.full[, list(compoundname = ChEBI, 
                                                                       description = DEFINITION,
                                                                       baseformula = FORMULA,
                                                                       identifier = ID, 
                                                                       charge = gsub(CHARGE,pattern = "$\\+\\d", replacement = "")
                                 )
                                 ]
                                 )
                                 # --- check formulae ---
                                 checked <- as.data.table(check.chemform.joanna(isotopes,
                                                                                db.formatted$baseformula))
                                 db.formatted$baseformula <- checked$new_formula
                                 keep <- checked[warning == FALSE, which = TRUE]
                                 db.formatted <- db.formatted[keep]
                                 # ----------------------
                                 RSQLite::dbWriteTable(conn, "base", db.formatted, append=TRUE)
                               }, 
                               wikipathways = function(dbname){
                                 chebi.loc <- file.path(options$db_dir, "chebi.full.db")
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
                                 chebi.join.table <- data.table(identifier = chebi.ids,
                                                                description = chebi$results$description,
                                                                widentifier = chebi$results$mb,
                                                                pathway = chebi$results$pathway)
                                 RSQLite::dbWriteTable(conn.chebi, "wikipathways", chebi.join.table, overwrite=TRUE)
                                 db.formatted <- RSQLite::dbGetQuery(conn.chebi, "SELECT DISTINCT  b.compoundname, 
                                                                                          w.description,
                                                                                          b.baseformula, 
                                                                                          w.widentifier as identifier, 
                                                                                          b.charge,
                                                                                          w.pathway 
                                                                         FROM base b
                                                                         JOIN wikipathways w
                                                                         ON b.identifier = w.identifier")
                                 # --- get pathway info ---
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
                                 base.loc <- file.path(options$db_dir, "smpdb_source")
                                 if(!dir.exists(base.loc)) dir.create(base.loc)
                                 zip.file <- file.path(base.loc, "SMPDB.zip")
                                 utils::download.file(file.url, zip.file)
                                 utils::unzip(zip.file, exdir = base.loc)
                                 # -------------------------------
                                 smpdb.tab <- fread(file.path(base.loc, "metabolites.csv"))
                                 # --- get charges ---
                                 charges <- c(gsub(str_match(smpdb.tab$`InChI`, pattern = "q[+-](\\d*)")[,1], pattern = "q|\\+", replacement = ""))
                                 charges[is.na(charges)] <- "0" # set all missing to zero
                                 db.cpds <- data.table(compoundname = smpdb.tab$`Metabolite Name`,
                                                       baseformula = smpdb.tab$Formula,
                                                       identifier = smpdb.tab$`Metabolite ID`,
                                                       widentifier = gsub(smpdb.tab$`HMDB ID`, pattern = "(HMDB0*)", replacement = ""),
                                                       charge = charges,
                                                       pathway = smpdb.tab$`SMPDB ID`)
                                 db.pathways <- data.table(name = smpdb.tab$`Pathway Name`,
                                                           identifier = smpdb.tab$`SMPDB ID`)
                                 checked <- as.data.table(check.chemform.joanna(isotopes,
                                                                                db.cpds$baseformula))
                                 db.cpds$baseformula <- checked$new_formula
                                 keep <- checked[warning == FALSE, which = TRUE]
                                 db.cpds <- db.cpds[keep]
                                 # --- get descriptions from hmdb!!! ---
                                 hmdb.loc <- file.path(options$db_dir, "hmdb.full.db")
                                 conn.hmdb <- RSQLite::dbConnect(RSQLite::SQLite(), hmdb.loc)
                                 RSQLite::dbWriteTable(conn.hmdb, "smpdb", db.cpds, overwrite=TRUE)
                                 db.formatted <- RSQLite::dbGetQuery(conn.hmdb, "SELECT DISTINCT  
                                                                        s.compoundname,
                                                                        b.description,
                                                                        s.baseformula, 
                                                                        s.charge as charge,
                                                                        s.identifier as identifier, 
                                                                        s.pathway 
                                                                        FROM smpdb s
                                                                        LEFT JOIN base b
                                                                        ON b.identifier = s.widentifier")
                                 RSQLite::dbRemoveTable(conn.hmdb, "smpdb")
                                 RSQLite::dbDisconnect(conn.hmdb)
                                 # --- create ---
                                 RSQLite::dbWriteTable(conn, "base", db.formatted, overwrite=TRUE)
                                 RSQLite::dbWriteTable(conn, "pathways", db.pathways, overwrite=TRUE)
                               },
                               kegg = function(dbname){
                                 # ---------------
                                 batches <- split(0:2000, ceiling(seq_along(0:2000)/100))
                                 cpds <- pbapply::pblapply(batches, cl=cl, FUN=function(batch){
                                   names(keggFind("compound", batch, "mol_weight"))
                                 })
                                 cpd.ids <- Reduce(c, cpds)
                                 id.batches <- split(cpd.ids, ceiling(seq_along(cpd.ids)/10))
                                 # ---------------
                                 # keggGet("C00032") # cpd
                                 # keggGet("map00860") # pathway
                                 # keggGet("M00121") # module
                                 # keggGet("H00201") # disease
                                 # --- GET COMPOUNDS ---
                                 kegg.cpd.list <- pbapply::pblapply(id.batches, cl=cl, FUN=function(batch){
                                   rest.result <- KEGGREST::keggGet(batch)
                                   # ---------------------------
                                   base.list <- lapply(rest.result, FUN=function(cpd){
                                     cpd$NAME_FILT <- gsub(cpd$NAME, pattern = ";", replacement = "")
                                     data.table(compoundname = c(paste(cpd$NAME_FILT, collapse=", ")),
                                                description = c(if("BRITE" %in% names(cpd)) paste(cpd$BRITE, collapse=", ") else{NA}),
                                                baseformula = c(cpd$FORMULA),
                                                identifier = c(cpd$ENTRY),
                                                charge = c(kegg.charge(cpd$ATOM)),
                                                pathway = if("PATHWAY" %in% names(cpd)) names(cpd$PATHWAY) else{NA}
                                     )
                                   })
                                   rbindlist(base.list)
                                 })
                                 db.cpds <- rbindlist(kegg.cpd.list)
                                 db.cpds$baseformula <- gsub(pattern = " |\\.", replacement = "", db.cpds$baseformula) # fix salt and spacing issues                                 
                                 # -------------------------------
                                 checked <- as.data.table(check.chemform.joanna(isotopes,
                                                                                db.cpds$baseformula))
                                 db.cpds$baseformula <- checked$new_formula
                                 keep <- checked[warning == FALSE, which = TRUE]
                                 db.cpds <- db.cpds[keep]
                                 # --- GET PATHWAYS ---
                                 pathways <- unique(db.cpds$pathway)
                                 pw.batches <- split(pathways, ceiling(seq_along(pathways)/10))
                                 kegg.pw.list <- pbapply::pblapply(pw.batches, cl=cl, FUN=function(batch){
                                   rest.result <- KEGGREST::keggGet(batch)
                                   # ---------------------------
                                   base.list <- lapply(rest.result, FUN=function(pw){
                                     data.table(identifier = as.character(pw$ENTRY),
                                                name = as.character(pw$NAME),
                                                class = if("CLASS" %in% names(pw)) pw$CLASS else{NA},
                                                module = if("MODULE" %in% names(pw)) pw$MODULE else{NA},
                                                disease = if("DISEASE" %in% names(pw)) as.character(pw$DISEASE) else{NA}
                                     )
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
                               pubchem = function(dbname, ...){
                                 # --- create working space ---
                                 baseLoc <- file.path(options$db_dir, "pubchem_source")
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
                                 pbsapply(file.urls, cl=cl, FUN=function(url){
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
                                 pbsapply(cl=cl, sdf.files, FUN=function(sdf.file){
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
                                                      db.formatted <- data.table(
                                                        compoundname = blockmatrix[, if("PUBCHEM_IUPAC_TRADITIONAL_NAME" %not in% colnames(blockmatrix)) "PUBCHEM_MOLECULAR_FORMULA" else("PUBCHEM_IUPAC_TRADITIONAL_NAME")],
                                                        description = c("PubChem"),
                                                        baseformula = gsub(x = blockmatrix[, "PUBCHEM_MOLECULAR_FORMULA"], 
                                                                           pattern="[\\+\\-]\\d*$", 
                                                                           replacement=""),
                                                        identifier = as.numeric(blockmatrix[, "PUBCHEM_COMPOUND_CID"]),
                                                        charge = blockmatrix[,"PUBCHEM_TOTAL_CHARGE"]
                                                      )
                                                      # --- check formulae ---
                                                      checked <- as.data.table(check.chemform.joanna(isotopes,
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
                                   read.csv.sql(file, sep="\t", sql = c(fn$paste("insert into base select compoundname, description, baseformula, identifier, charge from file")), dbname = db)
                                 })
                                 print("Done!")
                                 # --- cleanup --- downloading takes forever, would not recommend unless compressing... ---
                                 # remove.residuals <- readline("Done! Remove downloaded SDF / XLS files? (yes/no): ")
                                 # if(remove.residuals == "yes" | remove.residuals == "y") unlink(sdf.loc, recursive=TRUE) else return("Alrighty, files kept.")
                               })
  # --- execute ;) ---
  function.of.choice(dbname)
  RSQLite::dbDisconnect(conn)
}

#' @export
build.extended.db <- function(dbname, 
                              outfolder, 
                              adduct.table, 
                              continue=F, 
                              cl=FALSE, 
                              fetch.limit=-1,
                              cpd.limit=-1){
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
  limit.query <- if(cpd.limit == -1) "" else fn$paste("LIMIT $cpd.limit")
  if(continue){
    RSQLite::dbExecute(base.conn, fn$paste("ATTACH '$full.db' as db"))
    continue.query <- strwrap("SELECT DISTINCT baseformula, charge FROM base b
                              WHERE NOT EXISTS(SELECT DISTINCT baseformula, charge 
                              FROM db.done d
                              WHERE b.baseformula = d.baseformula
                              AND b.charge = d.basecharge)", width=10000, simplify=TRUE)
    total.formulae <- RSQLite::dbGetQuery(base.conn, fn$paste("SELECT Count(*)
                                                     FROM ($continue.query)"))
    formula.count <- total.formulae[1,]
    results <- RSQLite::dbSendQuery(base.conn, continue.query)
  } else{
    # --- add base db to the new one ---
    print("Attaching base...")
    RSQLite::dbExecute(full.conn, fn$paste("ATTACH '$base.db' AS tmp"))
    RSQLite::dbExecute(full.conn, fn$paste("CREATE TABLE IF NOT EXISTS done(baseformula text, basecharge text)"))
    RSQLite::dbExecute(full.conn, fn$paste("CREATE TABLE IF NOT EXISTS base AS SELECT * FROM tmp.base"))
    if(RSQLite::dbExistsTable(base.conn, "pathways")){
      RSQLite::dbExecute(full.conn, fn$paste("CREATE TABLE IF NOT EXISTS pathways AS SELECT * FROM tmp.pathways"))
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
    get.query <- fn$paste("SELECT DISTINCT b.baseformula, b.charge FROM base b $limit.query")
    total.formulae <- RSQLite::dbGetQuery(base.conn, fn$paste("SELECT Count(*)
                                                     FROM ($get.query)"))
    formula.count <- total.formulae[1,]
    results <- RSQLite::dbSendQuery(base.conn, get.query)
  }
  # --- start pb ---
  pb <- startpb(0, formula.count)
  print("Starting DB generation.")
  print(paste("Approximate batches:", formula.count / fetch.limit ))
  # --- waow, my first while in R ---
  while(!dbHasCompleted(results)){
    # --- fetch part of results ---
    partial.results <- as.data.table(RSQLite::dbFetch(results, fetch.limit))
    if(length(partial.results$baseformula) == 0) next
    # -----------------------
    print(paste(RSQLite::dbGetRowCount(results), formula.count, sep=" / "))
    # -----------------------
    checked.formulae <- as.data.table(check.chemform.joanna(isotopes, 
                                                            partial.results$baseformula))
    # -----------------------------
    # keep.rows <- checked.formulae[warning == FALSE & monoisotopic_mass %between% c(60, 600), which=TRUE]
    # -----------------------------
    keep.rows <- checked.formulae[warning == FALSE, which=TRUE]
    # -----------------------------
    if(length(keep.rows) == 0) next
    backtrack <- data.table(baseformula = checked.formulae[keep.rows, new_formula],
                            basemz = checked.formulae[keep.rows, monoisotopic_mass],
                            charge = c(partial.results$charge)[keep.rows])
    # -- go through each adduct --
    do.calc <- function(x){
      row <- adduct.table[x,]
      name <- row$Name
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
      isotopes <- enviPat::isopattern(
        isotopes,
        backtrack$final,
        threshold = 0.1,
        plotit = FALSE,
        charge = backtrack$final.charge,
        algo = 2,
        verbose = FALSE
      )
      # -------------------------
      isolist <- lapply(isotopes, function(isotable){
        if(isotable[[1]] == "error") return(NA)
        iso.dt <- data.table(isotable, fill=TRUE)
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
      meta.table <- data.table(baseformula = rep(backtrack.final$baseformula, repeat.times),
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
    tab.list <- parallel::parLapply(cl=cl, 1:nrow(adduct.table), fun=function(x) do.calc(x))
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
  ppm = 2 ### PPM HARD CODED HERE FOR NOW
  # ------------------------
  mzvals <- data.table(mzmed = total.table$fullmz,
                       foundinmode = total.table$foundinmode)
  mzranges <- data.table(mzmin = pbsapply(total.table$fullmz,cl=cl,
                                        FUN=function(mz, ppm){
                                          mz - mz * (ppm / 1E6)}, ppm=ppm),
                         mzmax = pbsapply(total.table$fullmz,cl=cl,
                                        FUN=function(mz, ppm){
                                          mz + mz * (ppm / 1E6)}, ppm=ppm))
  sql.make.meta <- strwrap("CREATE TABLE mzvals(
                           ID INTEGER PRIMARY KEY AUTOINCREMENT,
                           mzmed decimal(30,13),
                           foundinmode text)", width=10000, simplify=TRUE)
  RSQLite::dbExecute(full.conn, sql.make.meta)
  RSQLite::dbExecute(full.conn, "create index mzfind on mzvals(mzmed, foundinmode);")
  # --- write vals to table ---
  RSQLite::dbWriteTable(full.conn, "mzvals", mzvals, append=TRUE) # insert into
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
  RSQLite::dbExecute(full.conn, sql.make.rtree)
  # --- write ranges to table ---
  RSQLite::dbWriteTable(full.conn, "mzranges", mzranges, append=TRUE) # insert into
  # --- cleanup ---
  RSQLite::dbExecute(full.conn, "VACUUM")
  # ==========================
  print("Disconnecting...")
  RSQLite::dbDisconnect(base.conn)
  RSQLite::dbDisconnect(full.conn)
  on.exit(closepb(pb))
  print("Done! :-)")
}

#' @export
build.pat.db <- function(db.name, 
                         poslist, 
                         neglist, 
                         overwrite=FALSE,
                         rtree=TRUE,
                         ppm=2,
                         rmv.cols=c("mzmin", 
                                    "mzmax", 
                                    "npeaks",
                                    "fq.worst",
                                    "fq.best",
                                    "avg.int")){

  # ------------------------
  ppm = as.numeric(ppm)
  # ------------------------
  print(length(unique(poslist$mzmed)))
  mzvals <- data.table(mzmed = c(poslist$mzmed, neglist$mzmed),
                       foundinmode = c(rep("positive", nrow(poslist)), rep("negative", nrow(neglist))))
  mzranges <- data.table(mzmin = sapply(c(poslist$mzmed, neglist$mzmed), 
                                        FUN=function(mz, ppm){
                                          mz - mz * (ppm / 1E6)}, ppm=ppm),
                         mzmax = sapply(c(poslist$mzmed, neglist$mzmed), 
                                        FUN=function(mz, ppm){
                                          mz + mz * (ppm / 1E6)}, ppm=ppm))
  mzintensities <- reshape2::melt(as.data.table(rbind(poslist, neglist))[,(rmv.cols) := NULL], 
                        id="mzmed", 
                        variable="filename",
                        value="intensity")
  
  mzvals$foundinmode <- trimws(mzvals$foundinmode)
  mzintensities$filename <- trimws(mzintensities$filename)
  
  #filenames <- gsub(x=as.character(mzintensities$filename), "", pattern="[-_]\\d*$", perl=T)
  #replicates <- gsub(x=as.character(mzintensities$filename), "", pattern=".*[-_](?=\\d*$)", perl=T)
  # ------------------------
  #mzintensities$filename <- filenames
  #mzintensities$replicate <- replicates
  # ------------------------
  if(overwrite==TRUE & file.exists(db.name)) file.remove(db.name)
  # --- reconnect / remake ---
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), db.name)
  # ------------------------
  sql.make.int <- strwrap("CREATE TABLE mzintensities(
                          ID INTEGER PRIMARY KEY AUTOINCREMENT,
                          mzmed decimal(30,13),
                          filename text,
                          intensity float)", width=10000, simplify=TRUE)
  RSQLite::dbExecute(conn, sql.make.int)
  # --- write intensities to table and index ---
  RSQLite::dbWriteTable(conn, "mzintensities", mzintensities, append=TRUE) # insert into
  RSQLite::dbExecute(conn, "CREATE INDEX intindex ON mzintensities(filename,'mzmed',intensity)")
  # ------------------------
  sql.make.meta <- strwrap("CREATE TABLE mzvals(
                           ID INTEGER PRIMARY KEY AUTOINCREMENT,
                           mzmed decimal(30,13),
                           foundinmode text)", width=10000, simplify=TRUE)
  RSQLite::dbExecute(conn, sql.make.meta)
  RSQLite::dbExecute(conn, "create index mzfind on mzvals(mzmed, foundinmode);")
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
  # --- write ranges to table ---
  RSQLite::dbWriteTable(conn, "mzranges", mzranges, append=TRUE) # insert into
  # --- cleanup ---
  RSQLite::dbExecute(conn, "VACUUM")
  # ----------------
  RSQLite::dbDisconnect(conn)
  print("Made!")}

#' @export
load.excel <- function(path.to.xlsx, 
                       path.to.patdb = patdb,
                       tabs.to.read = c("General",
                                        "Setup",
                                        "Individual Data",
                                        "Pen Data",
                                        "Admin")){
  # --- connect to sqlite db ---
  # path.to.xlsx = "/Users/jwolthuis/Documents/umc/metaboshiny_all/ExcelSheetsDSM/TR/DSM_TR1.xlsx"
  # path.to.patdb = "/Users/jwolthuis/Documents/umc/turkey.db"
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), path.to.patdb)
  # -------------------------------
  data.store <- pbapply::pblapply(tabs.to.read, FUN=function(tab.name){
    tab <- as.data.table(xlsx::read.xlsx(path.to.xlsx, sheetName = tab.name))
    # --- reformat colnames ---
    colnames(tab) <- tolower(gsub(x=colnames(tab), pattern = "\\.$|\\.\\.$", replacement = ""))
    colnames(tab) <- gsub(x=colnames(tab), pattern = "\\.|\\.\\.", replacement = "_")
    colnames(tab) <- gsub(x=colnames(tab), pattern= "sampling_date_dd_mon_yy", "sampling_date")
    # -----------------------------------------------------------------------------------------
    as.data.table(tab,keep.rownames=F)
  })
  # --- convert to data table --- ## make this nicer loooking in the future
  general <- data.store[[1]]
  setup <- data.store[[2]]
  individual.data <- data.store[[3]]
  
  individual.data$sampling_date <- as.factor(as.Date(individual.data$sampling_date, 
                                                     format = "%d-%m-%y"))
  if(is.na(individual.data$sampling_date[1])) levels(individual.data$sampling_date) <- factor(1) 
  print(head(individual.data))
  #names(individual.data) <- make.unique(names(individual.data))
  #individual.data
  individual.data$sampling_date <- as.numeric(as.factor(as.Date(individual.data$sampling_date, format = "%d-%m-%y")))
  pen.data <- data.store[[4]]
  admin <- data.store[[5]]
  
  setup <- as.data.table(apply(setup, MARGIN=2, trimws))
  individual.data <- as.data.table(apply(individual.data, MARGIN=2, trimws))
  general <- as.data.table(apply(general, MARGIN=2, trimws))
  pen.data <- as.data.table(apply(pen.data, MARGIN=2, trimws))
  admin <- as.data.table(apply(admin, MARGIN=2, trimws))
  
  # --- import to patient sql file ---
  #RSQLite::dbWriteTable(conn, "general", general, overwrite=TRUE) # insert into BUGGED FIX LATER
  RSQLite::dbWriteTable(conn, "setup", setup, overwrite=TRUE) # insert into
  RSQLite::dbWriteTable(conn, "individual_data", individual.data, overwrite=TRUE) # insert into
  RSQLite::dbWriteTable(conn, "pen_data", pen.data, overwrite=TRUE) # insert into
  RSQLite::dbWriteTable(conn, "admin", admin, overwrite=TRUE) # insert into
  # --- disconnect ---
  RSQLite::dbDisconnect(conn)
}
