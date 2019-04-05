#' @export
build.base.db <- function(dbname=NA,
                          outfolder=getOptions()$db_dir,
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
  print(dbname)
  db.formatted <- switch(tolower(dbname),
                         maconda = function(dbname, ...){ #ok
                           file.url = "https://www.maconda.bham.ac.uk/downloads/MaConDa__v1_0__csv.zip"

                           base.loc <- file.path(getOptions()$db_dir, "maconda_source")

                           if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
                           zip.file <- file.path(base.loc, "maconda.zip")
                           utils::download.file(file.url, zip.file,mode = "w")
                           utils::unzip(zip.file,files = "MaConDa__v1_0__extensive.csv",exdir = base.loc)

                           base.table <- data.table::fread(file = file.path(base.loc, "MaConDa__v1_0__extensive.csv"))

                           mysterious = which(base.table$name == "unknown")

                           base.table$formula[mysterious] <- paste0("IDK", 1:length(mysterious))

                           db.base <- data.table::data.table(compoundname = base.table$name,
                                                             description = paste(base.table$type_of_contaminant,
                                                                                 base.table$ion_source_type,
                                                                                 base.table$ion_mode),
                                                             baseformula = base.table$formula,
                                                             identifier=base.table$id,
                                                             charge=c(0),
                                                             structure=base.table$std_inchi)

                           db.base <- unique(db.base)
                           desc.orig <- db.base$description
                           desc.gsub.1 <- gsub(desc.orig, pattern = "^ |^  ", replacement = "")
                           desc.gsub.2 <- gsub(desc.gsub.1, pattern = "  ", replacement = " ")
                           desc.gsub.3 <- gsub(desc.gsub.2, pattern = "POS", replacement = "found in positive mode.")
                           desc.gsub.4 <- gsub(desc.gsub.3, pattern = "NEG", replacement = "found in negative mode.")
                           db.base$description <- Hmisc::capitalize(desc.gsub.4)

                           #RSQLite::dbWriteTable(conn, "base", db.base, append=TRUE)

                           ### do extended table and write to additional db (exception...)
                           db.full <- file.path(outfolder, "extended.db")
                           full.conn <- RSQLite::dbConnect(RSQLite::SQLite(), db.full)
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
                                                    foundinmode text,
                                                    source text)", width=10000, simplify=TRUE)
                           RSQLite::dbExecute(full.conn, sql.make.meta)
                           # -- reformat noise table ---
                           db.formatted <-
                             data.table(
                               baseformula = base.table$formula,
                               fullformula = base.table$formula,
                               basemz = base.table$exact_mass,
                               fullmz = base.table$exact_adduct_mass,
                               adduct = gsub(x = base.table$ion_form, pattern = "^\\[|\\]\\+|\\]\\-", replacement = ""),
                               basecharge = 0,
                               totalcharge = sapply(base.table$ion_mode, function(x){ if(x == "POS") 1 else -1 }),
                               isoprevalence = 100,
                               foundinmode = sapply(base.table$ion_mode, function(x){ if(x == "POS") "positive" else "negative" }),
                               source = c('maconda')
                             )

                           missing <- db.formatted$basemz == 0
                           db.formatted.full <- db.formatted[!missing,]

                           # --- write to db ---

                           RSQLite::dbWriteTable(full.conn,
                                                 "extended",
                                                 db.formatted.full,
                                                 append=TRUE)
                           # --- index ---

                           db.base
                         },
                         hmdb = function(dbname){ #ok
                           print("Downloading XML database...")
                           file.url <- "http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip"
                           # ----
                           base.loc <- file.path(getOptions()$db_dir, "hmdb_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
                           zip.file <- file.path(base.loc, "HMDB.zip")
                           utils::download.file(file.url, zip.file,mode = "w")
                           utils::unzip(zip.file, exdir = base.loc)

                           # --- go through xml ---

                           print("Parsing XML...")

                           require(XML)

                           input = file.path(base.loc,"hmdb_metabolites.xml")

                           library(XML)
                           library(RCurl)
                           library(rlist)
                           theurl <- getURL("http://www.hmdb.ca/statistics",.opts = list(ssl.verifypeer = FALSE) )
                           tables <- readHTMLTable(theurl)
                           stats = data.table::rbindlist(tables)
                           n = as.numeric(as.character(gsub(x = stats[Description == "Total Number of Metabolites"]$Count, pattern = ",", replacement="")))

                           db.formatted <- data.frame(
                             compoundname = rep(NA, n),
                             baseformula = rep(NA, n),
                             identifier = rep(NA, n),
                             structure = rep(NA, n),
                             charge = rep(NA, n),
                             description = rep(NA, n)
                           )

                           pb <- pbapply::startpb(min = 0, max = n)

                           idx <<- 0

                           metabolite = function(currNode){

                             if(idx %% 1000 == 0){
                               pbapply::setpb(pb, idx)
                             }

                             idx <<- idx + 1

                             currNode <<- currNode

                             db.formatted[idx, "compoundname"] <<- xmlValue(currNode[['name']])
                             db.formatted[idx, "identifier"] <<- xmlValue(currNode[['accession']])
                             db.formatted[idx, "baseformula"] <<- xmlValue(currNode[['chemical_formula']])
                             db.formatted[idx, "structure"] <<- xmlValue(currNode[['smiles']])
                             db.formatted[idx, "description"] <<- paste("HMDB:",
                                                                        xmlValue(currNode[['cs_description']]),
                                                                        "CHEMSPIDER:",
                                                                        xmlValue(currNode[['description']])
                             )
                             x <- currNode[['predicted_properties']]
                             properties <- currNode[['predicted_properties']]
                             db.formatted[idx, "charge"] <<- str_match(xmlValue(properties),
                                                                       pattern = "formal_charge([+|\\-]\\d*|\\d*)")[,2]
                           }

                           xmlEventParse(input, branches =
                                           list(metabolite = metabolite))

                           # - - - - - - - - - - - -
                           db.formatted
                         },
                         chebi = function(dbname, ...){ #ok

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

                           # ----------------------
                           db.formatted
                         },
                         wikipathways = function(dbname){ #ok
                           chebi.loc <- file.path(getOptions()$db_dir, "chebi.base.db")
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
                                                               b.identifier,
                                                               w.description,
                                                               b.baseformula,
                                                               b.charge,
                                                               b.structure
                                                               FROM base b
                                                               JOIN wikipathways w
                                                               ON b.identifier = w.identifier")
                           # --- get pathway info ---
                           # -- R PACKAGE EXISTS: source("https://bioconductor.org/biocLite.R"); biocLite("rWikiPathways"); --
                           # sparql.pathways <- SPARQL::SPARQL(url="http://sparql.wikipathways.org/",
                           #                                   query='PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
                           #                                   PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
                           #                                   PREFIX dc:  <http://purl.org/dc/elements/1.1/>
                           #                                   PREFIX foaf: <http://xmlns.com/foaf/0.1/>
                           #                                   PREFIX schema: <http://schema.org/>
                           #                                   PREFIX wp:      <http://vocabularies.wikipathways.org/wp#>
                           #                                   PREFIX dcterms:  <http://purl.org/dc/terms/>
                           #
                           #                                   SELECT DISTINCT str(?titleLit) as ?name ?identifier ?ontology
                           #                                   WHERE {
                           #                                   ?pathway dc:title ?titleLit .
                           #                                   ?pathway dc:identifier ?identifier .
                           #                                   OPTIONAL {?pathway wp:ontologyTag ?ontology .}
                           #                                   } ')
                           # db.pathways <- sparql.pathways$results
                           db.formatted$identifier <- gsub(db.formatted$identifier, pattern = "<|>", replacement = "")
                           #db.formatted$pathway <- gsub(db.formatted$pathway, pattern = "<|_.*", replacement = "")
                           #db.pathways$identifier <- gsub(db.pathways$identifier, pattern = "<|>", replacement = "")
                           # ---------------------------------
                           #RSQLite::dbWriteTable(conn, "pathways", db.pathways, overwrite=TRUE)
                           RSQLite::dbRemoveTable(conn.chebi, "wikipathways")

                           db.formatted
                         },
                         bloodexposome = function(dbname){ #ok

                           file.url = "http://exposome.fiehnlab.ucdavis.edu/blood_exposome_database_v1_08_2018.xlsx"

                           base.loc <- file.path(getOptions()$db_dir, "bloodexposome_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc)
                           excel.file <- file.path(base.loc, "exposome.xlsx")
                           utils::download.file(file.url, excel.file)

                           db.full <- openxlsx::read.xlsx(excel.file, sheet = 2, colNames=T, startRow = 3)

                           db.formatted <- unique(data.table::data.table(compoundname = db.full$CompoundName,
                                                                         description = c("No description available."),
                                                                         baseformula = db.full$Formula,
                                                                         identifier = db.full$PubChemCID,
                                                                         charge = c(0),
                                                                         structure = db.full$SMILES))

                           # --- get charges from smiles ---

                           charges = pbapply::pbsapply(db.formatted$structure, cl=session_cl, FUN=function(smile){
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

                           db.formatted$charge <- charges

                           # --- check formulae ---
                           checked <- data.table::as.data.table(check.chemform.joanna(isotopes,
                                                                                      db.formatted$baseformula))
                           db.formatted$baseformula <- checked$new_formula
                           keep <- checked[warning == FALSE, which = TRUE]
                           db.formatted <- db.formatted[keep]

                           # ----------------------

                           db.formatted

                         },
                         smpdb = function(dbname){ #ok
                           # ---------------
                           file.url <- "http://smpdb.ca/downloads/smpdb_metabolites.csv.zip"
                           # ----
                           base.loc <- file.path(getOptions()$db_dir, "smpdb_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc)
                           zip.file <- file.path(base.loc, "SMPDB.zip")
                           utils::download.file(file.url, zip.file)
                           utils::unzip(zip.file, exdir = base.loc)
                           # -------------------------------

                           smpdb.paths <- list.files(path = base.loc, pattern = "\\.csv$", full.names = T)

                           subtables <- pbapply::pblapply(smpdb.paths, cl = session_cl, fread)
                           smpdb.tab <- rbindlist(subtables, fill = T)

                           db.formatted <- data.table::data.table(compoundname = smpdb.tab$`Metabolite Name`,
                                                                  description = paste0(smpdb.tab$`Pathway Name`,
                                                                                       " (",
                                                                                       smpdb.tab$`Pathway Subject`,
                                                                                       " pathway)"),
                                                                  baseformula = smpdb.tab$Formula,
                                                                  identifier = smpdb.tab$`Metabolite ID`,
                                                                  structure = smpdb.tab$SMILES)


                           db.formatted <- db.formatted[,.(description=paste0(unique(description),collapse=", ")),by=list(compoundname, baseformula, identifier, structure)]
                           db.formatted <- db.formatted[-1,]

                           igroups = split(1:nrow(db.formatted), ceiling(seq_along(1:nrow(db.formatted))/200))

                           db.rows <- pbapply::pblapply(igroups, cl = F, function(is, db.formatted){
                             db.frag = db.formatted[is,]
                             iatoms = rcdk::parse.smiles(db.frag$structure)

                             db.frag$charge <- sapply(iatoms, function(iat){
                               charge = 0 # set default charge to zero
                               try({
                                   charge = rcdk::get.total.charge(iat)
                                 })
                             })

                             db.frag
                           }, db.formatted = db.formatted)

                           db.formatted <- rbindlist(db.rows)

                           db.formatted
                         },
                         kegg = function(dbname){ #ok
                           # ---------------
                           batches <- split(0:2300, ceiling(seq_along(0:2300)/100))
                           cpds <- pbapply::pblapply(batches, cl=session_cl, FUN=function(batch){
                             names(KEGGREST::keggFind("compound", batch, "mol_weight"))
                           })
                           cpd.ids <- Reduce(c, cpds)
                           id.batches <- split(cpd.ids, ceiling(seq_along(cpd.ids)/10))

                           # --- GET COMPOUNDS ---
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
                                                      #widentifier = c(paste(cpd$DBLINKS, collapse=";")),
                                                      charge = c(kegg.charge(cpd$ATOM)),
                                                      structure = c(NA)
                                                      #,pathway = if("PATHWAY" %in% names(cpd)) names(cpd$PATHWAY) else{NA}
                               )
                             })
                             rbindlist(base.list)
                           })

                           db.formatted <- rbindlist(kegg.cpd.list)

                           db.formatted$baseformula <- gsub(pattern = " |\\.", replacement = "", db.formatted$baseformula) # fix salt and spacing issues

                           # checked <- data.table::as.data.table(check.chemform.joanna(isotopes,
                           #                                                            db.formatted$baseformula))
                           # db.formatted$baseformula <- checked$new_formula
                           # keep <- checked[warning == FALSE, which = TRUE]
                           # db.formatted <- db.formatted[keep]
                           # --- GET PATHWAYS ---
                           # pathways <- unique(db.formatted$pathway)
                           # pw.batches <- split(pathways, ceiling(seq_along(pathways)/10))
                           # kegg.pw.list <- pbapply::pblapply(pw.batches, cl=session_cl, FUN=function(batch){
                           #   rest.result <- KEGGREST::keggGet(batch)
                           #   # ---------------------------
                           #   base.list <- lapply(rest.result, FUN=function(pw){
                           #     data.table::data.table(identifier = as.character(pw$ENTRY),
                           #                            name = as.character(pw$NAME),
                           #                            class = if("CLASS" %in% names(pw)) pw$CLASS else{NA},
                           #                            module = if("MODULE" %in% names(pw)) pw$MODULE else{NA},
                           #                            disease = if("DISEASE" %in% names(pw)) as.character(pw$DISEASE) else{NA})
                           #   })
                           #   rbindlist(base.list)
                           # })
                           # db.pathways <- rbindlist(kegg.pw.list)
                           # -------------------------------------------------------

                           #RSQLite::dbWriteTable(conn, "pathways", db.pathways, overwrite=TRUE)

                           db.formatted
                         },
                         metacyc = function(dbname, ...){ #ok
                           # NOTE: Requires downloading this SmartTable as delimited file: https://metacyc.org/group?id=biocyc17-31223-3729417004
                           # May need to remake smartTable if anything on the website changes unfortunately
                           # TODO: download file directly from link, will need a javascript. Maybe Rselenium??

                           base.loc <- file.path(getOptions()$db_dir, "metacyc_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc)

                           source.file = file.path(base.loc, "All_compounds_of_MetaCyc.txt")
                           if(!file.exists(source.file)){
                             message("Please download SmartTable from 'https://metacyc.org/group?id=biocyc17-31223-3729417004' as 'All_compounds_of_MetaCyc.txt' and save in the backend/db/metacyc_source folder.")
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

                           # pathways <- pbapply::pbsapply(metacyc.raw$Pathways.of.compound, cl=session_cl, FUN=function(pw){
                           #   pw <- unlist(strsplit(pw, split = '\\"'))
                           #   pw <- pw[pw != " // "]
                           #   pw <- gsub(pw, pattern = "&", replacement="")
                           #   pw <- gsub(pw, pattern = ";", replacement="")
                           #   res <- gsub(pw, pattern = "<((i|\\/i)|sub)>|\\/|\\|", replacement = "",perl = T)
                           #   paste0(res, collapse=" --- ")
                           # })
                           #
                           # # give the pathwys some identifiers
                           # uniq.pws <- unique(unlist(pbapply::pbsapply(pathways, cl=session_cl, function(x){unlist(strsplit(x, split = " --- "))})))
                           #
                           # db.pathways <- data.table::data.table(name = uniq.pws,
                           #                                       identifier = paste0("METACYC_PW_", 1:length(uniq.pws)))

                           db.formatted <- data.table(compoundname = compounds,
                                                      description = metacyc.raw$Summary,
                                                      baseformula = metacyc.raw$Chemical.Formula,
                                                      identifier = paste0("METACYC_CP_", 1:length(charges)),
                                                      charge = charges,
                                                      structure = metacyc.raw$SMILES
                                                      #, pathway = pathways
                           )

                           # # --- fix the multiple pathway thingy ---
                           # db.fixed.rows <- pbapply::pblapply(1:nrow(db.formatted), cl = session_cl, FUN=function(i, db.cpd = db.cpd, db.pw = db.pw){
                           #   row <- db.cpd[i,]
                           #   pathways <- unlist(strsplit(row$pathway, split = " --- "))
                           #   if(length(pathways) == 0) pathways <- c(1)
                           #   pw.ids <- sapply(pathways, function(pw) db.pw[name == pw]$identifier)
                           #   res <- data.table(compoundname = rep(row$compoundname, length(pathways)),
                           #                     description = rep(row$description, length(pathways)),
                           #                     baseformula = rep(row$baseformula, length(pathways)),
                           #                     identifier = rep(row$identifier, length(pathways)),
                           #                     charge = rep(row$charge, length(pathways)),
                           #                     structure = rep(row$structure, length(pathways)),
                           #                     pathway = pw.ids)
                           #   # - - - - - - - - - -
                           #   res
                           # }, db.cpd=db.formatted, db.pw = db.pathways)

                           #db.formatted <- rbindlist(db.fixed.rows)

                           #db.formatted$pathway <- as.character(db.formatted$pathway)

                           #RSQLite::dbWriteTable(conn, "pathways", db.pathways, overwrite=TRUE)

                           db.formatted
                         },
                         lipidmaps = function(dbname, ...){ #ok

                           file.url = "http://www.lipidmaps.org/resources/downloads/LMSDFDownload3Jan19.zip"

                           # ----
                           base.loc <- file.path(getOptions()$db_dir, "lipidmaps_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc)
                           zip.file <- file.path(base.loc, "lipidmaps.zip")
                           utils::download.file(file.url, zip.file)
                           utils::unzip(zip.file, exdir = base.loc)
                           # -------------------------------

                           sdf.path <- list.files(base.loc,
                                                  pattern = "sdf",
                                                  full.names = T)

                           desc <- function(sdfset){

                             datablock = data.table::as.data.table(datablock2ma(datablocklist=datablock(sdfset)))

                             last_cn <<- colnames(datablock)

                             if(!("FORMULA" %in% last_cn)){
                               mat = as.matrix(data.table::data.table(
                                 identifier=datablock$PUBCHEM_COMPOUND_CID,
                                 compoundname = datablock$PUBCHEM_IUPAC_NAME,
                                 baseformula = datablock$PUBCHEM_MOLECULAR_FORMULA,
                                 structure = datablock$PUBCHEM_OPENEYE_CAN_SMILES,
                                 description = paste0("Alternative names: ",
                                                      apply( datablock[ , grep(colnames(datablock), pattern="NAME"),with=F ] , 1 , paste , collapse = "-" ))
                               ))
                             }else{
                               mat = as.matrix(data.table::data.table(
                                 identifier = as.character(datablock$LM_ID),
                                 compoundname = as.character(if("NAME" %in% colnames(datablock)) datablock$NAME else datablock$SYSTEMATIC_NAME),
                                 baseformula = as.character(datablock$FORMULA),
                                 structure = as.character(if("SMILES" %in% last_cn) datablock$SMILES else datablock$INCHI),
                                 description = as.character(paste0("Main class: ", datablock$MAIN_CLASS,
                                                                   ". Subclass: ", datablock$SUB_CLASS))
                               ))
                             }

                             return(mat)
                           }

                           sdfStream.joanna(input=sdf.path, output=file.path(base.loc, "lipidmaps_parsed.csv"), append=FALSE, fct=desc, silent = T)

                           db.base <- data.table::fread(file.path(base.loc, "lipidmaps_parsed.csv"), fill = T, header=T)

                           db.base$charge <- pbapply::pbsapply(db.base$structure, cl = 0, function(smi){
                             charge = 0 # set default charge to zero
                             try({
                               iatom <- rcdk::parse.smiles(smi)[[1]]
                               charge = rcdk::get.total.charge(iatom)
                             })
                             charge
                           })

                           db.formatted <- db.base[,-1,with=F]

                           db.formatted

                         },
                         metabolights = function(dbname, ...){ # times out...
                           require(webchem)
                           require(rcdk)
                           require(curl)

                           all_ids <- read.table("ftp://ftp.ebi.ac.uk/pub/databases/metabolights/eb-eye/unichem.tsv", sep="\t")
                           colnames(all_ids) <- c("identifier", "inchi", "inchikey")

                           metabs <- all_ids$identifier

                           metabs <- pbapply::pblapply(metabs, cl = F, FUN=function(id){
                             met_info = NA

                             print(id)

                             try({
                               url <- paste0("https://www.ebi.ac.uk/metabolights/webservice/beta/compound/", id)
                               tries = 4
                               while(!is.list(met_info) & tries > 0){
                                 #if(tries < 4) print(paste0("retrying ",id))
                                 met_info <- jsonlite::fromJSON(txt = url)
                                 tries = tries - 1
                               }
                             })
                             Sys.sleep(time = .1)
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

                           db.formatted

                         }, dimedb = function(dbname, ...){ #ok
                           files = c(#"structures.zip",
                             "dimedb_pathways.zip",
                             "dimedb_sources.zip",
                             "dimedb_pc_info.zip",
                             "dimedb_id_info.zip")
                           file.url <- "https://dimedb.ibers.aber.ac.uk/help/downloads/"
                           file.urls <- paste0(file.url, files)
                           # ----
                           print("Downloading files...")
                           base.loc <- file.path(getOptions()$db_dir, "dimedb_source")
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
                                                      structure = casted$SMILES
                                                      #pathway = c(NA)
                           )

                           # - - -

                           db.formatted[which(db.formatted$description == "-")] <- c("Unknown")
                           db.formatted$description <- gsub(db.formatted$description, pattern = "^-|-$", replacement = "")

                           # - - -

                           db.formatted

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
                                                                               charge = 0,
                                                                               baseformula = NA)
                                                              # - - - - -
                                                              try({
                                                                iatom <- rcdk::parse.smiles(smi)[[1]]
                                                                #charge = rcdk::get.total.charge(iatom)
                                                                formula = rcdk::get.mol2formula(iatom)
                                                                row = data.table(canonical_SMILES = smi,
                                                                                 charge = formula@charge,
                                                                                 baseformula = as.character(formula@string))
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
                             #print(x)
                             paste(x[!is.na(x)], collapse=", ")
                           })

                           db.formatted <- data.table::data.table(compoundname = as.character(db.3$chemical_compoundLabel),
                                                                   description = as.character(db.3$description),
                                                                   baseformula = as.character(db.3$baseformula),
                                                                   identifier= as.character(db.3$chemical_compound),
                                                                   charge= as.character(db.3$charge),
                                                                   structure = as.character(db.3$canonical_SMILES))

                           #from = "\u2080\u2081\u2082\u2083\u2084\u2085\u2086\u2087\u2088\u2089"
                           #to = "0123456789"
                           #db.formatted$baseformula <- chartr(from, to, db.formatted$baseformula)

                           db.formatted$description[db.formatted$description == ""] <- "Unknown"
                           db.formatted$baseformula <- as.character(db.formatted$baseformula)
                           db.formatted <- db.formatted[!is.na(db.formatted$baseformula),]

                           # - - write - -

                           db.formatted

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
                                                                  identifier= table_main$abbreviation,
                                                                  charge= table_main$charge,
                                                                  structure= table_main$smile,
                                                                  isHuman = table_main$isHuman,
                                                                  isMicrobe = table_main$isMicrobe)

                           #filter for weights that are way too much w/ weird formulas

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

                           # check integrity of formulae

                           checked <- data.table::as.data.table(check.chemform.joanna(isotopes,
                                                                                      db.formatted$baseformula))
                           db.formatted$baseformula <- checked$new_formula

                           keep <- checked[warning == FALSE, which = TRUE]

                           db.formatted <- db.formatted[keep,]

                           db.formatted

                         }, mona = function(dbname, ...){
                           json <- jsonlite::read_json("~/Downloads/MoNA-export-All_Spectra.json")

                         }, knapsack = function(dbname, ...){
                           # TOO SLOW
                           require(data.table)
                           url <- "http://kanaya.naist.jp/KNApSAcK/"
                           response <- XML::htmlParse(url)
                           tab <- as.data.table(XML::readHTMLTable(response)[[1]])
                           metab_count <- as.numeric(gsub(tab[V1 == "metabolite", 2], pattern = " entries", replacement = ""))

                           #url2 = "http://kanaya.naist.jp/knapsack_jsp/information.jsp?sname=C_ID&word=C00001001"
                           #response <- XML::htmlParse(url2)


                           url <- "http://kanaya.naist.jp/knapsack_jsp/information.jsp?sname=C_ID&word="
                           hrefs <- list()

                           max_digits = 8 #C00000001
                           ids <- sapply(1:metab_count, function(i){
                             paste0("C", stringr::str_pad(i, max_digits, pad = "0"))
                           })
                           #cl = makeCluster(3, "FORK")
                           responses = pbapply::pblapply(ids, cl=0, function(id){
                             print(id)
                             response <- XML::readHTMLTable(paste0(url,id))
                             # - - - return - - -
                             response
                           })



                         }, respect = function(dbname, ...){
                           # - - download reSpect database, phytochemicals - -

                           file.url <- "http://spectra.psc.riken.jp/menta.cgi/static/respect/respect.zip"

                           # ----
                           #base.loc <- getOptions()$db_dir
                           base.loc <- file.path(getOptions()$db_dir, "respect_source")
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

                           db.formatted
                         },
                         expoexplorer = function(dbname, ...){

                           file.url <- "http://exposome-explorer.iarc.fr/system/downloads/current/biomarkers.csv.zip"

                           base.loc <- file.path(getOptions()$db_dir, "exex_source")

                           if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
                           zip.file <- file.path(base.loc, "expoexpo_comp.zip")
                           utils::download.file(file.url, zip.file,mode = "w")
                           utils::unzip(zip.file, exdir = base.loc)

                           base.table <- data.table::fread(file = file.path(base.loc, "biomarkers.csv"))

                           db.formatted <- data.table::data.table(compoundname = base.table$Name,
                                                                  description = base.table$Description,
                                                                  baseformula = base.table$Formula,
                                                                  identifier= base.table$ID,
                                                                  charge= c(NA),
                                                                  structure= base.table$SMILES)

                           db.formatted <- unique(db.formatted)

                           # - - use correlations to get some custom descriptions :) - -

                           file.url <- "http://exposome-explorer.iarc.fr/system/downloads/current/correlation_values.csv.zip"

                           zip.file <- file.path(base.loc, "expoexpo_corr.zip")
                           utils::download.file(file.url, zip.file,mode = "w")
                           utils::unzip(zip.file, exdir = base.loc)

                           corr.table <- data.table::fread(file = file.path(base.loc, "correlation_values.csv"))

                           descriptions <- pbapply::pbsapply(1:nrow(corr.table), function(i){
                             row = corr.table[i,]
                             desc <- paste("Found in", R.utils::decapitalize(row$Biospecimen),
                                           "of", R.utils::decapitalize(row$`Subject group`),
                                           "under", R.utils::decapitalize(row$Population),
                                           "in", row$Country,
                                           "after taking in", R.utils::decapitalize(row$Intake),
                                           paste0("(", row$`Analytical method`, ", p ", row$`Correlation p-value`, ")."))
                           })

                           corr.table$`Pasted` <- descriptions

                           df <- corr.table[,c("Excretion ID", "Pasted")]
                           aggr = aggregate( Pasted ~ `Excretion ID`, df, function(x) toString(paste(unique(x),collapse = " ")))

                           final.table <- merge(db.formatted, aggr, by.x = "identifier", by.y = "Excretion ID", all.x=T)

                           final.table$description <- pbapply::pbsapply(1:nrow(final.table), function(i){
                             row = final.table[i,]
                             a = if(!is.na(row$description) & row$description != "NA") row$description else ""
                             b = row$Pasted
                             paste0(a,b)
                           })

                           db.formatted <- final.table[,-"Pasted"]

                           # - - -

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

                           db.formatted

                         },
                         foodb = function(dbname, ...){ #ok

                           file.url <- "http://www.foodb.ca/system/foodb_2017_06_29_csv.tar.gz"

                           base.loc <- file.path(getOptions()$db_dir, "foodb_source")

                           if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
                           zip.file <- file.path(base.loc, "foodb.zip")
                           utils::download.file(file.url, zip.file, mode = 'wb', method = 'libcurl')
                           utils::unzip(zip.file, exdir = base.loc)

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

                           db.formatted
                         },
                         massbank = function(dbname, ...){
                           file.url <- "https://github.com/MassBank/MassBank-data/archive/master.zip"
                           base.loc <- file.path(getOptions()$db_dir, "massbank_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
                           zip.file <- file.path(base.loc, "massbank.zip")
                           utils::download.file(file.url, zip.file,mode = "w", method='libcurl')
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

                           db.formatted

                         }, supernatural = function(dbname, ...){
                           
                           base.loc <- file.path(getOptions()$db_dir, "supernatural_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
                           
                           library(XML)
                           library(RCurl)
                           library(rlist)
                           theurl <- getURL("http://bioinf-applied.charite.de/supernatural_new/",.opts = list(ssl.verifypeer = FALSE) )
                           tables <- readHTMLTable(theurl)
                           stats = data.table::rbindlist(tables)
                           n = as.numeric(
                             gsub(str_match(theurl, pattern="contains (.*) natural compounds")[,2], pattern = ",", replacement = '')
                           )
                           
                           # http://bioinf-applied.charite.de/supernatural_new/src/download_mol.php?sn_id=SN00000001
                           base.url = "http://bioinf-applied.charite.de/supernatural_new/src/download_mol.php?sn_id="
                           id.nr = str_pad(i, 8, pad = "0")
                           id.str = paste0("SN", id.nr)
                           file.url = paste0(base.url, id.str)
                           
                           all.ids = paste0("SN", str_pad(1:n, 8, pad = "0"))
                           
                           print("Downloading molfiles...")
                           
                           pbapply::pbsapply(all.ids, cl = F, function(id, file.url, base.loc){
                             mol.file <- file.path(base.loc, paste0(id, ".mol"))
                             utils::download.file(file.url, mol.file, mode = "w",quiet = T)  
                           }, file.url = file.url, base.loc = base.loc)
                           
                           # this might be too big... mail the ppl if they want to upload the whole thing..
                           
                         })()


  if(dbname != "maconda"){

    formulae <- unique(db.formatted$baseformula)
    # - - check formulae - -
    checked <- data.table::as.data.table(check.chemform.joanna(isotopes,
                                                               formulae))
    conv.table <- data.table(old = formulae, new = checked$new_formula)
    keep <- checked[warning == FALSE, which = TRUE]
    conv.table.keep <- conv.table[keep,]

    order = match(db.formatted$baseformula, conv.table.keep$old)

    db.formatted <- as.data.table(merge(db.formatted, conv.table.keep, by.x = "baseformula", by.y = "old"))
    db.formatted$baseformula <- db.formatted$new

    db.formatted <- db.formatted[!is.na(db.formatted$baseformula),-"new", with=F]

  }

  # - - write - -

  RSQLite::dbWriteTable(conn, "base", db.formatted, append=TRUE)
  RSQLite::dbExecute(conn, "create index b_idx1 on base(baseformula, charge)")

  RSQLite::dbDisconnect(conn)

  print("did base db!")
}

#' @export
build.extended.db <- function(dbname,
                              outfolder=getOptions()$db_dir,
                              adduct.table=adducts,
                              cl = FALSE,
                              fetch.limit = -1){

  data(isotopes, package = "enviPat")
  base.db <- file.path(outfolder, paste0(dbname, ".base.db"))
  full.db <<- file.path(outfolder, paste0("extended.db"))

  base.conn <- RSQLite::dbConnect(RSQLite::SQLite(), base.db)

  # - - - - - - - - - - - - - - - -

  first.db = if(!file.exists(full.db)) T else F

  full.conn <- RSQLite::dbConnect(RSQLite::SQLite(), full.db)

  # --- add base db to the new one ---
  print("Attaching base...")
  RSQLite::dbExecute(full.conn, gsubfn::fn$paste("ATTACH '$base.db' AS tmp"))

  sql.make.meta <- strwrap("CREATE TABLE IF NOT EXISTS extended(
                           baseformula text,
                           fullformula text,
                           basemz decimal(30,13),
                           fullmz decimal(30,13),
                           adduct text,
                           basecharge int,
                           totalcharge int,
                           isoprevalence float,
                           foundinmode text,
                           source text)", width=10000, simplify=TRUE)
  RSQLite::dbExecute(full.conn, sql.make.meta)

  if(first.db){
    # do all of the formulae
    get.query = "SELECT DISTINCT a.baseformula, a.charge FROM tmp.base a"
  }else{
    # only do the ones not yet present in the extended database
    get.query = "SELECT DISTINCT a.baseformula, a.charge FROM tmp.base a LEFT JOIN extended b ON a.baseformula = b.baseformula AND a.charge = b.basecharge WHERE b.baseformula IS NULL AND b.basecharge IS NULL"
  }

  results <- RSQLite::dbGetQuery(full.conn, get.query)
  formula.count <- nrow(results)

  if(formula.count == 0){
    print("Everything already in current database, aborting! ")
    return(NULL)
  }

  # --- start pb ---

  chunks = split(1:nrow(results), ceiling(seq_along(1:nrow(results))/fetch.limit))

  pbapply::pblapply(chunks, function(selection, full.db){

    partial.results <- as.data.frame(results[selection,])


    checked.formulae <- data.table::as.data.table(check.chemform.joanna(isotopes,
                                                                        partial.results$baseformula))

    keep.rows <- checked.formulae[warning == FALSE, which=TRUE]

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

      isotables <- enviPat::isopattern(
        isotopes,
        backtrack$final,
        threshold = 0.1,
        plotit = FALSE,
        charge = backtrack$final.charge,
        algo = 2,
        verbose = FALSE
      )

      isolist <- lapply(isotables, function(isotable){
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
                                           fullmz = isotable$fullmz,
                                           adduct = c(name),
                                           basecharge = rep(backtrack.final$charge, repeat.times),
                                           totalcharge = rep(backtrack.final$final.charge, repeat.times),
                                           isoprevalence = isotable$isoprevalence,
                                           foundinmode = c(mode),
                                           source = c(dbname))

      # - - return - -

      core.conn <- DBI::dbConnect(RSQLite::SQLite(), dbname = full.db)
      res <- DBI::dbSendQuery(core.conn, "PRAGMA busy_timeout=5000;")
      DBI::dbClearResult(res)

      repeat {
        rv <- try(
          DBI::dbWriteTable(core.conn, "extended", meta.table, append= T )
        )
        if(!is(rv, "try-error")) break
      }
      DBI::dbDisconnect(core.conn)
    }

    if(length(cl)==1){
      if(cl == FALSE){
        cl = NULL
      }
    }
    # this stores everyone in MEMORY! use insert statements...
    if(is.null(cl)){
      sapply(1:nrow(adduct.table), FUN=function(x) do.calc(x))
    }else{
      parallel::clusterExport(cl, "full.db")
      parallel::parSapply(cl=cl, 1:nrow(adduct.table), FUN=do.calc)
    }
  }, full.db = full.db)

  # --- indexy ---
  print("Indexing extended table...")

  if(first.db){
    RSQLite::dbExecute(full.conn, "create index e_idx1 on extended(baseformula, basecharge)")
    RSQLite::dbExecute(full.conn, "create index e_idx2 on extended(fullmz, foundinmode)")
  }else{
    RSQLite::dbExecute(full.conn, "ANALYZE extended")
  }

  # --- cleanup ---
  RSQLite::dbExecute(full.conn, "VACUUM")
  # ==========================
  print("Disconnecting...")
  RSQLite::dbDisconnect(base.conn)
  RSQLite::dbDisconnect(full.conn)
  print("Done! :-)")
}

#' @export
build.pat.db <- function(db.name,
                         pospath,
                         negpath,
                         overwrite=FALSE,
                         rtree=TRUE,
                         make.full = TRUE,
                         ppm=2,
                         inshiny=F){


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
  if(inshiny) setProgress(.20)

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
    batch_info = NULL
  }

  gc()

  if(inshiny) setProgress(.30)

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
  if(inshiny) setProgress(.40)

  #poslist$filename <- trimws(poslist$filename)
  #neglist$filename <- trimws(neglist$filename)

  # COMPUTER CANT HANDLE THE RBIND, WRITE THEM SEPERATELY

  #mzintensities = rbind(poslist, neglist)

  # ------------------------

  if(overwrite==TRUE & file.exists(db.name)) file.remove(db.name)

  # --- reconnect / remake ---

  conn <- RSQLite::dbConnect(RSQLite::SQLite(), db.name)

  # ------------------------

  if(!is.null(batch_info)){
    RSQLite::dbWriteTable(conn, "batchinfo", batch_info, overwrite=T) # insert into
  }

  # ------------------------

  sql.make.int <- strwrap("CREATE TABLE mzintensities(
                          ID INTEGER PRIMARY KEY AUTOINCREMENT,
                          mzmed decimal(30,13),
                          filename text,
                          intensity float)", width=10000, simplify=TRUE)

  RSQLite::dbExecute(conn, sql.make.int)

  if(inshiny) setProgress(.60)

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

  if(inshiny) setProgress(.70)

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

  if(inshiny) setProgress(.80)

  # --- write ranges to table ---
  RSQLite::dbWriteTable(conn, "mzranges", mzranges, append=TRUE) # insert into
  # --- cleanup ---
  RSQLite::dbExecute(conn, "VACUUM")

  if(inshiny) setProgress(.90)

  # ----------------
  RSQLite::dbDisconnect(conn)
  print("Made!")}


#' @export
load.metadata.csv <- function(path.to.csv,
                              path.to.patdb){

  #path.to.csv = "~/Downloads/marc_data/meta.csv"
  #path.to.patdb = "~/Downloads/marc_data/marc.db"

  print(path.to.patdb)
  print(path.to.csv)

  # --- connect to sqlite db ---
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), path.to.patdb)

  csv <- data.table::fread(path.to.csv)
  setup <- data.table(group = as.character(unique(csv[,c("group")][[1]])))

  RSQLite::dbWriteTable(conn, "setup", setup, overwrite=TRUE) # insert into
  RSQLite::dbWriteTable(conn, "individual_data", csv, overwrite=TRUE) # insert into

  RSQLite::dbDisconnect(conn)
}

#' @export
load.metadata.excel <- function(path.to.xlsx,
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

db.build.custom <- function(db.name = "MyDb",
                            db.short = "mydb",
                            db.description = "Personal custom database.",
                            db.icon = "www/questionmark.png",
                            outfolder = getOptions()$db_dir,
                            csv){

  db.base = data.table::fread(csv)

  columns =  c("compoundname",
               "description",
               "baseformula",
               "identifier",
               "charge",
               "structure")

  keep.columns <- intersect(columns,colnames(db.base))

  db.formatted <- db.base[, ..keep.columns]

  if(all(is.na(db.formatted$identifier))){
    db.formatted$identifier <- c(1:nrow(db.formatted))
  }

  # check the formulas
  checked <- data.table::as.data.table(check.chemform.joanna(isotopes,
                                                             db.formatted$baseformula))
  db.formatted$baseformula <- checked$new_formula
  keep <- checked[warning == FALSE, which = TRUE]
  db.formatted <- db.formatted[keep]

  # open db
  outfolder <- file.path(outfolder, "custom")
  if(!dir.exists(outfolder)) dir.create(outfolder)

  db <- file.path(outfolder, paste0(db.short, ".base.db"))
  if(file.exists(db)) file.remove(db)
  print(db)
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), db)
  RSQLite::dbExecute(conn, statement = "create table base(compoundname text, description text, baseformula text, identifier text, charge text, structure text)")

  print(db.formatted)
  # create folder for db
  RSQLite::dbWriteTable(conn, "base", db.formatted, overwrite=TRUE)
  RSQLite::dbDisconnect(conn)

  # write metadata to file.. json?
  meta.dbpage =
    list(title = db.name,
         description = db.description,
         image_id = paste0(db.short, "_icon"))

  meta.img =
    list(name = paste0(db.short, "_icon"), path = db.icon, dimensions = c(200, 200))

  save(list = c("meta.img", "meta.dbpage"), file = file.path(outfolder, paste0(db.short, ".RData")))
}
