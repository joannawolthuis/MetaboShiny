#' @export
build.base.db <- function(dbname=NA,
                          outfolder,
                          optfile,
                          cl=FALSE){
    # --- check if user chose something ---
  if(is.na(dbname)) return("~ Please choose one of the following options: HMDB, ChEBI, PubChem, MetaCyc, internal, noise, KEGG! d(>..w.<)b ~")
  # --- make dir if not exists ---
  if(!dir.exists(outfolder)) dir.create(outfolder)
  # --- create .db file ---
  db <- file.path(normalizePath(outfolder), paste0(dbname, ".base.db"))
  print(db)
  if(file.exists(db)) file.remove(db)
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), db)
  RSQLite::dbExecute(conn, statement = "create table base(compoundname text, description text, baseformula text, identifier text, charge text, structure text)")
  # --- create base ---
  db.formatted <- switch(tolower(dbname),
                         maconda = function(dbname, ...){ #ok
                           file.url = "https://www.maconda.bham.ac.uk/downloads/MaConDa__v1_0__csv.zip"
                           
                           base.loc <- file.path(getOptions(optfile)$db_dir, "maconda_source")
                           
                           if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
                           zip.file <- file.path(base.loc, "maconda.zip")
                           utils::download.file(file.url, zip.file,mode = "w",extra = "-k")
                           utils::unzip(normalizePath(zip.file),files = "MaConDa__v1_0__extensive.csv",exdir = normalizePath(base.loc))
                           
                           base.table <- data.table::fread(file = file.path(base.loc, "MaConDa__v1_0__extensive.csv"))
                           
                           mysterious = which(base.table$name == "unknown")
                           
                           base.table$formula[mysterious] <- paste0("IDK", 1:length(mysterious))
                           
                           has.inchi <- which(base.table$std_inchi != "")
                           inchis <- base.table$std_inchi[has.inchi]
                           smiles = pbapply::pbsapply(inchis, webchem::cs_inchi_smiles)
                           
                           iatoms <- pbapply::pbsapply(smiles, function(x){
                             mol=NULL
                             try({
                               mol = rcdk::parse.smiles(x)[[1]]
                               rcdk::do.aromaticity(mol)
                               rcdk::do.typing(mol)
                               rcdk::do.isotopes(mol)
                             })
                             mol
                           })
                           
                           new.smiles <- pbapply::pbsapply(1:length(iatoms), function(i){
                             mol = iatoms[[i]]
                             smi <- if(is.null(mol)) smi = "" else rcdk::get.smiles(mol)#, flavor = rcdk::smiles.flavors("Canonical"))
                           })
                           
                           base.table$smiles = c("")
                           base.table$smiles[has.inchi] <- new.smiles
                           charges <- lapply(gsub(base.table$ion_form, pattern = ".*]", replacement=""), function(ch){
                             if(grepl("\\d", ch)){
                               print("multicharge")
                               mode = if(grepl("+", ch)) "pos" else "neg"
                               num = as.numeric(gsub("\\+|\\-", replacement="", ch))
                               ch = switch(mode,
                                           pos = abs(num),
                                           neg = -1*num)
                             }else{
                               ch = switch(ch,
                                           "+" = 1,
                                           "-" = -1)
                             }
                             ch
                           })
                           charges[sapply(charges, is.null)] <- NA
                           base.table$charge <- unlist(charges)
                           base.table$charge[is.na(base.table$charge)] <- sapply(base.table$ion_mode[is.na(base.table$charge)], function(mode) switch(mode, POS=1, NEG=-1))
                           
                           # fake formulas
                           missing.smi <- which(base.table$smiles=="")
                           fake.smi <- paste0("[", base.table$formula[missing.smi], "]", base.table$charge[missing.smi])
                           base.table$smiles[missing.smi] <- fake.smi
                           
                           db.base <- data.table::data.table(compoundname = base.table$name,
                                                             description = paste(base.table$type_of_contaminant,
                                                                                 base.table$ion_source_type,
                                                                                 base.table$ion_mode),
                                                             baseformula = base.table$formula,
                                                             identifier=base.table$id,
                                                             charge=c(base.table$charge),
                                                             structure=base.table$smiles)
                           
                           db.base <- unique(db.base)
                           desc.orig <- db.base$description
                           desc.gsub.1 <- gsub(desc.orig, pattern = "^ |^  ", replacement = "")
                           desc.gsub.2 <- gsub(desc.gsub.1, pattern = "  ", replacement = " ")
                           desc.gsub.3 <- gsub(desc.gsub.2, pattern = "POS", replacement = "found in positive mode.")
                           desc.gsub.4 <- gsub(desc.gsub.3, pattern = "NEG", replacement = "found in negative mode.")
                           db.base$description <- Hmisc::capitalize(desc.gsub.4)
                           
                           RSQLite::dbWriteTable(conn, "base", db.base, append=TRUE)
                           
                           
                           ### do extended table and write to additional db (exception...)
                           db.full <- file.path(outfolder, "extended.db")
                           full.conn <- RSQLite::dbConnect(RSQLite::SQLite(), db.full)
                           # --- create extended table ---
                           sql.make.meta <- strwrap("CREATE TABLE IF NOT EXISTS extended(
                           struct_id text,
                           fullformula text,
                           fullmz decimal(30,13),
                           adduct text,
                           isoprevalence float,
                           foundinmode text,
                           source text,
                           FOREIGN KEY(struct_id) REFERENCES structures(id))", width=10000, simplify=TRUE)
                           RSQLite::dbExecute(full.conn, sql.make.meta)
                           
                           # get inchi, etc..
                           # -- reformat noise table ---
                           db.formatted <-
                             data.table(
                               struct_id = c(""),
                               fullformula = base.table$formula,
                               fullmz = base.table$exact_adduct_mass,
                               adduct = gsub(x = base.table$ion_form, pattern = "^\\[|\\]\\+|\\]\\-", replacement = ""),
                               isoprevalence = 100,
                               foundinmode = sapply(base.table$ion_mode, function(x){ if(x == "POS") "positive" else "negative" }),
                               source = c('maconda')
                             )
                           
                           done.structures = RSQLite::dbGetQuery(full.conn, "SELECT COUNT(*) FROM structures")[,1]
                           start.id = done.structures + 1
                           db.formatted$struct_id <- seq(start.id, start.id + nrow(db.formatted) - 1, 1)
                           #missing <- db.formatted$basemz == 0
                           #db.formatted.full <- db.formatted[!missing,]
                           
                           partial.results <- data.table::data.table(id = db.formatted$struct_id,
                                                                     smiles = base.table$smiles,
                                                                     source = c("maconda"))
                           
                           DBI::dbAppendTable(full.conn, "structures", partial.results)
                           
                           # --- write to db ---
                           
                           RSQLite::dbWriteTable(full.conn,
                                                 "extended",
                                                 db.formatted,
                                                 append=TRUE)
                           
                         },
                         t3db = function(dbname, ...){
                           # t3db
                           file.url <- "http://www.t3db.ca/system/downloads/current/toxins.csv.zip"
                           # ----
                           # base.loc = "~/MetaboShiny/databases/t3db_source"
                           base.loc <- file.path(getOptions(optfile)$db_dir, "t3db_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
                           zip.file <- file.path(base.loc, "T3DB.zip")
                           utils::download.file(file.url, zip.file,mode = "w")
                           utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
                           db.base <- data.table::fread(file.path(base.loc, "toxins.csv"), fill=T)  
                           
                           db.formatted <- data.table::data.table(
                             compoundname = db.base$Name,
                             baseformula = db.base$`Chemical Formula`,
                             description = db.base$Description,
                             charge = c(0),
                             identifier = db.base$`T3DB ID`,
                             structure = db.base$SMILES 
                           )
                           
                           iatoms = rcdk::parse.smiles(db.formatted$structure)
                           
                           db.formatted$charge <- sapply(iatoms, function(iat){
                             charge = 0 # set default charge to zero
                             try({
                               charge = rcdk::get.total.charge(iat)
                             })
                             charge
                           })
                           
                           db.formatted
                         },
                         hsdb = function(dbname, ...){
                           file.url = "ftp://ftp.nlm.nih.gov/nlmdata/.hsdblease/hsdb.xml.20190528.zip"
                           base.loc <- file.path(getOptions(optfile)$db_dir, "hsdb_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
                           zip.file <- file.path(base.loc, "HSDB.zip")
                           utils::download.file(file.url, zip.file,mode = "w")
                           utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
                           
                           library(XML)
                           library(RCurl)
                           library(rlist)
                           
                           input = list.files(base.loc, pattern = "\\.xml",full.names = T)
                           
                           theurl <- getURL("https://toxnet.nlm.nih.gov/help/hsdbcasrn.html",.opts = list(ssl.verifypeer = FALSE) )
                           
                           n = str_count(as.character(theurl), pattern = "cgi-bin")
                           
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
                           
                           parsed <- XML::xmlTreeParse(input)
                           i=1
                           
                           db.formatted <- data.table::data.table(
                             compoundname = xmlApply(parsed[[1]]$children$hsdb, function(x) xmlValue(x['NameOfSubstance'][[1]])),
                             baseformula = xmlApply(parsed[[1]]$children$hsdb, function(x) xmlValue(x['mf'][[1]])),
                             identifier = xmlApply(parsed[[1]]$children$hsdb, function(x) xmlValue(x['DOCNO'][[1]]))
                           )
                           parsed[[1]]$children$hsdb[[i]]
                           parsed[[1]]$children$hsdb[[i]]['mf'][1]
                           parsed[[1]]$children$hsdb[[i]]['ocpp']
                           parsed[[1]]$children$hsdb[[i]]['sy'] 
                           parsed[[1]]$children$hsdb[[i]]['mf'] 
                           
                           identifier = xmlValue(parsed[[1]]$children$hsdb[[i]]['CASRegistryNumber']$CASRegistryNumber)
                           
                           require(webchem)
                           cas_ids <- xmlApply(parsed[[1]]$children$hsdb, function(x) xmlValue(x['CASRegistryNumber'][[1]]))
                           
                           smiles = pbapply::pbsapply(cas_ids, cl = F, function(id) webchem::cir_query(id, representation = "smiles", resolver = NULL,
                                                                                                       first = FALSE)[[1]])
                           
                           mol = rcdk::parse.smiles(smiles[[1]])
                           charge = rcdk::get.total.formal.charge(mol[[1]])
                           
                           # - - - - - - - - - - - -
                           db.formatted
                           
                           
                         },
                         hmdb = function(dbname){ #ok
                           file.url <- "http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip"
                           # ----
                           base.loc <- file.path(getOptions(optfile)$db_dir, "hmdb_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
                           zip.file <- file.path(base.loc, "HMDB.zip")
                           # utils::download.file(file.url, zip.file,mode = "w")
                           # utils::unzip(zip.file, exdir = base.loc)
                           
                           # --- go through xml ---
                           
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
                           
                           # BiocManager::install("rWikiPathways")
                           # require(rWikiPathways)
                           # orgs = listOrganisms()
                           # for(org in orgs){
                           #   pws = listPathways(org)
                           #   for(pw in pws){
                           #     mets = getXrefList(pw, 'Ce') 
                           #     
                           #   }
                           # }
                           
                           chebi.loc <- file.path(getOptions(optfile)$db_dir, "chebi.base.db")
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
                         bloodexposome = function(dbname){ #currently ded, pls check later
                           
                           file.url = "https://exposome1.fiehnlab.ucdavis.edu/download/BloodExpsomeDatabase_version_1.0.xlsx"
                           
                           base.loc <- file.path(getOptions(optfile)$db_dir, "bloodexposome_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc)
                           excel.file <- file.path(base.loc, "exposome.xlsx")
                           utils::download.file(file.url, excel.file)
                           
                           db.full <- openxlsx::read.xlsx(excel.file, sheet = 1, colNames=T, startRow = 3)
                           
                           db.formatted <- unique(data.table::data.table(compoundname = db.full$Compound.Name,
                                                                         description = paste0("Found in ", db.full$BloodPaperCount, " papers related to blood exposome."),
                                                                         baseformula = db.full$Molecular.Formula,
                                                                         identifier = db.full$PubChem.CID,
                                                                         charge = db.full$Charge,
                                                                         structure = db.full$CanonicalSMILES))
                           
                           db.formatted
                           
                         },
                         smpdb = function(dbname){ #ok
                           # ---------------
                           file.url <- "http://smpdb.ca/downloads/smpdb_metabolites.csv.zip"
                           # ----
                           base.loc <- file.path(getOptions(optfile)$db_dir, "smpdb_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc)
                           zip.file <- file.path(base.loc, "SMPDB.zip")
                           utils::download.file(file.url, zip.file)
                           utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
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
                           
                           kegg.cpd.list <- pbapply::pblapply(id.batches, cl=0, FUN=function(batch){
                             rest.result <- KEGGREST::keggGet(batch)
                             # ---------------------------
                             base.list <- lapply(rest.result, FUN=function(cpd){
                               cpd$NAME_FILT <- gsub(cpd$NAME, pattern = ";", replacement = "")
                               data.table::data.table(compoundname = c(paste(cpd$NAME_FILT, collapse=", ")),
                                                      description = paste0("Involved in pathways: ", 
                                                                           paste0(cpd$PATHWAY, collapse = ", "),
                                                                           ". More specifically: ", 
                                                                           paste0(cpd$MODULE, collapse = ", "),
                                                                           ". Also associated with compound classes: ",
                                                                           paste0(
                                                                             unique(trimws(
                                                                                      gsub(cpd$BRITE, pattern = "\\[.*\\]|  D\\d* |\\(.*\\)|\\d*", replacement= "")
                                                                                      )
                                                                                    ), collapse = ", ")
                                                                           ),
                                                      baseformula = c(cpd$FORMULA),
                                                      identifier = c(cpd$ENTRY),
                                                      charge = 0,
                                                      structure = NA
                                                      ,pathway = if("PATHWAY" %in% names(cpd)) names(cpd$PATHWAY) else{NA}
                               )
                             })
                             res = rbindlist(base.list)
                             res
                           })
                           
                           # - - - 
                           db.formatted <- data.table::rbindlist(kegg.cpd.list)
                           
                           base.loc <- file.path(getOptions(optfile)$db_dir, "kegg_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc)
                           
                           kegg.mol.paths = pbapply::pblapply(id.batches, cl=0, FUN=function(batch){
                              #mols = Rcpi::getMolFromKEGG(batch, parallel = 1)
                              bigmol = KEGGREST::keggGet(batch, "mol")
                              mols = strsplit(x = paste0("\n \n \n",bigmol), split = "\\$\\$\\$\\$\n")[[1]]
                              fps = normalizePath(file.path(base.loc, paste0(str_match(mols, pattern = "<ENTRY>\ncpd:(.*)\n")[,2], ".mol")),mustWork = F)
                              sapply(1:length(mols), function(i) writeLines(text = mols[[i]], 
                                                                              con = fps[i]
                                                                              )
                                     )
                              fps
                           })
                           
                           fns = unlist(flattenlist(kegg.mol.paths))

                           rJava::.jcall("java/lang/System","V","gc")
                           gc()
                           
                           smiles.rows = pbapply::pblapply(fns, function(fn){
                             smiles=NA
                             try({
                               iatom = rcdk::load.molecules(molfiles = fn)[[1]]
                               smiles = rcdk::get.smiles(molecule = iatom)
                               })
                             id = gsub(basename(fn), pattern = "\\.mol", replacement="")
                             data.table::data.table(identifier = id, smiles = smiles)
                             })
                 
                           smitable <- data.table::rbindlist(smiles.rows)
                           
                           db.merged <- merge(db.formatted, smitable, by = "identifier")
                           
                           db.formatted <- data.table(compoundname = db.merged$compoundname,
                                                      description = db.merged$description,
                                                      baseformula = db.merged$baseformula,
                                                      identifier = db.merged$identifier,
                                                      charge = c(0),
                                                      structure = db.merged$smiles
                                                      #, pathway = pathways
                           )
                         },
                         metacyc = function(dbname, ...){ #ok
                           # NOTE: Requires downloading this SmartTable as delimited file: https://metacyc.org/group?id=biocyc17-31223-3729417004
                           # May need to remake smartTable if anything on the website changes unfortunately
                           # TODO: download file directly from link, will need a javascript. Maybe Rselenium??
                           
                           base.loc <- file.path(getOptions(optfile)$db_dir, "metacyc_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc)
                           
                           source.file = file.path(base.loc, "All_compounds_of_MetaCyc.txt")
                           if(!file.exists(source.file)){
                             message("Please download SmartTable from 'https://metacyc.org/group?id=biocyc17-31223-3729417004' as 'All_compounds_of_MetaCyc.txt' and save in the backend/db/metacyc_source folder.")
                             return(NULL)
                           }
                           metacyc.raw = read.table(source.file,header = T, sep = "\t",
                                                    check.names = F,fill=T,quote="")
                           
                           colnames(metacyc.raw) <- gsub(x=as.character(colnames(metacyc.raw)), pattern = '\\"', replacement="")
                           metacyc.raw[] <- lapply(metacyc.raw, gsub, pattern = '\\"', replacement = "")
                           
                           charges = pbapply::pbsapply(metacyc.raw$SMILES, cl=session_cl, FUN=function(smile){
                             m <- rcdk::parse.smiles(smile)
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
                             pw <- iconv(pw, "latin1", "UTF-8",sub='')
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
                                                      baseformula = metacyc.raw$`Chemical Formula`,
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
                         drugbank = function(dbname, ...){
                           
                           # ----
                           
                           base.loc <- file.path(getOptions(optfile)$db_dir, "drugbank_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc)
                           
                           zip.file <- file.path(base.loc, "drugbank.zip")
                           if(!file.exists(zip.file)){
                             file.rename(file.path(base.loc, "drugbank_all_full_database.xml.zip"), zip.file)
                           }
                           utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
                           
                           theurl <- getURL("https://www.drugbank.ca/stats",.opts = list(ssl.verifypeer = FALSE) )
                           tables <- readHTMLTable(theurl,header = F)
                           stats = as.data.table(tables[[1]])
                           colnames(stats) <- c("Description", "Count")
                           n = as.numeric(as.character(gsub(x = stats[Description == "Total Number of Drugs"]$Count, pattern = ",", replacement="")))
                           
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
                             
                             if(idx %% 10 == 0){
                               pbapply::setpb(pb, idx)
                             }
                             
                             idx <<- idx + 1
                             
                             currNode <<- currNode

                             properties <- currNode[['calculated-properties']]
                             
                             if(is.null(properties)){
                               properties <- currNode[['experimental-properties']]
                             }
                             
                             proplist <- xmlToList(properties)
                             if(length(proplist) == 0){
                               return(NULL)
                             }
                             
                             # find formula
                             which.form <- which(sapply(proplist, function(x){
                               if("kind" %in% names(x)){
                                 res = x[['kind']] == "Molecular Formula"  
                               }else{
                                 res = FALSE
                               }
                               res
                               }))
                             
                             which.struc <- which(sapply(proplist, function(x){
                               if("kind" %in% names(x)){
                                 res = x[['kind']] == "SMILES"  
                               }else{
                                 res = FALSE
                               }
                               res
                               }))
                             
                             which.charge <- which(sapply(proplist, function(x){
                               if("kind" %in% names(x)){
                                 res = x[['kind']] == "Physiological Charge"  
                               }else{
                                 res = FALSE
                               }
                               res
                               }))
                             
                             # find structure
                             
                             if(length(which.form) == 0 & length(which.struc) == 0){
                               return(NULL)
                             }
                             
                             db.formatted[idx, "compoundname"] <<- xmlValue(currNode[['name']])
                             db.formatted[idx, "identifier"] <<- xmlValue(currNode[['drugbank-id']])
                             db.formatted[idx, "baseformula"] <<- proplist[[which.form]][['value']]
                             db.formatted[idx, "structure"] <<- if(length(which.struc) > 0){
                               proplist[[which.struc]][['value']]
                             }else{
                               ""
                             }
                             db.formatted[idx, "description"] <<- xmlValue(currNode[['description']])
                             db.formatted[idx, "charge"] <<- if(length(which.charge) > 0){
                               proplist[[which.charge]][['value']]
                               }else{
                                 0
                               }
                           }
                           
                           res = xmlEventParse(file = file.path(base.loc, "full database.xml"), branches =
                                           list("drug" = metabolite, "drugbank-metabolite-id-value" = print))
                           
                           # - - - - - - - - - - - -
                           db.formatted
                           
                         },
                         lipidmaps = function(dbname, ...){ #ok
                           
                           file.url = "http://www.lipidmaps.org/resources/downloads/LMSDFDownload3Jan19.zip"
                           
                           # ----
                           base.loc <- file.path(getOptions(optfile)$db_dir, "lipidmaps_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc)
                           zip.file <- file.path(base.loc, "lipidmaps.zip")
                           utils::download.file(file.url, zip.file)
                           utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
                           # -------------------------------
                           
                           sdf.path <- list.files(base.loc,
                                                  pattern = "sdf",
                                                  full.names = T,
                                                  recursive = T)
                           
                           desc <- function(sdfset){
                             mat <- NULL
                             try({
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
                             })
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
                           
                           # require(RCurl)
                           # 
                           uris <- pbapply::pblapply(metabs, cl = F, FUN=function(id){
                             met_info = NA
                             url <- paste0("https://www.ebi.ac.uk/metabolights/webservice/beta/compound/", id)
                           })
                           # 
                           # all_info <- jsonlite::read_json("~/Downloads/mapping.json")
                           # View(all_info$studyMapping$MTBLS59)
                           # View(all_info$compoundMapping$`CSID 16569894`)
                           # 
                           # for(item in all_info$compoundMapping){
                           #   for(info in item){
                           #     print(info$study)
                           #   }
                           # }
                           # 
                           # metab_rows <- pbapply::pblapply(all_info$compoundMapping, cl = F, FUN=function(cpd){
                           #   met_info <- cpd[[1]]$mafEntry
                           #   data.table(compoundname = met_info$metaboliteIdentification,
                           #              description = if(!is.null(met_info$description)) met_info$description else "unknown",
                           #              baseformula = met_info$chemicalFormula,
                           #              identifier = if(!is.null(met_info$identifier)) met_info$identifier else "unknown",
                           #              structure = met_info$smiles,
                           #              charge = met_info$charge)
                           #   
                           # })
                           #metab_tbl <- rbindlist(metab_rows)
                           
                           
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
                           base.loc <- file.path(getOptions(optfile)$db_dir, "dimedb_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc)
                           pbapply::pbsapply(file.urls, function(url){
                             zip.file <- file.path(base.loc, basename(url))
                             utils::download.file(url, zip.file)
                             utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
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
                             response <- XML::readHTMLTable(paste0(url,id))
                             # - - - return - - -
                             response
                           })
                           
                           
                           
                         },phenolexplorer = function(dbname, ...){
                           
                           file.urls = paste0(
                             "http://phenol-explorer.eu/system/downloads/current/",
                             c("composition-data.xlsx.zip",
                               "compounds.csv.zip",
                               "compounds-structures.csv.zip",
                               "metabolites.csv.zip",
                               "metabolites-structures.csv.zip"))
                           
                           # ---
                           base.loc <- file.path(getOptions(optfile)$db_dir, "phenolexplorer_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc)
                           
                           for(url in file.urls){
                             zip.file <- file.path(base.loc, basename(url))
                             utils::download.file(url, destfile = zip.file)
                             utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
                           }
                             
                           # start from the composition table
                           compo_tbl <- openxlsx::read.xlsx(file.path(base.loc, "composition-data.xlsx"), sheet = 1)
                           compo_tbl$description <- paste0("Present in ", tolower(compo_tbl$food),
                                                           "(", tolower(compo_tbl$food_group), ", ", 
                                                           tolower(compo_tbl$food_sub_group), "). ", 
                                                           "Belongs to the compound class of ",
                                                           tolower(compo_tbl$compound_group),
                                                           " (", tolower(compo_tbl$compound_sub_group), "). ",
                                                           "PMIDS: ", compo_tbl$pubmed_ids)
                           
                           db.base <- unique(compo_tbl[, c("compound", "description")])
                           
                           merged.cpds <- merge(db.base, struct_tbl_keep, by.x = "compound", by.y = "name")
                           
                           # unique structure info
                           struct_tbl <- fread(file.path(base.loc,"compounds-structures.csv"))
                           struct_tbl_keep <- struct_tbl[, c("id", "smiles", "name")]
                           
                           # metabolites
                           met_tbl <- fread(file.path(base.loc,"metabolites.csv"))
                           met_struct <- fread(file.path(base.loc,"metabolites-structures.csv"))
                           met_struct_keep <- met_struct[, c("id", "smiles", "name")]
                           
                           missing = (!(met_tbl$name %in% compo_tbl$compound))
                           mis_mets <- met_tbl[missing,]
                            
                           merged.mets <- merge(mis_mets, met_struct_keep, by.x = "name", by.y = "name")
                           merged.mets <- merged.mets[, c("name", "synonyms", "id.x", "smiles")]
                           colnames(merged.mets) <- c("compound", "description", "id", "smiles")
                           merged.mets$description <- paste0("Metabolite of compound in food.", merged.mets$description)
                           
                           merged <- unique(rbind(merged.cpds, merged.mets))
                           
                           # mtbls
                           iatoms = rcdk::parse.smiles(merged$smiles)
                           
                           merged$charge <- pbapply::pbsapply(iatoms, function(iat){
                             charge = 0 # set default charge to zero
                             try({
                               charge = rcdk::get.total.charge(iat)
                             })
                             charge
                           })
                           
                           merged$baseformula <- pbapply::pbsapply(iatoms, function(iat){
                             charge = 0 # set default charge to zero
                             try({
                               charge = rcdk::get.mol2formula(iat)@string
                             })
                             charge
                           })
                           
                           db.formatted <- data.table::data.table(compoundname = merged$compound,
                                                                  description = merged$description,
                                                                  baseformula = merged$baseformula,
                                                                  identifier= merged$id,
                                                                  charge= merged$charge,
                                                                  structure= merged$smiles)
                           
                           db.formatted <- db.formatted[ , .(description = paste(description, collapse=". ")), by = c("compoundname", "baseformula", "identifier", "charge", "structure")]
                           
                           db.formatted
                         }, 
                         respect = function(dbname, ...){
                           # - - download reSpect database, phytochemicals - -
                           
                           file.url <- "http://spectra.psc.riken.jp/menta.cgi/static/respect/respect.zip"
                           
                           # ----
                           #base.loc <- getOptions(optfile)$db_dir
                           base.loc <- file.path(getOptions(optfile)$db_dir, "respect_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
                           zip.file <- file.path(base.loc, "respect.zip")
                           utils::download.file(file.url, zip.file,mode = "w")
                           utils::unzip(normalizePath(zip.file), exdir = (base.loc))
                           
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
                           
                           base.loc <- file.path(getOptions(optfile)$db_dir, "exex_source")
                           
                           if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
                           zip.file <- file.path(base.loc, "expoexpo_comp.zip")
                           utils::download.file(file.url, zip.file,mode = "w")
                           utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
                           
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
                           utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
                           
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
                           base.loc <- file.path(getOptions(optfile)$db_dir, "foodb_source")
                           
                           if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
                           zip.file <- file.path(base.loc, "foodb.tar.gz")
                           utils::download.file(file.url, zip.file, mode = 'wb', method = 'libcurl')
                           utils::untar(normalizePath(zip.file), exdir = normalizePath(base.loc))                       
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
                           base.loc <- file.path(getOptions(optfile)$db_dir, "massbank_source")
                           if(!dir.exists(base.loc)) dir.create(base.loc,recursive = T)
                           zip.file <- file.path(base.loc, "massbank.zip")
                           utils::download.file(file.url, zip.file,mode = "w", method='libcurl')
                           utils::unzip(normalizePath(zip.file), exdir = normalizePath(base.loc))
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
                           
                           base.loc <- file.path(getOptions(optfile)$db_dir, "supernatural_source")
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
                           
                           pbapply::pbsapply(all.ids, cl = F, function(id, file.url, base.loc){
                             mol.file <- file.path(base.loc, paste0(id, ".mol"))
                             utils::download.file(file.url, mol.file, mode = "w",quiet = T)
                           }, file.url = file.url, base.loc = base.loc)
                           
                           # this might be too big... mail the ppl if they want to upload the whole thing..
                           
                         },
                         toxnetdb=function(dbname, ...){
                           
                         })()
  
  if(dbname != "maconda"){
    db.formatted.bu <<- db.formatted
    uniques =  unique(db.formatted$baseformula)
    uniq.nona <- uniques[!is.na(uniques)]
    checked <- data.table::as.data.table(check.chemform.joanna(isotopes,
                                                               uniq.nona))
    
    conv.table <- data.table(old = uniq.nona, 
                             new = checked$new_formula,
                             warning = checked$warning)
    order = match(db.formatted$baseformula, conv.table$old)
    
    db.formatted <- as.data.table(merge(db.formatted, conv.table, by.x = "baseformula", by.y = "old"))
    db.formatted$baseformula <- db.formatted$new
    
    
    # load in all available SMILES
    iatoms <- pbapply::pbsapply(db.formatted$structure, function(x){
      mol=NULL
      try({
        mol = rcdk::parse.smiles(x)[[1]]
        rcdk::do.aromaticity(mol)
        rcdk::do.typing(mol)
        rcdk::do.isotopes(mol)
      })
      mol
    })
    
    rJava::.jcall("java/lang/System","V","gc")
    gc()
    new.smiles <- pbapply::pbsapply(1:length(iatoms), function(i){
      mol = iatoms[[i]]
      smi <- if(is.null(mol)) smi = "" else rcdk::get.smiles(mol)#, flavor = rcdk::smiles.flavors("Canonical"))
    })
    
    db.formatted$structure <- as.character(new.smiles)
    rJava::.jcall("java/lang/System","V","gc")
    gc()
    
    ok.form <- !db.formatted$warning
    ok.smi <- grepl(db.formatted$structure, pattern="[A-Z]") #db.formatted$structure != "" & db.formatted$structure != "*=*"
    
    # which have formula but ddont have structure?
    yes.form.no.struct <- ok.form & !ok.smi
    
    # which have structure but dont have formula?
    no.form.yes.struct <- !ok.form & ok.smi
    
    # - - generate formula from structure - - 
    if(any(no.form.yes.struct)){
      rJava::.jcall("java/lang/System","V","gc")
      gc()
      generated.formulae <- sapply(which(no.form.yes.struct), function(i){
        mol = iatoms[[i]]
        form = NULL
        try({
          form = rcdk::get.mol2formula(mol)@string
        })
        form
      })
      rJava::.jcall("java/lang/System","V","gc")
      gc()  
      db.formatted$baseformula[no.form.yes.struct] <- unlist(lapply(generated.formulae, function(x) if(is.null(x)) NA else as.character(x)))
    }
    
    ok.form <- ok.form & !is.na(db.formatted$baseformula)
    
    # generate fake SMILES from ones that don't have structure but do have formula and charge
    fake.smiles <- paste0("FORMONLY[", db.formatted$baseformula[yes.form.no.struct], "]",db.formatted$charge[yes.form.no.struct])
    db.formatted$structure[yes.form.no.struct] <- as.character(fake.smiles)
    
    keep = !is.na(db.formatted$baseformula)
    
    # # adduct rules
    # chonks <- split(which(ok.smi), ceiling(seq_along(which(ok.smi))/200))
    # rule.rows <- pbapply::pblapply(chonks, function(chunk){
    #   add.cols <- lapply(1:nrow(adduct_rules), function(i){
    #     curr = adduct_rules$SHORT[i]
    #     query <- adduct_rules$SMARTS[i]
    #     if(query != "charge"){
    #       matches <- rcdk::matches(query = query, target = iatoms[chunk], T)
    #     }        
    #     vals = as.numeric(unlist(lapply(matches, function(x) length(x$mapping))))
    #     vals
    #   })
    #   bound = do.call("cbind", add.cols)
    #   cbind()
    # })
    
    # remove rows with neither
    db.formatted <- db.formatted[keep, -c("new","warning")]
    
    db.formatted$charge <- sapply(db.formatted$structure, function(smi){
      if(grepl("FORMONLY", x = smi)){
        ch = str_match(smi, "\\](.*$)")[[2]]
      }else{
        mol = rcdk::parse.smiles(smi)[[1]]
        ch = rcdk::get.total.formal.charge(molecule = mol)
      }
      ch
    })
    rJava::.jcall("java/lang/System","V","gc")
    gc()  
    
    # reorder just in case
    # compoundname text, description text, baseformula text, 
    # identifier text, charge text, structure text)
    
    db.formatted <- db.formatted[,c("compoundname", "description", "baseformula", 
                                    "identifier", "charge", "structure")]
    
    # - - write - -
    
    RSQLite::dbWriteTable(conn, "base", db.formatted, append=TRUE)
    #RSQLite::dbExecute(conn, "create index b_idx1 on base(baseformula, charge, structure)")
    
    RSQLite::dbDisconnect(conn)
  }
}

#' @export
build.extended.db <- function(dbname,
                              outfolder = getOptions(lcl$paths$opt.loc)$db_dir,
                              adduct.table = adducts,
                              cl = FALSE,
                              fetch.limit = -1,
                              mzrange = c(60, 600)){
  
  outfolder <- normalizePath(outfolder)
  data(isotopes, package = "enviPat")
  base.db <- file.path(outfolder, paste0(dbname, ".base.db"))
  full.db <<- file.path(outfolder, paste0("extended.db"))
  
  first.db = if(!file.exists(full.db)) T else F
  full.conn <- RSQLite::dbConnect(RSQLite::SQLite(), full.db)
  
  # --- add base db to the new one ---
  
  RSQLite::dbExecute(full.conn, gsubfn::fn$paste("PRAGMA foreign_keys = ON"))
  RSQLite::dbExecute(full.conn, gsubfn::fn$paste("ATTACH '$base.db' AS tmp"))
  
  sql.make.meta <- strwrap("CREATE TABLE IF NOT EXISTS extended(
                           struct_id text,
                           fullformula text,
                           fullmz decimal(30,13),
                           adduct text,
                           isoprevalence float,
                           foundinmode text,
                           source text,
                           FOREIGN KEY(struct_id) REFERENCES structures(id))", width=10000, simplify=TRUE)
  
  RSQLite::dbExecute(full.conn, sql.make.meta)
  
  #print("a")
  
  # create ids for the new structures...
  if(first.db){
    RSQLite::dbExecute(full.conn, gsubfn::fn$paste("PRAGMA auto_vacuum = 1"))
    RSQLite::dbExecute(full.conn, "CREATE TABLE IF NOT EXISTS structures(id INT PRIMARY KEY, smiles TEXT, source TEXT)")
    new.structures = RSQLite::dbGetQuery(full.conn, "SELECT structure FROM tmp.base")
  }else{
    new.structures = RSQLite::dbGetQuery(full.conn, "SELECT structure FROM tmp.base WHERE structure NOT IN(SELECT smiles FROM structures)")
  }
  
  #print("b")
  
  if(nrow(new.structures) == 0){
    print("all already done")
    return(NULL)
  }
  
  done.structures = RSQLite::dbGetQuery(full.conn, "SELECT COUNT(*) FROM structures")[,1]
  start.id = done.structures + 1
  
  #print("c")
  
  mapper <- data.table::data.table(id = seq(start.id, start.id + nrow(new.structures) - 1, 1),
                                   smiles = new.structures$structure,
                                   source = c(dbname))
  
  # --- start pb ---
  
  #print("d")
  
  chunks = split(1:nrow(mapper), ceiling(seq_along(1:nrow(mapper))/fetch.limit))
  
  print(cl)
  pbapply::pblapply(chunks, function(selection, cl=cl, full.db){
    
    #print("e")
    
    partial.results <- as.data.frame(mapper[selection,])
    
    smiles <- partial.results$smiles #glucose
    
    valid.struct = !grepl(smiles, pattern = "FORMONLY")
    
    mols <- lapply(1:nrow(partial.results), function(i) return(NA))
    mols[valid.struct] <- rcdk::parse.smiles(smiles[valid.struct])#sapply(smiles, rcdk::parse.smiles)
    
    #print(
    backtrack <- data.table::data.table(structure = smiles,
                                        parse_smile = valid.struct)
    
    # get formula from smiles (or substitute)
    backtrack_molinfo <- lapply(1:nrow(backtrack), function(i){
      row = backtrack[i,]
      mol = mols[[i]]
      if(row$parse_smile){
        # get molecular formula
        mf = rcdk::get.mol2formula(mol)@string
        # get charge
        ch = rcdk::get.total.formal.charge(mol)
      }else{
        # generate.. "FORMONLY[C6H12O6]0"
        mf = str_match(row$structure, "\\[(.*)\\]")[[2]]
        ch = str_match(row$structure, "\\](.*$)")[[2]]
      }
      data.table(baseformula = mf, charge = ch)
    })
    
    #print("g")
    
    backtrack <- cbind(backtrack, rbindlist(backtrack_molinfo))
    
    backtrack$baseformula <- check.chemform.joanna(isotopes, backtrack$baseformula)$new
    
    # adduct rules
    for(i in 1:nrow(adduct_rules)){
      #print("g1")
      
      curr = adduct_rules$SHORT[i]
      query <- adduct_rules$SMARTS[i]
      
      #print("g2")
      
      if(curr == "Nch"){
        backtrack[, (curr) := backtrack$charge]
      }else{
        #print("g3")
        all_matches <- lapply(1:nrow(backtrack), function(i){
          row = backtrack[i,]
          if(row$parse_smile){
            rcdk::matches(query = query, target = mols[[i]], T)
          }
          else{list(mapping=list(a=c(1:10)))}
        })
        #print("g4")
        vals = as.numeric(sapply(all_matches, function(x) length(x[[1]]$mapping)))
        backtrack[, (curr) := vals]
      }
    }
    
    # -- go through each adduct --
    do.calc <- function(x, smimap){
      
      #print("i")
      
      row <- adduct.table[x,]
      name <- row$Name
      mode <- row$Ion_mode
      # - - CHECK IF ADDUCT CAN EVEN BE FORMED BASED ON STRUCTURE - - - 
      rules_raw = row$Rule
      rules_split = strsplit(rules_raw, "AND| AND ",)[[1]]
      qualified_per_rule = as.data.table(sapply(rules_split, function(rule){
        middle = if(grepl("<", rule)) "below" else if(grepl(">", rule)) "above" else "equals"
        leftright = strsplit(rule, ">|<|=")[[1]]
        left_val = leftright[1]
        right = as.numeric(leftright[2])
        left = as.numeric(backtrack[,..left_val][[1]])
        qualifies = switch(middle,
                           below = left < right,
                           equals = left == right,
                           above = left > right)
      }))
      
      qualified = apply(qualified_per_rule, MARGIN = 1, all)
      
      #print("j")
      
      backtrack_filt <- backtrack[qualified,]
      
      unique_formulas <- unique(backtrack_filt[,c("baseformula", "charge")])
      
      if(nrow(backtrack_filt) == 0) return(NULL)
      
      # --- adduct before multiplication ---
      adduct_before <- row$AddAt
      deduct_before <- row$RemAt
      adduct_before <- if(is.na(adduct_before)) FALSE else adduct_before
      deduct_before <- if(is.na(deduct_before)) FALSE else deduct_before
      
      if(adduct_before != FALSE){
        formulae.add <- mergeform.joanna(formula1 = unique_formulas$baseformula,
                                         formula2 = adduct_before)
        unique_formulas$adducted <- formulae.add
      }else{
        unique_formulas$adducted <- unique_formulas$baseformula
      }
      #print("k")
      
      # --- placeholder ---
      unique_formulas$final <- unique_formulas$adducted
      # --- is deduction necessary? ---
      
      if(deduct_before != FALSE){
        formulae <- unique_formulas$final
        
        #print("k1")
        
        can.deduct <- which(!check.ded.joanna(formulas = unique_formulas$adducted,
                                              deduct = deduct_before))
        
        #print("k2")
        
        if(length(can.deduct) == 0) return(NA)
        deductibles <- formulae[can.deduct]
        
        #print("k3")
        
        formulae.ded <- subform.joanna(deductibles, deduct_before)
        unique_formulas$final[can.deduct] <- formulae.ded
        unique_formulas <- unique_formulas[can.deduct]
      }
      
      #print("l")
      
      # --- multiplication ---
      multiplier <- as.numeric(row$xM)
      
      if(multiplier > 1){
        formulae.mult <- multiform.joanna(unique_formulas$adducted, multiplier)
      }else{
        formulae.mult <- unique_formulas$adducted
      }
      
      unique_formulas$multiform <- formulae.mult
      unique_formulas$adducted <- unique_formulas$multiform
      
      #print("l1")
      
      # --- adduct after multiplication ---
      adduct_after <- row$AddEx
      deduct_after <- row$RemEx
      
      adduct_after <- if(is.na(adduct_after)) FALSE else adduct_after
      deduct_after <- if(is.na(deduct_after)) FALSE else deduct_after
      
      #print("l2")
      
      if(adduct_after != FALSE){
        formulae.add <- mergeform.joanna(formula1 = unique_formulas$adducted,
                                         formula2 = adduct_after)
        unique_formulas$adducted <- formulae.add
      }
      
      #print("l3")
      
      # --- placeholder ---
      unique_formulas$final <- unique_formulas$adducted
      
      # --- is deduction necessary? ---
      if(deduct_after != FALSE){
        formulae <- unique_formulas$final
        can.deduct <- which(!check.ded.joanna(formulas = unique_formulas$adducted,
                                              deduct = deduct_after))
        #print("l4")
        
        if(length(can.deduct) == 0) return(NA)
        deductibles <- formulae[can.deduct]
        formulae.ded <- subform.joanna(deductibles, deduct_after)
        
        #print("l5")
        
        unique_formulas$final[can.deduct] <- formulae.ded
        unique_formulas <- unique_formulas[can.deduct]
      }
      # - - - fixing... 
      if(nrow(backtrack_filt) == 0) return(NA)
      unique_formulas$final.charge <- c(as.numeric(unique_formulas$charge)) + c(as.numeric(row$Charge))
      unique_formulas <- unique_formulas[final.charge != 0]
      if(nrow(backtrack) == 0) return(NA)
      
      # --- get isotopes ---
      
      testformulas <<- unique_formulas
      
      # === RCDK METHOD ===
      # isotables <- lapply(1:nrow(unique_formulas), 
      #                     function(i){
      #                       isos = NULL
      #                       #try({
      #                       print("- - - - - ")
      #                         row <- unique_formulas[i,]
      #                         print(row)
      #                         form = rcdk::get.formula(row$final,charge = row$final.charge)
      #                         print(form)
      #                         if(form@mass %between% c(mzrange[1], mzrange[2])){
      #                           isos = rcdk::get.isotopes.pattern(form) 
      #                         }  
      #                       #})
      #                       # - - - 
      #                       isos
      #                     })
      # 
      # isolist <- lapply(isotables, function(isotable){
      #   if(is.null(isotable)) return(NA)
      #   iso.dt <- data.table::data.table(isotable, fill=TRUE)
      #   result <- iso.dt[,1:2]
      #   names(result) <- c("fullmz", "isoprevalence")
      #   # --- return ---
      #   result
      # })
      # 
      # names(isolist) <- unique_formulas$final
      
      # === ENVIPAT METHOD ===
      
      # check mass
      checked <- check.chemform.joanna(isotopes, unique_formulas$final)
      checked$mz <- checked$monoisotopic_mass/abs(unique_formulas$final.charge)
      keep <- which(checked$mz %between% mzrange)
      
      if(length(keep) == 0){
        return(NULL)
      }else{
        unique_formulas <- unique_formulas[keep,]
        
       #x print(unique_formulas)
         
        isotables <- enviPat::isopattern(
          isotopes,
          unique_formulas$final,
          threshold = 0.1,
          plotit = FALSE,
          charge = unique_formulas$final.charge,
          #algo = 2,
          verbose = FALSE,
        )
        isolist <- lapply(isotables, function(isotable){
          if(isotable[[1]] == "error"){
            return(NA)
          }
          iso.dt <- data.table::data.table(isotable, fill=TRUE)
          result <- iso.dt[,1:2]
          names(result) <- c("fullmz", "isoprevalence")
          # --- return ---
          result
        })
        
        #print("m")
        
        isolist.nonas <- isolist[!is.na(isolist)]
        isotable <- rbindlist(isolist.nonas)
        keep.isos <- names(isolist.nonas)
        
        charges <- unique_formulas$final.charge[which(!is.na(isolist))]
        
        repeat.times <- c(unlist(lapply(isolist.nonas, FUN=function(list) nrow(list))))
        isotable$final <- rep(keep.isos, repeat.times)
        isotable$final.charge <- rep(charges, repeat.times)
        
        #print("n")
        
        # --- remove 'backtrack' rows that couldn't be calculated ---
        unique_formulas <- unique_formulas[final %in% keep.isos]
        
        formula_plus_iso <- merge(unique_formulas, isotable, by = c("final", "final.charge"))
        
        backtrack.final <- merge(backtrack_filt, formula_plus_iso, 
                                 by=c("baseformula", "charge"), allow.cartesian = T)
        
        try({
          meta.table <- data.table::data.table(structure = backtrack.final$structure,
                                               fullformula = backtrack.final$final,
                                               fullmz = backtrack.final$fullmz,
                                               adduct = c(name),
                                               isoprevalence = backtrack.final$isoprevalence,
                                               foundinmode = c(mode),
                                               source = c(dbname))
          
          #print("o")
          
          # map SMILES to smile_id
          ids <- smimap$id[match(meta.table$structure, smimap$smiles)]
          meta.table$struct_id <- ids
          
          meta.table <- data.table::data.table(struct_id = meta.table$struct_id,
                                               fullformula = meta.table$fullformula,
                                               fullmz = meta.table$fullmz,
                                               adduct = meta.table$adduct,
                                               isoprevalence = meta.table$isoprevalence,
                                               foundinmode = meta.table$foundinmode,
                                               source = meta.table$source)
          # - - return - -
          
          #print("p")
          
          core.conn <- DBI::dbConnect(RSQLite::SQLite(), dbname = full.db)
          res <- DBI::dbSendQuery(core.conn, "PRAGMA busy_timeout=5000;")
          DBI::dbClearResult(res)
          
          #print("q")
          repeat{
            rv <- try({
              DBI::dbAppendTable(core.conn, "extended", meta.table)
            })
            if(!is(rv, "try-error")) break
            DBI::dbDisconnect(core.conn)
          }
        })
      }
    }
    
    #print("s")
    # if(length(cl)==1){
    #   if(cl == FALSE){
    #     cl = NULL
    #   }
    # }
    # this stores everyone in MEMORY! use insert statements...
    #if(is.null(cl)){
      sapply(1:nrow(adduct.table), FUN=function(x) do.calc(x, mapper))
    #}else{
    #  parallel::clusterExport(cl, "full.db")
    #  parallel::parSapply(cl=cl, 1:nrow(adduct.table), FUN=function(x) do.calc(x, mapper))
    #}
    
    #print("t")
    
    repeat{
      rv <- try({
        DBI::dbAppendTable(full.conn, "structures", partial.results)
      })
      if(!is(rv, "try-error")) break
      DBI::dbDisconnect(full.conn)
    }
    #print("u")
    
  }, full.db = full.db)
  
  # --- indexy ---
  #RSQLite::dbExecute(full.conn, "CREATE INDEX IF NOT EXISTS e_idx1 on extended(structure)") # slightly smaller index...
  RSQLite::dbExecute(full.conn, "CREATE INDEX IF NOT EXISTS e_idx2 on extended(fullmz, foundinmode)")
  #RSQLite::dbExecute(full.conn, "create index if not exists e_idx3 on extended(baseformula, adduct)") # new one for isotope backtracking??
  
  # --- cleanup ---
  #RSQLite::dbExecute(full.conn, "VACUUM")
  RSQLite::dbDisconnect(full.conn)
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
    new.qc.name <- paste0("QC", qc.i)
    new.qc.name <- gsub(colnames(poslist)[qc], pattern = "(^QC[\\d|\\d\\d])", replacement = new.qc.name,perl = T)
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
}


#' @export
load.metadata.csv <- function(path.to.csv,
                              path.to.patdb){
  
  #path.to.csv = "~/Downloads/maria_meta.csv"
  #path.to.patdb = "~/Downloads/maria_3ppm.db"
  
  
  # --- connect to sqlite db ---
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), path.to.patdb)
  
  csv <- data.table::fread(path.to.csv)
  colnames(csv) <- tolower(colnames(csv))
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
                            outfolder = getOptions(lcl$paths$opt.loc)$db_dir,
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
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), db)
  RSQLite::dbExecute(conn, statement = "create table base(compoundname text, description text, baseformula text, identifier text, charge text, structure text)")
  
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
