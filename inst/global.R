# === GENERAL OPTIONS ===

options(stringsAsFactors = FALSE,"java.parameters" = c("-Xmx8G")) # give java enough memory for smiles parsing
if(Sys.getenv('SHINY_PORT') == "") options(shiny.maxRequestSize=10000*1024^2)

library(ggplot2)
library(data.table)
library(plotly)
library(MetaboShiny)

# set the home path
home = normalizePath("~")

# load default adduct table
#TODO: add option to put user custom tables in user directory
data(adducts, package = "MetaDBparse")
adducts <- data.table::as.data.table(adducts)
data(adduct_rules, package = "MetaDBparse")
adduct_rules <- data.table::as.data.table(adduct_rules)

caret.mdls <- caret::getModelInfo()

# === THE BELOW LIST CONTAINS ALL GLOBAL VARIABLES THAT METABOSHINY CALLS UPON LATER ===
gbl <- list(constants = list(ppm = 2, # TODO: re-add ppm as option for people importing their data through csv
                             ml.twoonly = c("adaboost","logicBag","bartMachine","binda",
                                            "ada","gamboost","glmboost","chaid",
                                            "C5.0Cost","rpartCost","deepboost",
                                            "dwdPoly","dwdRadial","glm","glmnet",
                                            "glmStepAIC","glmnet_h2o","svmLinearWeights2",
                                            "dwdLinear","svmLinearWeights","logreg","mlpKerasDropoutCost",
                                            "mlpKerasDecayCost","ORFlog","ORFpls","ORFridge","ORFsvm",
                                            "plsRglm","rotationForest","rotationForestCp",
                                            "svmRadialWeights","nodeHarvest"),
                             # get all caret models that can do classification and have some kind of importance metric
                             ml.models = names(caret.mdls)[sapply(1:length(caret.mdls), function(i){
                               curr.mdl = caret.mdls[[i]]
                               can.classify = if("Classification" %in% curr.mdl$type) TRUE else FALSE
                               has.importance = if("varImp" %in% names(curr.mdl)) TRUE else FALSE
                               can.classify & has.importance
                             })],
                             max.cols = 20,
                             images = list(list(name = 'load_icon', path = 'www/cute.png', dimensions = c(100, 100)),
                                           list(name = 'empty', path = 'www/empty.png', dimensions = c("100%", 1)),
                                           list(name = 'cute_package', path = 'www/cat.png', dimensions = c(80, 80)),
                                           list(name = 'internal_logo', path = 'www/umcinternal.png', dimensions = c(120, 120)),
                                           list(name = 'login_header', path = 'www/login_icon.png', dimensions = c(300,200)),
                                           list(name = 'noise_logo', path = 'www/umcnoise.png', dimensions = c(120, 120)),
                                           list(name = 'hmdb_logo', path = 'www/hmdblogo.png', dimensions = c(150, 100)),
                                           list(name = 'metacyc_logo', path = 'www/metacyc.png', dimensions = c(300, 80)),
                                           list(name = 'chebi_logo', path = 'www/chebilogo.png', dimensions = c(140, 140)),
                                           list(name = 'wikipathways_logo', path = 'www/wikipathways.png', dimensions = c(130, 150)),
                                           list(name = 'kegg_logo', path = 'www/kegglogo.gif', dimensions = c(200, 150)),
                                           list(name = 'smpdb_logo', path = 'www/smpdb_logo_adj.png', dimensions = c(200, 160)),
                                           list(name = 'dimedb_logo', path = 'www/dimedb_logo.png', dimensions = c(310, 120)),
                                           list(name = 'wikidata_logo', path = 'www/wikidata.png', dimensions = c(250, 200)),
                                           list(name = 'respect_logo', path = 'www/respect_logo.png', dimensions = c(250, 100)),
                                           list(name = 'massbank_logo', path = 'www/massbank_logo.jpg', dimensions = c(250, 100)),
                                           list(name = 'metabolights_logo', path = 'www/metabolights_logo.png', dimensions = c(200, 200)),
                                           list(name = 'vmh_logo', path = 'www/vmh.jpg', dimensions = c(200, 200)),
                                           list(name = 'foodb_logo', path = 'www/foodb_logo.png', dimensions = c(250, 90)),
                                           list(name = 'pos_icon', path = 'www/handpos.png', dimensions = c(120, 120)),
                                           list(name = 'neg_icon', path = 'www/handneg.png', dimensions = c(120, 120)),
                                           list(name = 'excel_icon', path = 'www/excel.png', dimensions = c(120, 120)),
                                           list(name = 'excel_icon_2', path = 'www/excel.png', dimensions = c(120, 120)),
                                           list(name = 'db_icon', path = 'www/servers.png', dimensions = c(150, 150)),
                                           list(name = 'lipidmaps_logo', path = 'www/lipidmaps.png', dimensions = c(200, 150)),
                                           list(name = 'bloodexposome_logo', path = 'www/bloodexposome.png', dimensions = c(250, 200)),
                                           list(name = 'csv_icon', path = 'www/office.png', dimensions = c(100, 100)),
                                           list(name= 'magicball', path = 'www/magic-ball2.png', dimensions = c(200,200)),
                                           list(name = 'dataset_icon', path = 'www/office.png', dimensions = c(100, 100)),
                                           list(name = 'plus', path = 'www/add.png', dimensions = c(150, 150)),
                                           list(name= 'maconda_logo', path = 'www/maconda.png', dimensions = c(250,100)),
                                           list(name= 'expoexplorer_logo', path = 'www/exposome_explorer.png', dimensions = c(250,100)),
                                           list(name= 't3db_logo', path = 'www/t3db_logo.png', dimensions = c(200,80)),
                                           list(name= 'drugbank_logo', path = 'www/drugbank_logo_2.png', dimensions = c(230,80)),
                                           list(name= 'cmmmediator_logo', path = 'www/ceummlogo.jpeg', dimensions = c(230,80)),
                                           list(name= 'pubchem_logo', path = 'www/pubchem_logo.png', dimensions = c(260,70)),
                                           list(name= 'chemspider_logo', path = 'www/chemspider_logo.png', dimensions = c(130,150)),
                                           list(name= 'phenolexplorer_logo', path = 'www/phenolexplorer_logo.jpg', dimensions = c(200,60)),
                                           list(name = 'sidebar_icon', path = 'www/detective.png', dimensions = c(60, 60)),
                                           list(name = 'merge_icon', path = 'www/merge.png', dimensions = c(150, 150)),
                                           list(name = 'db_icon', path = 'www/database.png', dimensions = c(150, 150)),
                                           list(name = 'laptop_icon', path = 'www/laptop.png', dimensions = c(150, 150))
                                           
                             ),# all image paths, if you add an image you can add it here
                             default.text = list(list(name='curr_definition', text="No m/z selected"),
                                                 list(name="curr_cpd", text="..."),# default text options at startup
                                                 list(name="ml_train_ss", text="all"),
                                                 list(name="ml_test_ss", text="all")
                             ),
                             db.build.info = list(
                               hmdb = list(title = "HMDB",
                                           description = "Metabolites commonly found in human biological samples.",
                                           image_id = "hmdb_logo"),
                               metacyc = list(title = "MetaCyc",
                                              description = "Large pathway database with over 10 000 available compounds. Spans several organisms!",
                                              image_id = "metacyc_logo"),
                               chebi = list(title = "ChEBI",
                                            description = "A broad database with known chemicals of biological interest.",
                                            image_id = "chebi_logo"),
                               wikipathways = list(title = "WikiPathways",
                                                   description = "Open source biological pathway database. Currently only partially available. Requires CHEBI to be built.",
                                                   image_id = "wikipathways_logo"),
                               kegg = list(title = "KEGG",
                                           description = "Large pathway database with info on pathways in various organisms, involved enzymes, and connected disease phenotypes.",
                                           image_id = "kegg_logo"),
                               smpdb = list(title = "SMPDB",
                                            description = "Small molecule pathway database. Compounds overlap with HMDB.",
                                            image_id = "smpdb_logo"),
                               dimedb = list(title = "DimeDB",
                                             description = "A direct infusion database of biologically relevant metabolite structures and annotations.",
                                             image_id = "dimedb_logo"),
                               wikidata = list(title = "Wikidata",
                                               description = "Central storage for the data of its Wikimedia sister projects including Wikipedia, Wikivoyage, Wikisource, and others.",
                                               image_id = "wikidata_logo"),
                               respect = list(title = "ReSpect",
                                              description = "RIKEN MSn spectral database for phytochemicals (ReSpect) is a collection of literature and in-house MSn spectra data for research on plant metabolomics.",
                                              image_id = "respect_logo"),
                               massbank = list(title = "MassBank",
                                               description = "This site presents the database of comprehensive, high-resolution mass spectra of metabolites. Supported by the JST-BIRD project, it offers various query methods for standard spectra from Keio Univ., RIKEN PSC, and others.
                                                    In 2008, MassBank was authorized as the official mass spectral database of The Mass Spectrometry Society of Japan.",
                                               image_id = "massbank_logo"),
                               metabolights = list(title = "MetaboLights",
                                                   description = "MetaboLights is a database for Metabolomics experiments and derived information. The database is cross-species, cross-technique and covers metabolite structures and their reference spectra as well as their biological roles, locations and concentrations, and experimental data from metabolic experiments.",
                                                   image_id = "metabolights_logo"),
                               vmh = list(title = "VMH",
                                          description = "Virtual Metabolic Human (VMH) hosts ReconMap, an extensive network of human metabolism, and bacterial metabolites.",
                                          image_id = "vmh_logo"),
                               foodb = list(title = "FooDB",
                                            description = "FooDB is the world’s largest and most comprehensive resource on food constituents, chemistry and biology. It provides information on both macronutrients and micronutrients, including many of the constituents that give foods their flavor, color, taste, texture and aroma.",
                                            image_id = "foodb_logo"),
                               bloodexposome = list(title = "Blood Exposome DB",
                                                    description = " This new blood exposome database can be applied to prioritize literature-based chemical reviews, developing target assays in exposome research, identifying compounds in untargeted mass spectrometry and biological interpretation in metabolomics data.",
                                                    image_id = "bloodexposome_logo"),
                               maconda = list(title = "MaConDa",
                                              description = "MaConDa currently contains ca. 200 contaminant records detected across several MS platforms. The majority of records include theoretical as well as experimental MS data.",
                                              image_id = "maconda_logo"),
                               lipidmaps = list(title = "LIPID MAPS",
                                                description = "The LIPID MAPS Structure Database (LMSD) is a relational database encompassing structures and annotations of biologically relevant lipids.",
                                                image_id = "lipidmaps_logo"),
                               expoexplorer = list(title = "Exposome-Explorer",
                                                   description = "Exposome-Explorer is the first database dedicated to biomarkers of exposure to environmental risk factors for diseases.",
                                                   image_id = "expoexplorer_logo"),
                               t3db = list(title = "T3DB",
                                           description = "The Toxin and Toxin Target Database (T3DB), or, soon to be referred as, the Toxic Exposome Database, is a unique bioinformatics resource that combines detailed toxin data with comprehensive toxin target information.",
                                           image_id = "t3db_logo"),
                               drugbank = list(title = "DrugBank",
                                               description = "The DrugBank database is a comprehensive, freely accessible, online database containing information on drugs and drug targets.",
                                               image_id = "drugbank_logo"),
                               phenolexplorer = list(title = "Phenol-Explorer",
                                                     description = "Phenol-Explorer is the first comprehensive database on polyphenol content in foods. The database contains more than 35,000 content values for 500 different polyphenols in over 400 foods.",
                                                     image_id = "phenolexplorer_logo"),
                               # - - leave magicball last - -
                               cmmmediator = list(title = "CEU Mass Mediator",
                                              description = "(ONLINE ONLY) CEU Mass Mediator is a tool for searching metabolites in different databases (Kegg, HMDB, LipidMaps, Metlin, MINE and an in-house library).",
                                              image_id = "cmmmediator_logo"),
                               magicball = list(title = "MagicBall",
                                                description = "Algorithm to predict molecular formula from m/z value",
                                                image_id = "magicball"),
                               pubchem = list(title = "PubChem",
                                                description = "(VIA MAGICBALL) PubChem is the world's largest collection of freely accessible chemical information.",
                                                image_id = "pubchem_logo"),
                               chemspider = list(title = "ChemSpider",
                                                description = "(VIA MAGICBALL) 	A chemical structure database providing fast access to over 77 million structures, properties and associated information.",
                                                image_id = "chemspider_logo"),
                               custom = list(title = "Custom",
                                             description = "Please select to start your own database creation!",
                                             image_id = "plus")
                             )
),
functions = list(# default color functions at startup, will be re-loaded from options
  cf = rainbow,
  color.function = rainbow,
  color.vec = rainbow,
  # available plot themes for ggplot. Can add more,also user-defined ones,
  # but put them in shiny_general.R first so they are sourced properly.
  plot.themes = list(bw=ggplot2::theme_bw,
                     classic=ggplot2::theme_classic,
                     gray=ggplot2::theme_gray,
                     min=ggplot2::theme_minimal,
                     dark=ggplot2::theme_dark,
                     light=ggplot2::theme_light,
                     line=ggplot2::theme_linedraw),
  color.functions = {
    # available colorbrewer themes to load into ggplot. These are the standard brew names used in their functions color.brewer etc.
    brew.cols <- c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", # - - sequential - -
                   "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples",
                   "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd",
                   "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", # - - diverging - -
                   "Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3" # - - qualitative - -
    )
    
    # generate direct functions from the brewer colours
    brew.opts <- lapply(brew.cols, function(opt) grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[opt,]$maxcolors, opt))))
    names(brew.opts) <- brew.cols
    
    # add the more general color scale functions as options
    base.opts <- list("rb"=rainbow,
                      "y2b"=colorRamps::ygobb,
                      "ml1"=colorRamps::matlab.like2,
                      "ml2"=colorRamps::matlab.like,
                      "m2g"=colorRamps::magenta2green,
                      "c2y"=colorRamps::cyan2yellow,
                      "b2y"=colorRamps::blue2yellow,
                      "g2r"=colorRamps::green2red,
                      "b2g"=colorRamps::blue2green,
                      "b2r"=colorRamps::blue2red,
                      "b2p"=grDevices::cm.colors,
                      "bgy"=grDevices::topo.colors,
                      "gyw"=grDevices::terrain.colors,
                      "ryw"=grDevices::heat.colors,
                      "bw"=MetaboShiny::blackwhite.colors)
    
    # add into a single list for use in interface
    append(base.opts, brew.opts)
  }),
# set default paths
paths = list(
  # available paths when selecting a new file or folder
  volumes =  c('MetaboShiny' = getwd(),
               'Home'=home,
               '~' = normalizePath("~"),
               'Documents'=file.path(home, "Documents"),
               'Downloads'=file.path(home, "Downloads"),
               'R Installation'=R.home(),
               'Desktop'=file.path(home, "Desktop"))
),
# default vectors to go through in metaboshiny
vectors = list(
  hide_match_cols = c("structure", "baseformula","isocat", "fullformula", "finalcharge", "query_mz"),#, "source"), #"query_mz", "identifier"), #c("description","structure", "baseformula", "dppm", "source"),
  # default indices of chosen adducts
  pos_selected_add = c(2),
  neg_selected_add = c(2),
  # list of available databases!!
  db_no_build = c("cmmmediator",
                  "chemspider",
                  "magicball",
                  "custom",
                  "pubchem"),
  db_list = c( # this determines the show order of dbs in the app
    "hmdb",
    "chebi",
    "kegg",
    "metacyc",
    #"dsstox",
    #"wikipathways",
    "smpdb",
    "dimedb",
    "wikidata",
    "vmh",
    "respect",
    "massbank",
    "metabolights",
    "foodb",
    "maconda",
    "bloodexposome",
    "expoexplorer",
    "lipidmaps",
    't3db',
    'drugbank',
    'phenolexplorer',
    "magicball",
    'cmmmediator',
    'pubchem',
    'chemspider',
    "custom"
  ),
  # list of positive adducts
  pos_adducts = adducts[Ion_mode == "positive",
                        c("Name")],
  # list of negative adducts
  neg_adducts = adducts[Ion_mode == "negative",
                        c("Name")],
  wordcloud = list(top = 20),
  calc_adducts = c("M+H", "M-H")
)
)

gbl$vectors$wordcloud$skip <- unique(c( # manual curation(
  "on", "in", "and", "at", "sglt", "association",
  "an", "by", "is", "it", "that",
  "as", "be", "like", "can", "a", "of",
  "to", "but", "not", "mainly", "the",
  "", "which", "from", "changes", "weeks",
  "found", "its", "two", "one", "if", "no",
  "yes", "any", "were", "observed", "also",
  "why", "other", "only", "known", "so", "do",
  "with", "resulting", "reaction", "via",
  "mtblc", "group", "groups",
  "metabolism", "mh", "ms", "position", "positions",
  "produced", "this", "ce", "mtblc",
  "has", "ko", "predicted", "are", "been", "isolated",
  "occurs", "form", "generated", "obtained", "et", "insource",
  "or", "nu", "substituted", "exhibits", "for", "ev",
  "attached", "constituent", "negative", "ph", "pmid",
  "compound", "residue", "unknown", "residues", "through",
  "lcesiitft", "lcesitof","lcesiqtof","have", "derived",
  "compounds", "having", "lcesiqq", "different", "more",
  "ec", "activity", "metabolite", "biotransformer","biotransformer¹",
  # pubmed words https://www.ncbi.nlm.nih.gov/books/NBK3827/table/pubmedhelp.T.stopwords/
  c("a", "about", "again", "all", "almost", "also", "although",
    "always", "among", "an", "and", "another", "any", "are", "as",
    "at", "be", "because", "been", "before", "being", "between",
    "both", "but", "by", "can", "could", "did", "do", "does", "done",
    "due", "during", "each", "either", "enough", "especially", "etc",
    "for", "found", "from", "further", "had", "has", "have", "having",
    "here", "how", "however", "i", "if", "in", "into", "is", "it",
    "its", "itself", "just", "kg", "km", "made", "mainly", "make",
    "may", "mg", "might", "ml", "mm", "most", "mostly", "must", "nearly",
    "neither", "no", "nor", "obtained", "of", "often", "on", "our",
    "overall", "perhaps", "pmid", "quite", "rather", "really", "regarding",
    "seem", "seen", "several", "should", "show", "showed", "shown",
    "shows", "significantly", "since", "so", "some", "such", "than",
    "that", "the", "their", "theirs", "them", "then", "there", "therefore",
    "these", "they", "this", "those", "through", "thus", "to", "upon",
    "use", "used", "using", "various", "very", "was", "we", "were",
    "what", "when", "which", "while", "with", "within", "without",
    "would", LETTERS, letters, "acid", "cell", "cells", "human", "humans",
    "practically", "containing", "belongs", "class", "chemspider", "considered",
    "primarily", "pathway", "novo", "tgi", "tgii", "acids", "molecule", "enzyme",
    "available", "description", "neutral", "logp", "dgi", "located", "relatively",
    "family", "least", "common", "four", "species", "skeleton", "total", "contain", "arising",
    "substituents", "bond", "oxo", "alpha", "formal", "brbr", "effects", "nomenclature", "metacyccp",
    "additional", "range", "events", "principle", "involved", "product", "catalyzed", 
    "qqq", "qtof", "nmr", "eib", "cid", "exp", "fragv", "lcesiqft", "mhho", "major", 
    "constit", "lpls", "formed", "lpes", "main", "subclass",
    "analysis", "approach", "area", "assessment", "assume", "authority", 
    "available", "benefit", "concept", "consistent", "constitutional", 
    "context", "contract", "create", "data", "definition", "derived", 
    "distribution", "economic", "environment", "established", "estimate", 
    "evidence", "export", "factors", "financial", "formula", "function", 
    "identified", "income", "indicate", "individual", "interpretation", 
    "involved", "issues", "labour", "legal", "legislation", "major", 
    "method", "occur", "percent", "period", "policy", "principle", 
    "procedure", "process", "required", "research", "response", "role", 
    "section", "sector", "significant", "similar", "source", "specific", 
    "structure", "theory", "variables", "achieve", "acquisition", 
    "administration", "affect", "appropriate", "aspects", "assistance", 
    "categories", "chapter", "commission", "community", "complex", 
    "computer", "conclusion", "conduct", "consequences", "construction", 
    "consumer", "credit", "cultural", "design", "distinction", "elements", 
    "equation", "evaluation", "features", "final", "focus", "impact", 
    "injury", "institute", "investment", "items", "journal", "maintenance", 
    "normal", "obtained", "participation", "perceived", "positive", 
    "potential", "previous", "primary", "purchase", "range", "region", 
    "regulations", "relevant", "resident", "resources", "restricted", 
    "security", "sought", "select", "site", "strategies", "survey", 
    "text", "traditional", "transfer", "alternative", "circumstances", 
    "comments", "compensation", "components", "consent", "considerable", 
    "constant", "constraints", "contribution", "convention", "coordination", 
    "core", "corporate", "corresponding", "criteria", "deduction", 
    "demonstrate", "document", "dominant", "emphasis", "ensure", 
    "excluded", "framework", "funds", "illustrated", "immigration", 
    "implies", "initial", "instance", "interaction", "justification", 
    "layer", "link", "location", "maximum\t", "minorities", "negative", 
    "outcomes", "partnership", "philosophy", "physical", "proportion", 
    "published", "reaction", "registered", "reliance", "removed", 
    "scheme", "sequence", "sex", "shift", "specified", "sufficient", 
    "task", "technical", "techniques", "technology", "validity", 
    "volume", "access", "adequate", "annual", "apparent", "approximated", 
    "attitudes", "attributed", "civil", "code", "commitment", "communication", 
    "concentration", "conference", "contrast", "cycle", "debate", 
    "despite", "dimensions", "domestic", "emerged", "error", "ethnic", 
    "goals", "granted", "hence", "hypothesis", "implementation", 
    "implications", "imposed", "integration", "internal", "investigation", 
    "job", "label", "mechanism", "obvious", "occupational", "option", 
    "output", "overall", "parallel", "parameters", "phase", "predicted", 
    "principal", "prior", "professional", "project", "promote", "regime", 
    "resolution", "retained", "series", "statistics", "status", "stress", 
    "subsequent", "sum", "summary", "undertaken", "academic", "adjustment", 
    "alter", "amendment", "aware", "capacity", "challenge", "clause", 
    "compounds", "conflict", "consultation", "contact", "decline", 
    "discretion", "draft", "enable", "energy", "enforcement", "entities", 
    "equivalent", "evolution", "expansion", "exposure", "external", 
    "facilitate", "fundamental", "generated", "generation", "image", 
    "liberal", "licence", "logic", "marginal", "medical", "mental", 
    "modified", "monitoring", "network", "notion", "objective", "orientation", 
    "perspective", "precise", "prime", "psychology", "pursue", "ratio", 
    "rejected", "revenue", "stability", "styles", "substitution", 
    "sustainable", "symbolic", "target", "transition", "trend", "version", 
    "welfare", "whereas", "abstract", "accurate", "acknowledged", 
    "aggregate", "allocation", "assigned", "attached", "author", 
    "bond", "brief", "capable", "cited", "cooperative", "discrimination", 
    "display", "diversity", "domain", "edition", "enhanced", "estate", 
    "exceed", "expert", "explicit", "federal", "fees", "flexibility", 
    "furthermore", "gender", "ignored", "incentive", "incidence", 
    "incorporated", "index", "inhibition", "initiatives", "input", 
    "instructions", "intelligence", "interval", "lecture", "migration", 
    "minimum", "ministry", "motivation", "neutral", "nevertheless", 
    "overseas", "preceding", "presumption", "rational", "recovery", 
    "revealed", "scope", "subsidiary", "tapes", "trace", "transformation", 
    "transport", "underlying", "utility", "adaptation", "adults", 
    "advocate", "aid", "channel", "chemical", "classical", "comprehensive", 
    "comprise", "confirmed", "contrary", "converted", "couple", "decades", 
    "definite", "deny", "differentiation", "disposal", "dynamic", 
    "eliminate", "empirical", "equipment", "extract", "file", "finite", 
    "foundation", "global", "grade", "guarantee", "hierarchical", 
    "identical", "ideology", "inferred", "innovation", "insert", 
    "intervention", "isolated", "media", "mode", "paradigm", "phenomenon", 
    "priority", "prohibited", "publication", "quotation", "release", 
    "reverse", "simulation", "solely", "somewhat", "submitted", "successive", 
    "survive", "thesis", "topic", "transmission", "ultimately", "unique", 
    "visible", "voluntary", "abandon", "accompanied", "accumulation", 
    "ambiguous", "appendix", "appreciation", "arbitrary", "automatically", 
    "bias", "chart", "clarity", "conformity", "commodity", "complement", 
    "contemporary", "contradiction", "crucial", "currency", "denote", 
    "detected", "deviation", "displacement", "dramatic", "eventually", 
    "exhibit", "exploitation", "fluctuations", "guidelines", "highlighted", 
    "implicit", "induced", "inevitably", "infrastructure", "inspection", 
    "intensity", "manipulation", "minimised", "nuclear", "offset", 
    "paragraph", "plus", "practitioners", "predominantly", "prospect", 
    "radical", "random", "reinforced", "restore", "revision", "schedule", 
    "tension", "termination", "theme", "thereby", "uniform", "vehicle", 
    "via", "virtually", "widespread", "visual", "accommodation", 
    "analogous", "anticipated", "assurance", "attained", "behalf", 
    "bulk", "ceases", "coherence", "coincide", "commenced", "incompatible", 
    "concurrent", "confined", "controversy", "conversely", "device", 
    "devoted", "diminished", "distorted", "distortion", "duration", 
    "erosion", "ethical", "format", "founded", "inherent", "insights", 
    "integral", "intermediate", "manual", "mature", "mediation", 
    "medium", "military", "minimal", "mutual", "norms", "overlap", 
    "passive", "portion", "preliminary", "protocol", "qualitative", 
    "refine", "relaxed", "restraints", "revolution", "rigid", "route", 
    "scenario", "sphere", "subordinate", "supplementary", "suspended", 
    "team", "temporary", "trigger", "unified", "violation", "vision", 
    "adjacent", "albeit", "assembly", "collapse", "colleagues", "compiled", 
    "conceived", "convinced", "depression", "encountered", "enormous", 
    "forthcoming", "inclination", "integrity", "intrinsic", "invoked", 
    "levy", "likewise", "nonetheless", "notwithstanding", "odd", 
    "ongoing", "panel", "persistent", "posed", "reluctant", "socalled", 
    "straightforward", "undergo", "whereby","makes", "organic", 
    "moiety", "based", "pca", "ring", "acidic", "lmpk", "phosphate", 
    "disease", "atom", "cyclic", "body", "action", "enzymes", "enzyme", "methoxy", 
    "exists", "solid", "products", "amino", "derivatives", "backbone", 
    "units", "slightly", "member", "moderately", "orange", "aromatic", 
    "fused", "pos", "neg", "methyl", "protein", "linked", "agents", 
    "include", "classification", "replaced", "pka", "outside", "derivative",
    "weakly", "lcesiq", "hydroxyl", "structures", "structure", "esa", 
    "respectively", "consists", "metabolic", "atoms", "dimethyl", 
    "ether", "multiple", "examples", "example", "dm", "level", 
    "clinical", "higher", "control", "increased", "associated", 
    "mode", "decreased", "studies", "study", "patient", "patients", 
    "years", "type", "gene", "compared", "including", "lower", 
    "ci", "concentration", "concentrations", "reduced", "increase", 
    "results", "result", "levels", "It", "lt", "model", "effect", 
    "high", "low", "expression", "test", "production", "aim",
    "tested","studied","ng", "se"),
  qdapDictionaries::Top200Words,
  tm::stopwords("english"),
  gbl$vectors$db_list
))

#' Squishes HTML elements close together.
data(isotopes, package = "enviPat")

# interleave for sorting later ...
add_idx <- order(c(seq_along(gbl$vectors$pos_adducts$Name), seq_along(gbl$vectors$neg_adducts$Name)))
sort_order <<- unlist(c(gbl$vectors$pos_adducts$Name, gbl$vectors$neg_adducts$Name))[add_idx]

session_cl <- NULL
debug_mSet <- NULL
debug_lcl <- NULL
debug_input <- NULL

try({
  orca_server <- plotly::orca_serve()
})

msg.vec <- c()
