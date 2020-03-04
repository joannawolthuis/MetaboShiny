# === GENERAL OPTIONS ===

options(stringsAsFactors = FALSE,"java.parameters" = c("-Xmx8G")) # give java enough memory for smiles parsing

if(Sys.getenv('SHINY_PORT') == "") options(shiny.maxRequestSize=10000*1024^2)

library(ggplot2)
library(data.table)
library(plotly)
library(shinyBS)
library(shinyjs)
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
                                           list(name = 'empty2', path = 'www/empty.png', dimensions = c("100%", 1)),
                                           list(name = 'empty3', path = 'www/empty.png', dimensions = c("100%", 1)),
                                           list(name = 'empty4', path = 'www/empty.png', dimensions = c("100%", 1)),
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
                                           list(name= 'drugbank_logo', path = 'www/drugbank_logo_2.png', dimensions = c(280,80)),
                                           list(name= 'cmmmediator_logo', path = 'www/ceummlogo.jpeg', dimensions = c(230,80)),
                                           list(name= 'pubchem_logo', path = 'www/pubchem_logo.png', dimensions = c(260,70)),
                                           list(name= 'chemspider_logo', path = 'www/chemspider_logo.png', dimensions = c(100,150)),
                                           list(name= 'phenolexplorer_logo', path = 'www/phenolexplorer_logo.jpg', dimensions = c(200,60)),
                                           list(name = 'sidebar_icon', path = 'www/detective.png', dimensions = c(60, 60)),
                                           list(name = 'merge_icon', path = 'www/merge.png', dimensions = c(150, 150)),
                                           list(name = 'db_icon', path = 'www/database.png', dimensions = c(150, 150)),
                                           list(name = 'ecmdb_logo', path = 'www/ecmdb_logo.png', dimensions = c(250, 80)),
                                           list(name = 'bmdb_logo', path = 'www/hmp_logo.png', dimensions = c(150, 150)),
                                           list(name = 'rmdb_logo', path = 'www/hmp_logo.png', dimensions = c(150, 150)),
                                           list(name = 'mvoc_logo', path = 'www/mvoc_logo.png', dimensions = c(230, 60)),
                                           list(name = 'mcdb_logo', path = 'www/mcdb_logo.png', dimensions = c(250, 80)),
                                           list(name = 'pamdb_logo', path = 'www/pamdb_logo.png', dimensions = c(150, 150)),
                                           list(name = 'ymdb_logo', path = 'www/ymdb_logo.png', dimensions = c(330, 90)),
                                           list(name = 'nanpdb_logo', path = 'www/nanpdb_logo.png', dimensions = c(200, 70)),
                                           list(name = 'lmdb_logo', path = 'www/lmdb_logo.png', dimensions = c(230, 100)),
                                           list(name = 'supernatural2_logo', path = 'www/supernatural2_logo.png', dimensions = c(300, 60)),
                                           list(name = 'stoff_logo', path = 'www/stoff_logo.png', dimensions = c(250, 130)),
                                           list(name = 'knapsack_logo', path = 'www/knapsack_logo.gif', dimensions = c(200, 100)),
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
                               lmdb = list(title = "LMDB",
                                           description = "The Livestock Metabolome Database (LMDB) is a freely available electronic database containing detailed information about small molecule metabolites found in different livestock species. It is intended to be used for applications in metabolomics, clinical chemistry, biomarker discovery and general education.",
                                           image_id = "lmdb_logo"),
                               ymdb = list(title = "YMDB",
                                           description = "The Yeast Metabolome Database (YMDB) is a manually curated database of small molecule metabolites found in or produced by Saccharomyces cerevisiae (also known as Baker’s yeast and Brewer’s yeast).",
                                           image_id = "ymdb_logo"),
                               ecmdb = list(title = "ECMDB",
                                            description = "The ECMDB is an expertly curated database containing extensive metabolomic data and metabolic pathway diagrams about Escherichia coli (strain K12, MG1655).",
                                            image_id = "ecmdb_logo"),
                               rmdb = list(title = "RMDB",
                                           description = "The Bovine Rumen Metabolome Database (RMDB) makes available tables containing the set of 246 ruminal fluid metabolites or metabolite species from the bovine ruminal fluid metabolome, along with their concentrations, related literature reference and links to their known diet associations.",
                                           image_id = "rmdb_logo"),
                               bmdb = list(title = "BMDB",
                                           description = "The Bovine Metabolome Database (BMDB) The Bovine Metabolome Database (BMDB) is a freely available electronic database containing detailed information about small molecule metabolites found in beef and dairy cattle.",
                                           image_id = "bmdb_logo"),
                               mcdb = list(title = "MCDB",
                                           description = "The Milk Composition Database (MCDB) is a freely available electronic database containing detailed information about small molecule metabolites found in cow milk.",
                                           image_id = "mcdb_logo"),
                               nanpdb = list(title = "NANPDB",
                                             description = "To the best of our knowledge this is the largest database of natural products isolated from native organisms of Northern Africa. In a recent survey of the African flora, including sources from Northern Africa, it was shown that this part of the world could be a huge repository of bioactive NPs with diverse scaffolds and activities. Consequently, many NPs are used in traditional Northern African medicine.",
                                             image_id = "nanpdb_logo"),
                               stoff = list(title = "STOFF-IDENT",
                                            description = "STOFF-IDENT is a database of water relevant substances collated from various sources within the STOFF-IDENT and FOR-IDENT projects, hosted by the Bavarian Environment Agency (Bayerisches Landesamt für Umwelt, LfU), the University of Applied Sciences Weihenstephan-Triesdorf (HSWT) and the Technical University of Munich (TUM).",
                                            image_id = "stoff_logo"),
                               pamdb = list(title = "PAMDB",
                                            description = "The PAMDB is an expertly curated database containing extensive metabolomic data and metabolic pathway diagrams about Pseudomonas aeruginosa (reference strain PAO1).",
                                            image_id = "pamdb_logo"),
                               mvoc = list(title = "mVOC",
                                           description = "The mVOC 2.0 Database is based on extensive literature search for microbial volatile organic compounds (mVOCs)",
                                           image_id = "mvoc_logo"),
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
                               supernatural2 = list(title = "SUPER NATURAL II",
                                                    description = "(VIA MAGICBALL) A database of natural products. It contains 325,508 natural compounds (NCs), including information about the corresponding 2d structures, physicochemical properties, predicted toxicity class and potential vendors.",
                                                    image_id = "supernatural2_logo"),
                               chemspider = list(title = "ChemSpider",
                                                 description = "(VIA MAGICBALL) A chemical structure database providing fast access to over 77 million structures, properties and associated information.",
                                                 image_id = "chemspider_logo"),
                               knapsack = list(title = "KNApSAcK",
                                               description = "(VIA MAGICBALL) The purpose of the KNApSAcK Metabolomics is to search metabolites from MS peak, molecular weight and molecular formula, and species.",
                                               image_id = "knapsack_logo"),
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
  volumes =  list('MetaboShiny' = getwd(),
                  'Home' = home,
                  '~' = normalizePath("~"),
                  'Documents' = file.path(home, "Documents"),
                  'Downloads' = file.path(home, "Downloads"),
                  'R Installation' = R.home(),
                  'Desktop' = file.path(home, "Desktop"))
),
# default vectors to go through in metaboshiny
vectors = list(
  hide_match_cols = c("structure", "identifier","baseformula",
                      "isocat", "fullformula", "finalcharge", "query_mz"),
  # list of available databases!!
  db_no_build = c("cmmmediator",
                  "chemspider",
                  "magicball",
                  "knapsack",
                  'supernatural2',
                  "custom",
                  "pubchem"),
  db_categories = list(versatile = c("wikidata", "dimedb", "metacyc", "chebi", "massbank", "cmmediator"),
                       verbose = c("hmdb", "chebi", "t3db", "metabolights", "ymdb", "ecmdb", "pamdb"),
                       livestock = c("lmdb", "rmdb", "bmdb", "metacyc", "mcdb"),
                       human = c("hmdb", "metacyc", "expoexplorer", "t3db", "bloodexposome"),
                       microbial = c("ymdb", "ecmdb", "pamdb", "vmh", "mvoc"),
                       pathway = c("vmh", "smpdb", "kegg"),
                       food = c("foodb", "phenolexplorer"),
                       plant = c("nanpdb", "respect", "metacyc", "supernatural2", "knapsack"),
                       chemical = c("chebi", "massbank", "maconda", "stoff", "lipidmaps", "chemspider"),
                       massspec = c("massbank", "respect", "maconda", "dimedb"),
                       online = c("cmmmediator", "chemspider","pubchem","knapsack"),
                       custom = c(),
                       predictive = c("magicball", "chemspider", "pubchem", "supernatural2", "knapsack")
  ),
  db_list = c( # this determines the show order of dbs in the app
    "hmdb",
    "chebi",
    "kegg",
    "metacyc",
    "mcdb",
    "pamdb",
    "mvoc",
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
    'lmdb',
    'ymdb',
    'ecmdb',
    'rmdb',
    'bmdb',
    'stoff',
    'nanpdb',
    'phenolexplorer',
    "magicball",
    'cmmmediator',
    'pubchem',
    'supernatural2',
    'chemspider',
    'knapsack',
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

gbl$vectors$db_categories$all <- gbl$vectors$db_list

gbl$vectors$wordcloud$filters <- list(
  stopwords = unique(c(tidytext::stop_words$word, 
                       qdapDictionaries::Top200Words,
                       tm::stopwords("english"))),
  metabolomics = c("metabolism", "metabolic",
                   "metabolomic", "metabolomics",
                   "biochemical", "mass", "spectrometry", 
                   "nmr", "direct", "infusion","exposome","papers",
                   "compounds","compound"),
  default = c(gbl$vectors$db_list, "exposome",
              "Synonyms", "synonyms"))

data(isotopes, package = "enviPat")

radioTooltip <- function(id, choice, title, placement = "bottom", trigger = "hover", options = NULL){
  options = shinyBS:::buildTooltipOrPopoverOptionsList(title, placement, trigger, options)
  options = paste0("{'", paste(names(options), options, sep = "': '", collapse = "', '"), "'}")
  bsTag <- shiny::tags$script(shiny::HTML(paste0("
    $(document).ready(function() {
      setTimeout(function() {
        $('input', $('#", id, "')).each(function(){
          if(this.getAttribute('value') == '", choice, "') {
            opts = $.extend(", options, ", {html: true});
            $(this.parentElement).tooltip('destroy');
            $(this.parentElement).tooltip(opts);
          }
        })
      }, 500)
    });
  ")))
  htmltools::attachDependencies(bsTag, shinyBS:::shinyBSDep)
}

# interleave for sorting later ...
add_idx <- order(c(seq_along(gbl$vectors$pos_adducts$Name), seq_along(gbl$vectors$neg_adducts$Name)))
sort_order <<- unlist(c(gbl$vectors$pos_adducts$Name, gbl$vectors$neg_adducts$Name))[add_idx]

session_cl <- NULL
debug_mSet <- NULL
debug_lcl <- NULL
debug_input <- NULL

try({
  success=F
  orca_server <- plotly::orca_serve()
  success=T
},silent=T)
if(!success) print("Orca isn't working, please check your installation. If on Mac, please try starting Rstudio from the command line with the command 'open -a Rstudio'")

msg.vec <- c()
