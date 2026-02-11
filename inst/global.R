# === GENERAL OPTIONS ===

options(stringsAsFactors = FALSE, "java.parameters" = c("-Xmx8G")) # give java enough memory for smiles parsing

library(MetaboShiny)
library(shinyBS)

if (Sys.getenv("SHINY_PORT") == "") options(shiny.maxRequestSize = 10000 * 1024^2)

# set the home path
home <- path.expand("~")

# TODO: add option to put user custom tables in user directory
adducts <- MetaboShiny:::load_adducts_from_admin(home)
adduct_rules <- MetaboShiny:::load_adduct_rules_from_admin(home)

MetaboShiny:::init_global_state(environment())
