# === GENERAL OPTIONS ===

options(stringsAsFactors = FALSE, "java.parameters" = c("-Xmx8G")) # give java enough memory for smiles parsing

library(MetaboShiny)
library(shinyBS)

if (Sys.getenv("SHINY_PORT") == "") options(shiny.maxRequestSize = 10000 * 1024^2)

# set the home path
home <- path.expand("~")

# TODO: add option to put user custom tables in user directory
if (!exists("load_adducts_from_admin", envir = asNamespace("MetaboShiny"), inherits = FALSE)) {
  stop(
    "Installed package 'MetaboShiny' does not contain `load_adducts_from_admin()`.\n",
    "This usually means your installed package is older than this repo.\n",
    "Reinstall/update MetaboShiny from this source tree, then restart R and try again."
  )
}
if (!exists("load_adduct_rules_from_admin", envir = asNamespace("MetaboShiny"), inherits = FALSE)) {
  stop(
    "Installed package 'MetaboShiny' does not contain `load_adduct_rules_from_admin()`.\n",
    "This usually means your installed package is older than this repo.\n",
    "Reinstall/update MetaboShiny from this source tree, then restart R and try again."
  )
}
adducts <- MetaboShiny:::load_adducts_from_admin(home)
adduct_rules <- MetaboShiny:::load_adduct_rules_from_admin(home)

MetaboShiny:::init_global_state(environment())
