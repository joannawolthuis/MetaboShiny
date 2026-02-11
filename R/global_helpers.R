load_adducts_from_admin <- function(home){
  admin_dir <- file.path(home, "MetaboShiny", "saves", "admin")
  if("adducts.csv" %in% list.files(admin_dir)){
    data.table::fread(file.path(admin_dir, "adducts.csv"))
  }else{
    data(adducts, package = "MetaDBparse")
    data.table::as.data.table(adducts)
  }
}

load_adduct_rules_from_admin <- function(home){
  admin_dir <- file.path(home, "MetaboShiny", "saves", "admin")
  if("adduct_rules.csv" %in% list.files(admin_dir)){
    data.table::fread(file.path(admin_dir, "adduct_rules.csv"))
  }else{
    data(adduct_rules, package = "MetaDBparse")
    data.table::as.data.table(adduct_rules)
  }
}
