#' @export

get_ref_vars <- function(fac="Label"){
  req(csv_loc)
  csv <- fread(csv_loc, sep="\t", header = T)
  # --- return ---
  unique(csv[,fac])
}

