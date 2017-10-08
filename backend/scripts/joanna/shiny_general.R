#' @export
get_ref_vars <- function(fac="Label"){
  req(csv_loc)
  csv <- fread(csv_loc, sep="\t", header = T)
  # --- return ---
  unique(csv[,fac])
}

getProfile <-function (varName, title=varName, mode="stat") {
  # ---------------
  require(ggplot2)
  # ---------------
  varInx <- colnames(dataSet$norm) == varName;
  var <- as.data.table(dataSet$norm, keep.rownames = T)[,varInx, with=FALSE];
  samp.names <- rownames(dataSet$norm)
  exp.fac <<- dataSet$filt.cls
  # ---------------
  if(mainmode == "time"){
    time.fac <<- dataSet$time.fac;
    translator <- data.table(
      index = 1:length(samp.names),
      Sample = gsub(x = samp.names, pattern = "T\\d$", replacement=""),
      Group = exp.fac,
      Time = time.fac,
      Abundance = dataSet$norm[,varInx]
    )
  }else if(mainmode == "stat"){
    translator <- data.table(
      index = 1:length(samp.names),
      Sample = gsub(x = samp.names, pattern = "T\\d$", replacement=""),
      Group = exp.fac,
      Abundance = dataSet$norm[,varInx]
    )
  }
  # ---------------
  return(translator)
}

kegg.charge <- function(atomlist){
  atom.str <- paste(atomlist, collapse = " ") 
  charges  <- str_match(atom.str, pattern = "#[+-]|#\\d*[+-]")
  formal_charge = 0
  for(ch in charges[!is.na(charges)]){
    ch.base <- gsub(ch, pattern = "#", replacement = "")
    ch.base <- if(ch.base == "-" | ch.base == "+") paste0(ch.base, 1) else(ch.base)
    ch.base <- gsub(ch.base, pattern = "\\+", replacement = "")
    ch.base <- as.numeric(sub("([0-9.]+)-$", "-\\1", ch.base))
    # -------
    formal_charge = formal_charge + ch.base
  }
  formal_charge
}
