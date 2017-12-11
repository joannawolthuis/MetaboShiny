#' @export
get_ref_vars <- function(fac="Label"){
  req(csv_loc)
  # csv_loc = "backend/appdata/euronutrition/Test.csv"
  csv <- fread(csv_loc, sep="\t", header = T)[,1:5]
  # --- return ---
  c(unique(csv[,..fac]))[[fac]]
}

getProfile <- function(varName, title=varName, mode="stat"){
  # ---------------
  require(ggplot2)
  print(varName)
  print(colnames(mSet$dataSet$norm))
  # ---------------
  varInx <- colnames(mSet$dataSet$norm) == varName;
  var <- as.data.table(mSet$dataSet$norm, 
                       keep.rownames = T)[,varInx, with=FALSE];
  samp.names <- rownames(mSet$dataSet$norm)
  # ---------------
  if(mode == "time"){
    time.fac <<- mSet$dataSet$time.fac;
    translator <- data.table(
      index = 1:length(samp.names),
      Sample = gsub(x = samp.names, pattern = "_T\\d$", replacement=""),
      Group = mSet$dataSet$facB,
      Time = time.fac,
      Abundance = mSet$dataSet$norm[,varInx]
    )
  }else if(mode == "stat"){
    translator <- data.table(
      index = 1:length(samp.names),
      Sample = gsub(x = samp.names, pattern = "T\\d$", replacement=""),
      Group = mSet$dataSet$filt.cls,
      Abundance = mSet$dataSet$norm[,varInx]
    )
  }
  # ---------------
  return(translator)
}

kegg.charge <- function(atomlist){
  charges  <-regmatches(
    atomlist,
    regexpr(atomlist, pattern = "#[+-]|#\\d*[+-]",perl = T)
  ) 
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
