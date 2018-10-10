#' @export
get_ref_vars <- function(fac="Label"){
  req(csv_loc)
  # csv_loc = "backend/appdata/euronutrition/Test.csv"
  csv <- fread(csv_loc, sep="\t", header = T)[,1:5]
  # --- return ---
  c(unique(csv[,..fac]))[[fac]]
}

get_ref_cpds <- function(){
  req(csv_loc)
  # csv_loc = "backend/appdata/euronutrition/Test.csv"
  csv <- fread(csv_loc, sep="\t", header = T)[1,]
  keep.colnames <- colnames(csv)[!colnames(csv) %in% c("Sample", "Label", "Time")]
  # --- return ---
  c(keep.colnames)
}

getProfile <- function(varName, title=varName, sourceTable = mSet$dataSet$norm, mode="stat"){
  # ---------------
  varInx <- colnames(sourceTable) == varName;
  var <- as.data.table(sourceTable, 
                       keep.rownames = T)[,varInx, with=FALSE];
  samp.names <- rownames(sourceTable)
  # ---------------
  if(mode == "time"){
    time.fac <<- mSet$dataSet$time.fac;
    translator <- data.table(
      index = 1:length(samp.names),
      Sample = gsub(x = samp.names, pattern = "_T\\d$", replacement=""),
      Group = mSet$dataSet$facB,
      Time = time.fac,
      Abundance = sourceTable[,varInx]
    )
  }else if(mode == "stat"){
    translator <- data.table(
      index = 1:length(samp.names),
      Sample = gsub(x = samp.names, pattern = "T\\d$", replacement=""),
      Group = mSet$dataSet$cls,
      Abundance = sourceTable[,varInx]
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

mape <- function(actual,pred){
  mape <- mean(abs((actual - pred)/actual))*100
  return (mape)
}

flattenlist <- function(x){  
  morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
  out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE, use.names = T))
  if(sum(morelists)){ 
    Recall(out)
  }else{
    return(out)
  }
}

get.col.map <- function(optionfile){
  options <- getOptions(optionfile)
  unparsed.cols <- options$gcols
  col.items <- strsplit(unparsed.cols, split = "&")
  keys <- gsub(col.items[[1]], pattern = ":.*$", replacement="")
  items <- gsub(col.items[[1]], pattern = "^.*:", replacement="")
  names(items) <- keys
  # - - - -
  items
}

set.col.map <- function(optionfile, colmap){
  joined <- paste0(sapply(names(colmap), function(gr){
    paste0(gr, ":", colmap[gr])
  }), collapse="&")
  # - - - -
  setOption(optionfile, "gcols", joined)
}
