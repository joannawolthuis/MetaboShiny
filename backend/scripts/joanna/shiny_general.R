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
  # ---------------
  if(mainmode == "time"){
    time.fac <- dataSet$time.fac;
    translator <- data.table(
      index = 1:length(samp.names),
      Sample = gsub(x = samp.names, pattern = "T\\d$", replacement=""),
      Group = exp.fac,
      Time = time.fac,
      Abundance = dataSet$norm[,varInx]
    )
  }else if(mainmode == "stat"){
    exp.fac <- dataSet$filt.cls
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