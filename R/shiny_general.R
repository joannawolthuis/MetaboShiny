getProfile <- function(mSet, varName, title=varName, mode="stat"){
  sourceTable = mSet$dataSet$norm
  # ---------------
  varInx <- colnames(sourceTable) == varName;
  var <- as.data.table(sourceTable,
                       keep.rownames = T)[,varInx, with=FALSE];
  samp.names <- rownames(sourceTable)
  # ---------------
  if(mode == "multi"){
    translator <- data.table(
      index = 1:length(samp.names),
      Sample = gsub(x = samp.names, pattern = "_T|_t\\d$", replacement=""),
      GroupA = mSet$dataSet$facA,
      GroupB = mSet$dataSet$facB,
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
  charges <- regmatches(
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
  col.items <- strsplit(unparsed.cols, split = "&")[[1]]
  # - - - -
  col.items
}

set.col.map <- function(optionfile, colmap){
  joined <- paste0(
    colmap, collapse="&")
  # - - - -
  setOption(optionfile, key="gcols", value=joined)
}

abstracts2wordcloud <-function(abstracts, queryword,  top=20){
  abstracts1 <- data.frame('Abstract' = RISmed::AbstractText(abstracts))#, 'Year'=YearPubmed(fetch))
  abstractsOnly <- as.character(abstracts1$Abstract)
  abstractsOnly <- paste(abstractsOnly, 
                       sep="", 
                       collapse="")
  abstractsOnly <- as.vector(abstractsOnly)
  abstractsOnly <- qdap::strip(abstractsOnly)
  stsp <- qdap::rm_stopwords(abstractsOnly, 
                           stopwords = c(gbl$vectors$wordcloud$skip, queryword))
  ord <- as.data.frame(table(stsp))
  ord <- ord[order(ord$Freq, decreasing=TRUE),]
  head(ord, top)
}

p2stars = function(pval){
  if(length(pval) == 0){
    stars <- ""
  }else{
    if(pval > 0.05) stars <- "n.s."
    else if(pval < 0.05 & pval > 0.01) stars <- "*"
    else if(pval < 0.01 & pval > 0.001) stars <- "***"
    else stars <- "****"
  }
}


havingIP <- function(){
  if(.Platform$OS.type == "windows"){
    ipmessage = system("ipconfig", intern=T)
  }else{
    ipmessage = system("ifconfig", intern=T)
  }
  print(ipmessage)
  validIP = "((25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)[.]){3}(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)"
  any(grep(validIP, ipmessage))
}

internetWorks <- function(testsite = "http://www.google.com"){
  works = FALSE
  try({
    GET(testsite)
    works=TRUE
  },silent = T)
  works
}
