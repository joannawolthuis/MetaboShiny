#' @title Get m/z intensity profile
#' @description Given a m/z, return intensities for all samples for this m/z.
#' @param mSet mSet object
#' @param varName m/z
#' @param title Title of plot, Default: varName
#' @param mode stat or multi mode? , Default: 'stat'
#' @return Data table with m/z intensities per sample.
#' @rdname getProfile
#' @export 
getProfile <- function(mSet, varName, title=varName, mode="stat"){
  sourceTable = mSet$dataSet$norm
  # ---------------
  varInx <- colnames(sourceTable) == varName;
  var <- data.table::as.data.table(sourceTable,
                       keep.rownames = T)[,varInx, with=FALSE];
  samp.names <- rownames(sourceTable)
  # ---------------
  if(mode == "multi"){
    if(mSet$dataSet$exp.type == "t"){
      translator <- data.table(
        index = 1:length(samp.names),
        Sample = gsub(x = samp.names, pattern = "_T|_t\\d$", replacement=""),
        GroupA = mSet$dataSet$time.fac,
        GroupB = mSet$dataSet$time.fac,
        Abundance = sourceTable[,varInx]
      )
    }else{
      translator <- data.table(
        index = 1:length(samp.names),
        Sample = gsub(x = samp.names, pattern = "_T|_t\\d$", replacement=""),
        GroupA = mSet$dataSet$facA,
        GroupB = mSet$dataSet$facB,
        Abundance = sourceTable[,varInx]
      )  
    }
    
  }else if(mode == "stat"){
    translator <- data.table(
      index = 1:length(samp.names),
      Sample = gsub(x = samp.names, pattern = "T\\d$", replacement=""),
      Group = mSet$dataSet$cls,
      Abundance = sourceTable[,varInx]
    )
  }
#  translator = translator[complete.cases(translator),]
  # ---------------
  return(translator)
}

#' @title Mean absolute percentage error calculation
#' @description Calculates MAPE (Mean absolute percentage error). Will be used as one of the options in isotope pattern scoring.
#' @param actual Data values
#' @param pred Predicted values (for example, isotope pattern)
#' @return MAPE score
#' @rdname mape
#' @export 
mape <- function(actual,pred){
  mape <- mean(abs((actual - pred)/actual))*100
  return (mape)
}

#' @title Un-nest a list
#' @description Flatten a list of lists into one long list
#' @param x list object
#' @return list object
#' @rdname flattenlist
#' @export 
flattenlist <- function(x){
  morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
  out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE, use.names = T))
  if(sum(morelists)){
    Recall(out)
  }else{
    return(out)
  }
}

#' @title Get color options
#' @description Retrieves color options from user options file.
#' @param optionfile Path to options file.
#' @return Color options (hex)
#' @rdname get.col.map
#' @export 
get.col.map <- function(optionfile){
  options <- getOptions(optionfile)
  unparsed.cols <- options$gcols
  col.items <- strsplit(unparsed.cols, split = "&")[[1]]
  # - - - -
  col.items
}

#' @title Set color options
#' @description Sets color options in user options file.
#' @param optionfile Path to options file.
#' @param colmap Vector of hex colors.
#' @rdname set.col.map
#' @export 
set.col.map <- function(optionfile, colmap){
  joined <- paste0(
    colmap, collapse="&")
  # - - - -
  setOption(optionfile, key="gcols", value=joined)
}

#' @title Generate word cloud from list of abstracts
#' @description Wrapper function to create a wordcloud from a literature search.
#' @param abstracts List of abstracts (RISmed format)
#' @param queryword Searched word (will be added to words not to include in the cloud)
#' @param top Top x words used in plot, Default: 20
#' @return Data table that can be used to generate word cloud (word vs. frequency)
#' @seealso 
#'  \code{\link[RISmed]{AbstractText}}
#'  \code{\link[qdap]{strip}},\code{\link[qdap]{rm_stopwords}}
#' @rdname abstracts2wordcloud
#' @export 
#' @importFrom RISmed AbstractText
#' @importFrom qdap strip rm_stopwords
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

#' @title P value to 'stars'
#' @description Get one to four star significance from a given p value.
#' @param pval P value.
#' @return Character vector of one to four stars (or n.s. non significant)
#' @rdname p2stars
#' @export 
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

#' @title Are you currently online?
#' @description Function to check if the user can connect to the internet.
#' @param testsite Website to connect to for the test, Default: 'http://www.google.com'
#' @return TRUE/FALSE
#' @rdname internetWorks
#' @export 
internetWorks <- function(testsite = "http://www.google.com"){
  works = FALSE
  try({
    GET(testsite)
    works=TRUE
  },silent = T)
  works
}

#' @title Show alert in MetaboShiny
#' @description Function to create a SweetAlert in MetaboShiny with user message.
#' @param message User message
#' @param title Title of alert
#' @param myImg Path to image in /www folder
#' @param session Shiny session, Default: shiny::getDefaultReactiveDomain()
#' @return SweetAlert object to render in shiny
#' @seealso 
#'  \code{\link[shiny]{domains}},\code{\link[shiny]{builder}}
#'  \code{\link[shinyWidgets]{sendSweetAlert}}
#' @rdname metshiAlert
#' @export 
#' @importFrom shiny getDefaultReactiveDomain img
#' @importFrom shinyWidgets sendSweetAlert
metshiAlert <- function(content,
                        session = shiny::getDefaultReactiveDomain(),
                        title = "Error",
                        myImg = "gemmy_rainbow.png"){
  if(typeof(content) == "character"){
    content = h3(content)
  }
  
  shinyWidgets::sendSweetAlert(
    session = session,
    title = title,
    text = tags$div(
      shiny::img(#class = "rotategem", 
                 src = myImg, 
                 #width = "30px", 
                 height = "30px"),
      br(),
      content
    ),
    html = TRUE
  )
}

#' @title Reformat metadata
#' @description Clean matadata to be imported in MetShi. Often used to import 'new' metadata into an existing dataset.
#' @param metadata Metadata data table
#' @return Cleaned data table
#' @rdname reformat.metadata
#' @export 
reformat.metadata <- function(metadata){
  colnames(metadata) <- tolower(colnames(metadata))
  colnames(metadata) <- tolower(gsub(x=colnames(metadata), pattern = "\\.$|\\.\\.$|\\]", replacement = ""))
  colnames(metadata) <- gsub(x=colnames(metadata), pattern = "\\[|\\.|\\.\\.| ", replacement = "_")
  colnames(metadata) <- gsub(colnames(metadata), pattern = "characteristics_|factor_value_", replacement="")
  metadata[metadata == "" | is.na(metadata)] <- c("unknown")
  setnames(metadata, old = c("sample_name", "source_name"), new = c("sample", "individual"), skip_absent = T)
  
  # - - - 
  return(metadata)
}
