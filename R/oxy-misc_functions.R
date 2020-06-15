#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname trim
#' @export 
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname trim.trailing
#' @export 
trim.trailing <- function (x) sub("\\s+$", "", x)

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param lst PARAM_DESCRIPTION
#' @param type PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[data.table]{rbindlist}}
#' @rdname list.to.df
#' @export 
#' @importFrom pbapply pblapply
#' @importFrom data.table rbindlist
list.to.df <- function(lst,
                       type){
  dfStorage <- pbapply::pblapply(seq_along(lst), FUN=function(y, n, i){
    matches <- as.character(y[[i]])
    if((length(matches) == 0) && (typeof(matches) == "character")){ matches = NA }
    temp.df <- data.frame(stringsAsFactors = FALSE, 
                          "mz"=rep(n[[i]], length(matches)), 
                          "CompoundName"=matches[[1]], 
                          "Adduct"=matches[[2]], 
                          "Isotope"=matches[[3]], 
                          "Source"=rep(type, length(matches)))
    temp.df
  },y=lst, n=names(lst))
  df <- data.table::rbindlist(dfStorage)
  # --- return ---
  df
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param vector PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname factorize
#' @export 
factorize <- function(vector){
  items <- as.factor(vector)
  new.levels <- c(0:(length(levels(items)) - 1 ))
  levels(items) <- as.numeric(new.levels)
  # --- return ---
  items
}

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mz PARAM_DESCRIPTION
#' @param ppm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ppm_range
#' @export 
ppm_range <- function(mz, ppm) c((mz - (ppm/1000000 * mz)), (mz + (ppm/1000000 * mz)))

# plot molecules in R plot window instead of separate Java window
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param smi PARAM_DESCRIPTION
#' @param width PARAM_DESCRIPTION, Default: 500
#' @param height PARAM_DESCRIPTION, Default: 500
#' @param marg PARAM_DESCRIPTION, Default: 0
#' @param main PARAM_DESCRIPTION, Default: ''
#' @param style PARAM_DESCRIPTION, Default: 'bow'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[rcdk]{parse.smiles}},\code{\link[rcdk]{c("get.depictor", "get.depictor")}},\code{\link[rcdk]{c("view.image.2d", "view.image.2d")}}
#' @rdname plot.mol
#' @export 
#' @importFrom rcdk parse.smiles get.depictor view.image.2d
plot.mol = function(smi,
                    width=500,
                    height=500,
                    marg=0,
                    main='',
                    style="bow") {

  if(is.null(smi)){
    return(NULL)
  }
  molecule = rcdk::parse.smiles(smi)[[1]]
  par(mar=c(marg,marg,marg,marg)) # set margins to zero since this isn't a real plot
  dept = rcdk::get.depictor(width = width, height = height, zoom = 3, style = style,
                      annotate = "off", abbr = "on", suppressh = TRUE,
                      showTitle = FALSE, smaLimit = 100, sma = NULL)
  temp1 = rcdk::view.image.2d(molecule, dept) # get Java representation into an image matrix. set number of pixels you want horiz and vertical

  # - - return - -
  # A temp file to save the output. It will be deleted after renderImage
  # sends it, because deleteFile=TRUE.
  a <- tempfile(fileext='.png')
 
  plot(NA,NA,xlim=c(1,10),ylim=c(1,10),xaxt='n',yaxt='n',xlab='',ylab='',main=main,bg=NA) # create an empty plot
  rasterImage(temp1,1,1,10,10) # boundaries of raster: xmin, ymin, xmax, ymax. here i set them equal to plot boundaries

  # Return a list
  list(src = a)
}

# @export
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mzvals PARAM_DESCRIPTION
#' @param cl PARAM_DESCRIPTION, Default: FALSE
#' @param ppm PARAM_DESCRIPTION, Default: 3
#' @param charge PARAM_DESCRIPTION, Default: 1
#' @param element.counts PARAM_DESCRIPTION, Default: list(c("C", 0, 50), c("H", 0, 50), c("N", 0, 50), c("O", 0, 50), 
#'    c("S", 0, 50), c("Na", 0, 5), c("Cl", 0, 5), c("P", 0, 5))
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[pbapply]{pbapply}}
#'  \code{\link[rcdk]{generate.formula}}
#'  \code{\link[data.table]{rbindlist}}
#' @rdname find.formulas
#' @export 
#' @importFrom pbapply pblapply
#' @importFrom rcdk generate.formula
#' @importFrom data.table data.table rbindlist
find.formulas <- function(mzvals, cl=FALSE, ppm=3, charge=1, element.counts = list(c("C",0,50),c("H",0,50),
                                                                                   c("N",0,50),c("O",0,50),
                                                                                   c("S",0,50),c("Na", 0, 5),
                                                                                   c("Cl", 0, 5), c("P", 0,5))){

  found.rows <- pbapply::pblapply(mzvals,cl=cl, function(mz){
    window = mz * (ppm / 1e6)
    # --- generate molecular formulae ---
    found.mfs <- rcdk::generate.formula(mz, window = window, 
                                  element.counts, 
                                  validation=TRUE, charge=charge)
    rows <- if(length(found.mfs) == 0) NA else(
      rows <- lapply(found.mfs, function(formula){
        # --- check for ppm range ---
        mz.found <- formula@mass
        within.ppm <- abs(mz - mz.found) < window
        if(within.ppm){
          data.table::data.table(origMZ = mz,
                     genMZ = mz.found,
                     BaseFormula = formula@string)
        } else(data.table::data.table(origMZ=mz,
                          genMZ=NA,
                          BaseFormula=NA))
      })
    )
    # --- return ---
    unique(data.table::rbindlist(rows[!is.na(rows)]))
  })
  data.table::rbindlist(found.rows[!is.na(found.rows)])
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[shiny]{isolate}},\code{\link[shiny]{reactiveValuesToList}}
#' @rdname joanna_debugger
#' @export 
#' @importFrom shiny isolate reactiveValuesToList
joanna_debugger <- function(){
  lcl <<- debug_lcl
  mSet <<- debug_mSet
  input <<- debug_input
  enrich <<- debug_enrich 
  shown_matches <<- shiny::isolate({shiny::reactiveValuesToList(debug_matches)})
  my_selection <<- shiny::isolate({shiny::reactiveValuesToList(debug_selection)})
  browse_content <<- shiny::isolate({shiny::reactiveValuesToList(debug_browse_content)})
  pieinfo <<- shiny::isolate({shiny::reactiveValuesToList(debug_pieinfo)})
  result_filters <<- shiny::isolate({shiny::reactiveValuesToList(debug_result_filters)})
}

# @export
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param testsite PARAM_DESCRIPTION, Default: 'http://www.google.com'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param file.loc PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname getOptions
#' @export 
getOptions <- function(file.loc){
  opt_conn <- file(file.loc)
  # ----------------
  options_raw <- readLines(opt_conn)
  close(opt_conn)
  # --- list-ify ---
  options <- list()
  for(line in options_raw){
    split  <- (strsplit(line, ' = '))[[1]]
    options[[split[[1]]]] = split[[2]]
  }
  # --- return ---
  options
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param file.loc PARAM_DESCRIPTION
#' @param key PARAM_DESCRIPTION
#' @param value PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname setOption
#' @export 
setOption <- function(file.loc, key, value){
  opt_conn <- file(file.loc)
  # -------------------------
  options <- getOptions(file.loc)
  # --- add new or change ---
  options[[key]] = value
  # --- list-ify ---
  new_options <- lapply(seq_along(options), FUN=function(i){
    line <- paste(names(options)[i], options[i], sep=" = ")
    line
  })
  writeLines(opt_conn, text = unlist(new_options))
  
  close(opt_conn)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname get_os
#' @export 
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

# Customised TRUE-FALSE switch button for Rshiny
# Only sing CSS3 code (No javascript)
#
# Sébastien Rochette
# http://statnmap.com/en/
# April 2016
#
# CSS3 code was found on https://proto.io/freebies/onoff/
# For CSS3 customisation, refer to this website.


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param inputId PARAM_DESCRIPTION
#' @param label PARAM_DESCRIPTION
#' @param value PARAM_DESCRIPTION, Default: FALSE
#' @param col PARAM_DESCRIPTION, Default: 'GB'
#' @param type PARAM_DESCRIPTION, Default: 'TF'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname switchButton
#' @export 
switchButton <- function(inputId, label, value=FALSE, col = "GB", type="TF") {
  
  # # color class
  # if (col != "RG" & col != "GB" & col != "BW") {
  #   stop("Please choose a color between \"RG\" (Red-Green) 
  #        and \"GB\" (Grey-Blue) and Black-white(added).")
  # }
  # if (!type %in% c("OO", "TF", "YN", "TTFC", "2d3d")){
  #   warning("No known text type (\"OO\", \"TF\" or \"YN\") have been specified, 
  #           button will be empty of text") 
  # }
  if(col == "RG"){colclass <- "RedGreen"}
  if(col == "GB"){colclass <- "GreyBlue"}
  if(col == "BW"){colclass <- "BlackWhite"}
  if(type == "OO"){colclass <- paste(colclass,"OnOff")}
  if(type == "TF"){colclass <- paste(colclass,"TrueFalse")}
  if(type == "YN"){colclass <- paste(colclass,"YesNo")}
  if(type == "TTFC"){colclass <- paste(colclass,"TtFc")}
  if(type == "ASMB"){colclass <- paste(colclass,"AsMb")}
  if(type == "CLBR"){colclass <- paste(colclass,"ClBr")}
  
  if(type == "2d3d"){
    colclass <- paste(colclass,"twodthreed")
  }
  
  # No javascript button - total CSS3
  # As there is no javascript, the "checked" value implies to
  # duplicate code for giving the possibility to choose default value
  
  if(value){
    tagList(
      tags$div(class = "form-group shiny-input-container",
               tags$div(class = colclass,
                        tags$label(label, class = "control-label"),
                        tags$div(
                          class = "onoffswitch",
                          tags$input(type = "checkbox", name = "onoffswitch", class = "onoffswitch-checkbox",
                                     id = inputId, checked = ""
                          ),
                          tags$label(class = "onoffswitch-label", `for` = inputId,
                                     tags$span(class = "onoffswitch-inner"),
                                     tags$span(class = "onoffswitch-switch")
                          )
                        )
               )
      )
    )
  } else {
    tagList(
      tags$div(class = "form-group shiny-input-container",
               tags$div(class = colclass,
                        tags$label(label, class = "control-label"),
                        tags$div(class = "onoffswitch",
                                 tags$input(type = "checkbox", name = "onoffswitch", class = "onoffswitch-checkbox",
                                            id = inputId
                                 ),
                                 tags$label(class = "onoffswitch-label", `for` = inputId,
                                            tags$span(class = "onoffswitch-inner"),
                                            tags$span(class = "onoffswitch-switch")
                                 )
                        )
               )
      )
    ) 
  }
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param content PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[shiny]{builder}}
#' @rdname sardine
#' @export 
#' @importFrom shiny div
sardine <- function(content) shiny::div(style="display: inline-block;vertical-align:top;", content)

# loading screen
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param failed PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[shiny]{modalDialog}},\code{\link[shiny]{fluidPage}},\code{\link[shiny]{builder}}
#' @rdname loadModal
#' @export 
#' @importFrom shiny modalDialog fluidRow br h3 img
loadModal <- function(failed = FALSE) {
  shiny::modalDialog(
    shiny::fluidRow(align="center",
             shiny::br(),shiny::br(),
             shiny::h3("Starting MetaboShiny..."),
             shiny::br(),shiny::br(),
             #shiny::helpText("ଘ(੭ˊᵕˋ)੭* ੈ✩‧₊˚"), 
             #shiny::br(),shiny::br(),
             shiny::img(class="rotategem", src="gemmy_rainbow.png", width="70px", height="70px"),
             shiny::br(),shiny::br(),shiny::br()
    )
  )
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param y PARAM_DESCRIPTION
#' @param include.equals PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname expand.grid.unique
#' @export 
expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  
  y <- unique(y)
  
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param searchTerms PARAM_DESCRIPTION
#' @param retmax PARAM_DESCRIPTION, Default: 500
#' @param mindate PARAM_DESCRIPTION, Default: 2000
#' @param maxdate PARAM_DESCRIPTION, Default: 2019
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[RISmed]{EUtilsSummary}},\code{\link[RISmed]{QueryId}},\code{\link[RISmed]{EUtilsGet}}
#' @rdname getAbstracts
#' @export 
#' @importFrom RISmed EUtilsSummary QueryId EUtilsGet
getAbstracts <- function(searchTerms, retmax=500, mindate=2000, maxdate=2019){
  searchTerms = strsplit(searchTerms, " ")[[1]]
  # ==== SEARCH A METABOLITE TERM =====
  SUMMARYx <- RISmed::EUtilsSummary(paste0(searchTerms, collapse="+"),type = "esearch", db = "pubmed",
                            datetype = "edat",retmax = 500, 
                            mindate = 2000, maxdate = 2019)
  Idsx <- RISmed::QueryId(SUMMARYx)
  # Get Download(Fetch)
  MyDownloadsx <- RISmed::EUtilsGet(SUMMARYx, type = "efetch", db = "pubmed")
  #Make Data.frame of MyDownloads
  abstractsx <- data.frame(title = MyDownloadsx@ArticleTitle,
                           abstract = MyDownloadsx@AbstractText,
                           journal = MyDownloadsx@Title,
                           DOI = MyDownloadsx@PMID,
                           year = MyDownloadsx@YearPubmed)
  #Constract a charaterized variable for abstract of abstracts data.frame
  abstractsx <- abstractsx %>% mutate(abstract = as.character(abstract))
  #abstractsx$abstract <- as.character(abstractsx$abstract) # alternative for above line
  return(abstractsx)
}

#wordfrequency of whole abstract #frequencies = getWordFrequency(abstracts)#it should be # from here
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param abstractsx PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname getWordFrequency
#' @export 
getWordFrequency <- function(abstractsx){
  require(dplyr)
  abstractsx <- data.table(abstract = abstractsx)
  #Split a column into tokens using the tokenizers package
  CorpusofMyCloudx <- abstractsx %>% unnest_tokens(word, abstract) %>% count(word, sort = TRUE)
  CorpusofMyCloudx$word <- gsub("^\\d+$", "", CorpusofMyCloudx$word)
  return(CorpusofMyCloudx)
}

# frequencies <- getWordFrequency(abstractsx)#? or abstracts# why is it problem???????????????????
#filterList are medicalwords
#filterWords are summerized and uniqued of filterList
#getFilteredWordFreqency is filterd from medicalwords
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param frequencies PARAM_DESCRIPTION
#' @param filterList PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname getFilteredWordFreqency
#' @export 
getFilteredWordFreqency <- function(frequencies, filterList){
  filterWords <- unique(filterList)#$word)
  filteredWords <- frequencies[!(frequencies$word %in% filterWords),]
  filteredWords <- filteredWords[]
  return(filteredWords)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param str PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname escape
#' @export 
escape = function(str){
  str_adj = gsub(str, pattern = "( )", replacement = "\\\\\\1")
  str_adj = gsub(str_adj, pattern = "(\\))", replacement = "\\\\\\1")
  str_adj = gsub(str_adj, pattern = "(\\()", replacement = "\\\\\\1")
  str_adj = gsub(str_adj, pattern = "(&)", replacement = "\\\\\\1")
  # - - - - - 
  cat(str_adj)
}


