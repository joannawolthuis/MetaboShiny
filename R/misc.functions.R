trim <- function (x) gsub("^\\s+|\\s+$", "", x)

trim.trailing <- function (x) sub("\\s+$", "", x)

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

factorize <- function(vector){
  items <- as.factor(vector)
  new.levels <- c(0:(length(levels(items)) - 1 ))
  levels(items) <- as.numeric(new.levels)
  # --- return ---
  items
}

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

ppm_range <- function(mz, ppm) c((mz - (ppm/1000000 * mz)), (mz + (ppm/1000000 * mz)))

# plot molecules in R plot window instead of separate Java window
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

#'@export
joanna_debugger <- function(){
  lcl <<- debug_lcl
  mSet <<- debug_mSet
  input <<- debug_input
  shown_matches <<- shiny::isolate({shiny::reactiveValuesToList(debug_matches)})
  my_selection <<- shiny::isolate({shiny::reactiveValuesToList(debug_selection)})
}

# @export
internetWorks <- function(testsite = "http://www.google.com"){
  works = FALSE
  try({
    GET(testsite)
    works=TRUE
  },silent = T)
  works
}  

#' Load user options saved in file.
#'
#' \code{getOptions} returns all current user options defined in the given options file.
#'
#' @param file.loc Path to user options file to read in.
#' @return R list with keys as option types and values as option values.
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


#' Changes an option in the given MetaboShiny user options file.
#'
#' \code{setOption} changes an option in the options file. Can also be used to add new options to the file.
#'
#' @param file.loc Location of user options file. Usually .txt but any format is fine.
#' @param key Name of the new option / to change option
#' @param value Value of the option to change or add 
setOption <- function(file.loc, key, value){
  opt_conn <- file(file.loc)
  # -------------------------
  options <- MetaboShiny::getOptions(file.loc)
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

#' Gets the current used operating system. Important for parallel/multithreaded functions if using makeCluster("FORK")
#' \code{get_os} finds the name of the OS the user is running this function on.
#'
#' @return osx, windows/win or linux
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

#' A function to change the Original checkbox of rshiny
#' into a nice true/false or on/off switch button
#' No javascript involved. Only CSS code.
#' 
#' To be used with CSS script 'button.css' stored in a 'www' folder in your Shiny app folder
#' 
#' @param inputId The input slot that will be used to access the value.
#' @param label Display label for the control, or NULL for no label.
#' @param value Initial value (TRUE or FALSE).
#' @param col Color set of the switch button. Choose between "GB" (Grey-Blue) and "RG" (Red-Green)
#' @param type Text type of the button. Choose between "TF" (TRUE - FALSE), "OO" (ON - OFF) or leave empty for no text.

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

#' Squishes HTML elements close together.
sardine <- function(content) shiny::div(style="display: inline-block;vertical-align:top;", content)

# loading screen
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
