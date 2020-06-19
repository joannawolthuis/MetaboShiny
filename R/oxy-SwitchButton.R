# Customised TRUE-FALSE switch button for Rshiny
# Only sing CSS3 code (No javascript)
#
# Sébastien Rochette
# http://statnmap.com/en/
# April 2016
#
# CSS3 code was found on https://proto.io/freebies/onoff/
# For CSS3 customisation, refer to this website.


#' @title Generate custom switch button
#' @description Generates css for custom switch button.
#' @param inputId Shiny ID to use (behaves as checkbox)
#' @param label HTML label
#' @param value Default on/off? TRUE/FALSE, Default: FALSE
#' @param col Color binary used (GB green blue, GB grey blue), Default: 'GB'
#' @param type Pick display option TF true-false, OO on-off, YN yes-no, Default: 'TF'
#' @return HTML tag to display in shiny app
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  switchButton("doThing", "Do thing?", TRUE, "GB", "YN")
#'  }
#' }
#' @seealso 
#'  \code{\link[shiny]{tag}},\code{\link[shiny]{builder}}
#' @rdname switchButton
#' @export 
#' @importFrom shiny tagList tags
switchButton <- function(inputId, label, value=FALSE, col = "GB", type="TF") {
  
  # color class
  if (col != "RG" & col != "GB") {
    stop("Please choose a color between \"RG\" (Red-Green) 
      and \"GB\" (Grey-Blue).")
  }
  if (!type %in% c("OO", "TF", "YN")){
   warning("No known text type (\"OO\", \"TF\" or \"YN\") have been specified, 
     button will be empty of text") 
  }
  if(col == "RG"){colclass <- "RedGreen"}
  if(col == "GB"){colclass <- "GreyBlue"}
  if(type == "OO"){colclass <- paste(colclass,"OnOff")}
  if(type == "TF"){colclass <- paste(colclass,"TrueFalse")}
  if(type == "YN"){colclass <- paste(colclass,"YesNo")}
  
  # No javascript button - total CSS3
  # As there is no javascript, the "checked" value implies to
  # duplicate code for giving the possibility to choose default value
  
  if(value){
    shiny::tagList(
      shiny::tags$div(class = "form-group shiny-input-container",
        shiny::tags$div(class = colclass,
          shiny::tags$label(label, class = "control-label"),
          shiny::tags$div(class = "onoffswitch",
            shiny::tags$input(type = "checkbox", name = "onoffswitch", class = "onoffswitch-checkbox",
              id = inputId, checked = ""
            ),
            shiny::tags$label(class = "onoffswitch-label", `for` = inputId,
              shiny::tags$span(class = "onoffswitch-inner"),
              shiny::tags$span(class = "onoffswitch-switch")
            )
          )
        )
      )
    )
  } else {
    shiny::tagList(
      shiny::tags$div(class = "form-group shiny-input-container",
        shiny::tags$div(class = colclass,
          shiny::tags$label(label, class = "control-label"),
          shiny::tags$div(class = "onoffswitch",
            shiny::tags$input(type = "checkbox", name = "onoffswitch", class = "onoffswitch-checkbox",
              id = inputId
            ),
            shiny::tags$label(class = "onoffswitch-label", `for` = inputId,
              shiny::tags$span(class = "onoffswitch-inner"),
              shiny::tags$span(class = "onoffswitch-switch")
            )
          )
        )
      )
    ) 
  }
}


#' @title Small image checkbox that fades upon deselection
#' @description Used in MetaboShiny database selection when searching or pre-matching
#' @param inputId Shiny input ID to connect to
#' @param img.path Path to image to use, Default: NULL
#' @param value Default value (T/F), Default: FALSE
#' @return Shiny HTML tag
#' @seealso 
#'  \code{\link[shiny]{tag}},\code{\link[shiny]{builder}}
#' @rdname fadeImageButton
#' @export 
#' @importFrom shiny tagList tags
fadeImageButton <- function(inputId, img.path=NULL,value=FALSE) {
  # ---------------
  #if(is.null(img)) stop("Please enter an image name (place in www folder please)")
  # ---------------
  shiny::tagList(
    shiny::tags$div(class = "form-group shiny-input-container",
             if(value){
               shiny::tags$input(type="checkbox", class="fadebox", id=inputId)
             } else{
               shiny::tags$input(type="checkbox", class="fadebox", id=inputId, checked="")
             },
    shiny::tags$label(class="btn", "for"=inputId,
               shiny::tags$img(src=img.path, id="btnLeft"))
    )
  )
}
