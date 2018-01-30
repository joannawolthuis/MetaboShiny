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

switchButton <- function(inputId, label, value=FALSE) {
  tagList(
    tags$div(class = "form-group shiny-input-container",
             tags$div(class = "plus-minus",
                      tags$label(class = "css-label", `for` = "a"),
                      tags$div(class = "fa fa-plus",
                               tags$input(type = "checkbox", class = "fa fa-minus",
                                          id = inputId
                               )
                      )
             )
    )
  )
}

# '<input type="checkbox" id="btnControl"/>
#                                                        <label class="btn" for="btnControl">
# <img src="http://placekitten.com.s3.amazonaws.com/homepage-samples/200/140.jpg" id="btnLeft"/>
# </label>'
# 
# img = "hereisyourpicture"
# tags$img(src="...")
# tagList(
#   tags$input(type="checkbox", id="btnControl"),
#   tags$div(class="form-group shiny-input-container",
#            tags$label(class="btn", "for"="btnControl", tags$img(src="http://placekitten.com.s3.amazonaws.com/homepage-samples/200/140.jpg", id="btnLeft"))
#   )
# )


fadeImageButton <- function(inputId, img.path=NA,value=FALSE) {
  # ---------------
  if(is.na(img)) stop("Please enter an image name (place in www folder please)")
  # ---------------
  tagList(
    tags$div(class = "form-group shiny-input-container",
             if(value){
               tags$input(type="checkbox", class="fadebox", id=inputId)
             } else{
               tags$input(type="checkbox", class="fadebox", id=inputId, checked="")
             },
    tags$label(class="btn", "for"=inputId,
               tags$img(src=img.path, id="btnLeft"))
    )
  )
}
