#' @title Generate css for MetShi footer
#' @description Takes user aesthetic settings and generates footer CSS and hTML
#' @param RGB.bg Hex background color, Default: '#ffc8c2'
#' @param font.col Font color, Default: 'white'
#' @return CSS
#' @seealso 
#'  \code{\link[gsubfn]{fn}}
#' @rdname footer.css
#' @export 
#' @importFrom gsubfn fn
footer.css <- function(RGB.bg = "#ffc8c2", font.col="white"){
  base <- ".footer {position:fixed;
                    bottom:0;
                    width:100%;
                    height:40px;
                    z-index:1005;
                    color: white;
                    background-color:$RGB.bg;
                    border-style:solid; 
                    border-color:black;
                    border-width:1px; 
                    border-bottom: black;
                    border-left: black;
                    border-right: black;}
            .footer .btn {margin-top:2px; 
                          background-color:$RGB.bg; 
                          color:$font.col; 
                          border-style:none;}"
  no.na <- gsub(x = base, pattern = "\\n", replacement = "")
  filled <- gsubfn::fn$paste(no.na)
  # --- return ---
  filled
}

#' @title Generate navigation bar CSS for use in MetaboShiny
#' @param RGB.bg Hex background color, Default: '#ffc8c2'
#' @param RGB.active.bg Hex active background color, Default: '#e1897f'
#' @param RGB.col Hex text color, Default: '#ffffff'
#' @param RGB.active.col Hex active text color, Default: '#ffffff'
#' @return CSS string
#' @seealso 
#'  \code{\link[gsubfn]{fn}}
#' @rdname nav.bar.css
#' @export 
#' @importFrom gsubfn fn
nav.bar.css <- function(RGB.bg = "#ffc8c2",
                        RGB.active.bg = "#e1897f",
                        RGB.col = "#ffffff",
                        RGB.active.col = "#ffffff"){

  base <- ".navbar { background-color: $RGB.bg !important;
  border: 0px !important;
  margin-bottom: 5px;}
  .nav > li[class=active] > a {
  background-color: $RGB.active.bg !important;
  color: $RGB.active.col !important; }
  .nav > li > a {
  background-color: $RGB.bg !important;
  color: $RGB.col !important;
  }
  .navbar:after{
      height:30px;
      width: 100%;
      background: -webkit-gradient(radial, 50% 0, 10, 50% 0, 40, from($RGB.bg), color-stop(0.49, $RGB.bg), color-stop(0.51, #fff), to(white));
    -webkit-background-size: 29px 100%;
  }
  .focus_b 
  {
    background-color: $RGB.active.bg !important;
    color: $RGB.active.col !important; 
    font-weight: bolder;
    /*padding: 3px; add this if you want a bit of space around the text*/
  }
  
  .focus_a
  {
    background-color: $RGB.bg !important;
    color: $RGB.col !important;
    font-weight: bold;
    /*padding: 3px; add this if you want a bit of space around the text*/
  }"

  no.na <- gsub(x = base, pattern = "\\n", replacement = "")
  filled <- gsubfn::fn$paste(no.na)
  # --- return ---
  filled
}

#' @title Generate font CSS
#' @description Takes user options and generates CSS to set font correctly on start.
#' @param font.h1 h1 font, Default: 'Press Start 2P'
#' @param font.h2 h2 font, Default: 'Raleway bold'
#' @param font.h3 h3 font, Default: 'Raleway'
#' @param font.body body font, Default: 'Raleway serif'
#' @param size.h1 h1 size, Default: 33
#' @param size.h2 h3 size, Default: 22
#' @param size.h3 h3 size, Default: 15
#' @param size.body body size, Default: 13
#' @param font.col Font color, Default: 'white'
#' @param online Are we online? (needed for proper font setting), Default: T
#' @return CSS string
#' @seealso 
#'  \code{\link[stringr]{str_extract}}
#'  \code{\link[gsubfn]{fn}}
#' @rdname app.font.css
#' @export 
#' @importFrom stringr str_extract
#' @importFrom gsubfn fn
app.font.css <- function(font.h1 = "Press Start 2P",
                        font.h2 = "Raleway bold",
                        font.h3 = "Raleway",
                        font.body = "Raleway serif",
                        size.h1 = 33,
                        size.h2 = 22,
                        size.h3 = 15,
                        size.body = 13,
                        font.col = "white",
                        online=T){

  size.h1 <- paste0(size.h1, "pt")
  size.h2 <- paste0(size.h2, "pt")
  size.h3 <- paste0(size.h3, "pt")
  size.body <- paste0(size.body, "pt")

  font.h1.w <- stringr::str_extract(font.h1, "(serif)|(bold)|(italic)")
  font.h2.w <- stringr::str_extract(font.h2, "(serif)|(bold)|(italic)")
  font.h3.w <- stringr::str_extract(font.h3, "(serif)|(bold)|(italic)")
  font.body.w <- stringr::str_extract(font.body, "(serif)|(bold)|(italic)")

  if(is.na(font.h1.w)) font.h1.w = ""
  if(is.na(font.h2.w)) font.h2.w = ""
  if(is.na(font.h3.w)) font.h3.w = ""
  if(is.na(font.body.w)) font.body.w = ""

  font.h1 <- gsub(font.h1, pattern = "( serif)|( bold)|( italic)", replacement = "")
  font.h2 <- gsub(font.h2, pattern = "( serif)|( bold)|( italic)", replacement = "")
  font.h3 <- gsub(font.h3, pattern = "( serif)|( bold)|( italic)", replacement = "")
  font.body <- gsub(font.body, pattern = "( serif)|( bold)|( italic)", replacement = "")

  fonts <- c(font.h1, font.h2, font.h3, font.body)

  font_weights <- list(font.h1.w,
                       font.h2.w,
                       font.h3.w,
                       font.body.w)

  names(font_weights) <- fonts

  font_joined <- sapply(unique(names(font_weights)), function(x) unname(unlist(font_weights[names(font_weights)==x])), simplify=FALSE)

  font.url <- sapply(1:length(font_joined), function(i){
    font.types <- font_joined[[i]]
    font.name <- names(font_joined)[[i]]
    font.name <- gsub(x = font.name, pattern=" ", replacement = "+")
    font.types <- font.types[font.types != ""]
    if(length(font.types) == 0){
      font.name
    }else{
      paste0(font.name, ":", paste0(font.types, collapse=","))
    }
  })

  font.url <- paste0(font.url, collapse = "|")

  font.h1.w = if(font.h1.w == "") "" else{paste0("font-weight: ", font.h1.w, ";")}
  font.h2.w = if(font.h2.w == "") "" else{paste0("font-weight: ", font.h2.w, ";")}
  font.h3.w = if(font.h3.w == "") "" else{paste0("font-weight: ", font.h3.w, ";")}
  font.body.w = if(font.body.w == "") "" else{paste0("font-weight: ", font.body.w, ";")}

  base.import <- if(online) gsubfn::fn$paste('@import url("https://fonts.googleapis.com/css?family=$font.url&display=swap");') else ""

  base <- '$base.import

  #appHeader {
  margin: 2px;
  font-family: "$font.h1";
  color: $font.col;
  $font.h1.w
  font-size: $size.h1;
  line-height: 0.5;
  margin-top: 8px;
  z-index:99998 !important;
  height: 200px;
  line-height: 60px;
  padding-right: 20px;
  padding-left: 10px;
  }
  
  h1 {
  font-family: "$font.h1";
  }

  h2 {
  margin: 1px;
  color: black;
  font-family: "$font.h2";
  $font.h2.w
  font-size: $size.h2;
  line-height: 1;
  }

  h3 {
  margin: 1px;
  color: black;
  font-family: "$font.h3";
  $font.h3.w
  font-size: $size.h3;
  line-height: 1;
  }

  body {
  color: black;
  font-family: "$font.body";
  $font.body.w
  font-size: $size.body;
  background-color: black;
  }'

  no.na <- gsub(x = base, pattern = "\\n", replacement = "")
  filled <- gsubfn::fn$paste(no.na)
  # -----
  filled
}

#https://stackoverflow.com/questions/31425841/css-for-each-page-in-r-shiny
## Modify the CSS style of a given selector
modifyStyle <- function(selector, ...) {
  
  values <- as.list(substitute(list(...)))[-1L]
  parameters <- names(values)
  
  args <- Map(function(p, v) paste0("'", p,"': '", v,"'"), parameters, values)
  jsc <- paste0("$('",selector,"').css({", paste(args, collapse = ", "),"});")
  
  shinyjs::runjs(code = jsc)
  
}
