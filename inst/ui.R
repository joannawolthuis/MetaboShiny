shiny::fluidPage(
  class = "hidden", id = "metshi",
  ECharts2Shiny::loadEChartsLibrary(),
  shinyjs::useShinyjs(),
  shinyDarkmode::use_darkmode(),
  shiny::includeCSS("www/metaboshiny.css"),
  tags$head(tags$script(src = "cursor.js")),
  shiny::div(
    style = "position: absolute;
    left: 79%;
    top: 1%;",
    class = "plus",
    shiny::img(
      src = "metshi_gemmo.png",
      id = "metshiGem",
      style = "position: relative;
      height: 100px;
      top: 7%;
      left: 88%;
      z-index: 1005;"
    ),
    shiny::div(
      id = "heartHolder",
      style = "left: 22%;
      bottom: 4%;
      position: relative;",
      MetaboShiny::fadeImageButton("fancy",
                                   img.path = "metshi_heart.png",
                                   value = F
      )
    )
  ),
  shinyjs::extendShinyjs(
    text = "shinyjs.closeWindow = function() { window.close(); }",
    functions = c("closeWindow")
  ),
  shiny::div(
    id = "loading-page",
    style = "position: fixed;
    width: 200%;
    height: 200%;
    z-index: 4000;
    background-color: black;
    opacity: 0.6;
    margin-left: -20px;",
    div(
      id = "load-img-holder",
      class = "imagetop",
      style = "left: 25%;
      position: absolute;
      top: 25%;
      height:100px;
      width:100px;
      margin-top: -60px;
      margin-left: -60px;",
      div(
        id = "loading-bg",
        style = "background-image: url(metshi_heart_bezel.png);
                                    width:120px;
                                    height:100px;"
      )
    )
  ),
  shiny::div(
    shiny::navbarPage(
      windowTitle = "MetaboShiny",
      # use this for title
      # https://codepen.io/maxspeicher/pen/zrVKLE
      title = shiny::div(
        id = "appHeader",
        class = "outlined",
        "MetaboShiny"
      ), # make it use the sparkle.js for unnecessary sparkle effects ;)
      id = "nav_general",
      # this tab shows the available databases, if they are installed, and buttons to install them. generated as output$db_build_ui in 'server'
      MetaboShiny:::ui_tab_database(gbl, adducts),
      MetaboShiny:::ui_tab_data_import(gbl, adducts),
      MetaboShiny:::ui_tab_normalize(gbl, adducts),
      MetaboShiny:::ui_tab_prematch(gbl, adducts),
      MetaboShiny:::ui_tab_analyse(gbl, adducts),
      MetaboShiny:::ui_tab_settings(gbl, adducts),
      MetaboShiny:::ui_tab_help(gbl, adducts),
      # prompt user on opening the quit tab.
      # TODO: add 'save project?' dialog
      shiny::div(class = "scallop-down"),
      shiny::div(class = "cursorHolder"),
      shiny::div(class = "line"),
      footer = MetaboShiny:::ui_footer()
    ),
    style = "margin-bottom:100px;"
  )
)
