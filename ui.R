shinyUI(
  fluidPage(theme = "metaboshiny.css",
            tags$head(tags$script(src = "sparkle.js")),
            ECharts2Shiny::loadEChartsLibrary(),
            shinyalert::useShinyalert(), # enable shinyalert (necessary for on close message)
            shinyjqui::includeJqueryUI(),
            windowTitle = "MetaboShiny",
    uiOutput("currUI")
  )
)