# for future reference: https://www.r-bloggers.com/deploying-desktop-apps-with-r/ :-)

function(input, output, session) {
  library(data.table)
  library(plotly)
  library(ggplot2)
  library(MetaboShiny)
  library(MetaDBparse)
  library(plyr)
  library(bslib)
  library(dplyr)
  library(pathview)

  Sys.setenv(VROOM_CONNECTION_SIZE = "500000000")

  options(shiny.maxRequestSize = 50000 * 1024^2)
  setTimeLimit(cpu = Inf)

  # FIX FOR NORMALIZATION
  OFFtoJSON <- function(obj, ...) {
    print("Disabled in MetaboShiny!")
  }

  # `rlang::env_unlock()` is defunct (rlang >= 1.1.5). We only need to unlock the
  # specific binding we overwrite.
  try({
    rjsonio_ns <- asNamespace("RJSONIO")
    if (exists("toJSON", envir = rjsonio_ns, inherits = FALSE)) {
      was_locked <- bindingIsLocked("toJSON", rjsonio_ns)
      if (was_locked) {
        unlockBinding("toJSON", rjsonio_ns)
      }
      assign("toJSON", OFFtoJSON, envir = rjsonio_ns)
      if (was_locked) {
        lockBinding("toJSON", rjsonio_ns)
      }
    }
  })

  AddErrMsg <- function(msg) {
    print(msg)
    try({
      shiny::showNotification(msg)
    })
  }

  localize_fn <- function(fn) {
    environment(fn) <- environment()
    fn
  }

  ensure_dir <- localize_fn(MetaboShiny:::ensure_dir)
  show_tabs <- localize_fn(MetaboShiny:::show_tabs)
  hide_tabs <- localize_fn(MetaboShiny:::hide_tabs)
  read_table_if_in_dir <- localize_fn(MetaboShiny:::read_table_if_in_dir)
  render_combi_picker <- localize_fn(MetaboShiny:::render_combi_picker)
  format_options_lines <- localize_fn(MetaboShiny:::format_options_lines)
  default_options_template <- localize_fn(MetaboShiny:::default_options_template)
  set_option_default <- localize_fn(MetaboShiny:::set_option_default)
  update_theme_inputs <- localize_fn(MetaboShiny:::update_theme_inputs)
  set_option_value <- localize_fn(MetaboShiny:::set_option_value)
  render_text_outputs <- localize_fn(MetaboShiny:::render_text_outputs)
  update_ml_method <- localize_fn(MetaboShiny:::update_ml_method)
  stop_session_cluster <- localize_fn(MetaboShiny:::stop_session_cluster)
  start_session_cluster <- localize_fn(MetaboShiny:::start_session_cluster)
  render_image_outputs <- localize_fn(MetaboShiny:::render_image_outputs)
  update_theme_controls <- localize_fn(MetaboShiny:::update_theme_controls)
  show_stat_panels <- localize_fn(MetaboShiny:::show_stat_panels)
  set_stat_collapse <- localize_fn(MetaboShiny:::set_stat_collapse)
  make_sel_adducts <- localize_fn(MetaboShiny:::make_sel_adducts)
  update_adducts_from_filter <- localize_fn(MetaboShiny:::update_adducts_from_filter)

  shinyDarkmode::darkmode_toggle(inputid = "night_mode")

  shiny::observe({
    if (input$night_mode) {
      shinyjs::addClass(
        selector = ".dbimg > .shiny-image-output > img",
        class = "antirevert"
      )
    } else {
      shinyjs::removeClass(
        selector = ".dbimg > .shiny-image-output > img",
        class = "antirevert"
      )
    }
  })

  assignInNamespace("AddErrMsg", AddErrMsg,
    ns = "MetaboAnalystR",
    envir = as.environment("package:MetaboAnalystR")
  )

  shiny::showNotification("Starting server process...")

  # detach("package:MetaboShiny", unload=T)
  # used to be in startshiny.R
  options("download.file.method" = "libcurl")
  options(expressions = 5e5)
  online <- MetaboShiny::internetWorks()

  # make MetaboShiny_storage dir in home first..
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/userfiles/:cached --rm -it metaboshiny/master /bin/bash
  # with autorun
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/userfiles/:cached --rm metaboshiny/master Rscript startShiny.R
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny /bin/bash
  # docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny Rscript startShiny.R
  # current instructions
  # Rshiny app to analyse untargeted metabolomics data! BASH INSTRUCTIONS: STEP 1: mydir=~"/MetaboShiny" #or another of your choice | STEP 2: mkdir $mydir | STEP 3: docker run -p 8080:8080 -v $mydir:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny /start.sh

  # rjava.so error.. or rdb corrupt.. 'sudo R CMD javareconf'

  runmode <- if (file.exists(".dockerenv")) "docker" else "local"

  mSet <- NULL
  opts <- list()
  showtext::showtext_auto(enable = T)
  # showtext::showtext_opts(dpi=72)

  lcl <- list(
    proj_name = "",
    last_mset = "",
    hasChanged = FALSE,
    load_ui = FALSE,
    last_dir = c(),
    lists = list(),
    prev_mz = "",
    prev_struct = "",
    beep = F,
    tables = list(
      last_matches = data.table::data.table(query_mz = "none"),
      prev_pie = data.table::data.table()
    ),
    functions = list(),
    aes = list(
      font = list(),
      mycols = c(),
      spectrum = "rb",
      theme = "min"
    ),
    vectors = list(
      proj_names = c(),
      built_dbs = gbl$vectors$db_list,
      custom_db_list = c()
    ),
    paths = list(
      opt.loc = "",
      patdb = "",
      work_dir = ""
    )
  )

  reinstall_restart <- function() {
    devtools::install(reload = T, upgrade = F)
    # require(MetaboShiny)
    MetaboShiny::start_metshi(inBrowser = T)
  }

  shiny::showModal(loadModal())

  # == REACTIVE VALUES ==

  MetaboShiny:::init_reactive_values(environment())


  MetaboShiny:::init_user_settings_observer(environment())

  # ================================= DEFAULTS ===================================

  # set progress bar style to 'old' (otherwise it's not movable with CSS)
  shiny::shinyOptions(progress.style = "old")

  options(
    digits = 22,
    spinner.size = 0.5,
    spinner.type = 6,
    spinner.color = "black",
    spinner.color.background = "white"
  )

  # create default text objects in UI
  render_text_outputs(gbl$constants$default.text)

  showtext::showtext_auto() ## Automatically use showtext to render text for future devices

  # this toggles when 'interface' values change (for example from 'bivar' to 'multivar' etc.)
  MetaboShiny:::init_interface_tabs_observer(environment())

  MetaboShiny:::init_core_observers(environment())

  render_combi_picker("combi_anal1", "combi_anal1_picker", "A result column:")
  render_combi_picker("combi_anal2", "combi_anal2_picker", "B result column:")

  MetaboShiny:::init_debug_observer(environment())

  MetaboShiny:::init_export_observer(environment())

  # triggered when user enters the statistics tab
  shinyjs::runjs('$("#mainPanel").resizable({
                                              handles: "e",
                                              resize: function() {
                                                $("#sidePanel").outerWidth($("#panelContainer").innerWidth() - $("#mainPanel").outerWidth());
                                              }
                                            });')

  MetaboShiny:::init_statistics_observer(environment())

  MetaboShiny:::init_analysis_constants(environment())

  MetaboShiny:::init_analysis_observers(environment())

  # ==== LOAD LOGIN UI ====

  # init all observer
  for (fp in list.files("reactive", full.names = T)) {
    source(fp, local = T)
  }

  MetaboShiny:::init_quit_observers(environment())

  MetaboShiny:::init_zoom_observers(environment())

  MetaboShiny:::init_update_check(environment())

  MetaboShiny:::init_fancy_observer(environment())

  onStop(function() {
    print("Closing MetaboShiny...")
    if (!is.null(session_cl)) {
      parallel::stopCluster(session_cl)
    }
    session_cl <<- NULL
    gc()
    rmv <- list.files(".", pattern = ".csv|.log", full.names = T)
    if (all(file.remove(rmv))) NULL
  })
}
