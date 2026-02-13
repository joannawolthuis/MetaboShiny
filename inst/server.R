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

  # FIX FOR OPTIONS (robust parsing + allow empty values like `apikey =`)
  # The app calls `MetaboShiny::getOptions()`/`MetaboShiny::setOption()` from the
  # installed package; patch them in-place so `runApp()` works without reinstall.
  try({
    metshi_ns <- asNamespace("MetaboShiny")

    replace_ns_fun <- function(ns, name, fun) {
      if (!exists(name, envir = ns, inherits = FALSE)) {
        return(invisible(NULL))
      }
      was_locked <- bindingIsLocked(name, ns)
      if (was_locked) {
        unlockBinding(name, ns)
      }
      assign(name, fun, envir = ns)
      if (was_locked) {
        lockBinding(name, ns)
      }
      invisible(NULL)
    }

    getOptions_safe <- function(file.loc) {
      lines <- readLines(file.loc, warn = FALSE)
      options <- list()

      for (line in lines) {
        if (!nzchar(trimws(line))) {
          next
        }

        eq_pos <- regexpr("=", line, fixed = TRUE)[[1]]
        if (eq_pos < 1) {
          next
        }

        key <- trimws(substr(line, 1, eq_pos - 1))
        value <- trimws(substr(line, eq_pos + 1, nchar(line)))

        if (nzchar(key)) {
          options[[key]] <- value
        }
      }

      options
    }

    setOption_safe <- function(file.loc, key, value) {
      options <- getOptions_safe(file.loc)
      options[[key]] <- value

      out_lines <- vapply(names(options), FUN.VALUE = character(1), FUN = function(option_name) {
        option_value <- options[[option_name]]
        if (is.null(option_value) || length(option_value) == 0) {
          option_value <- ""
        }
        paste0(option_name, " = ", option_value)
      })

      writeLines(text = out_lines, con = file.loc)
      invisible(NULL)
    }

    replace_ns_fun(metshi_ns, "getOptions", getOptions_safe)
    replace_ns_fun(metshi_ns, "setOption", setOption_safe)
  })

  # DEV-ONLY: override installed MetaboShiny functions with local `R/` sources.
  # Opt-in via either:
  # - `Sys.setenv(METABOSHINY_USE_LOCAL = "1")` OR
  # - `options(metaboshiny.use_local = TRUE)`
  #
  # This is intended for quick testing without reinstalling the package.
  try({
    use_local <- isTRUE(getOption("metaboshiny.use_local", FALSE)) ||
      identical(Sys.getenv("METABOSHINY_USE_LOCAL"), "1")

    if (use_local) {
      metshi_ns <- asNamespace("MetaboShiny")

      replace_ns_fun_safe <- function(ns, name, fun) {
        if (exists("replace_ns_fun", inherits = TRUE)) {
          replace_ns_fun(ns, name, fun)
        } else {
          was_locked <- bindingIsLocked(name, ns)
          if (was_locked) {
            unlockBinding(name, ns)
          }
          assign(name, fun, envir = ns)
          if (was_locked) {
            lockBinding(name, ns)
          }
          invisible(NULL)
        }
      }

      r_dir_candidates <- c("R", file.path("..", "R"), file.path("..", "..", "R"))
      r_dir <- r_dir_candidates[dir.exists(r_dir_candidates)][1]

      if (is.na(r_dir) || !nzchar(r_dir)) {
        message("DEV mode: requested local overrides, but no `R/` directory found from wd: ", getwd())
        invisible(NULL)
      } else {
        message("DEV mode: loading local overrides from ", normalizePath(r_dir), " (wd: ", getwd(), ").")

        # Use `globalenv()` as parent so sourced code can resolve functions from the
        # normal search path (e.g. `utils::globalVariables` via `package:utils`).
        local_env <- new.env(parent = globalenv())
        r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
        r_files <- r_files[!grepl(paste0("\\", .Platform$file.sep, "Rserve\\", .Platform$file.sep), r_files)]

        for (fp in r_files) {
          sys.source(fp, envir = local_env)
        }

        overridden <- character(0)
        added <- character(0)
        for (nm in ls(local_env, all.names = TRUE)) {
          obj <- get(nm, envir = local_env, inherits = FALSE)
          if (!is.function(obj)) {
            next
          }

          if (exists(nm, envir = metshi_ns, inherits = FALSE)) {
            replace_ns_fun_safe(metshi_ns, nm, obj)
            overridden <- c(overridden, nm)
          } else {
            assign(nm, obj, envir = metshi_ns)
            added <- c(added, nm)
          }
        }

        overridden <- unique(overridden)
        added <- unique(added)
        message(
          "DEV mode: overridden ", length(overridden),
          " + added ", length(added),
          " functions from local `", r_dir, "`."
        )
        try(
          shiny::showNotification(
            paste0(
              "DEV mode: overridden ", length(overridden),
              " + added ", length(added),
              " functions from local `", r_dir, "`"
            ),
            type = "warning"
          ),
          silent = TRUE
        )
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

  diag_plots <- isTRUE(getOption("metaboshiny.diag_plots", FALSE)) ||
    identical(Sys.getenv("METABOSHINY_DIAG_PLOTS"), "1")
  if (diag_plots) {
    message("Diagnostic plot mode enabled (fixed heights, no resizable panel, no conditional plot containers).")
    try(shiny::showNotification("Diagnostic plot mode enabled", type = "message"), silent = TRUE)
  }

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

  enable_showtext <- function(enable = TRUE) {
    if (!requireNamespace("showtext", quietly = TRUE)) {
      return(invisible(FALSE))
    }

    ok <- TRUE
    tryCatch(
      {
        showtext::showtext_auto(enable = enable)
      },
      error = function(e) {
        ok <<- FALSE
        try(
          shiny::showNotification(
            paste("showtext disabled:", conditionMessage(e)),
            type = "warning"
          ),
          silent = TRUE
        )
      }
    )

    invisible(ok)
  }

  enable_showtext(enable = TRUE)
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

  enable_showtext(enable = TRUE) ## Automatically use showtext to render text for future devices

  # this toggles when 'interface' values change (for example from 'bivar' to 'multivar' etc.)
  MetaboShiny:::init_interface_tabs_observer(environment())

  MetaboShiny:::init_core_observers(environment())

  render_combi_picker("combi_anal1", "combi_anal1_picker", "A result column:")
  render_combi_picker("combi_anal2", "combi_anal2_picker", "B result column:")

  MetaboShiny:::init_debug_observer(environment())

  MetaboShiny:::init_export_observer(environment())

  # Re-enable resizable main panel and trigger widget resize during interaction.
  shinyjs::runjs('
    try {
      if ($("#mainPanel").hasClass("ui-resizable")) {
        $("#mainPanel").resizable("destroy");
      }
      $("#mainPanel .ui-resizable-handle").remove();
      $("#mainPanel").removeClass("ui-resizable ui-resizable-autohide");
      $("#mainPanel").resizable({
        handles: "se",
        minHeight: 500,
        minWidth: 700,
        start: function() {
          try {
            var extraH = 16;
            var extraW = 16;
            var contentMinH = Math.ceil(this.scrollHeight + extraH);
            var contentMinW = Math.ceil(this.scrollWidth + extraW);
            $(this).resizable("option", "minHeight", Math.max(500, contentMinH));
            $(this).resizable("option", "minWidth", Math.max(700, contentMinW));
          } catch (e) {}
        },
        resize: function() {
          try { $(window).trigger("resize"); } catch (e) {}
        },
        stop: function() {
          try { $(window).trigger("resize"); } catch (e) {}
        }
      });
    } catch (e) {}
  ')

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
