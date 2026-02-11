ui_tab_database <- function(gbl, adducts) {
  shiny::tabPanel("database",
    icon = shiny::icon("database",
      class = "outlined"
    ),
    value = "database",
    shiny::fluidRow(
      align = "center",
      # ok
      shinyBS::bsCollapse(shinyBS::bsCollapsePanel(
        title = shiny::h2("Settings"), value = "panel4",
        shiny::sliderInput(
          inputId = "db_mz_range", label = "What mass range can your machine detect?",
          min = 0, max = 3000, value = c(60, 600),
          step = 1, post = "m/z", dragRange = TRUE
        ),
        shiny::tags$i("Warning: increasing the upper m/z boundary may drastically increase database build times
                                                                                     - calculating isotopes is more time-consuming for larger molecules!"),
        shiny::radioButtons("db_build_mode",
          label = "Build base database, extended (isotopes+adducts) or both?",
          choices = c("base", "extended", "both"), selected = "both"
        ),
        shiny::tags$i("For example: if you only defined a new adduct, pick 'extended' as the source database doesn't change."),
        shiny::helpText("Include minor isotopes?"),
        shinyWidgets::switchInput(
          inputId = "db_all_iso",
          size = "mini",
          onLabel = "Yes",
          offLabel = "No",
          value = T
        ),
        shiny::helpText("Include heavy C, N, O count info for each isotope?"),
        shinyWidgets::switchInput(
          inputId = "db_count_iso",
          size = "mini",
          onLabel = "Yes",
          offLabel = "No",
          value = F
        ),
        shiny::helpText("Queue multiple DB builds?"),
        shinyWidgets::switchInput(
          inputId = "db_build_multi",
          size = "mini",
          onLabel = "Yes",
          offLabel = "No",
          value = F
        ),
        shiny::conditionalPanel(
          "input.db_build_multi == true",
          shiny::actionButton("db_build_sel_all",
            label = "Select all"
          ),
          shiny::actionButton("db_build_multi_all",
            label = "Build selected DBs",
            icon = shiny::icon("wrench")
          )
        )
      ))
    ),
    shiny::uiOutput("db_build_ui")
  )
}

ui_tab_data_import <- function(gbl, adducts) {
  shiny::tabPanel("data import",
    icon = shiny::icon("upload", class = "outlined"),
    # === NEW ===
    shiny::fluidRow(
      align = "center",
      shiny::column(
        4,
        shiny::icon("signature", class = "fa-8x"),
        shiny::textInput("proj_name_new", label = "STEP 1: What is your project name?", value = ""),
        shiny::icon("object-group", class = "fa-8x"),
        br(),
        br(),
        shiny::tags$b("STEP 2: Click buttons to select files"),
        br(),
        br(),
        div(shiny::uiOutput("wipe_regex_ui"), style = "font-size:70%"),
        shiny::fileInput(inputId = "metadata", "Select metadata", buttonLabel = "Browse", accept = c(".csv", ".tsv")),
        # shinyFiles::shinyFilesButton('metadata', 'metadata', 'Select metadata', FALSE),
        br(),
        shinyWidgets::checkboxGroupButtons(
          inputId = "ms_modes",
          label = "Mass spec modes:",
          choices = c(
            `<b class='fa fa-plus'> positive</b>` = "pos",
            `<b class='fa fa-minus'> negative</b>` = "neg"
          ),
          selected = c("pos", "neg"),
          justified = F
        ),
        shiny::uiOutput("outlist_pickers"),
        shiny::div(shiny::imageOutput("proj_merge_check"), style = "height:70px")
      ),
      shiny::column(
        4,
        shiny::icon("bullseye", class = "fa-8x"),
        shiny::numericInput("ppm", "STEP 3: What level accuracy does your mass spectrometer have?", min = 0.01, max = 50, value = 3),
        shiny::helpText("Round m/z values based on ppm error?"),
        shinyWidgets::switchInput(
          inputId = "roundMz",
          size = "mini",
          onLabel = "Yes",
          offLabel = "No",
          value = TRUE
        ),
        shiny::icon("star-half-alt", class = "fa-8x"),
        shinyWidgets::sliderTextInput("perc_limit_samp", "STEP 4: What percentage missing m/z for a sample is allowed?",
          choices = c(0, 0.0001, 0.001, 0.01, 0.1, seq(1, 100, 1)),
          selected = 100, grid = T
        )
      ),
      shiny::column(
        4,
        shiny::icon("user-check", class = "fa-8x"),
        shinyWidgets::sliderTextInput("perc_limit_mz", "STEP 5: What percentage missing samples for a m/z value is allowed?",
          choices = c(0, 0.0001, 0.001, 0.01, 0.1, seq(1, 100, 1)),
          selected = 80, grid = T
        ),
        shinyWidgets::actionBttn("checkMiss", label = "Check distribution (recommended)", size = "xs"),
        shiny::plotOutput("missMzPlot"),
        shiny::uiOutput("missMzRatio")
      )
    ),
    shiny::hr(),
    shiny::fluidRow(
      align = "center",
      shiny::tags$b("STEP 6: Convert to input-ready format"),
      shiny::br(), shiny::br(),
      shinyWidgets::circleButton("create_csv", icon = shiny::icon("long-arrow-alt-right", class = "fa-2x"), size = "lg")
    ),
    shiny::fluidRow(
      align = "center",
      shiny::imageOutput("proj_csv_check", inline = T), shiny::br(), shiny::br(),
      shiny::tags$b("STEP 7: If "), shiny::icon("check-circle"), shiny::tags$b(" continue to normalization"),
    )
  )
}

ui_tab_normalize <- function(gbl, adducts) {
  tabPanel("normalize",
    icon = shiny::icon("shower", class = "outlined"), value = "filter",
    shiny::fluidRow(align = "center", shiny::column(
      3,
      shiny::selectizeInput(width = "80%", "samp_var", "Which variable represents sample amount/concentration?", choices = c(" ")), # TODO: only show this when normalize by sample specific factor (specnorm) is selected
      shiny::selectizeInput("batch_var", "What are your batch variables?", choices = c("batch"), multiple = TRUE, options = list(maxItems = 2L)),
      shiny::selectizeInput(
        inputId = "batch_method_a",
        label = "Batch correction method if machine batch and injection order are present:",
        choices = list(
          "WaveICA - WaveICA" = "waveica",
          "BatchCorrMetabolomics - batchCorr" = "batchCorr",
          "limma - removeBatchEffect" = "limma",
          "ComBat" = "combat",
          "CovBat" = "covbat"
        ),
        selected = "waveica"
      ),
      shiny::selectizeInput(
        inputId = "batch_method_b",
        label = "Batch correction method for other variables",
        choices = list(
          "limma - removeBatchEffect" = "limma",
          "ComBat" = "combat",
          "CovBat" = "covbat"
        ),
        selected = "combat"
      ),
      shiny::conditionalPanel(
        "input.batch_method_b == 'combat'",
        shiny::selectizeInput("batch_covar", "What variable(s) covariate(s) with batch?",
          choices = c(), multiple = TRUE
        ),
      ),
      shiny::actionButton("check_csv", "Get options", icon = shiny::icon("sync")),
      shiny::hr(),
      shiny::selectizeInput(
        width = "80%", "filt_type", "How will you filter your m/z values?", choices = list(
          "Interquantile range" = "iqr",
          "Relative stdev" = "rsd",
          "Non-parametric relative stdev" = "nrsd",
          "Mean" = "mean",
          "Standard deviation" = "sd",
          "Median absolute deviation" = "mad",
          "Median" = "median",
          "None" = "none"
        ),
        selected = "none"
      ),
      shiny::conditionalPanel(
        "input.filt_type != 'none'",
        shiny::numericInput(
          inputId = "maxMz",
          label = "Max. m/z to keep",
          value = 5000
        )
      ),
      shiny::selectizeInput(
        width = "80%", "norm_type", "What type of normalization do you want to do?", choices = list(
          "Quantile normalization" = "QuantileNorm",
          "QCs per batch normalization" = "QcNorm",
          "By reference feature" = "ProbNorm",
          "By reference compound" = "CompNorm",
          "By sample specific factor" = "SpecNorm",
          "Sum" = "SumNorm",
          "Median" = "MedianNorm",
          "None" = "NULL"
        ),
        selected = "QuantileNorm"
      ),
      # ref 	289.18874-
      shiny::conditionalPanel(
        "input.norm_type == 'ProbNorm'",
        div(shinyWidgets::pickerInput(
          inputId = "ref_mz",
          label = div(icon("search"), style = "font-size: xx-large;margin-top: -30px;color: black;-webkit-text-fill-color: white;-webkit-text-stroke-width: 1.5px;-webkit-text-stroke-color: #DFDCDC"),
          choices = "click below button",
          choicesOpt = list(
            subtext = "",
            style = "text-align:center;"
          ),
          options = list(
            `live-search` = TRUE,
            size = 10
          )
        ), class = "mzpicker"),
        shiny::actionButton("check_ref_mzs", "Get options", icon = shiny::icon("sync"))
      ),
      shiny::selectizeInput(
        width = "80%", "trans_type", "How will you transform your data?", choices = list(
          "Log transform" = "LogNorm",
          "Cubic root transform" = "CrNorm",
          "None" = "NULL"
        ),
        selected = "LogNorm"
      ),
      shiny::selectizeInput(
        width = "80%", "scale_type", "How will you scale your data?", choices = list(
          "Autoscale/Z-transform" = "AutoNorm",
          "Mean-center" = "MeanCenter",
          "Pareto Scaling" = "ParetoNorm",
          "Range scaling" = "RangeNorm",
          "None" = "NULL"
        ),
        selected = "AutoNorm"
      ),
      shinyWidgets::sliderTextInput("miss_perc_2", "What percentage missing samples for a m/z is allowed?",
        choices = c(0, 0.0001, 0.001, 0.01, 0.1, seq(1, 100, 1)),
        selected = 20, grid = T
      ),
      shinyWidgets::sliderTextInput("miss_perc_samp", "What percentage missing m/z for a sample is allowed?",
        choices = c(0, 0.0001, 0.001, 0.01, 0.1, seq(1, 100, 1)),
        selected = 100, grid = T
      ),
      shiny::selectizeInput(
        width = "80%", "miss_type", "How to deal with missing values?", choices = list(
          "Half feature minimum" = "colmin",
          "Half sample minimum" = "rowmin",
          "Total minimum" = "min",
          "Random forest" = "rf",
          # "Impute w/ regression" = "regr",
          "KNN imputation" = "knn",
          "SVD imputation" = "svdImpute",
          "BPCA imputation" = "bpca",
          "PPCA imputation" = "ppca",
          "Median" = "median",
          "Mean" = "mean",
          "Leave them out" = "exclude",
          "Leave them alone" = "none"
        ),
        selected = "rowmin"
      ),
      shinyWidgets::switchInput(
        inputId = "miss_minority_filter",
        value = FALSE,
        label = "Only check missing values on the minority class?",
        onLabel = "yes",
        offLabel = "no",
        size = "small" # "<div class=\"fa-flip-vertical\"><i class=\"fas fa-chart-bar fa-rotate-90\"></i></div>"
      ),
      shiny::conditionalPanel(
        "input.miss_type == 'rf'",
        shinyWidgets::switchInput("rf_norm_method",
          label = "Method:",
          value = TRUE, onLabel = "missRanger",
          offLabel = "missForest", size = "small"
        ),
        shiny::sliderInput("rf_norm_ntree", label = "Trees built per variable", value = 10, min = 1, max = 50, step = 1),
        # numericInput("rf_norm_mtry", label = "Trees built per variable", value = 10, min = 1, max = 50)
        shiny::radioButtons("rf_norm_parallel",
          label = "Parallelize?", choices = list(
            "no",
            "forests",
            "variables"
          ),
          selected = "variables"
        )
      ),
      shinyWidgets::switchInput(
        inputId = "redo_upon_change",
        value = F,
        label = "Renormalize after subsetting?:",
        onLabel = "yes",
        offLabel = "no",
        size = "small" # "<div class=\"fa-flip-vertical\"><i class=\"fas fa-chart-bar fa-rotate-90\"></i></div>"
      ),
      shinyWidgets::switchInput(
        inputId = "miss_upon_subset",
        value = F,
        label = "Refilter by missing values after subsetting?:",
        onLabel = "yes",
        offLabel = "no",
        size = "small" # "<div class=\"fa-flip-vertical\"><i class=\"fas fa-chart-bar fa-rotate-90\"></i></div>"
      ),
      # - - - - - -
      shinyWidgets::switchInput(
        inputId = "repl_merge",
        value = FALSE,
        label = "Merge tech. replicates?",
        onLabel = "yes",
        offLabel = "no",
        size = "small" # "<div class=\"fa-flip-vertical\"><i class=\"fas fa-chart-bar fa-rotate-90\"></i></div>"
      ),
      shinyWidgets::switchInput(
        inputId = "miss_group_replicates",
        value = FALSE,
        label = "Group missing values by replicates?",
        onLabel = "yes",
        offLabel = "no",
        size = "small" # "<div class=\"fa-flip-vertical\"><i class=\"fas fa-chart-bar fa-rotate-90\"></i></div>"
      ),
      shiny::conditionalPanel(
        "input.repl_merge == true",
        shiny::selectInput("repl_merge_fun",
          "How to merge m/z intensity over replicates?",
          choices = c(
            "mean",
            "sum",
            "min",
            "max",
            "median"
          )
        )
      ),
      # - - - - - -
      shinyWidgets::switchInput(
        inputId = "pca_corr",
        value = FALSE,
        label = "PCA correction",
        onLabel = "yes",
        offLabel = "no",
        size = "small" # "<div class=\"fa-flip-vertical\"><i class=\"fas fa-chart-bar fa-rotate-90\"></i></div>"
      ),
      shiny::conditionalPanel(
        "input.pca_corr == true",
        sliderInput("keep_pcs",
          label = h3("Keep variance in PCs:"), min = 0,
          max = 300, value = c(10, 50)
        )
      ),
      MetaboShiny::switchButton(
        inputId = "remove_outliers",
        label = "Exclude outliers?",
        value = FALSE, col = "BW", type = "YN"
      ),
      shiny::actionButton("initialize", "Go", icon = shiny::icon("hand-point-right")),
      br(), br(), br(), br()
    ), shiny::column(
      9,
      # show the summary plots post-normalization
      shiny::tabsetPanel(
        tabPanel(
          "m/z values", # icon=shiny::icon("braille"),
          shiny::fluidRow(
            shiny::column(
              6,
              shiny::imageOutput("empty2", width = "100%", height = "1px"),
              shinyjqui::jqui_resizable(shiny::uiOutput("var1_wrap"))
            ),
            shiny::column(
              6,
              shinyjqui::jqui_resizable(shiny::uiOutput("var3_wrap"))
            )
          ),
          shiny::fluidRow(
            shiny::column(
              6,
              shinyjqui::jqui_resizable(shiny::uiOutput("var2_wrap"))
            ),
            shiny::column(
              6,
              shinyjqui::jqui_resizable(shiny::uiOutput("var4_wrap"))
            )
          )
        ),
        shiny::tabPanel(
          "samples", # icon=shiny::icon("tint"),
          shiny::fluidRow(
            shiny::column(
              6,
              shinyjqui::jqui_resizable(shiny::uiOutput("samp1_wrap"))
            ),
            shiny::column(
              6,
              shinyjqui::jqui_resizable(shiny::uiOutput("samp3_wrap"))
            )
          ),
          shiny::fluidRow(
            shiny::column(
              6,
              shinyjqui::jqui_resizable(shiny::uiOutput("samp2_wrap"))
            ),
            shiny::column(
              6,
              shinyjqui::jqui_resizable(shiny::uiOutput("samp4_wrap"))
            )
          )
        )
      )
    ))
  )
}

ui_tab_prematch <- function(gbl, adducts) {
  shiny::tabPanel("prematch",
    icon = shiny::icon("search", class = "outlined"), value = "prematch",
    # - - - pre-matching part - - -
    shiny::fluidRow(
      align = "center",
      MetaboShiny::switchButton(
        inputId = "do_prematch",
        label = "Do matching beforehand?",
        col = "BW",
        type = "YN"
      ),
      shiny::tags$i("All m/z values will be searched in the databases of choice and the results will be saved to your save file for fast access."), shiny::br(),
      shiny::tags$i("Search results can still be overridden by manual searching. Don't forget to save after!")
    ),
    shiny::br(),
    shiny::fluidRow(
      align = "center",
      shiny::column(2),
      shiny::column(
        8, shiny::conditionalPanel(
          "input.do_prematch == true",
          shiny::h2("Included databases:"),
          shiny::uiOutput("db_prematch_select"),
          shinyWidgets::circleButton("select_db_prematch_all",
            icon = shiny::icon("shopping-cart"),
            size = "default"
          )
        ),
        shiny::hr(),
        shiny::fluidRow(
          shiny::column(
            6, shiny::h2("Find matches"),
            shinyWidgets::circleButton(
              inputId = "prematch",
              icon = shiny::icon("searchengin")
            )
          ),
          shiny::column(
            6, shiny::h2("Clear matches"),
            shinyWidgets::circleButton(
              inputId = "clear_prematch",
              icon = shiny::icon("trash")
            )
          )
        )
      ),
      shiny::column(2)
    )
  )
  # this tab is the main analysis tab. all tabs for all analyses are listed here, but the visibility is changed depending on the current experiment
}

ui_tab_analyse <- function(gbl, adducts) {
  shiny::tabPanel("analyse",
    icon = shiny::icon("chart-bar", class = "outlined"), value = "analysis",
    div(id = "panelContainer", sidebarLayout(
      position = "right",
      mainPanel = # shinyjqui::jqui_resizable(
        mainPanel(
          width = 9, id = "mainPanel",
          shiny::fluidRow(
            align = "center",
            shiny::imageOutput("empty3", width = "100%", height = "1px"),
            # nokay
            shinyBS::bsCollapse(
              id = "collapse_summary", multiple = F,
              shinyBS::bsCollapsePanel(
                title = tags$b(
                  icon("caret-down"),
                  "show/hide selected m/z abundance",
                  icon("caret-down")
                ),
                style = "info",
                shinyjqui::jqui_resizable(
                  shiny::uiOutput("summary_plot_wrap")
                )
              )
            )
          ),
          shiny::tabsetPanel(
            id = "statistics", selected = "inf",
            # TODO: T-SNE
            # this tab shows general information, mostly a message with 'please give me some data' :-)
            shiny::tabPanel("start",
              icon = shiny::icon("star", "fa-2x"), value = "inf",
              shiny::fluidRow(align = "center", shiny::column(
                width = 12,
                shiny::br(), shiny::br(), shiny::br(), shiny::br(),
                # shiny::hr(),
                # shiny::icon("arrow-right","fa-lg"), shiny::icon("arrow-right","fa-lg"), shiny::icon("arrow-right","fa-lg"),
                shiny::br(), shiny::br(),
                shiny::h2("Please select a variable of interest in the sidebar!"), shiny::br(),
                shiny::icon("exchange-alt", "fa-4x"),
                shiny::h2("Alternatively, load an existing dataset from the bottom toolbar."), shiny::br(),
                shiny::icon("folder-open", "fa-4x"),
                shiny::br(), shiny::br(), shiny::br()
                # shiny::hr()
                # shiny::icon("arrow-right","fa-lg"), shiny::icon("arrow-right","fa-lg"), shiny::icon("arrow-right","fa-lg")
              ))
            ),
            shiny::navbarMenu("dimension reduction",
              icon = icon("cube", "fa-2x"), menuName = "dimred",
              MetaboShiny:::ui_pca(),
              MetaboShiny:::ui_plsda(),
              MetaboShiny:::ui_tsne(),
              MetaboShiny:::ui_ica(),
              MetaboShiny:::ui_umap()
            ),
            shiny::navbarMenu("per m/z",
              icon = shiny::icon("fingerprint", "fa-2x"), menuName = "permz",
              MetaboShiny:::ui_t_test(),
              MetaboShiny:::ui_proda(),
              MetaboShiny:::ui_anova(),
              MetaboShiny:::ui_fold_change(),
              MetaboShiny:::ui_cliffs_delta(),
              MetaboShiny:::ui_meba(),
              MetaboShiny:::ui_asca(),
              MetaboShiny:::ui_pattern()
            ),
            shiny::navbarMenu("overview analyses",
              icon = shiny::icon("globe", "fa-2x"), menuName = "overview",
              MetaboShiny:::ui_intersection_plot(),
              MetaboShiny:::ui_heatmap(),
              MetaboShiny:::ui_network(),
              MetaboShiny:::ui_multirank(),
              MetaboShiny:::ui_venn(),
              MetaboShiny:::ui_power(),
              MetaboShiny:::ui_enrichment(),
              MetaboShiny:::ui_feature_selection()
            ),
            shiny::tabPanel("machine learning",
              value = "ml", icon = shiny::icon("signature", "fa-2x"),
              shiny::br(),
              shiny::tabsetPanel(
                id = "ml2",
                MetaboShiny:::ui_ml_configure(),
                MetaboShiny:::ui_ml_results()
              )
            ),
          ),
        ),
        # this is the sidebar that shows in the analysis tab. contains a lot of settings on the current variable of interest, plot themes and colours, and venn diagrams.
        sidebarPanel =
          shinyjqui::jqui_resizable(
            sidebarPanel(
              align = "center", width = 3, id = "sidePanel",
              shiny::fluidRow(
                align = "center",
                  shiny::tabsetPanel(
                    id = "anal_sidebar", selected = "start",
                    shiny::tabPanel(
                      title = shiny::icon("star"), value = "start",
                      shiny::br(), shiny::br(), shiny::br(), shiny::br(), shiny::br(),
                      tags$i("Load a dataset to get started!"),
                      br(),
                      icon("caret-down", "fa-3x"),
                      br(),
                      icon("box-open", "fa-5x"),
                      shiny::br(), shiny::br(), shiny::br(),
                      shiny::br(), shiny::br()
                    ),
                    shiny::tabPanel(
                      value = "switch",
                      title = shinyBS::tipify(shiny::icon("exchange-alt"),
                        title = "switch to different variable and/or subset based on variable",
                        trigger = "hover", options = list(container = "body")
                      ),
                      shiny::h2("Current experiment:"),
                      shiny::div(
                        shiny::div(style = "display: inline-block;vertical-align:top;", shiny::uiOutput("curr_name")),
                        style = "background-color:white;
                              width:100%;
                              margin: 0 auto;
                              border-top: 1px solid #DFDCDC;
                              border-bottom: 1px solid #DFDCDC;"
                      ),
                      shiny::hr(),
                      shiny::h2("Change variable of interest"),
                      shiny::selectizeInput("stats_type",
                        label = "Normal, two-factor or time series?",
                        choices = list(
                          "one factor" = "1f",
                          "two factors" = "2f",
                          "time series" = "t",
                          "time series + one factor" = "t1f"
                        ),
                        selected = "1f", width = "80%"
                      ),
                      fluidRow(align = "center", shiny::uiOutput("stats_picker")),
                      shiny::conditionalPanel(
                        "input.stats_type == '1f'",
                        shiny::checkboxInput("paired",
                          "Matched samples only?",
                          value = F
                        )
                      ),
                      shiny::conditionalPanel(
                        "input.stats_type == 't' || input.stats_type == 't1f'",
                        shiny::selectizeInput("time_var",
                          label = "Time series variable:",
                          choices = c(" "), width = "80%"
                        )
                      ),
                      shinyWidgets::actionBttn(
                        inputId = "change_cls",
                        label = "do stats on selected",
                        style = "bordered",
                        icon = icon("magic"),
                        size = "sm"
                      )
                    ),
                    shiny::tabPanel(
                      value = "subset",
                      title = shinyBS::tipify(shiny::icon("filter"),
                        title = "Subset based on variable or m/z(s)",
                        trigger = "hover", options = list(container = "body")
                      ),
                      shiny::hr(),
                      shiny::h2("Subset data"),
                      helpText("Current sample count:"),
                      h2(shiny::textOutput("samp_count")),
                      helpText("Current m/z count:"),
                      h2(shiny::textOutput("mz_count")),
                      br(),
                      shiny::selectizeInput("subset_var", label = "Subset data based on:", choices = c(" "), width = "80%"),
                      shinyWidgets::pickerInput(
                        inputId = "subset_group",
                        label = "Group(s) in subset:",
                        choices = c(" "),
                        choicesOpt = list(
                          subtext = c(" ")
                        ), multiple = T
                      ),
                      shinyWidgets::actionBttn(
                        inputId = "change_subset",
                        label = "click to subset by metadata",
                        style = "bordered",
                        icon = icon("filter"),
                        size = "sm"
                      ), br(), br(),
                      shinyWidgets::pickerInput(
                        inputId = "storage_choice",
                        label = "Load existing meta-dataset:",
                        choices = c(" "),
                        choicesOpt = list(
                          subtext = c(" ")
                        )
                      ),
                      shiny::span(
                        shinyWidgets::actionBttn(
                          inputId = "load_storage",
                          label = "load",
                          style = "bordered",
                          icon = icon("folder-open"),
                          size = "sm",
                          color = "success"
                        ),
                        shinyWidgets::actionBttn(
                          inputId = "remove_storage",
                          label = "delete",
                          style = "bordered",
                          icon = icon("eraser"),
                          size = "sm",
                          color = "warning"
                        )
                      ),
                      br(), br(),
                      helpText("Reset to non-subsetted or paired dataset"),
                      shinyWidgets::actionBttn(
                        inputId = "reset_subset",
                        label = "click to un-subset data",
                        style = "bordered",
                        icon = icon("undo"),
                        size = "sm"
                      )
                    ),
                    shiny::tabPanel(
                      value = "search", title = shinyBS::tipify(shiny::icon("search"),
                        title = "search m/z > compound or compound > mz",
                        trigger = "hover", options = list(container = "body")
                      ),
                      shiny::br(),
                      shinyBS::bsCollapse(shinyBS::bsCollapsePanel(
                        title = shiny::h2("Options"),
                        value = "search_panel",
                        shiny::tabsetPanel(
                          id = "tab_iden_1", selected = "pick_databases",
                          # forward searching
                          shiny::tabPanel(
                            title = shiny::icon("database"),
                            # title="DBs",
                            value = "pick_databases",
                            shiny::uiOutput("db_search_select"),
                            br(),
                            shiny::div(
                              id = "curly-brace", shiny::div(id = "left", class = "brace"),
                              shiny::div(id = "right", class = "brace")
                            ),
                            br(), br(),
                            shinyWidgets::circleButton("select_db_all",
                              icon = shiny::icon("shopping-cart"),
                              size = "default"
                            ) # shiny::icon("fingerprint"), size = "sm")
                          ),
                          shiny::tabPanel(
                            title = icon("medal"),
                            shiny::fluidRow(
                              shiny::helpText("Only consider highly correlated m/z values?"),
                              shinyWidgets::switchInput(
                                inputId = "score_use_corr",
                                size = "mini",
                                onLabel = "Yes",
                                offLabel = "No",
                                value = TRUE
                              ),
                              shiny::conditionalPanel(
                                "input.score_use_corr == true",
                                shiny::selectizeInput(
                                  width = "80%", "score_corr_method",
                                  "Which correlation method used?",
                                  selected = "pearson",
                                  choices = list(
                                    "pearson",
                                    "spearman",
                                    "kendall"
                                  )
                                ),
                                shiny::sliderInput("score_corr_min",
                                  label = "Min required correlation: ",
                                  min = 0, max = 1,
                                  value = 0.9,
                                  step = 0.05
                                )
                              )
                            ),
                            shiny::tabsetPanel(
                              shiny::tabPanel(
                                title = icon("percentage"), # title="isoScore",
                                shiny::fluidRow(
                                  align = "center",
                                  shiny::br(),
                                  shiny::h2("isoMatcher"),
                                  shiny::h3("Score matches by isotope pattern"),
                                  shiny::helpText("Per match, investigates if other same-formula isotopes are present."),
                                  shiny::helpText("Score pattern intensity simularity? (if not, only check for isotope presence)"),
                                  shinyWidgets::switchInput(
                                    inputId = "iso_use_int",
                                    size = "mini",
                                    onLabel = "Yes",
                                    offLabel = "No",
                                    value = TRUE
                                  ),
                                  shiny::conditionalPanel(
                                    "input.iso_use_int == true",
                                    shiny::selectizeInput(
                                      width = "80%", "iso_score_method",
                                      "Which method used to score compounds of same weight?",
                                      selected = "mscore",
                                      choices = list(
                                        "M-score" = "mscore"
                                        # "Chi-square" = "chisq",
                                        # "Mean absolute percentage error" = "mape",
                                        # "SIRIUS" = "sirius",
                                        # "Network-based" = "network"
                                      )
                                    ),
                                    shiny::sliderInput("int_prec",
                                      label = "Intensity error margin",
                                      min = 1, max = 100,
                                      value = 2, post = "%"
                                    )
                                  ),
                                  shiny::uiOutput("iso_rt_ui"),
                                  shinyWidgets::circleButton("score_iso", icon = shiny::icon("medal"), size = "sm")
                                )
                              ),
                              shiny::tabPanel(
                                title = shiny::icon("clock"),
                                shiny::fluidRow(
                                  align = "center",
                                  shiny::br(),
                                  shiny::h2("addMatcher"),
                                  shiny::h3("Score matches by fellow adducts"),
                                  shiny::helpText("Per match, investigates if other main isotope adducts are present."),
                                  shinyWidgets::pickerInput(
                                    inputId = "score_adducts",
                                    choices = adducts$Name,
                                    label = "Adducts to look for:",
                                    selected = c("[M+H]1+", "[M+Na]1+", "[M+2H]2+", "[M+K]1+"),
                                    multiple = T,
                                    options = list(
                                      `actions-box` = TRUE
                                    )
                                  ),
                                  shiny::uiOutput("score_adducts_cats"),
                                  shiny::uiOutput("add_rt_ui"),
                                  shinyWidgets::circleButton("score_add",
                                    icon = shiny::icon("medal"),
                                    size = "sm"
                                  )
                                )
                              )
                            )
                          ),
                          shiny::tabPanel(
                            title = shiny::icon("cloud"), # title="wordcloud",
                            tabsetPanel(
                              id = "wordclouds", selected = "settings",
                              tabPanel(
                                title = "settings",
                                shiny::br(),
                                shinyWidgets::switchInput(
                                  inputId = "wordcloud_manual", value = TRUE,
                                  onLabel = "own word", #<i class=\"fas fa-cloud\"></i>",
                                  offLabel = "from matches" # "<div class=\"fa-flip-vertical\"><i class=\"fas fa-chart-bar fa-rotate-90\"></i></div>"
                                ),
                                shiny::tabPanel(
                                  "input.wordcloud_manual == true",
                                  helpText("Type a metabolite below and search PubMed to find documents that contain that word in the text."),
                                  textInput("wordcloud_searchTerm", label = h3("Enter your search terms"), placeholder = "enter your search terms"),
                                  helpText("You can specify the start and end years of your search, use the format YYYY"),
                                  sliderInput(inputId = "wordcloud_dateRange", label = "Date range input:", min = 2000, max = 2025, value = c(2010, 2020), sep = ""),
                                  sliderInput("wordcloud_absFreq", "Amount of abstracts to use:", min = 1, max = 4000, value = 500)
                                ),
                                actionButton("do_wordcloud", "Plot")
                              ),
                              tabPanel(
                                title = "filters",
                                helpText("Type a search term: a word filter will be created with the resulting top words."),
                                textInput("wordcloud_filterTerm", label = h3("Enter your search terms"), placeholder = "enter your search terms"),
                                helpText("You can specify the start and end years of your search, use the format YYYY"),
                                sliderInput(inputId = "wordcloud_filterDateRange", label = "Date range input:", min = 2000, max = 2025, value = c(2010, 2020), sep = ""),
                                sliderInput("wordcloud_filterAbsFreq", "Amount of abstracts to use:", min = 10, max = 500000, value = 5000),
                                sliderInput("wordcloud_filterTopN", "Top words used for filter:", min = 1, max = 5000, value = 500),
                                actionButton("wordcloud_make_filter", "Go! :)")
                              ),
                              tabPanel(
                                title = "plot",
                                shiny::conditionalPanel(
                                  "input.wordbar == true",
                                  wordcloud2::wordcloud2Output(outputId = "wordcloud"),
                                  tags$script(HTML(
                                    "$(document).on('click', '#canvas', function() {",
                                    'word = document.getElementById("wcSpan").innerHTML;',
                                    "Shiny.onInputChange('wordcloud_selected_word', word);",
                                    "});"
                                  ))
                                ),
                                conditionalPanel(
                                  "input.wordbar == false",
                                  shinyjqui::jqui_resizable(shiny::uiOutput("wordbar_wrap"))
                                ),
                                shinyWidgets::switchInput(
                                  inputId = "wordbar", value = TRUE,
                                  onLabel = "cloud", #<i class=\"fas fa-cloud\"></i>",
                                  offLabel = "barchart" # "<div class=\"fa-flip-vertical\"><i class=\"fas fa-chart-bar fa-rotate-90\"></i></div>"
                                ),
                                br(),
                                shiny::sliderInput("wordcloud_topWords", "Top words used for plot:", min = 1, max = 10000, value = 100),
                                shiny::selectizeInput("wordcloud_filter",
                                  "Filter",
                                  choices = c("stopwords"),
                                  selected = "stopwords",
                                  multiple = T
                                ),
                                shiny::wellPanel(
                                  id = "def", style = "overflow-y:scroll; max-height: 400px; font-size: 80%",
                                  helpText("Click a word to see abstracts!"),
                                  shiny::uiOutput("wordcloud_abstracts")
                                )
                              )
                            )
                          )
                        )
                      )),
                      br(),
                      shiny::uiOutput("manual_search"),
                      br(),
                      div(
                        shiny::fluidRow(
                          align = "center",
                          div(
                            br(), br(),
                            div(shiny::icon("paw", "fa-xs fa-rotate-90"), style = "margin-top: 20px; display: flex;"),
                            div(shiny::icon("paw", "fa-xs fa-rotate-90"), style = "display: flex;"),
                            div(shiny::icon("paw", "fa-xs fa-rotate-90"), style = "margin-top: 20px; display: flex;"),
                            div(shinyWidgets::pickerInput(
                              inputId = "curr_mz",
                              label = div(icon("search"), style = "font-size: xx-large;margin-top: -30px;color: black;-webkit-text-fill-color: white;-webkit-text-stroke-width: 1.5px;-webkit-text-stroke-color: #DFDCDC"),
                              choices = "fa-cat",
                              choicesOpt = list(
                                subtext = "",
                                style = "text-align:center;"
                              ),
                              options = list(
                                `live-search` = TRUE,
                                size = 10
                              )
                            ), class = "mzpicker"),
                            div(shiny::icon("paw", "fa-xs fa-rotate-90"), style = "margin-top: 20px; display: flex;"),
                            div(shiny::icon("paw", "fa-xs fa-rotate-90"), style = "display: flex;"),
                            div(shiny::icon("paw", "fa-xs fa-rotate-90"), style = "margin-top: 20px; display: flex;"),
                            br(),
                            style = "align-items: center; display: inline-flex;"
                          ),
                        ),
                        shinyWidgets::prettyToggle(
                          status_off = "default",
                          status_on = "warning",
                          bigger = T,
                          animation = "pulse",
                          inputId = "star_mz",
                          label_on = "",
                          label_off = "",
                          outline = TRUE,
                          plain = TRUE,
                          value = FALSE,
                          icon_on = icon("star"),
                          icon_off = htmltools::browsable(tags$i(class = "far fa-star"))
                        ),
                        style = "border-top: 1px solid #DFDCDC; border-bottom: 1px solid #DFDCDC; margin-left: -15px; margin-right: -15px; background-color:white;"
                      ),
                      br(),
                      shiny::tabsetPanel(
                        id = "tab_search", shiny::tabPanel(
                          title = shinyBS::tipify(
                            shiny::icon("paw"),
                            "mz <> molecule matches"
                          ),
                          shinyBS::bsCollapse(shinyBS::bsCollapsePanel(
                            title = shiny::h2("Compound info"),
                            value = "panel2",
                            shiny::tabsetPanel(
                              id = "tab_iden_3",
                              shiny::tabPanel(
                                title = shiny::icon("atlas"),
                                wellPanel(
                                  id = "def", style = "overflow-y:scroll; max-height: 200px",
                                  shiny::uiOutput("desc_ui")
                                )
                              ), shiny::tabPanel(
                                title = shiny::icon("atom"),
                                br(),
                                shiny::imageOutput("empty4", width = "100%", height = "1px"),
                                shiny::uiOutput("curr_formula"),
                                shiny::plotOutput("curr_struct", inline = T)
                              )
                            )
                          )),
                          shiny::tabsetPanel(
                            id = "tab_iden_2",
                            # forward searching
                            shiny::tabPanel(
                              title = "mz > molecule", value = "mzmol",
                              br(),
                              shiny::tags$i("Found matches for this m/z value:"),
                              br(),
                              shinyWidgets::searchInput(
                                inputId = "match_search_query",
                                label = "Search compound descriptions",
                                placeholder = "",
                                btnSearch = icon("search"),
                                btnReset = icon("eraser"),
                                width = "100%"
                              ),
                              shiny::fluidRow(
                                align = "left",
                                shiny::div(DT::dataTableOutput("match_tab"),
                                  style = "font-size:80%;"
                                )
                              )
                            ),
                            # reverse searching
                            shiny::tabPanel(
                              title = "molecule > mz", value = "molmz",
                              br(),
                              shiny::tags$i("Press the below button to browse compounds of the selected databases."), br(),
                              shiny::actionButton("browse_db", "Browse", icon = shiny::icon("eye")),
                              br(),
                              shiny::div(DT::dataTableOutput("browse_tab"), style = "font-size:80%"),
                              br(),
                              shiny::tags$i("M/z values matching adducts or isotopes of this compound:"), br(),
                              shiny::div(DT::dataTableOutput("hits_tab"), style = "font-size:80%")
                            )
                          ),
                          shiny::fluidRow(
                            shiny::div(shiny::helpText("Show heavy C, N, O count if available?"),
                              style = "font-size:80%;"
                            ),
                            shinyWidgets::switchInput(
                              inputId = "show_iso_labels",
                              size = "mini",
                              onLabel = "show",
                              offLabel = "hide",
                              value = FALSE
                            )
                          ),
                          br(),
                          shiny::helpText("Current match filtering:"),
                          shiny::fluidRow(
                            align = "center",
                            shiny::imageOutput("empty", width = "100%", height = "1px"),
                            shiny::div(style = "display: inline-block;vertical-align:top;", shiny::icon("plus-circle")),
                            shiny::div(style = "display: inline-block;vertical-align:top;", shiny::div(verbatimTextOutput("curr_add"), style = "font-size:60%;padding:2px")),
                            shiny::div(style = "display: inline-block;vertical-align:top;", shiny::icon("percentage")),
                            shiny::div(style = "display: inline-block;vertical-align:top;", shiny::div(verbatimTextOutput("curr_iso"), style = "font-size:60%;padding:2px")),
                            shiny::div(style = "display: inline-block;vertical-align:top;", shiny::icon("database")),
                            shiny::div(style = "display: inline-block;vertical-align:top;", shiny::div(verbatimTextOutput("curr_db"), style = "font-size:60%;padding:2px"))
                          ),
                          shiny::helpText("Undo match filtering"),
                          shinyWidgets::circleButton("undo_match_filt", icon = shiny::icon("undo-alt")),
                          br()
                        ),
                        shiny::tabPanel(
                          title = shinyBS::tipify(shiny::icon("filter"),
                            title = "filter matches on adduct, isotope or database"
                          ),
                          value = "match_filters_tab",
                          shinyWidgets::verticalTabsetPanel(
                            id = "match_filters",
                            contentWidth = 10,
                            menuSide = "right",
                            shinyWidgets::verticalTabPanel(
                              box_height = "80px",
                              title = NULL, icon = icon("plus", "fa-2x"), value = "pie_add",
                              plotly::plotlyOutput("match_pie_add"),
                              shiny::uiOutput("filter_adducts_cats")
                            ),
                            shinyWidgets::verticalTabPanel(
                              box_height = "80px",
                              title = NULL, icon = icon("percentage", "fa-2x"), value = "pie_iso",
                              plotly::plotlyOutput("match_pie_iso")
                            ),
                            shinyWidgets::verticalTabPanel(
                              box_height = "80px",
                              title = NULL, icon = icon("database", "fa-2x"), value = "pie_db",
                              plotly::plotlyOutput("match_pie_db")
                            )
                          ), br()
                        )
                      )
                    ),
                    # this tab is used to select user plot theme and user colours (discrete and continuous)
                    shiny::tabPanel(
                      value = "plot_aes", title = shinyBS::tipify(shiny::icon("paint-brush"),
                        title = "change plot style and colours",
                        trigger = "hover", options = list(container = "body")
                      ),
                      shiny::h2("Plot style"), br(),
                      helpText("Click to reload plot:"),
                      shinyWidgets::circleButton("reload_plots",
                        icon = shiny::icon("paint-brush")
                      ),
                      br(),
                      helpText("Use ggplot or plotly for plots? (ggplot -> faster, but no interactivity)"),
                      shinyWidgets::switchInput(
                        inputId = "ggplotly",
                        onLabel = "plotly",
                        offLabel = "ggplot",
                        offStatus = "warning",
                        onStatus = "info",
                        value = F
                      ),
                      helpText("Export plots as .png or .svg?"),
                      shinyWidgets::switchInput(
                        inputId = "plotsvg",
                        onLabel = "svg",
                        offLabel = "png",
                        offStatus = "info",
                        onStatus = "warning",
                        value = T
                      ),
                      # helpText("Plot resolution (ggplot):"),
                      # shiny::numericInput("plot_resolution", "PPI", min = 5, max = 1200, value = 300),
                      shiny::h2("Legend"), br(),
                      shinyWidgets::switchInput(
                        inputId = "legend",
                        size = "mini",
                        onLabel = "show",
                        offLabel = "hide",
                        value = FALSE
                      ),
                      shiny::h2("M/z labels"),
                      shiny::tags$i("only ggplot mode"),
                      shinyWidgets::switchInput("plot_mzlabels", value = T, onLabel = "show", offLabel = "hide", size = "small"),
                      shiny::numericInput("ggplot_font_size", label = "Plot font base size:", value = 19),
                      shiny::numericInput("ggplot_font_size_spec", label = "Plot other font size:", value = 8),
                      plotly::plotlyOutput("ggplot_font_size_example"),
                      shiny::selectizeInput("ggplot_sum_style",
                        multiple = T, label = "Style(s)", choices = list(
                          "Box" = "box",
                          "Violin" = "violin",
                          "Beeswarm" = "beeswarm",
                          "Scatterplot" = "scatter",
                          "Sina" = "sina"
                        ),
                        selected = c("beeswarm"), width = "80%"
                      ),
                      shiny::selectizeInput("ggplot_sum_stats", label = "Stats shown", choices = list("median", "mean", "none"), width = "80%"),
                      shiny::h2("Shape"),
                      shiny::selectizeInput("shape_var", label = "Marker shape based on:", choices = c(" "), width = "80%"),
                      shiny::h2("Color"),
                      shiny::selectizeInput("fill_var", label = "Marker fill color based on:", choices = c(" "), width = "80%"),
                      shiny::selectizeInput("col_var", label = "Marker outline color based on:", choices = c(" "), width = "80%"),
                      shiny::h2("Hover text"),
                      shiny::selectizeInput("txt_var", label = "Marker hover text based on:", choices = c(" "), width = "80%"),
                      # input$plot_mzlabels
                      shiny::h2("Plot theme"),
                      shiny::selectizeInput("ggplot_theme",
                        label = "Theme", choices = list(
                          "Grid, white bg" = "bw",
                          "No grid, white bg" = "classic",
                          "Grid, gray bg" = "gray",
                          "Minimal" = "min",
                          "Grid, black bg" = "dark",
                          "Grid, white bg, gray axes" = "light",
                          "Line drawing" = "line"
                        ),
                        selected = "min", width = "80%"
                      ),
                      plotly::plotlyOutput("ggplot_theme_example", inline = F),
                      shiny::h2("Continuous data"),
                      # the below options need to match with the corresponding function storage in 'global'. if you want to add more it'll go here!
                      shiny::selectizeInput("color_ramp", label = "Color scheme", choices = list(
                        "RAINBOW!" = "rb",
                        "Yellow - blue" = "y2b",
                        "Matlab 1" = "ml1",
                        "Matlab 2 " = "ml2",
                        "Magenta - Green" = "m2g",
                        "Cyan - yellow" = "c2y",
                        "Blue - yellow" = "b2y",
                        "Green - red" = "g2r",
                        "Blue - green" = "b2g",
                        "Blue - red" = "b2r",
                        "Blue - pink (pastel)" = "b2p",
                        "Blue - green - yellow" = "bgy",
                        "Green - yellow - white" = "gyw",
                        "Red - yellow - white" = "ryw",
                        "Grayscale" = "bw",
                        "Blues (brew)" = "Blues",
                        "Blue - green (brew)" = "BuGn",
                        "Blue - purple (brew)" = "BuPu",
                        "Green - blue (brew)" = "GnBu",
                        "Greens (brew)" = "Greens",
                        "Grayscale (brew)" = "Greys",
                        "Oranges (brew)" = "Oranges",
                        "Orange - red (brew)" = "OrRd",
                        "Purple - blue (brew)" = "PuBu",
                        "Purple - blue - green (brew)" = "PuBuGn",
                        "Purple - red (brew)" = "PuRd",
                        "Purples (brew)" = "Purples",
                        "Red - purple (brew)" = "RdPu",
                        "Reds (brew)" = "Reds",
                        "Yellow - green (brew)" = "YlGn",
                        "Yellow - green - blue (brew)" = "YlGnBu",
                        "Yellow - orange - brown (brew)" = "YlOrBr",
                        "Yellow - orange - red (brew)" = "YlOrRd",
                        "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", # TODO: add descriptions (or remove all?)
                        "RdGy", "RdYlBu", "RdYlGn", "Spectral",
                        "Accent", "Dark2", "Paired", "Pastel1",
                        "Pastel2", "Set1", "Set2", "Set3"
                      ), selected = "rainbow", width = "80%"),
                      # preview plot
                      shinyjqui::jqui_resizable(plotly::plotlyOutput("ramp_plot")),
                      shiny::h2("Discrete data"),
                      shiny::uiOutput("colorPickers")
                    ),
                    shiny::tabPanel(
                      value = "metadata",
                      title = shinyBS::tipify(shiny::icon("plus"),
                        trigger = "hover",
                        title = "upload new metadata file", options = list(container = "body")
                      ),
                      br(), br(),
                      shiny::h2("New metadata"),
                      shiny::helpText("If you want to update your metadata, please use the below upload screen!
                                                                                                                                                   Any samples missing in your uploaded csv will be labeled as 'unknown' for new metadata columns."),
                      shiny::fileInput(inputId = "metadata_new", "Select metadata", buttonLabel = "Browse"),
                      # shinyFiles::shinyFilesButton('metadata_new', 'upload', 'Select metadata', FALSE),
                      shinyWidgets::circleButton("metadata_new_add", icon = shiny::icon("arrow-right"))
                    )
                  )
                )
              )
            )
    ))
  )
}

ui_tab_settings <- function(gbl, adducts) {
  shiny::tabPanel("settings",
    icon = shiny::icon("cog", class = "outlined"), value = "options",
    shiny::tabsetPanel(
      id = "tab_settings", selected = "Global",
      shiny::tabPanel("Global",
        icon = shiny::icon("box-open"),
        # MetaboShiny::switchButton(inputId = "db_only", label = "Run in database-only mode?",
        #                           value = F,
        #                           col = "BW", type = "YN")
        shiny::fluidRow(
          align = "center",
          # work folder
          shinyFiles::shinyDirButton("get_work_dir", label = "Change projects folder", title = "Pick projects directory"),
          shiny::textOutput("curr_exp_dir"),
          # db folder
          shinyFiles::shinyDirButton("get_db_dir", label = "Change database folder", title = "Pick database directory"),
          shiny::textOutput("curr_db_dir"), br(),
          shiny::sliderInput("ncores", "How many cores can MetShi use?",
            value = 1, min = 1,
            max = 2, step = 1, post = " cores"
          ),
          shiny::uiOutput("ncore_ui"),
          shiny::br(),
          shiny::numericInput("seed", "Which seed do you want to use for randomization?", value = 1337),
          shiny::actionButton("set_seed", label = "apply"),
          shiny::br(),
          shiny::helpText("How to handle samples with missing metadata? (re-checked each subset/switch step)?"),
          shinyWidgets::switchInput(
            inputId = "omit_unknown",
            onLabel = "omit",
            offLabel = "keep",
            offStatus = "success",
            onStatus = "danger",
            value = T
          ),
          shiny::br(),
          shiny::selectizeInput("beep",
            label = "Noise to make after task is done:",
            choices = list(
              "none",
              "ping",
              "coin",
              "fanfare",
              "complete",
              "treasure",
              "ready",
              "shotgun",
              "mario",
              "wilhelm",
              "facebook",
              "sword"
            ), selected = 1
          )
        )
      ),
      shiny::tabPanel("Project",
        icon = shiny::icon("gift"),
        # shiny::textInput(inputId="proj_name", label="Project name", value = ''),
        shiny::fluidRow(
          align = "center",
          shiny::selectizeInput(
            inputId = "proj_name",
            label = "Project name",
            choices = c("..."), # existing projects in user folder (generated in 'global')
            selected = "...",
            options = list(create = TRUE)
          ), # let users add new names
          shiny::actionButton("set_proj_name", label = "Apply"),
          shiny::hr(),
          shiny::helpText("Clone project (to skip file input w/ different normalization)"),
          shiny::textInput(inputId = "proj_clone_name", label = "Cloned project name:"),
          shiny::actionButton("clone_proj", label = "Clone")
        )
      ),
      shiny::tabPanel("Search",
        icon = shiny::icon("paw"),
        br(),
        shiny::fluidRow(
          align = "center",
          shinyWidgets::materialSwitch(
            inputId = "autocopy",
            label = "Copy compound info to clipboard?",
            value = TRUE,
            status = "success"
          ),
          shiny::conditionalPanel(
            "input.autocopy == true",
            shinyWidgets::radioGroupButtons(
              inputId = "matchAutocopy",
              label = "What info will you copy?",
              choices = c(
                "SMILES",
                "name",
                "formula"
              ),
              checkIcon = list(
                yes = icon("ok",
                  lib = "glyphicon"
                )
              )
            )
          )
        )
      ),
      # change list of adducts used, or add your own
      # TODO: fix, i think this is currently non-functional
      shiny::tabPanel("Adducts",
        icon = shiny::icon("plus-square"),
        shiny::tabsetPanel(
          shiny::tabPanel(
            "Current adducts",
            shiny::tabsetPanel(
              shiny::tabPanel(
                "Definitions",
                shiny::fluidRow(align = "center", shiny::h3("Current adduct table:")),
                rhandsontable::rHandsontableOutput("adduct_tab",
                  width = "100%",
                  height = "70%"
                ),
                shiny::fluidRow(
                  align = "center", shiny::fileInput(
                    inputId = "add_tab",
                    "Import adduct table",
                    multiple = F,
                    accept = c(".RData", ".csv")
                  ),
                  shiny::div(style = "display: inline-block;vertical-align:top;", shiny::actionButton("import_adducts", "Import definitions"))
                )
              ),
              shiny::tabPanel(
                "Rules",
                shiny::fluidRow(align = "center", shiny::h3("Current adduct rules:")),
                rhandsontable::rHandsontableOutput("adduct_rules_tab",
                  width = "100%",
                  height = "70%"
                ),
                shiny::fileInput(
                  inputId = "add_rule_tab",
                  "Import adduct rule table",
                  multiple = F,
                  accept = c(".RData", ".csv")
                ),
                shiny::fluidRow(
                  align = "center",
                  shiny::div(style = "display: inline-block;vertical-align:top;", shiny::actionButton("import_adduct_rules", "Import rules"))
                )
              )
            ),
            shiny::fluidRow(
              align = "center",
              shiny::actionButton("save_adducts", "Save changes"),
              shiny::hr()
            )
          )
        )
      ),
      shiny::tabPanel("Formula prediction & lookup",
        icon = shiny::icon("searchengin"),
        shiny::fluidRow(
          align = "center",
          br(),
          shiny::passwordInput("apikey",
            label = "If you want to use ChemSpider, please register on their website and enter your API key below.",
            value = ""
          ),
          shiny::div(style = "display: inline-block;vertical-align:top;", shinyWidgets::circleButton("set_api", icon = icon("save"))), shiny::div(style = "display: inline-block;vertical-align:top;", shiny::textOutput("api_set")),
          selectizeInput("predict_rules",
            label = "Use which chemical formula rules?",
            choices = c("senior", "lewis", "hc", "chnops", "nops"),
            multiple = T,
            selected = c()
          ),
          selectizeInput("predict_elements",
            label = "Consider which atoms? (including adducts!)",
            choices = c("C", "H", "N", "O", "P", "S"), selected = c("C", "H", "N", "O", "P", "S"), multiple = T
          ),
          MetaboShiny::switchButton(
            inputId = "predict_details",
            label = "Get detailed PubChem matches? (SLOW, be warned!)",
            col = "BW", type = "YN", value = F
          ),
          MetaboShiny::switchButton(
            inputId = "predict_structure_check",
            label = "Uniformize ChemSpider/PubChem found structures?",
            col = "BW", type = "YN", value = T
          ),
          MetaboShiny::switchButton(
            inputId = "predict_adduct_rules",
            label = "Filter ChemSpider/PubChem hits with available structures using your adduct rules?",
            col = "BW", type = "YN", value = F
          )
        )
      ),
      # change toolbar colour, text font and size
      shiny::tabPanel("Aesthetic",
        icon = shiny::icon("child"),
        shiny::fluidRow(
          align = "center",
          shiny::h3("Change app settings"),
          shiny::hr(),
          shiny::h2("Navigation bar colours"),
          colourpicker::colourInput(
            inputId = "bar.col.1",
            label = paste("Inactive background"),
            # value = opts$col1,
            allowTransparent = FALSE
          ),
          colourpicker::colourInput(
            inputId = "bar.col.2",
            label = paste("Active background"),
            # value = opts$col2,
            allowTransparent = FALSE
          ),
          colourpicker::colourInput(
            inputId = "bar.col.3",
            label = paste("Inactive tab text"),
            # value = opts$col3,
            allowTransparent = FALSE
          ),
          colourpicker::colourInput(
            inputId = "bar.col.4",
            label = paste("Active tab text"),
            # value = opts$col4,
            allowTransparent = FALSE
          ),
          shiny::br(),
          shiny::a(shiny::h2("Fonts (Google fonts)"), href = "https://fonts.google.com/", target = "_blank"),
          shiny::textInput(inputId = "font.1", label = "h1", value = "Leckerli One"),
          shiny::textInput(inputId = "font.2", label = "h2", value = "Leckerli One"),
          shiny::textInput(inputId = "font.3", label = "h3", value = "Open Sans"),
          shiny::textInput(inputId = "font.4", label = "body", value = "Open Sans"),
          shiny::br(), # TODO: font size modifier slider
          shiny::h2("Font size"),
          shiny::sliderInput("size.1", label = "h1", value = 30, min = 5, max = 50),
          shiny::sliderInput("size.2", label = "h2", value = 20, min = 5, max = 50),
          shiny::sliderInput("size.3", label = "h3", value = 13, min = 5, max = 50),
          shiny::sliderInput("size.4", label = "body", value = 10, min = 5, max = 50),
          shiny::hr(),
          shiny::actionButton("change_css", "Save settings (restart to apply)")
        )
      )
    )
  )
}

ui_tab_help <- function(gbl, adducts) {
  shiny::tabPanel("help",
    icon = shiny::icon("question", class = "outlined"),
    shiny::fluidRow(
      align = "center", style = "li { display:inline-table; }",
      shiny::uiOutput("help")
    )
  )
}

ui_footer <- function() {
  shiny::div(
    class = "footer",
    shiny::fluidRow(
      align = "center",
      shiny::span(
        shiny::actionButton("zoom_in",
          label = "",
          icon = shiny::icon("plus")
        ),
        shiny::actionButton("zoom_out",
          label = "",
          icon = shiny::icon("minus")
        ),
        shinyWidgets::prettyToggle(
          inputId = "night_mode",
          label_off = "day",
          label_on = "night",
          outline = TRUE,
          plain = TRUE,
          value = T,
          status_off = "default",
          status_on = "default",
          icon_off = icon("sun"),
          icon_on = icon("moon"),
          inline = T
        ), shiny::tags$i(icon("map-marker"),
          shiny::textOutput("proj_name", inline = T),
          style = "vertical-align: middle;"
        )
      ),
      shiny::actionButton("load_mset",
        label = "load",
        icon = shiny::icon("folder-open")
      ),
      shiny::actionButton("save_mset",
        label = "save",
        icon = shiny::icon("save")
      ),
      shiny::actionButton("debug_metshi",
        label = "send to R",
        icon = shiny::icon("r-project")
      ),
      shiny::actionButton("quit_metshi", label = "quit", icon = shiny::icon("eraser")),
      shiny::span(shiny::tags$i(icon("stopwatch"), shiny::textOutput("last_saved", inline = T), style = "vertical-align: middle;")),
      shiny::span(shiny::uiOutput("has_unsaved_changes", inline = T), style = "margin-left: 10px;")
    )
  )
}
