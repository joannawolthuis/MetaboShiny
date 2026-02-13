ui_pca <- function() {
  shiny::tabPanel("pca",
    value = "pca",
    shiny::fluidRow(
      align = "center",
      shinyBS::bsCollapse(
        multiple = T, id = "collapse_pca",
        shinyBS::bsCollapsePanel(
          title = h2("settings"),
          value = "collapse_pca_settings",
          shinyWidgets::radioGroupButtons("pca_source", "Used data:",
            choices = c(
              "original",
              "pre-batch correction",
              "normalized"
            ),
            selected = "normalized"
          ),
          shinyWidgets::actionBttn(
            inputId = "do_pca",
            label = "click to start PCA",
            style = "bordered",
            icon = icon("terminal"),
            size = "sm"
          )
        ),
        shinyBS::bsCollapsePanel(
          title = h2("plots"), value = "collapse_pca_plots",
          shiny::tabsetPanel(
            shiny::tabPanel(
              "samples",
              shiny::fluidRow(align = "center", shiny::column(
                12,
                (shinyjqui::jqui_resizable(shiny::uiOutput("plot_pca_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700)))
              ))
            ),
            shiny::tabPanel(
              "loadings",
              shiny::fluidRow(align = "center", shiny::column(
                12,
                (shinyjqui::jqui_resizable(shiny::uiOutput("plot_pca_loadings_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700)))
              ))
            )
          ),
          fluidRow(
            align = "center",
            br(),
            shinyWidgets::switchInput(
              inputId = "pca_2d3d",
              size = "mini",
              onLabel = "2D",
              offLabel = "3D",
              value = TRUE
            ),
            helpText("show ellipses"),
            shinyWidgets::switchInput(
              inputId = "pca_ellipse",
              size = "mini",
              onLabel = "Yes",
              offLabel = "No",
              value = TRUE
            ),
            fluidRow(
              column(4, shiny::selectizeInput("pca_x",
                label = "X axis:", choices = {
                  n <- paste0("PC", 1:500)
                  l <- as.list(c(1:500))
                  names(l) <- n
                  l
                },
                selected = 1, width = "80%"
              )),
              column(4, shiny::selectizeInput("pca_y", label = "Y axis:", choices = {
                n <- paste0("PC", 1:500)
                l <- as.list(c(1:500))
                names(l) <- n
                l
              }, selected = 2, width = "80%")),
              column(4, shiny::selectizeInput("pca_z", label = "Z axis:", choices = {
                n <- paste0("PC", 1:500)
                l <- as.list(c(1:500))
                names(l) <- n
                l
              }, selected = 3, width = "80%"))
            )
          )
        ),
        shinyBS::bsCollapsePanel(title = h2("tables"), value = "collapse_pca_tables", shiny::tabsetPanel(
          id = "pca_2",
          shiny::tabPanel(
            title = "table",
            shiny::div(DT::dataTableOutput("pca_tab", width = "100%"), style = "font-size:80%")
          ),
          shiny::tabPanel(
            title = "scree",
            (shinyjqui::jqui_resizable(shiny::uiOutput("pca_scree_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700)))
          ),
          shiny::tabPanel(
            title = "loadings",
            shiny::div(DT::dataTableOutput("pca_load_tab", width = "100%"), style = "font-size:80%")
          )
        )), open = 1
      )
    )
  )
}

ui_plsda <- function() {
  shiny::tabPanel("pls-da",
    value = "plsda",
    shiny::fluidRow(
      align = "center",
      shinyBS::bsCollapse(
        multiple = T, id = "collapse_plsda", shinyBS::bsCollapsePanel(
          title = h2("settings"), value = "collapse_plsda_settings",
          shiny::fluidRow(
            align = "center",
            shiny::div(
              style = "display:inline-block",
              shiny::selectizeInput("plsda_type",
                label = "Type:",
                choices = list("Normal" = "normal")
                # ,
                #             "Orthogonal" = "ortho",
                #             "Sparse" = "sparse")
                , width = "100px",
                selected = 1
              )
            ), br(),
            shinyWidgets::actionBttn(
              inputId = "do_plsda",
              label = "click to start PLS-DA",
              style = "bordered",
              icon = icon("terminal"),
              size = "sm"
            )
          )
        ),
        shinyBS::bsCollapsePanel(
          title = h2("plots"), value = "collapse_plsda_plots",
          shiny::tabsetPanel(
            shiny::tabPanel(
              "samples",
              shiny::fluidRow(align = "center", shiny::column(
                12,
                (shinyjqui::jqui_resizable(shiny::uiOutput("plot_plsda_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700)))
              ))
            ),
            shiny::tabPanel(
              "loadings",
              shiny::fluidRow(align = "center", shiny::column(
                12,
                (shinyjqui::jqui_resizable(shiny::uiOutput("plot_plsda_loadings_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700)))
              ))
            ),
            shiny::tabPanel(
              "performance",
              shiny::tabsetPanel(
                id = "plsda_2",
                shiny::tabPanel(
                  title = "cross-validation",
                  (shinyjqui::jqui_resizable(shiny::uiOutput("plsda_cv_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700)))
                ),
                shiny::tabPanel(
                  title = "permutation",
                  (shinyjqui::jqui_resizable(shiny::uiOutput("plsda_perm_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700)))
                )
              )
            )
          ),
          fluidRow(
            align = "center",
            br(),
            shinyWidgets::switchInput(
              inputId = "plsda_2d3d",
              size = "mini",
              onLabel = "2D",
              offLabel = "3D",
              value = TRUE
            ),
            helpText("show ellipses"),
            shinyWidgets::switchInput(
              inputId = "plsda_ellipse",
              size = "mini",
              onLabel = "Yes",
              offLabel = "No",
              value = TRUE
            ),
            fluidRow(
              column(4, shiny::selectizeInput("plsda_x", label = "X axis:", choices = {
                n <- paste("Component", 1:20)
                l <- as.list(c(1:20))
                names(l) <- n
                l
              }, selected = 1, width = "80%")),
              column(4, shiny::selectizeInput("plsda_y", label = "Y axis:", choices = {
                n <- paste("Component", 1:20)
                l <- as.list(c(1:20))
                names(l) <- n
                l
              }, selected = 2, width = "80%")),
              column(4, shiny::selectizeInput("plsda_z", label = "Z axis:", choices = {
                n <- paste("Component", 1:20)
                l <- as.list(c(1:20))
                names(l) <- n
                l
              }, selected = 3, width = "80%"))
            )
          )
        ),
        shinyBS::bsCollapsePanel(
          title = h2("tables"), value = "collapse_plsda_tables",
          shiny::tabPanel(
            title = "scree",
            shiny::div(DT::dataTableOutput("plsda_tab", width = "100%"), style = "font-size:80%")
          ),
          shiny::tabPanel(
            title = "loadings",
            shiny::div(DT::dataTableOutput("plsda_load_tab", width = "100%"), style = "font-size:80%")
          )
        )
      )
    )
  )
}

ui_tsne <- function() {
  shiny::tabPanel("t-sne",
    value = "tsne",
    shiny::fluidRow(
      align = "center",
      shinyBS::bsCollapse(
        id = "collapse_tsne",
        multiple = T, shinyBS::bsCollapsePanel(
          title = h2("settings"), value = "collapse_tsne_settings",
          shinyWidgets::radioGroupButtons("tsne_source", "Used data:",
            choices = c(
              "original",
              "pre-batch correction",
              "normalized"
            ),
            selected = "normalized"
          ),
          numericInput("tsne_dims", "Initial dimensions:", min = 5, step = 1, value = 30),
          numericInput("tsne_perplex", "Perplexity:", min = 5, value = 30),
          numericInput("tsne_maxiter", "Max iterations:", min = 10, value = 100),
          shinyWidgets::actionBttn(
            inputId = "do_tsne",
            label = "click to start t-SNE",
            style = "bordered",
            icon = icon("terminal"),
            size = "sm"
          )
        ),
        shinyBS::bsCollapsePanel(
          title = h2("plots"), value = "collapse_tsne_plots",
          (shinyjqui::jqui_resizable(shiny::uiOutput("tsne_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700))),
          br(),
          shiny::conditionalPanel(
            "input.ggplotly == true",
            shinyWidgets::switchInput(
              inputId = "tsne_2d3d",
              size = "mini",
              onLabel = "2D",
              offLabel = "3D",
              value = TRUE
            )
          ),
          helpText("show ellipses"),
          shinyWidgets::switchInput(
            inputId = "tsne_ellipse",
            size = "mini",
            onLabel = "Yes",
            offLabel = "No",
            value = TRUE
          ),
          fluidRow(
            column(4, shiny::selectizeInput("tsne_x", label = "X axis:", choices = {
              n <- paste("t-sne component", 1:20)
              l <- as.list(c(1:20))
              names(l) <- n
              l
            }, selected = 1, width = "80%")),
            column(4, shiny::selectizeInput("tsne_y", label = "Y axis:", choices = {
              n <- paste("t-sne component", 1:20)
              l <- as.list(c(1:20))
              names(l) <- n
              l
            }, selected = 2, width = "80%")),
            column(4, shiny::selectizeInput("tsne_z", label = "Z axis:", choices = {
              n <- paste("t-sne component", 1:20)
              l <- as.list(c(1:20))
              names(l) <- n
              l
            }, selected = 3, width = "80%"))
          )
        )
      )
    )
  )
}

ui_ica <- function() {
  shiny::tabPanel("ica",
    value = "ica",
    shiny::fluidRow(
      align = "center",
      shinyBS::bsCollapse(
        id = "collapse_ica",
        multiple = T, shinyBS::bsCollapsePanel(
          title = h2("settings"), value = "collapse_ica_settings",
          shinyWidgets::radioGroupButtons("ica_source", "Used data:",
            choices = c(
              "original",
              "pre-batch correction",
              "normalized"
            ),
            selected = "normalized"
          ),
          numericInput("ica_ncomp", "Components calculated:", min = 2, step = 1, value = 3),
          numericInput("ica_maxiter", "Max iterations:", min = 10, value = 100),
          shinyWidgets::radioGroupButtons(
            inputId = "ica_method",
            label = "Method:",
            choices = c(
              "fast",
              "imax", "jade"
            ),
            justified = TRUE, selected = "fast",
            checkIcon = list(
              yes = icon("ok",
                lib = "glyphicon"
              )
            )
          ),
          shinyWidgets::actionBttn(
            inputId = "do_ica",
            label = "click to start ICA (independant component analysis)",
            style = "bordered",
            icon = icon("terminal"),
            size = "sm"
          )
        ),
        shinyBS::bsCollapsePanel(
          title = h2("plots"), value = "collapse_ica_plots",
          (shinyjqui::jqui_resizable(shiny::uiOutput("ica_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700))),
          br(),
          shiny::conditionalPanel(
            "input.ggplotly == true",
            shinyWidgets::switchInput(
              inputId = "ica_2d3d",
              size = "mini",
              onLabel = "2D",
              offLabel = "3D",
              value = TRUE
            )
          ),
          helpText("show ellipses"),
          shinyWidgets::switchInput(
            inputId = "ica_ellipse",
            size = "mini",
            onLabel = "Yes",
            offLabel = "No",
            value = TRUE
          ),
          fluidRow(
            column(4, shiny::selectizeInput("ica_x", label = "X axis:", choices = {
              n <- paste0("IC", 1:20)
              l <- as.list(c(1:20))
              names(l) <- n
              l
            }, selected = 1, width = "80%")),
            column(4, shiny::selectizeInput("ica_y", label = "Y axis:", choices = {
              n <- paste0("IC", 1:20)
              l <- as.list(c(1:20))
              names(l) <- n
              l
            }, selected = 2, width = "80%")),
            column(4, shiny::selectizeInput("ica_z", label = "Z axis:", choices = {
              n <- paste0("IC", 1:20)
              l <- as.list(c(1:20))
              names(l) <- n
              l
            }, selected = 3, width = "80%"))
          )
        )
      )
    )
  )
}

ui_umap <- function() {
  shiny::tabPanel("umap",
    value = "umap",
    shiny::fluidRow(
      align = "center",
      shinyBS::bsCollapse(
        id = "collapse_umap",
        multiple = T, shinyBS::bsCollapsePanel(
          title = h2("settings"), value = "collapse_umap_settings",
          shinyWidgets::radioGroupButtons("umap_source", "Used data:",
            choices = c(
              "original",
              "pre-batch correction",
              "normalized"
            ),
            selected = "normalized"
          ),
          sliderInput("umap_ncomp", "Components calculated:", min = 2, max = 20, step = 1, value = 3),
          numericInput("umap_maxiter", "Max iterations:", min = 10, value = 100),
          numericInput("umap_neighbors", "Neighbors used", value = 15),
          shinyWidgets::actionBttn(
            inputId = "do_umap",
            label = "click to start UMAP",
            style = "bordered",
            icon = icon("terminal"),
            size = "sm"
          )
        ),
        shinyBS::bsCollapsePanel(
          title = h2("plots"), value = "collapse_umap_plots",
          (shinyjqui::jqui_resizable(shiny::uiOutput("umap_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700))),
          br(),
          shiny::conditionalPanel(
            "input.ggplotly == true",
            shinyWidgets::switchInput(
              inputId = "umap_2d3d",
              size = "mini",
              onLabel = "2D",
              offLabel = "3D",
              value = TRUE
            )
          ),
          helpText("show ellipses"),
          shinyWidgets::switchInput(
            inputId = "umap_ellipse",
            size = "mini",
            onLabel = "Yes",
            offLabel = "No",
            value = TRUE
          ),
          fluidRow(
            column(4, shiny::selectizeInput("umap_x", label = "X axis:", choices = {
              n <- paste("UMAP component", 1:20)
              l <- as.list(c(1:20))
              names(l) <- n
              l
            }, selected = 1, width = "80%")),
            column(4, shiny::selectizeInput("umap_y", label = "Y axis:", choices = {
              n <- paste("UMAP component", 1:20)
              l <- as.list(c(1:20))
              names(l) <- n
              l
            }, selected = 2, width = "80%")),
            column(4, shiny::selectizeInput("umap_z", label = "Z axis:", choices = {
              n <- paste("UMAP component", 1:20)
              l <- as.list(c(1:20))
              names(l) <- n
              l
            }, selected = 3, width = "80%"))
          )
        )
      )
    )
  )
}

ui_t_test <- function() {
  shiny::tabPanel("t-test",
    value = "tt",
    shiny::fluidRow(
      align = "center",
      shiny::fluidRow(
        align = "center",
        shinyBS::bsCollapse(
          multiple = T, id = "collapse_tt",
          shinyBS::bsCollapsePanel(
            title = h2("settings"), value = "collapse_tt_settings",
            MetaboShiny::switchButton("tt_nonpar", "Non-parametric?", col = "BW", type = "YN", value = T),
            MetaboShiny::switchButton("tt_eqvar", "Equal variance?", col = "BW", type = "YN", value = T),
            MetaboShiny::switchButton("tt_paired", "Paired analysis?", col = "BW", type = "YN", value = F),
            shiny::selectInput("tt_multi_test", "Multiple testing correction method:",
              choices = list(
                "Holm" = "holm",
                "Hochmerg" = "hochberg",
                "Hommel" = "hommel",
                "Bonferroni" = "bonferroni",
                "Benjamini & Hochberg" = "fdr",
                "Benjamini & Yekutieli" = "BY",
                "none"
              ),
              selected = "fdr"
            ),
            shiny::numericInput("tt_p_thresh",
              label = "Maximum p-value after multiple testing:",
              value = 0.05, max = 1, min = 0
            ),
            shinyWidgets::actionBttn(
              inputId = "do_tt",
              label = "click to start t-test",
              style = "bordered",
              icon = icon("terminal"),
              size = "sm"
            )
          ),
          shinyBS::bsCollapsePanel(
            title = h2("plots"), value = "collapse_tt_plots",
            (shinyjqui::jqui_resizable(shiny::uiOutput("tt_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700))),
            shiny::sliderInput("tt_topn", label = "Show top:", min = 5, max = 5000, value = 100)
          ),
          shinyBS::bsCollapsePanel(
            title = h2("tables"), value = "collapse_tt_tables",
            shiny::div(DT::dataTableOutput("tt_tab", width = "100%"), style = "font-size:80%")
          )
        )
      )
    )
  )
}

ui_logiscore <- function() {
  shiny::tabPanel("logistic scores",
    value = "logiscore",
    shiny::fluidRow(
      align = "center",
      shiny::fluidRow(
        align = "center",
        shinyBS::bsCollapse(
          multiple = T, id = "collapse_logiscore",
          shinyBS::bsCollapsePanel(
            title = h2("settings"), value = "collapse_logiscore_settings",
            shiny::selectInput("logiscore_multi_test", "Multiple testing correction method:",
              choices = list(
                "Holm" = "holm",
                "Hochmerg" = "hochberg",
                "Hommel" = "hommel",
                "Bonferroni" = "bonferroni",
                "Benjamini & Hochberg" = "fdr",
                "Benjamini & Yekutieli" = "BY",
                "none"
              ),
              selected = "fdr"
            ),
            shiny::numericInput("logiscore_p_thresh",
              label = "Maximum p-value after multiple testing:",
              value = 0.05, max = 1, min = 0
            ),
            shinyWidgets::actionBttn(
              inputId = "do_logiscore",
              label = "click to start logistic scores",
              style = "bordered",
              icon = icon("terminal"),
              size = "sm"
            )
          ),
          shinyBS::bsCollapsePanel(
            title = h2("plots"), value = "collapse_logiscore_plots",
            (shinyjqui::jqui_resizable(shiny::uiOutput("logiscore_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700))),
            shiny::sliderInput("logiscore_topn", label = "Show top:", min = 5, max = 5000, value = 100)
          ),
          shinyBS::bsCollapsePanel(
            title = h2("tables"), value = "collapse_logiscore_tables",
            shiny::div(DT::dataTableOutput("logiscore_tab", width = "100%"), style = "font-size:80%")
          )
        )
      )
    )
  )
}

ui_proda <- function() {
  shiny::tabPanel("proDA",
    value = "proda",
    shiny::helpText("T-test equivalent for data with missing values intact (uses non-imputed data)."),
    shiny::fluidRow(
      align = "center",
      shiny::fluidRow(
        align = "center",
        # ok
        shinyBS::bsCollapse(
          multiple = T, id = "collapse_proda",
          shinyBS::bsCollapsePanel(
            title = h2("settings"), value = "collapse_proda_settings",
            shinyWidgets::switchInput("proda_add_batch",
              "Add batch in model?",
              value = TRUE,
              onLabel = "yes",
              offLabel = "no"
            ),
            shiny::selectInput("proda_multi_test", "Multiple testing correction method:",
              choices = list(
                "Holm" = "holm",
                "Hochmerg" = "hochberg",
                "Hommel" = "hommel",
                "Bonferroni" = "bonferroni",
                "Benjamini & Hochberg" = "fdr",
                "Benjamini & Yekutieli" = "BY",
                "none"
              ),
              selected = "fdr"
            ),
            shinyWidgets::actionBttn(
              inputId = "do_proda",
              label = "click to start proDA t-test",
              style = "bordered",
              icon = icon("terminal"),
              size = "sm"
            )
          ),
          shinyBS::bsCollapsePanel(
            title = h2("plots"), value = "collapse_proda_plots",
            (shinyjqui::jqui_resizable(shiny::uiOutput("proda_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700))),
            shiny::sliderInput("proda_topn", label = "Show top:", min = 5, max = 5000, value = 100)
          ),
          shinyBS::bsCollapsePanel(
            title = h2("tables"), value = "collapse_proda_tables",
            shiny::div(DT::dataTableOutput("proda_tab", width = "100%"), style = "font-size:80%")
          )
        )
      )
    )
  )
}

ui_anova <- function() {
  shiny::tabPanel("anova",
    value = "aov",
    shiny::fluidRow(
      align = "center",
      # ok
      shinyBS::bsCollapse(
        id = "collapse_aov",
        multiple = T,
        shinyBS::bsCollapsePanel(
          title = h2("settings"), value = "collapse_aov_settings",
          shinyWidgets::actionBttn(
            inputId = "do_aov",
            label = "click to start ANOVA",
            style = "bordered",
            icon = icon("terminal"),
            size = "sm"
          )
        ),
        shinyBS::bsCollapsePanel(
          title = h2("plots"), value = "collapse_aov_plots",
          (shinyjqui::jqui_resizable(shiny::uiOutput("aov_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700))),
          shiny::sliderInput("aov_topn", label = "Show top:", min = 5, max = 5000, value = 100)
        ),
        shinyBS::bsCollapsePanel(
          title = h2("tables"), value = "collapse_aov_tables",
          shiny::div(DT::dataTableOutput("aov_tab", width = "100%"), style = "font-size:80%")
        )
      )
    )
  )
}

ui_fold_change <- function() {
  shiny::tabPanel("fold-change",
    value = "fc",
    shiny::fluidRow(
      align = "center",
      # ok
      shinyBS::bsCollapse(
        id = "collapse_fc", multiple = T,
        shinyBS::bsCollapsePanel(
          title = h2("settings"), value = "collapse_fc_settings",
          shiny::sliderInput("fc_thresh", "Fold-change threshold:",
            min = 0,
            max = 10,
            value = 1.2,
            step = 0.1
          ),
          MetaboShiny::switchButton("fc_paired", "Paired analysis?", col = "BW", type = "YN", value = F),
          shinyWidgets::actionBttn(
            inputId = "do_fc",
            label = "click to start fold-change analysis",
            style = "bordered",
            icon = icon("terminal"),
            size = "sm"
          )
        ),
        shinyBS::bsCollapsePanel(
          title = h2("plots"), value = "collapse_fc_plots",
          (shinyjqui::jqui_resizable(shiny::uiOutput("fc_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700))),
          shiny::sliderInput("fc_topn", label = "Show top:", min = 5, max = 5000, value = 100)
        ),
        shinyBS::bsCollapsePanel(
          title = h2("tables"), value = "collapse_fc_tables",
          shiny::div(DT::dataTableOutput("fc_tab", width = "100%"), style = "font-size:80%")
        )
      )
    )
  )
}

ui_cliffs_delta <- function() {
  shiny::tabPanel("cliff's delta",
    value = "cliffd",
    shiny::fluidRow(
      align = "center",
      # ok
      shinyBS::bsCollapse(
        id = "collapse_cliffd", multiple = T,
        shinyBS::bsCollapsePanel(
          title = h2("settings"), value = "collapse_cliffd_settings",
          shinyWidgets::actionBttn(
            inputId = "do_cliffd",
            label = "click to start cliff's delta analysis",
            style = "bordered",
            icon = icon("terminal"),
            size = "sm"
          )
        ),
        shinyBS::bsCollapsePanel(
          title = h2("plots"), value = "collapse_cliffd_plots",
          (shinyjqui::jqui_resizable(shiny::uiOutput("cliffd_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700))),
          shiny::sliderInput("cliffd_topn", label = "Show top:", min = 5, max = 5000, value = 100)
        ),
        shinyBS::bsCollapsePanel(
          title = h2("tables"), value = "collapse_cliffd_tables",
          shiny::div(DT::dataTableOutput("cliffd_tab", width = "100%"), style = "font-size:80%")
        )
      )
    )
  )
}

ui_meba <- function() {
  shiny::tabPanel("meba",
    value = "meba",
    shiny::fluidRow(
      align = "center",
      # ok
      shinyBS::bsCollapse(
        id = "collapse_meba", multiple = T,
        shinyBS::bsCollapsePanel(
          title = h2("settings"), value = "collapse_meba_settings",
          shinyWidgets::actionBttn(
            inputId = "do_meba",
            label = "click to start MEBA analysis",
            style = "bordered",
            icon = icon("terminal"),
            size = "sm"
          )
        ),
        shinyBS::bsCollapsePanel(
          title = h2("plots"), value = "collapse_meba_plots",
          (shinyjqui::jqui_resizable(shiny::uiOutput("meba_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700))),
          shiny::sliderInput("meba_topn", label = "Show top:", min = 5, max = 5000, value = 100)
        ),
        shinyBS::bsCollapsePanel(
          title = h2("tables"), value = "collapse_meba_tables",
          shiny::div(DT::dataTableOutput("meba_tab", width = "100%"), style = "font-size:80%")
        )
      )
    )
  )
}

ui_asca <- function() {
  shiny::tabPanel("asca",
    value = "asca",
    shiny::fluidRow(
      align = "center",
      # ok
      shinyBS::bsCollapse(
        id = "collapse_asca", multiple = T,
        shinyBS::bsCollapsePanel(
          title = h2("settings"), value = "collapse_asca_settings",
          shinyWidgets::actionBttn(
            inputId = "do_asca",
            label = "click to start ASCA analysis",
            style = "bordered",
            icon = icon("terminal"),
            size = "sm"
          )
        ),
        shinyBS::bsCollapsePanel(
          title = h2("plots"), value = "collapse_asca_plots",
          (shinyjqui::jqui_resizable(shiny::uiOutput("asca_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700)))
        ),
        shinyBS::bsCollapsePanel(
          title = h2("tables"), value = "collapse_asca_tables",
          shiny::div(DT::dataTableOutput("asca_tab", width = "100%"), style = "font-size:80%")
        )
      )
    )
  )
}

ui_pattern <- function() {
  shiny::tabPanel("pattern",
    value = "corr",
    shiny::fluidRow(
      align = "center",
      # ok
      shinyBS::bsCollapse(
        id = "collapse_corr", shinyBS::bsCollapsePanel(
          title = h2("settings"), value = "collapse_corr_settings",
          uiOutput("jqui_ui"), br(),
          shiny::selectizeInput("corr_corr",
            "Correlation metric:",
            choices = c("pearson", "spearman", "kendall"),
            selected = "spearman"
          ),
          # corr_p_thresh corr_r_thresh
          shiny::numericInput("corr_p_thresh",
            label = "Maximum p-value:",
            value = 0.1
          ),
          shiny::numericInput("corr_r_thresh",
            label = "Minimum correlation (absolute):",
            value = 0.1
          ),
          shinyWidgets::actionBttn(
            inputId = "do_corr",
            label = "click to start pattern finding",
            style = "bordered",
            icon = icon("terminal"),
            size = "sm"
          )
        ),
        shinyBS::bsCollapsePanel(
          title = h2("plots"), value = "collapse_corr_plots",
          (shinyjqui::jqui_resizable(shiny::uiOutput("corr_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700))),
          shiny::sliderInput("corr_topn", label = "Show top:", min = 5, max = 200, value = 20)
        ),
        shinyBS::bsCollapsePanel(
          title = h2("tables"), value = "collapse_corr_tables",
          shiny::div(DT::dataTableOutput("corr_tab", width = "100%"), style = "font-size:80%")
        )
      )
    )
  )
}

ui_intersection_plot <- function() {
  shiny::tabPanel("intersection plot",
    value = "combi",
    shiny::fluidRow(
      align = "center",
      shinyBS::bsCollapse(
        multiple = T, id = "collapse_combi",
        shinyBS::bsCollapsePanel(
          title = h2("settings"), value = "collapse_combi_settings",
          shiny::fluidRow(
            align = "center",
            column(
              6, shiny::selectizeInput("combi_anal1", label = "X-axis:", choices = c()),
              shiny::uiOutput("combi_anal1_picker"),
              shiny::selectizeInput("combi_anal1_trans", label = "Transformation:", choices = c(
                "none",
                "log10",
                "-log10",
                "abs"
              ), selected = "none")
            ),
            column(
              6, shiny::selectizeInput("combi_anal2", label = "Y-axis:", choices = c()),
              shiny::uiOutput("combi_anal2_picker"),
              shiny::selectizeInput("combi_anal2_trans", label = "Transformation:", choices = c(
                "none",
                "log10",
                "-log10",
                "abs"
              ), selected = "none")
            )
          ),
          shiny::fluidRow(
            align = "center",
            shiny::selectInput("combi_dist_metric",
              "Calculate x vs. y distance metrics",
              choices = list(
                "multiplication",
                "euclidean",
                "maximum",
                "manhattan",
                "minkowski"
              )
            )
          ),
          shinyWidgets::switchInput("combi_highlight_top",
            "Highlight top m/z values?",
            value = FALSE,
            onLabel = "yes",
            offLabel = "no"
          ),
          shiny::conditionalPanel(
            "input.combi_highlight_top == true",
            shiny::sliderInput("combi_highlight_top_n",
              label = "Highlight top m/z:",
              min = 1,
              max = 100,
              step = 1,
              value = 1
            )
          ),
          shinyWidgets::actionBttn(
            inputId = "do_combi",
            label = "click to make combi plot",
            style = "bordered",
            icon = icon("terminal"),
            size = "sm"
          )
        ),
        shinyBS::bsCollapsePanel(
          title = h2("plots"), value = "collapse_combi_plots",
          (shinyjqui::jqui_resizable(shiny::uiOutput("combi_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700)))
        ),
        shinyBS::bsCollapsePanel(
          title = h2("tables"), value = "collapse_combi_tables",
          shiny::div(DT::dataTableOutput("combi_tab", width = "100%"), style = "font-size:80%")
        )
      )
    )
  )
}

ui_heatmap <- function() {
  shiny::tabPanel("heatmap",
    value = "heatmap",
    shiny::fluidRow(
      align = "center",
      shinyBS::bsCollapse(
        multiple = T,
        shinyBS::bsCollapsePanel(
          title = h2("settings"), value = "collapse_heatmap_settings",
          shiny::selectizeInput("heattable",
            "Use which analysis?",
            choices = c("none available"),
            selected = 1,
            multiple = F
          ),
          MetaboShiny::switchButton("heatsign", label = "Only significant hits?", col = "GB", type = "YN"),
          MetaboShiny::switchButton("heatlimits", label = "Color based on -all- metabolites?", col = "GB", type = "YN"),
          shinyWidgets::radioGroupButtons("heatmap_source", "Used data:",
            choices = c(
              "original",
              "pre-batch correction",
              "normalized"
            ),
            selected = "normalized"
          ),
          shinyWidgets::actionBttn(
            inputId = "do_heatmap",
            label = "click to make heatmap",
            style = "bordered",
            icon = icon("terminal"),
            size = "sm"
          )
        ),
        shinyBS::bsCollapsePanel(
          title = h2("plots"), value = "collapse_heatmap_plots",
          shiny::verbatimTextOutput("heatmap_now", placeholder = F),
          (shinyjqui::jqui_resizable(shiny::uiOutput("heatmap_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700))),
          shiny::sliderInput("heatmap_topn",
            label = "Show top:",
            min = 10,
            max = 5000,
            step = 10,
            value = 20
          )
        )
      )
    )
  )
}

ui_network <- function() {
  shiny::tabPanel("network",
    value = "network",
    shiny::fluidRow(
      align = "center",
      shinyBS::bsCollapse(
        multiple = T,
        shinyBS::bsCollapsePanel(
          title = h2("settings"),
          value = "collapse_network_settings",
          shinyWidgets::switchInput("network_sel",
            "Include",
            value = TRUE,
            onLabel = "all",
            offLabel = "subset"
          ),
          shiny::conditionalPanel(
            "input.network_sel == false",
            shiny::selectizeInput("network_table",
              "Use which analysis?",
              choices = c("none available"),
              selected = 1,
              multiple = F
            ),
            shiny::sliderInput("network_topn",
              label = "Use top m/z:",
              min = 10,
              max = 5000,
              step = 10,
              value = 20
            )
          ),
          shinyWidgets::switchInput("network_highlight_top",
            "Highlight top m/z values?",
            value = FALSE,
            onLabel = "yes",
            offLabel = "no"
          ),
          shiny::conditionalPanel(
            "input.network_highlight_top == true",
            shiny::sliderInput("network_highlight_top_n",
              label = "Highlight top m/z:",
              min = 1,
              max = 100,
              step = 1,
              value = 1
            )
          ),
          shiny::numericInput("network_sign",
            label = "p-value threshhold",
            value = 0.05,
            min = 0.0000001,
            max = 1
          ),
          shiny::numericInput("network_minr",
            label = "min absolute edge correlation",
            value = 0.7,
            min = 0.1,
            max = 1
          ),
          shinyWidgets::switchInput("network_auto",
            "Network style",
            value = T,
            onLabel = "auto",
            offLabel = "choose"
          ),
          shiny::conditionalPanel(
            "input.network_auto == false",
            shiny::selectizeInput("network_style", "Style:",
              choices = list(
                "Hierarchical" = "hierarchical",
                "Circular" = "layout_in_circle",
                "Nicely" = "layout_nicely",
                "Sugiyama (slow)" = "layout_with_sugiyama",
                "Star" = "layout_as_star",
                "Tree" = "layout_as_tree",
                "Grid" = "layout_on_grid",
                "Sphere" = "layout_on_sphere",
                "Randomly" = "layout_randomly",
                "DH (slow)" = "layout_with_dh",
                "FR" = "layout_with_fr",
                "Gem" = "layout_with_gem",
                "Graphopt" = "layout_with_graphopt",
                "KK" = "layout_with_kk",
                "LGL" = "layout_with_lgl",
                "MDS (slow)" = "layout_with_mds"
              )
            )
          ),
          shiny::selectizeInput("network_corr",
            label = "Correlation method:",
            choices = c("pearson", "kendall", "spearman"),
            multiple = F,
            selected = "pearson"
          ),
          shinyWidgets::actionBttn(
            inputId = "do_network",
            label = "click to make network",
            style = "bordered",
            icon = icon("terminal"),
            size = "sm"
          )
        ),
        shinyBS::bsCollapsePanel(
          title = h2("plots"), value = "collapse_network_plots",
          shiny::verbatimTextOutput("network_now", placeholder = F),
          shiny::tabsetPanel(
            id = "network_results",
            shiny::tabPanel(
              title = "network", value = "network",
              icon = shiny::icon("project-diagram"),
              shiny::fluidRow(align = "center", (shinyjqui::jqui_resizable(shiny::uiOutput("network_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700)))),
            ),
            shiny::tabPanel(
              title = "heatmap", value = "heatmap",
              icon = shiny::icon("th"),
              shiny::fluidRow(align = "center", (shinyjqui::jqui_resizable(shiny::uiOutput("network_heatmap_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700)))),
            )
          )
        )
      )
    )
  )
}

ui_multirank <- function() {
  shiny::tabPanel("multirank",
    value = "multirank", # icon=shiny::icon("comments"),
    br(),
    sidebarLayout(
      position = "left",
      sidebarPanel = shiny::sidebarPanel(
        width = 4,
        shiny::fluidRow(shiny::div(DT::dataTableOutput("multirank_unselected"), style = "font-size:80%"), align = "center"),
        shiny::fluidRow(shinyWidgets::circleButton("multirank_add", icon = shiny::icon("arrow-down"), size = "sm"),
          shinyWidgets::circleButton("multirank_remove", icon = shiny::icon("arrow-up"), size = "sm"),
          align = "center"
        ),
        shiny::fluidRow(shiny::div(DT::dataTableOutput("multirank_selected"), style = "font-size:80%"), align = "center"),
        shiny::hr(),
        shiny::fluidRow(
          shiny::sliderInput("multirank_topn", label = "Only plot top:", min = 1, max = 200, post = " hits", value = 20),
          align = "center"
        ),
        shiny::fluidRow(
          shinyWidgets::actionBttn(
            inputId = "do_multirank",
            label = "click to find combined ranking",
            style = "bordered",
            icon = icon("terminal"),
            size = "sm"
          ),
          align = "center"
        )
      ),
      mainPanel = mainPanel(
        shiny::hr(),
        (shinyjqui::jqui_resizable(shiny::uiOutput("multirank_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700))),
        shiny::fluidRow(shiny::div(DT::dataTableOutput("multirank_tab"), style = "font-size:80%"),
          align = "center"
        )
      )
    )
  )
}

ui_venn <- function() {
  shiny::tabPanel("venn",
    value = "venn", # icon=shiny::icon("comments"),
    br(),
    sidebarLayout(
      position = "left",
      sidebarPanel = shiny::sidebarPanel(
        width = 4,
        shiny::fluidRow(shiny::div(DT::dataTableOutput("venn_unselected"), style = "font-size:80%"), align = "center"),
        shiny::fluidRow(shinyWidgets::circleButton("venn_add", icon = shiny::icon("arrow-down"), size = "sm"),
          shinyWidgets::circleButton("venn_remove", icon = shiny::icon("arrow-up"), size = "sm"),
          align = "center"
        ),
        shiny::fluidRow(shiny::div(DT::dataTableOutput("venn_selected"), style = "font-size:80%"), align = "center"),
        shiny::hr(),
        shiny::selectizeInput("venn_filter_mode",
          label = "Use top hits or an absolute threshold?",
          choices = c("top", "threshold", "top_and_threshold"), selected = "top"
        ),
        shiny::conditionalPanel(
          "input.venn_filter_mode == 'threshold' || input.venn_filter_mode == 'top_and_threshold'",
          shiny::uiOutput("venn_threshold_ui")
        ),
        shiny::conditionalPanel(
          "input.venn_filter_mode == 'top' || input.venn_filter_mode == 'top_and_threshold'",
          shiny::fluidRow(
            shiny::numericInput("venn_tophits", label = "Only include top:", min = 1, value = 20)
            # shiny::sliderInput("venn_tophits", label = "Only include top:", min = 1, max = 200, post = " hits", value=20)
            ,
            align = "center"
          )
        ),
        shinyWidgets::switchInput(
          inputId = "venn_plot_mode",
          value = FALSE,
          label = "Plot type",
          onLabel = "upset",
          offLabel = "venn",
          size = "small" # "<div class=\"fa-flip-vertical\"><i class=\"fas fa-chart-bar fa-rotate-90\"></i></div>"
        ),
        shiny::fluidRow(
          shinyWidgets::actionBttn(
            inputId = "venn_build",
            label = "click to make venn diagram",
            style = "bordered",
            icon = icon("terminal"),
            size = "sm"
          ),
          align = "center"
        )
      ),
      mainPanel = mainPanel(
        shiny::hr(),
        (shinyjqui::jqui_resizable(shiny::uiOutput("venn_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700))),
        # find the overlapping compounds between the groups you want to compare (user select)
        # TODO: enable this with clicking the numbers/areas
        shiny::fluidRow(shiny::selectizeInput(width = "80%", "intersect_venn", label = "Show hits from (only):", selected = 1, choices = "", multiple = T),
          align = "center"
        ),
        shiny::fluidRow(shiny::uiOutput("venn_pval"), align = "center"),
        shiny::br(),
        shiny::fluidRow(shiny::div(DT::dataTableOutput("venn_tab"), style = "font-size:80%"),
          align = "center"
        )
      )
    )
  )
}

ui_power <- function() {
  shiny::tabPanel("power calculation",
    value = "power",
    shiny::fluidRow(
      align = "center",
      shinyBS::bsCollapse(
        multiple = T, id = "collapse_power", shinyBS::bsCollapsePanel(
          title = h2("settings"), value = "collapse_power_settings",
          shiny::sliderInput("power_nsamp",
            label = "Up to how many samples (per group)?",
            value = 500,
            min = 5,
            max = 999,
            step = 1
          ),
          shiny::selectizeInput("power_comps", "Which comparisons do you want to make?",
            choices = c(" "),
            multiple = T
          ),
          shiny::numericInput("power_fdr", "False discovery rate:",
            value = 0.1,
            max = 1,
            min = 0.00001
          ),
          shinyWidgets::actionBttn(
            inputId = "do_power",
            label = "click to start power analysis",
            style = "bordered",
            icon = icon("terminal"),
            size = "sm"
          )
        ),
        shinyBS::bsCollapsePanel(
          title = h2("plots"), value = "collapse_power_plots",
          (shinyjqui::jqui_resizable(shiny::uiOutput("power_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700)))
        )
      )
    )
  )
}

ui_enrichment <- function() {
  shiny::tabPanel("enrichment",
    value = "enrich",
    sidebarLayout(
      position = "left",
      sidebarPanel = sidebarPanel(
        shiny::fluidRow(
          align = "center",
          h3("Enrichment parameters"),
          shiny::selectizeInput("mummi_lib",
            label = "Pathway database:",
            choices = list(
              "Homo sapiens (human) MFN" = "hsa_mfn",
              "Homo sapiens (human) KEGG" = "hsa01100",
              "Gutsy Full" = "custom_gutsy_full",
              "Gutsy Full (Genus level ID)" = "custom_gutsy_genus",
              "Gutsy Union (Genus level)" = "custom_gutsy_union",
              "Gutsy Intersect (Genus level)" = "custom_gutsy_intersect",
              "Gutsy Full (Pos Rho)" = "custom_gutsy_genus_poscor",
              "Gutsy Full (Neg Rho)" = "custom_gutsy_genus_negcor",
              "Sus scrofa (swine)" = "ssc01100",
              "Microbial metabolism in diverse environments" = "map01120",
              "Bos taurus (cow)" = "bsa01100",
              "Mus musculus (mouse)" = "mmu01100",
              "Rattus norvegicus (rat)" = "rno01100",
              "Danio rerio (zebrafish)" = "dre01100",
              "Caenorhabditis elegans (nematode)" = "cel01100",
              "Gallus gallus (chicken)" = "gga01100",
              "Escherichia coli K-12 MG1655" = "eco01100",
              # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4049927 clostridium, ruminococcus, lactobacillus, bacteroides
              "Clostridium perfringens 13" = "cpe01100",
              "Ruminococcus albus 7" = "ral01100",
              "Bacteroides salanitronis DSM 18170" = "bsa01100",
              "Lactobacillus salivarius UCC188" = "lsl01100",
              "Lipids - Main Chemical Class" = "main_lipid_class_mset",
              "Lipids - Sub Chemical Class" = "sub_lipid_class_mset",
              "Non-Lipids - Main Chemical Class" = "main_nolipid_class_mset",
              "Non-Lipids - Sub Chemical Class" = "sub_nolipid_class_mset",
              "Disease-associated Metabolite Sets (Blood)" = "blood_mset",
              "Disease-associated Metabolite Sets (CSF)" = "csf_mset",
              "Disease-associated Metabolite Sets (Urine)" = "urine_mset",
              "SNP-associated Metabolite Sets" = "snp_mset",
              "Location-based Metabolite Sets" = "location_mset",
              "Predicted Metabolite Sets" = "predicted_mset"
            ), selected = "hsa01120"
          ), # REGEX: [:|"|,][A-z]+?_[A-z]+?[:|"|,] ON the metaboanalyst page showing the table with options
          shinyWidgets::switchInput(
            inputId = "mummi_enr_method", value = TRUE,
            onLabel = "mummichog",
            offLabel = "gsea" # "<div class=\"fa-flip-vertical\"><i class=\"fas fa-chart-bar fa-rotate-90\"></i></div>"
          ),
          shinyWidgets::switchInput(
            inputId = "mummi_abbrev",
            value = FALSE,
            onLabel = "yes",
            label = "Abbreviate pathways in plot?",
            offLabel = "no" # "<div class=\"fa-flip-vertical\"><i class=\"fas fa-chart-bar fa-rotate-90\"></i></div>"
          ),
          shinyWidgets::switchInput(
            inputId = "mummi_rules",
            value = FALSE,
            onLabel = "yes",
            label = "Apply adduct rules?",
            offLabel = "no" # "<div class=\"fa-flip-vertical\"><i class=\"fas fa-chart-bar fa-rotate-90\"></i></div>"
          ),
          shiny::textInput("mummi_pval", label = "Required p-value (if applicable):", value = "0.005"),
          br(),
          h3("User data parameters"),
          shiny::selectizeInput("mummi_anal",
            "Use analysis:",
            choices = c("none available"),
            selected = 1,
            multiple = F
          ),
          shiny::selectizeInput("mummi_filter_mode",
            label = "Use top hits or top plus threshold?",
            choices = c("top", "top_and_threshold"), selected = "top"
          ),
          shiny::conditionalPanel(
            "input.mummi_filter_mode == 'threshold' || input.mummi_filter_mode == 'top_and_threshold'",
            shiny::fluidRow(
              column(width = 6, shiny::selectInput("mummi_threshold_sign",
                label = "sign:",
                choices = c(">", "<"),
                multiple = F,
                selected = ">"
              )),
              column(width = 6, shiny::numericInput("mummi_threshold_value",
                label = "value:",
                min = 0,
                max = 100,
                step = 0.001,
                value = 0
              ))
            )
          ),
          shiny::conditionalPanel(
            "input.mummi_filter_mode == 'top' || input.mummi_filter_mode == 'top_and_threshold'",
            shiny::fluidRow(
              shiny::numericInput("mummi_topn",
                "Top hits used:",
                min = 10,
                max = 20000,
                step = 1,
                value = 100
              ),
              align = "center"
            )
          ),
          shinyWidgets::pickerInput(
            inputId = "mummi_adducts",
            choices = adducts$Name,
            label = "Adducts to look for:",
            selected = c("[M+H]1+", "[M+Na]1+", "[M+2H]2+", "[M+K]1+"),
            multiple = T,
            options = list(
              `actions-box` = TRUE
            )
          ),
          shiny::uiOutput("mummi_adducts_cats"),
          shiny::hr(),
          shiny::helpText("Include correlated m/z values?"),
          shinyWidgets::prettyToggle("enrich_use_corr", label_on = "yes", label_off = "no", value = F),
          shiny::conditionalPanel(
            "input.enrich_use_corr == true",
            shiny::selectInput("enrich_corr_method", "Calculation method:",
              choices = c(
                "pearson",
                "kendall",
                "spearman"
              ),
              selected = "pearson"
            ),
            shiny::numericInput("enrich_min_corr",
              label = "Minimum correlation:", value = 0.9, min = 0, max = 1
            )
          ),
          shiny::hr(),
          shinyWidgets::actionBttn(
            inputId = "do_enrich",
            label = "click to start enrichment analysis",
            style = "bordered",
            icon = icon("terminal"),
            size = "sm"
          )
        )
      ),
      mainPanel = mainPanel(fluidRow(
        align = "center",
        shinyBS::bsCollapse(
          id = "collapse_enrich", multiple = T,
          shinyBS::bsCollapsePanel(
            title = h2("plots"), value = "collapse_enrich_plots",
            shiny::tabsetPanel(
              id = "enrich_results",
              shiny::tabPanel(
                "overview",
                (shinyjqui::jqui_resizable(shiny::uiOutput("enrich_plot_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700))),
                shinyWidgets::switchInput("enrich_plot_mode", "Plot type:",
                  onLabel = "point", offLabel = "bar",
                  value = T
                )
              ),
              shiny::tabPanel(
                "selected pathway",
                shinyWidgets::switchInput("enrich_pathway_plot_mode", "Style:",
                  onLabel = "kegg",
                  offLabel = "generic", value = T
                ),
                shiny::selectizeInput("enrich_pathway_projection",
                  "Project analysis results onto pathway:",
                  choices = c("none available"),
                  selected = 1,
                  multiple = F
                ),
                shinyWidgets::switchInput(
                  inputId = "enrich_summ_adds",
                  label = "Summarize multiple adducts matching one node?",
                  value = F,
                  onLabel = "yes",
                  offLabel = "no"
                ),
                shinyWidgets::switchInput(
                  inputId = "enrich_map_nonsig",
                  label = "Map nonsignificant m/z?",
                  value = F,
                  onLabel = "yes",
                  offLabel = "no"
                ),
                shiny::conditionalPanel(
                  "input.enrich_summ_adds == true",
                  shiny::selectInput("enrich_summ_adds_method", "Summary method:",
                    choices = c(
                      "sum",
                      "mean",
                      "median",
                      "absmax"
                    ),
                    selected = "sum"
                  )
                ),
                shiny::actionButton("enrich_plot_pathway", "Plot pathway"),
                shiny::wellPanel(
                  id = "pathway_panel",
                  style = "overflow-y:scroll; overflow-x:scroll; max-height: 600px",
                  imageOutput("enrich_pathway")
                ),
                DT::dataTableOutput("enrich_curr_direction_tbl", width = "100%"),
                shiny::verbatimTextOutput("enrich_curr_direction_txt")
              )
            )
          ),
          shinyBS::bsCollapsePanel(
            title = h2("tables"), value = "collapse_enrich_tables",
            shiny::div(
              DT::dataTableOutput("enrich_tab",
                width = "100%"
              ),
              style = "font-size:80%"
            ),
            shiny::div(
              DT::dataTableOutput("enrich_pw_tab",
                width = "100%"
              ),
              style = "font-size:80%"
            ) # ,
            # shiny::actionButton("enrich_export_excel")
          )
        )
      ))
    )
  )
}

ui_feature_selection <- function() {
  shiny::tabPanel("feature selection",
    value = "featsel",
    shiny::fluidRow(
      align = "center",
      shinyBS::bsCollapse(
        id = "collapse_featsel",
        multiple = T,
        shinyBS::bsCollapsePanel(
          title = h2("settings"),
          value = "collapse_featsel_settings",
          shinyWidgets::actionBttn(
            inputId = "do_featsel",
            label = "click to start Boruta feature selection",
            style = "bordered",
            icon = icon("terminal"),
            size = "sm"
          )
        ),
        shinyBS::bsCollapsePanel(
          title = h2("tables"), value = "collapse_featsel_tables",
          shiny::div(
            DT::dataTableOutput("featsel_tab",
              width = "100%"
            ),
            style = "font-size:80%"
          )
        )
      )
    )
  )
}


ui_ml_configure <- function() {
  shiny::tabPanel("configure",
    icon = icon("cog"),
    shiny::fluidRow(
      align = "center",
      shiny::column(
        width = 11,
        shiny::tags$h3("Job name"),
        shiny::icon("signature", "fa-3x"),
        shiny::textInput("ml_name",
          label = "",
          value = stringi::stri_rand_strings(1,
            10,
            pattern = "[A-Za-z0-9]"
          )
        )
      )
    ),
    shiny::hr(),
    shiny::fluidRow(
      align = "center",
      shiny::column(
        width = 3,
        shiny::tags$h3("Input data"),
        shiny::icon("table", "fa-3x"),
        hr(),
        shinyBS::bsCollapse(
          shinyBS::bsCollapsePanel(
            "Used data",
            shiny::selectizeInput("ml_used_table",
              label = "Used data:",
              choices = list("orig", "prebatch", "norm", "pca"),
              selected = "norm"
            ),
            shiny::selectizeInput("ml_include_covars",
              label = "Use which metadata for prediction?",
              choices = c(" "), multiple = TRUE
            )
          ),
          shinyBS::bsCollapsePanel(
            "Used samples",
            shiny::selectizeInput("ml_samp_distr",
              label = "Reuse train/test sample distribution from previous experiment?",
              choices = c(" "), multiple = FALSE
            ),
            shiny::conditionalPanel(
              "input.ml_samp_distr == ' '",
              shiny::selectInput("ml_tr_te_subset",
                "Train or test on a metadata-based subset of the data?",
                choices = c("yes", "no")
              ),
              shiny::conditionalPanel(
                "input.ml_tr_te_subset == 'yes'",
                shiny::helpText("Train on subset:"),
                shiny::fluidRow(
                  sardine(shiny::textOutput("ml_train_ss")),
                  shinyWidgets::circleButton(
                    inputId = "reset_ml_train",
                    icon = shiny::icon("undo"),
                    size = "xs",
                    status = "warning"
                  ), br(),
                  shinyWidgets::circleButton("ml_train_ss",
                    label = ":",
                    icon = shiny::icon("filter")
                  )
                ),
                shiny::helpText("Test on subset:"),
                shiny::fluidRow(
                  sardine(shiny::textOutput("ml_test_ss")),
                  shinyWidgets::circleButton(
                    inputId = "reset_ml_test",
                    icon = shiny::icon("undo"),
                    size = "xs",
                    status = "warning"
                  ), br(),
                  shinyWidgets::circleButton("ml_test_ss",
                    icon = shiny::icon("filter")
                  )
                )
              ),
              shiny::conditionalPanel(
                "input.ml_tr_te_subset == 'no'",
                shiny::sliderInput("ml_train_perc",
                  label = "Percentage of samples in training",
                  min = 1,
                  max = 100,
                  step = 1,
                  value = 80,
                  post = "%"
                )
              )
            )
          ),
          shinyBS::bsCollapsePanel(
            "Used m/z",
            shiny::selectizeInput("ml_specific_mzs", label = "Analysis-based top hits", choices = c(" "), multiple = FALSE),
            shiny::conditionalPanel(
              "input.ml_specific_mzs != 'no'",
              shiny::conditionalPanel(
                "input.ml_specific_mzs == 'manual'",
                div(shinyWidgets::pickerInput(
                  inputId = "ml_mzs",
                  label = div(icon("search"),
                    style = "font-size: xx-large;margin-top: -30px;color: black;-webkit-text-fill-color: white;-webkit-text-stroke-width: 1.5px;-webkit-text-stroke-color: #DFDCDC"
                  ),
                  choices = "fa-cat",
                  multiple = T,
                  choicesOpt = list(
                    subtext = "",
                    style = "text-align:center;"
                  ),
                  options = list(
                    `live-search` = TRUE,
                    size = 10
                  )
                ), class = "mzpicker")
              ),
              shiny::conditionalPanel(
                "input.ml_specific_mzs != 'manual'",
                shiny::sliderInput("ml_mzs_topn",
                  label = "Use top:", min = 1, max = 200, value = 20
                ),
                shiny::selectInput("ml_mzs_ordering", "Ordering method:",
                  choices = list(
                    "high ranking first" = "highfirst",
                    "low ranking first" = "lowfirst",
                    "randomize" = "random"
                  ), selected = "highfirst"
                ),
                shiny::selectInput("ml_mzs_handling", "Handling method:",
                  choices = list(
                    "remove", "permute", "keep"
                  ), selected = "keep"
                ),
                shiny::textOutput("ml_specific_mzs_sigcount")
              )
            )
          )
        )
      ),
      shiny::column(
        width = 1,
        shiny::icon("caret-right", "fa-5x vertical-center")
      ),
      shiny::column(
        width = 3,
        shiny::tags$h3("Data adjustment"),
        shiny::icon("shower", "fa-3x"),
        hr(),
        shinyBS::bsCollapse(shinyBS::bsCollapsePanel(
          title = "Batch correction",
          value = "ml_batch_panel",
          shiny::selectizeInput("ml_batch_covars",
            label = "Batch metadata column:",
            choices = c(" "),
            multiple = F
          ),
          shiny::helpText("Balance classes per batch factor?"),
          shinyWidgets::switchInput(
            inputId = "ml_batch_balance",
            value = FALSE,
            onLabel = "yes",
            offLabel = "no",
            size = "small" # "<div class=\"fa-flip-vertical\"><i class=\"fas fa-chart-bar fa-rotate-90\"></i></div>"
          ),
          shiny::helpText("Use seperate folds for batch members?"),
          shinyWidgets::switchInput(
            inputId = "ml_covar_fold_seperate",
            value = FALSE,
            onLabel = "yes",
            offLabel = "no",
            size = "small" # "<div class=\"fa-flip-vertical\"><i class=\"fas fa-chart-bar fa-rotate-90\"></i></div>"
          ),
          shiny::helpText("If batch present, make all batches same size in TOTAL? Otherwise each batch will have the right percentage train/test but one batch may be larger than another in TOTAL.:"),
          shinyWidgets::switchInput(
            inputId = "ml_batch_size_sampling",
            size = "mini",
            onLabel = "Yes",
            offLabel = "No",
            value = FALSE
          ),
          shiny::conditionalPanel(
            "input.ml_batch_size_sampling == true",
            shiny::numericInput("ml_groupsize",
              label = "Size per experimental group (auto: 0):",
              value = 0
            )
          )
        )),
        shinyBS::bsCollapse(shinyBS::bsCollapsePanel(
          title = "Class balancing",
          shinyWidgets::radioGroupButtons(
            inputId = "ml_sampling",
            label = "How to balance class labels?:",
            choices = c(
              `<i class='fa fa-arrow-down'></i> downsample` = "down",
              `ROSE` = "rose",
              `SMOTE` = "smote",
              `ADASYN` = "adasyn",
              `don't` = "none",
              `upsample <i class='fa fa-arrow-up'></i>` = "up"
            ),
            justified = FALSE,
            selected = "none"
          )
        )),
        shinyBS::bsCollapse(shinyBS::bsCollapsePanel(
          title = "Miscellanaous",
          shiny::selectizeInput("ml_preproc",
            label = "Data preprocessing",
            choices = c("center", "scale"),
            selected = c(), multiple = T
          ),
          shinyWidgets::switchInput(
            inputId = "ml_pca_corr",
            value = FALSE,
            label = "PCA correction",
            onLabel = "yes",
            offLabel = "no",
            size = "small" # "<div class=\"fa-flip-vertical\"><i class=\"fas fa-chart-bar fa-rotate-90\"></i></div>"
          ),
          shiny::conditionalPanel(
            "input.ml_pca_corr == true",
            sliderInput("ml_keep_pcs",
              label = "Keep variance in PCs:",
              min = 0,
              max = 300, value = c(10, 50)
            )
          )
        ))
      ),
      shiny::column(
        width = 1,
        shiny::icon("caret-right", "fa-5x vertical-center")
      ),
      shiny::column(
        width = 3,
        shiny::tags$h3("Model settings"),
        shiny::numericInput("ml_n_repls",
          "Number of regular models built:",
          value = 1, min = 1, max = NA
        ),
        shiny::icon("sliders-h", "fa-3x"),
        hr(),
        shinyBS::bsCollapse(shinyBS::bsCollapsePanel(
          title = "Algorithm",
          shiny::selectizeInput("ml_method",
            label = "Used algorithm",
            selected = "glmnet",
            choices = as.list(gbl$constants$ml.models),
            multiple = F
          ),
          shiny::div(shiny::uiOutput("ml_params"), style = "font-size:60%")
        )),
        shinyBS::bsCollapse(shinyBS::bsCollapsePanel(
          title = "Cross validation",
          shiny::selectizeInput("ml_perf_metr",
            label = "Performance metric",
            choices = c(
              "boot", "boot632", "optimism_boot",
              "boot_all", "cv", "repeatedcv",
              "LOOCV", "LGOCV", "none", "oob",
              "timeslice", "adaptive_cv",
              "adaptive_boot", "adaptive_LGOCV"
            ),
            multiple = F, selected = "cv"
          ),
          shiny::selectizeInput("ml_folds",
            label = "Fold CV", choices = c(
              "5",
              "10",
              "20",
              "50",
              "LOOCV"
            ),
            multiple = F
          )
        )),
        shinyBS::bsCollapse(shinyBS::bsCollapsePanel(
          title = "Significance",
          shiny::helpText("Calculate p-value using shuffled labels?"),
          shinyWidgets::switchInput(
            inputId = "ml_label_shuffle",
            size = "mini",
            onLabel = "Yes",
            offLabel = "No",
            value = FALSE
          ),
          shiny::conditionalPanel(
            "input.ml_label_shuffle == true",
            shiny::numericInput("ml_n_shufflings",
              "Number of randomized models:",
              value = 1, min = 1, max = NA
            ),
            shiny::helpText("Shuffle train or testing labels?"),
            shinyWidgets::switchInput(
              inputId = "ml_shuffle_mode",
              value = TRUE,
              onLabel = "train",
              offLabel = "test",
              size = "small"
            )
          )
        ))
      )
    ),
    shiny::hr(),
    shiny::fluidRow(
      align = "center",
      shiny::column(
        width = 11,
        shiny::uiOutput("ml_slurm_mem_ui"),
        shinyWidgets::actionBttn(
          inputId = "queue_ml",
          label = "add to queue",
          style = "bordered",
          icon = icon("plus"),
          size = "sm"
        ),
        shinyWidgets::actionBttn(
          inputId = "queue_ml_del",
          label = "remove from queue",
          style = "bordered",
          icon = icon("minus"),
          size = "sm"
        )
      )
    ),
    br(),
    shiny::fluidRow(
      align = "center",
      shiny::helpText("Recycle datasets with identical preprocessing settings?"),
      shinyWidgets::switchInput(
        inputId = "ml_resource_friendly",
        size = "mini",
        onLabel = "Yes",
        offLabel = "No",
        value = F
      ),
      shiny::column(
        width = 11,
        shinyWidgets::actionBttn(
          inputId = "do_ml",
          label = "start ML queue",
          style = "bordered",
          icon = icon("terminal"),
          size = "l"
        )
      )
    ),
    shiny::hr(),
    shiny::fluidRow(
      align = "center",
      shiny::column(
        width = 6,
        shiny::tags$h3("Job queue"),
        shiny::icon("list-ol", "fa-3x"),
        br(), br(),
        shinyBS::bsCollapse(shinyBS::bsCollapsePanel(
          title = "Advanced",
          value = "ml_queue_panel",
          shiny::textInput("queue_ml_name",
            label = "Queue (file) name:",
            value = "my_queue"
          ),
          shinyWidgets::actionBttn(
            inputId = "queue_ml_save",
            label = "save queue",
            style = "bordered",
            icon = icon("terminal"),
            size = "sm"
          ),
          shiny::fileInput(
            inputId = "queue_ml_load",
            label = "load queue"
          )
        )),
        DT::DTOutput("ml_queue_all")
      ),
      shiny::column(
        width = 5,
        shiny::tags$h3("Job settings"),
        shiny::icon("cogs", "fa-3x"),
        shiny::br(), shiny::br(),
        DT::DTOutput("ml_queue_sel")
      )
    )
  )
}

ui_ml_results <- function() {
  shiny::tabPanel("results",
    value = "res", icon = icon("chart-area"),
    shiny::fluidRow(
      align = "center",
      shiny::selectizeInput("show_which_ml",
        label = "Plot which model?",
        choices = c(" "),
        multiple = T
      ),
      shiny::actionButton("clear_ml_runs", "clear results", icon = icon("broom"))
    ),
    shiny::tabsetPanel(
      id = "ml_results",
      shiny::tabPanel(
        title = "curves", value = "roc", icon = shiny::icon("chart-area"),
        shiny::fluidRow(align = "center", (shinyjqui::jqui_resizable(shiny::uiOutput("ml_roc_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700)))),
        shiny::fluidRow(
          align = "center",
          shiny::column(
            4,
            shiny::selectizeInput("ml_plot_posclass",
              label = "Positive class",
              choices = c("placeholder"),
              selected = c("placeholder")
            )
          ),
          shiny::column(
            3,
            shiny::selectizeInput("ml_plot_x",
              "Metric for x axis:",
              choices = list(
                "False positive rate" = "fpr",
                "True positive rate" = "tpr",
                "Recall" = "recall",
                "Precision" = "precision",
                "Cutoff" = "cutoff"
              ),
              selected = "fpr"
            )
          ),
          shiny::column(
            3,
            shiny::selectizeInput("ml_plot_y",
              "Metric for y axis:",
              choices = list(
                "False positive rate" = "fpr",
                "True positive rate" = "tpr",
                "Recall" = "recall",
                "Precision" = "precision",
                "Cutoff" = "cutoff"
              ),
              selected = "tpr"
            )
          ),
          column(2, actionButton(inputId = "reload_ml_stats", "Go"))
        ),
        fluidRow(
          align = "center",
          shiny::selectInput("ml_plot_facet",
            label = "Facet ROC per metadata?",
            choices = c("don't facet"),
            multiple = F
          )
        ),
        fluidRow(
          shiny::column(6, shiny::div(
            DT::dataTableOutput("ml_overview_tab",
              width = "100%"
            ),
            style = "font-size:80%"
          )),
          shiny::column(6, shiny::plotOutput("conf_matr_plot"))
        )
      ),
      shiny::tabPanel("importance",
        value = "bar", icon = shiny::icon("star"),
        shiny::fluidRow(
          align = "center",
          (shinyjqui::jqui_resizable(shiny::uiOutput("ml_bar_wrap"), options = list(handles = "se", minHeight = 560, minWidth = 700))),
          shiny::sliderInput("ml_topn",
            label = "Show top:",
            min = 10,
            max = 200,
            step = 10,
            value = 20
          ),
          shiny::div(
            DT::dataTableOutput("ml_importance_tab",
              width = "100%"
            ),
            style = "font-size:80%"
          )
        )
      ),
      shiny::tabPanel("used settings",
        value = "param", icon = shiny::icon("cog"),
        shiny::div(DT::dataTableOutput("ml_param_tab",
          width = "50%"
        ))
      )
    )
  )
}
