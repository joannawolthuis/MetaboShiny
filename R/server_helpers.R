# Internal helpers used by inst/server.R
# Keep functions side-effect free where possible.

ensure_dir <- function(path){
  if(!dir.exists(path)) dir.create(path, recursive = TRUE)
}

show_tabs <- function(tabset_id, tabs){
  for(tab in tabs) shiny::showTab(tabset_id, tab)
}

hide_tabs <- function(tabset_id, tabs){
  for(tab in tabs) shiny::hideTab(tabset_id, tab)
}

read_table_if_in_dir <- function(dir_path, filename){
  if(filename %in% basename(list.files(dir_path))){
    data.table::fread(file.path(dir_path, filename))
  }else{
    NULL
  }
}

render_combi_picker <- function(input_id, output_id, label){
  output[[output_id]] <- shiny::renderUI({
    if(!is.null(input[[input_id]])){
      anal = mSet$analSet[[input[[input_id]]]]
      target.mat = grep("\\.mat", names(anal), value = TRUE)
      choices = colnames(anal[[target.mat]])
      shiny::selectInput(paste0(input_id, "_var"), label = label, choices = choices, selected = 1)
    }else{
      list()
    }
  })
}

format_options_lines <- function(opts){
  paste(paste0(names(opts), " = ", unlist(opts, use.names = FALSE)), collapse = "\n")
}

default_options_template <- function(){
  list(
    db_dir = "$home/MetaboShiny/databases",
    work_dir = "$home/MetaboShiny/saves/admin",
    proj_name = "MY_METSHI",
    ppm = "2",
    packages_installed = "Y",
    font1 = "Dosis",
    font2 = "Dosis",
    font3 = "Open Sans",
    font4 = "Open Sans",
    col1 = "#1961AB",
    col2 = "#FFFFFF",
    col3 = "#FFFFFF",
    col4 = "#000000",
    size1 = "40",
    size2 = "20",
    size3 = "15",
    size4 = "11",
    taskbar_image = "metshi_logo.png",
    gtheme = "classic",
    gcols = "#1C1400&#FFE552&#D49C1A&#EBC173&#8A00ED&#00E0C2&#95C200&#FF6BE4&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF&#FFFFFF",
    gspec = "RdBu",
    gfont = "15",
    mode = "complete",
    cores = "1",
    apikey = "",
    dbfavs = "",
    omit_unknown = "yes",
    beep = "no"
  )
}

set_option_default <- function(opts, name, default){
  if(!(name %in% names(opts))){
    opts[[name]] <- default
    MetaboShiny::setOption(lcl$paths$opt.loc, name, default)
  }
  opts
}

update_theme_inputs <- function(session, opts){
  shiny::updateSelectInput(session, "ggplot_theme", selected = opts$gtheme)
  shiny::updateSelectInput(session, "color_ramp", selected = opts$gspec)
}

set_option_value <- function(name, value){
  MetaboShiny::setOption(lcl$paths$opt.loc, name, value)
}

render_text_outputs <- function(items){
  lapply(items, function(item){
    output[[item$name]] <- shiny::renderText(item$text)
  })
}

update_ml_method <- function(session, choices){
  shiny::updateSelectInput(session,
                           "ml_method",
                           selected = "rf",
                           choices = as.list(choices))
}

stop_session_cluster <- function(){
  try({
    shiny::showNotification("Stopping threads...")
    parallel::stopCluster(session_cl)
  })
}

start_session_cluster <- function(cores){
  shiny::showNotification("Starting new threads...")
  logfile <<- tempfile()
  #print(paste("Log file:", logfile))
  #if(file.exists(logfile)) file.remove(logfile)
  session_cl <<- parallel::makeCluster(cores,
                                       outfile="")
  #,setup_strategy = "sequential") # leave 1 core for general use and 1 core for shiny session
  # send specific functions/packages to other threads
  parallel::clusterEvalQ(session_cl, {
    library(data.table)
    library(iterators)
    library(MetaboShiny)
    library(MetaDBparse)
  })
}

render_image_outputs <- function(images, tmp){
  lapply(images, function(image){
    output[[image$name]] <- shiny::renderImage({
      image$path <- if(grepl("tmp", image$path)) gsub("tmp", tmp, image$path) else image$path
      filename <- normalizePath(image$path)
      list(src = filename)
    }, deleteFile = FALSE)
  })
}

update_theme_controls <- function(session, opts, count){
  lapply(seq_len(count), function(i){
    colourpicker::updateColourInput(session=session,
                                    inputId = paste0("bar.col.", i),
                                    value = opts[[paste0("col", i)]])
    shiny::updateTextInput(session=session,
                           inputId = paste0("font.", i),
                           value = opts[[paste0("font", i)]])
    shiny::updateSliderInput(session=session,
                             inputId = paste0("size.", i),
                             value = opts[[paste0("size", i)]])
  })
}

show_stat_panels <- function(stat_id, show_plots){
  collapse_id <- paste0("collapse_", stat_id)
  shinyjs::show(selector = paste0("div.panel[value=", collapse_id, "_tables]"))
  if(show_plots){
    shinyjs::show(selector = paste0("div.panel[value=", collapse_id, "_plots]"))
  }else{
    shinyjs::hide(selector = paste0("div.panel[value=", collapse_id, "_plots]"))
  }
  collapse_id
}

set_stat_collapse <- function(session, stat_id, panes){
  collapse_id <- paste0("collapse_", stat_id)
  panel_values <- paste0(collapse_id, panes)
  for(panel_value in panel_values){
    shinyBS::updateCollapse(session, collapse_id, open = panel_value)
  }
}

make_sel_adducts <- function(adducts, selected_cats){
  sel_adducts = lapply(selected_cats, function(colu){
    col = adducts[[colu]] == "v"
    col[is.na(col)] <- FALSE
    col
  })
  if(length(sel_adducts) == 1){
    sel_adducts[[1]]
  }else{
    sel_adducts = lapply(sel_adducts, as.list)
    sel_adduct_table = data.table::rbindlist(sel_adducts)
    sel_adduct_sums = colSums(sel_adduct_table)
    sel_adduct_sums > 0
  }
}

update_adducts_from_filter <- function(sel_adducts){
  if(!is.null(pieinfo)){
    adds_in_search = as.character(pieinfo$add$Var.1)
    for(mzMode in c("positive", "negative")){
      sel_adducts_mode = sel_adducts & adducts$Ion_mode == mzMode
      keep_sel_adducts <- intersect(adducts$Name[sel_adducts_mode], adds_in_search)
      result_filters$add[[mzMode]] <- keep_sel_adducts
      search$go <- TRUE
    }
  }
}
