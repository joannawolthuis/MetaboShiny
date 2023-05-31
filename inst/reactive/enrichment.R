shiny::observeEvent(input$enrich_plot_pathway, {
  curr_pathway = rownames(enrich$overview)[input$enrich_tab_rows_selected]
  pws_all = gbl$vectors$kegg_pathways
  kegg_pathway_match = which(curr_pathway == mSet$analSet$enrich$path.all$name)
  if(length(kegg_pathway_match) > 0 & nrow(enrich$current) > 0){
    pw.code = mSet$analSet$enrich$path.all$code[kegg_pathway_match]

    dt = data.frame(id = enrich$current$identifier,
                    mz = enrich$current$rn)  
    
    uniq_ids = unique(dt$id)
    cpds = unique(dt$id)
    
    vals_for_projection <- getAllHits(mSet, 
                                      input$enrich_pathway_projection,
                                      randomize = F)
    vals_for_projection$significant <- NULL

    analysis = gsub(" \\(.*$", "", input$enrich_pathway_projection)
    if(analysis %in% c("tt","fc","combi")){
      updir = switch(analysis,
                     tt = levels(mSet$dataSet$cls)[1],
                     fc = levels(mSet$dataSet$cls)[2],
                     combi = levels(mSet$dataSet$cls)[2]
                     )
    }else{
      updir = NA
    }
    
    mzs_withvals = merge(enrich$current[,c("rn", "identifier", "significant")], 
                         vals_for_projection, by.x="rn", by.y="m.z")
    
    if(input$enrich_map_nonsig){
      mzs_withvals$significant <- c("yes")
    }
    
    if(!("statistic" %in% names(mzs_withvals))){
      mzs_withvals$statistic = c(1)
    }
    
    multi.state.table <- data.table::rbindlist(lapply(cpds, 
                                                      function(cpd_id){
      mz_rows = mzs_withvals[identifier == cpd_id]
      print(mz_rows)
      if(nrow(mz_rows) > 0){
        mz_rows$stastistic <- ifelse(mz_rows$significant == 'yes', mz_rows$statistic, 0)
        
        all.dat = data.table::data.table(id = cpd_id,
                                         mz = mz_rows$rn,
                                         value = mz_rows$stastistic)
        reshape2::dcast(all.dat,
                        formula = id ~ mz, value.var = "value")  
      }else{
        data.table::data.table()
      }
    }), fill = TRUE)
    
    colnames(multi.state.table) = c("id", paste0("add", 1:(ncol(multi.state.table)-1)))
    for(i in 1:nrow(multi.state.table)){
      row = multi.state.table[i, 2:ncol(multi.state.table)]
      new.order = c(1, order(abs(row),decreasing=T, na.last = T) + 1)
      multi.state.table[i,] <- multi.state.table[i, ..new.order]  
    }
    
    non.zero.cols = sapply(1:ncol(multi.state.table), 
                           function(i) any(!is.na(multi.state.table[[i]])))
    multi.state.table = multi.state.table[, ..non.zero.cols]
    multi.state.table <- as.data.frame(multi.state.table)
    rownames(multi.state.table) = multi.state.table$id
    
    species = substr(mSet$analSet$enrich$path.lib, 0, 3)
    pw.code = gsub("^[a-z]+", "", pw.code)
    
    data(bods, package = "pathview")
    
    if(species == "map"){
      species = "ko"
    }
    
    is.multi = ncol(multi.state.table) > 2
    if(!is.multi) colnames(multi.state.table)[2] <- "single"
      
    tmpdir = if(interactive()){
      getwd()
    }else{
      tempdir()
    }
    
    if(input$enrich_summ_adds){
      myfun = switch(
        input$enrich_summ_adds_method,
        sum = sum,
        mean = mean,
        median = median,
        absmax = function(x, ...) abs(max(x, ...))
      )
      multi.state.table <- data.table::rbindlist(lapply(1:nrow(multi.state.table), 
                                                        function(i){
                                                          row = multi.state.table[i,-1]
                                                          row = row[!is.na(row)]
                                                          summ.row = myfun(row, na.rm=T)
                                                          data.frame(id = multi.state.table$id[[i]], 
                                                                                add1 = summ.row)
                                                        }))
    }else{
      multi.state.table = data.table::rbindlist(lapply(1:nrow(multi.state.table), function(i){
        row = multi.state.table[i,]
        numvals = row[,2:ncol(row)]
        numvals.nona = unlist(numvals[, sapply(1:ncol(numvals), function(i) !is.na(numvals[[i]]))])
        numvals.filled = rep(numvals.nona, length.out = ncol(row)-1)
        row[2:ncol(row)] <- numvals.filled[order(numvals.filled)]
        row
      }))  
    }
    
    multi.state.table <- as.data.frame(multi.state.table)
    rownames(multi.state.table) <- multi.state.table$id
    
    library(pathview)
    if(!file.exists(paste0(species, pw.code, ".qs"))){
      #pathview::download.kegg(species = species,
      #                        pathway.id = pw.code)  
      download.kegg.jw(species = species,
                       pathway.id = pw.code)
    }
    
    # === UP/DOWN overview ===
    sig = enrich$current$significant == "yes"
    mapper = unique(enrich$current[sig, c("compoundname",
                                          "identifier")])
    names_mapper = merge(multi.state.table,
                         mapper,
                         by.x = "id", 
                         by.y = "identifier")
    names_mapper$val_summary <- sapply(1:nrow(names_mapper), function(i){
        vals = names_mapper[i, grepl("add", colnames(names_mapper))]
        vals = vals[!is.na(vals)]
        myfun(vals)
    })
    
    if(!is.na(updir)){
      upname = if(!is.na(updir)) paste("up in", mSet$settings$exp.var, "=", updir)
      downname = if(!is.na(updir)) paste("down in", mSet$settings$exp.var, "=", updir)
    }else{
      upname = "up"
      downname = "down"
    }
    
    names_mapper$direction <- ifelse(names_mapper$val_summary > 0, upname, downname)
    names_mapper <- names_mapper[,c("compoundname", "direction")]
    enrich$curr_direction <- names_mapper
    output$enrich_curr_direction_tbl <- DT::renderDataTable(metshiTable(names_mapper))
    output$enrich_curr_direction_txt <- shiny::renderText(paste0(Hmisc::capitalize(upname),":", 
                                                                 paste0(names_mapper[names_mapper$direction == upname,]$compoundname, collapse=", "),
                                                                 "\n\n",
                                                                 Hmisc::capitalize(downname), ":", 
                                                                 paste0(names_mapper[names_mapper$direction == downname,]$compoundname, collapse = ", ")))
    # ========================
    
    pv.out <- pathview::pathview(cpd.data = multi.state.table[,-1, drop=F],  
                                 cpd.idtype = "kegg",
                                 pathway.id = pw.code, 
                                 species = species,
                                 out.suffix = "metaboshiny", 
                                 keys.align = "y",
                                 match.data = TRUE,
                                 kegg.native = if(input$enrich_pathway_plot_mode) T else F,
                                 multi.state = is.multi,
                                 same.layer = TRUE,
                                 low = list(gene = "green", cpd = "blue"), 
                                 mid = list(gene = "gray", cpd = "gray"), 
                                 high = list(gene = "red", cpd = "red"),
                                 kegg.dir = tmpdir,
                                 min.nnodes = 1,
                                 res = 300 #ppi
                                 )   
    
    if(input$enrich_summ_adds){
      fn.partial = paste0(species, pw.code, ".metaboshiny")
    }else{
      fn.add = ifelse(is.multi, "multi", "single")
      fn.partial = paste0(species, pw.code, ".metaboshiny.", fn.add)
    }
    
    if(!input$enrich_pathway_plot_mode){
      # convert to png
      pathway_png = magick::image_convert(image = magick::image_read_pdf(paste0(fn.partial, ".pdf"),), 
                                          format = "png")
      magick::image_write(pathway_png, path = paste0(fn.partial, ".png"))
    }
    
    fn = paste0(fn.partial, ".png")
    
    output$enrich_pathway <- shiny::renderImage({
     filename <- normalizePath(file.path('.',
                                         fn))
      
      list(src = filename)
    },deleteFile = T)
  }
})