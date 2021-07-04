shiny::observeEvent(input$enrich_plot_pathway, {
  curr_pathway = rownames(enrich$overview)[input$enrich_tab_rows_selected]
  pws_all = gbl$vectors$kegg_pathways
  kegg_pathway_match = which(curr_pathway == mSet$analSet$enrich$path.all$name)
  if(length(kegg_pathway_match) > 0 & nrow(enrich$current) > 0){
    pw.code = mSet$analSet$enrich$path.all$code[kegg_pathway_match]

    dt = data.frame(id = enrich$current$identifier,
                    mz = enrich$current$rn)  
    
    uniq_ids = unique(dt$id)
    grey_nonsig = FALSE #any(enrich$current$significant == "yes")
    
    cpd2exp = mSet$analSet$enrich$cpd.value
    multi.state.table <- data.table::rbindlist(lapply(1:length(cpd2exp), function(i){
      cpd_id = names(cpd2exp)[i]
      value = cpd2exp[[i]]
      all.dat = data.table::data.table(id = cpd_id,
                                       value = value)
      reshape2::dcast(all.dat,formula = id ~ value)
    }), fill = TRUE)
    
    colnames(multi.state.table) = c("id", paste0("add", 1:(ncol(multi.state.table)-1)))
    for(i in 1:nrow(multi.state.table)){
      row = multi.state.table[i, 2:ncol(multi.state.table)]
      new.order = c(1, order(row, na.last = T)+1)
      multi.state.table[i,] <- multi.state.table[i, ..new.order]  
    }
    
    non.zero.cols = sapply(1:ncol(multi.state.table), function(i) any(!is.na(multi.state.table[[i]])))
    multi.state.table = multi.state.table[, ..non.zero.cols]
    multi.state.table <- as.data.frame(multi.state.table)
    rownames(multi.state.table) = names(cpd2exp)
    
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
    
    print(multi.state.table)
    pv.out <- pathview::pathview(cpd.data = multi.state.table, 
                                 #na.col = "red",
                                 pathway.id = pw.code, 
                                 species = species,
                                 out.suffix = "metaboshiny", 
                                 keys.align = "y",
                                 match.data = TRUE,
                                 multi.state = is.multi,
                                 same.layer = TRUE,
                                 low = list(gene = "green", cpd = "blue"), 
                                 mid = list(gene = "gray", cpd = "gray"), 
                                 high = list(gene = "red", cpd = "red"),
                                 kegg.dir = tmpdir,
                                 min.nnodes = 1)
    
    fn.add = ifelse(is.multi, "multi", "single")
    fn = paste0(species, pw.code, ".metaboshiny.", fn.add, ".png")
    
    
    output$enrich_pathway <- shiny::renderImage({
     filename <- normalizePath(file.path('.',
                                         fn))
      
      list(src = filename)
    },deleteFile = T)
    #setwd(old_dir)
  }
})