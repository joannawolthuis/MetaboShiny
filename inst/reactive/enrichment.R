shiny::observeEvent(input$enrich_plot_pathway, {
  curr_pathway = rownames(enrich$overview)[input$enrich_tab_rows_selected]
  pws_all = gbl$vectors$kegg_pathways
  kegg_pathway_match = which(curr_pathway == mSet$analSet$enrich$path.all$name)
  if(length(kegg_pathway_match) > 0 & nrow(enrich$current) > 0){
    pw.code = mSet$analSet$enrich$path.all$code[kegg_pathway_match]

    dt = data.frame(id = enrich$current$identifier,
                    mz = enrich$current$rn)  
    
    uniq_ids = unique(dt$id)
    grey_nonsig = any(enrich$current$significant == "yes") & input$enrich_sig_only
    cpds = unique(dt$id)
    
    mzs_withvals = merge(enrich$current[,c("rn", "identifier")], 
                         mSet$analSet$enrich$orig.input, by.x="rn", by.y="m.z")
    
    if(!("statistic" %in% names(mzs_withvals))){
      mzs_withvals$statistic = c(1)
    }
    
    multi.state.table <- data.table::rbindlist(lapply(cpds, 
                                                      function(cpd_id){
      mz_rows = mzs_withvals[identifier == cpd_id]
      if(nrow(mz_rows) > 0){
        if(grey_nonsig){
          mz_rows$stastistic <- ifelse(mz_rows$significant, 1, 0)
        }else{
          mz_rows$stastistic <- ifelse(mz_rows$significant, mz_rows$statistic, 0)
        }
        all.dat = data.table::data.table(id = cpd_id,
                                         value = mz_rows$stastistic)
        reshape2::dcast(unique(all.dat),
                        formula = id ~ value)  
      }else{
        data.table::data.table()
      }
    }), fill = TRUE)
    
    colnames(multi.state.table) = c("id", paste0("add", 1:(ncol(multi.state.table)-1)))
    for(i in 1:nrow(multi.state.table)){
      row = multi.state.table[i, 2:ncol(multi.state.table)]
      new.order = c(1, order(abs(row),decreasing=T, na.last = T)+1)
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
    
    multi.state.table = data.table::rbindlist(lapply(1:nrow(multi.state.table), function(i){
      row = multi.state.table[i,]
      numvals = row[,2:ncol(row)]
      numvals.nona = unlist(numvals[, sapply(1:ncol(numvals), function(i) !is.na(numvals[[i]]))])
      numvals.filled = rep(numvals.nona, length.out = ncol(row)-1)
      row[2:ncol(row)] <- numvals.filled[order(numvals.filled)]
      row
    }))
    
    multi.state.table <- as.data.frame(multi.state.table)
    rownames(multi.state.table) <- multi.state.table$id
    
    library(pathview)
    pv.out <- pathview::pathview(cpd.data = multi.state.table[,-1], 
                                 #na.col = "black",
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
                                 min.nnodes = 1)   
    
    fn.add = ifelse(is.multi, "multi", "single")
    fn.partial = paste0(species, pw.code, ".metaboshiny.", fn.add)
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
    #setwd(old_dir)
  }
})