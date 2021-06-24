shiny::observeEvent(input$enrich_plot_pathway, {
  print("hello")
  curr_pathway = rownames(enrich$overview)[input$enrich_tab_rows_selected]
  pws_all = gbl$vectors$kegg_pathways
  kegg_pathway_match = which(curr_pathway == mSet$analSet$enrich$path.all$name)
  if(length(kegg_pathway_match) > 0 & nrow(enrich$current) > 0){
    pw.code = mSet$analSet$enrich$path.all$code[kegg_pathway_match]
    #species = stringr::str_split(mSet$analSet$enrich$path.lib, "_")[[1]][1]
    dt = data.frame(id = enrich$current$identifier,
                    mz = enrich$current$rn,
                    mark = ifelse(enrich$current$significant == "yes", 1, 0))  
    
    # fc = sapply(dt$mz, function(mz){
    #   mSet$analSet$fc$fc.log[names(mSet$analSet$fc$fc.log) == mz]
    # })
    # 
    # dt$fc = fc
    
    uniq_ids = unique(dt$id)
    multi.state.rows <- lapply(uniq_ids, function(cpd_id){
      all.dat = dt[dt$id == cpd_id,]
      #cast = reshape2::dcast(all.dat, id ~ mz, value.var = "mark", fun.aggregate = c, fill = 0)
      #colnames(cast) = paste0("add", 1:ncol(cast))
      #as.data.frame(cast)
      if(any(all.dat$mark == 1)){
        cast = data.frame(id = cpd_id, add1 = 1)
      }else{
        cast = data.frame(id = cpd_id, add1 = 0)
      }
      cast
    })
    
    multi.state.table = as.data.frame(data.table::rbindlist(multi.state.rows, fill = T))
    rownames(multi.state.table) = uniq_ids
    
    for.pathway = dt$mark
    names(for.pathway) = dt$id
    
    #old_dir = getwd()
    #setwd(tempdir())
    
    species = substr(mSet$analSet$enrich$path.lib, 0, 3)
    pw.code = gsub("^[a-z]+", "", pw.code)
    
    pv.out <- pathview::pathview(cpd.data = multi.state.table, 
                                 pathway.id = pw.code, 
                                 # limit = c(min(for.pathway),
                                 #           max(for.pathway)),
                                 species = species, 
                                 out.suffix = "metaboshiny", 
                                 keys.align = "y", 
                                 kegg.native = TRUE, 
                                 match.data = FALSE,
                                 multi.state = FALSE,#TRUE,
                                 same.layer = TRUE,
                                 low = list(gene = "green", cpd = "blue"), 
                                 mid = list(gene = "gray", cpd = "gray"), 
                                 high = list(gene = "red", cpd = "red"))
    
    #fn = paste0(species, pw.code, ".metaboshiny.multi.png")
    fn = paste0(species, pw.code, ".metaboshiny.add1.png")
    output$enrich_pathway <- shiny::renderPlot({
     filename <- normalizePath(file.path('.',
                                         fn))
      
      img <- magick::image_read(fn)
      plot(img)
    })
    #setwd(old_dir)
  }
})