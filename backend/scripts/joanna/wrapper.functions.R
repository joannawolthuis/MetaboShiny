#' @export
get.csv <- 
  function(patdb, 
                    time.series = F, 
                    #exp.condition = "diet",
                    max.vals = -1,
                    group_adducts = T,
                    which_dbs = file.path(options$db_dir, "kegg.full.db"),
                    which_adducts = c("M+H", "M-H", "M"),
                    groupfac = "mz"
                    #,var_table = "setup",
                    #batches = NULL
                    ){
  
  library(data.table)
  
  
  # - - announce some stuff - -
  
  adducts <- paste(which_adducts, collapse = ", ")
  
  max.cols <- if(max.vals == -1) "unlimited" else max.vals + 3
  
  cat(gsubfn::fn$paste("Creating csv for metabolomics analysis with max $max.cols columns.
                - Grouping by $groupfac
                - Using adducts $adducts"))
  
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), patdb)
  
  # - - - - - - - - - - - - - - - - 
  
  if(group_adducts){
    z = get_all_matches(#exp.condition, 
                        pat.conn = conn,
                        which_dbs,
                        which_adducts,
                        groupfac
                        #,var_table,
                        #batches
                        )
  }else{
    query <- strwrap(gsubfn::fn$paste("select distinct d.*, s.*, b.*,
                                      i.mzmed as identifier,
                                      i.intensity
                                      from mzintensities i
                                      join individual_data d
                                      on i.filename = d.card_id
                                      join setup s on d.[Group] = s.[Group]
                                      join batchinfo b on b.sample = d.card_id
                                      group by d.card_id, b.batch, b.injection, i.mzmed, d.sampling_date"),
                     #group by d.card_id, 
                     #, 
                     width=10000,
                     simplify=TRUE)
    print(query)
    z = RSQLite::dbGetQuery(conn, query)
  }
  
  RSQLite::dbDisconnect(conn)
  
  cast.dt <- dcast.data.table(as.data.table(z), 
                              formula = ... ~ identifier,
                              fun.aggregate = sum, 
                              value.var = "intensity",verbose = T) # what to do w/ duplicates? 
  
  # - - - check the first 100 rows for variables - - -
  
  as.numi <- as.numeric(colnames(cast.dt)[1:100])
  
  exp.vars <- which(is.na(as.numi))
  
  # --- cast to right format ---
  small.set <- cast.dt[,..exp.vars,]
  
  #small.set.test <<- small.set
  #small.set <- small.set.test
  
  # --- name for metaboanalyst ---
  #colnames(small.set)[which(colnames(small.set) == "label")] <- "#"
  colnames(small.set)[which(colnames(small.set) == "card_id")] <- "Sample"
  colnames(small.set)[grep(x=colnames(small.set), pattern="sampling_date")] <- "Time"
  colnames(small.set) <- tolower(colnames(small.set))
  #colnames(small.set) <- gsub("(?<=\\b)([a-z])", "\\U\\1", colnames(small.set), perl=TRUE)
  small.set <- small.set[,which(unlist(lapply(small.set, function(x)!all(is.na(x))))),with=F]
  #small.set <- cbind(small.set, "$" = c(0))
  
  # --- rejoin w/ rest ---
  
  small.set <- cbind(small.set, cast.dt[,-exp.vars, with=FALSE])
  
  # - - - fix QCs - - - 
  
  qc.locs <- which(small.set$Sample == "QC")
  
  small.set$Sample[qc.locs] <- paste0("QC", 1:length(qc.locs))
  
  # - - make time series if necessary (this factorizes sampling date) - -
  
  if(time.series){
      print("time series")
      small.set$time <- as.numeric(as.factor(as.Date(small.set$time)))
      small.set$sample <- small.set$animal_internal_id
  }
  
  # - - measure file size - -
  
  size <- object.size(small.set)
  
  print(small.set[1:5,1:10])
  
  cat(paste("Resulting file will be approximately "))
  
  print(size, units = "MB")
  
  # - - - return - - -
  
  return(small.set)
}