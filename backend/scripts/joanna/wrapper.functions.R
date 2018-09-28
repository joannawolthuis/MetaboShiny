#' @export
get.csv <- function(patdb, 
                    time.series = T, 
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
  
  #batches <- batches[-which(batches == "run")]
  #if(exp.condition == "run") var_table <- "run"
  
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
                                      group by d.card_id, b.batch, b.injection, i.mzmed"),
                     #group by d.card_id, 
                     #d.sampling_date, 
                     width=10000,
                     simplify=TRUE)
    z = RSQLite::dbGetQuery(conn, query)
  }
  RSQLite::dbDisconnect(conn)
  
  #nvars = ncol(z) - 2
  
  cast.dt <- dcast.data.table(as.data.table(z), 
                              formula = ... ~ identifier,
                              fun.aggregate = sum, 
                              value.var = "intensity",verbose = T) # what to do w/ duplicates? 
  
  # - - - check the first 100 rows for variables (if need more..)
  
  as.numi <- as.numeric(colnames(cast.dt)[1:100])
  
  exp.vars <- which(is.na(as.numi))
  
  # --- cast to right format ---
  small.set <- cast.dt[,..exp.vars,]
  
  # --- name for metaboanalyst ---
  colnames(small.set)[which(colnames(small.set) == "label")] <- "#"
  colnames(small.set)[which(colnames(small.set) == "card_id")] <- "Sample"
  colnames(small.set)[which(colnames(small.set) == "sampling_date")] <- "Time"
  names(small.set) <- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(colnames(small.set)), perl=TRUE)
  small.set[,which(unlist(lapply(small.set, function(x)!all(is.na(x))))),with=F]
  #small.set <- cbind(small.set, "$" = c(0))
  
  # --- rejoin w/ rest ---
  print(head(small.set))
  
  small.set <- cbind(small.set, cast.dt[,-exp.vars, with=FALSE])
  
  # - - - fix QCs - - - 
  
  qc.locs <- which(small.set$Sample == "QC")
  
  small.set$Sample[qc.locs] <- paste0("QC", 1:length(qc.locs))
  
  # - - make time series if necessary (this factorizes sampling date) - -
  if(time.series){
      small.set$Time <- as.numeric(as.factor(as.Date(small.set$Time)))
      small.set$Sample <- paste(small.set$Sample, as.character(small.set$Time), sep="_T")
  }
  # - - measure file size - -
  
  size <- object.size(small.set)
  
  print(small.set[1:5,1:20])
  
  cat(paste("Resulting file will be approximately "))
  
  print(size, units = "MB")
  
  # - - - return - - -
  
  return(small.set)
}