#' @export
make.metshi.csv <-
  function(patdb,
           csv = gsub(patdb, 
                      pattern = "\\.db", 
                      replacement = "\\.csv"),
           max.vals = -1,
           group_adducts = F,
           which_dbs = c(),
           which_adducts = c(),
           groupfac = "mz"
  ){
    
    conn <- RSQLite::dbConnect(RSQLite::SQLite(), normalizePath(patdb))
    
    cat("Checking for mismatches between peak tables and metadata... \n")
    
    fn_meta <- RSQLite::dbGetQuery(conn, "SELECT DISTINCT card_id FROM individual_data")[,1]
    fn_int <- RSQLite::dbGetQuery(conn, "SELECT DISTINCT filename FROM mzintensities")[,1]
    
    cat(paste0("-- in peaklist, not in metadata: --- \n", 
               paste0(setdiff(fn_int,
                              fn_meta), 
                      collapse=", "), 
               "\n"))
    cat(paste0("-- in metadata, not in peaklist: --- \n", 
               paste0(setdiff(fn_meta,
                              fn_int), 
                      collapse=", "), 
               "\n\n"))
    
    if(DBI::dbExistsTable(conn, "batchinfo")){
      query <- strwrap(gsubfn::fn$paste("select distinct d.card_id as sample, d.sampling_date as time, d.*, b.injection
                                        from mzintensities i
                                        join individual_data d
                                        on i.filename = d.card_id
                                        join setup s on d.[Group] = s.[Group]
                                        join batchinfo b on b.sample = d.card_id"),
                       width=10000,
                       simplify=TRUE)
    }else{
      query <- strwrap(gsubfn::fn$paste("select distinct d.card_id as sample, d.sampling_date as time, d.*, s.*
                                        from mzintensities i
                                        join individual_data d
                                        on i.filename = d.card_id
                                        join setup s on d.[Group] = s.[Group]"),
                       width=10000,
                       simplify=TRUE)
    }
    
    RSQLite::dbExecute(conn, "PRAGMA journal_mode=WAL;")
    RSQLite::dbExecute(conn, "CREATE INDEX IF NOT EXISTS filenames ON mzintensities(filename)")
    all_mz = RSQLite::dbGetQuery(conn, "SELECT DISTINCT mzmed FROM mzvals")[,1]
    
    RSQLite::dbDisconnect(conn)
    
    # write rows to csv
    pbapply::pblapply(fn_meta, 
                      #cl = session_cl, 
                      function(filename){
      # connect
      conn <- RSQLite::dbConnect(RSQLite::SQLite(), normalizePath(patdb))
      
      # adjust query
      query_add = gsubfn::fn$paste(" WHERE i.filename = '$filename'")
      
      # get results for sample
      z.meta = data.table::as.data.table(RSQLite::dbGetQuery(conn, paste0(query, query_add)))
      
      if(nrow(z.meta)==0) return(NA)
      
      z.meta = z.meta[,-c("card_id", "sampling_date")]
      colnames(z.meta) <- tolower(colnames(z.meta))
      z.int = data.table::as.data.table(RSQLite::dbGetQuery(conn, paste0("SELECT DISTINCT 
                                                i.mzmed as identifier,
                                                i.intensity
                                                FROM mzintensities i", query_add)))
      
      if(nrow(z.int)==0) return(NA)
      
      missing_mz <- setdiff(all_mz, z.int$identifier)
      
      # cast to wide
      cast.dt <- data.table::dcast.data.table(z.int,
                                  formula = ... ~ identifier,
                                  fun.aggregate = sum,
                                  value.var = "intensity")
      
      complete = as.numeric(cast.dt[1,])
      names(complete) = colnames(cast.dt)
      
      missing = rep(NA, length(missing_mz))
      names(missing) <- missing_mz
      
      complete.row = c(complete[-1], missing)
      reordered <- order(as.numeric(names(complete.row)))
      
      RSQLite::dbDisconnect(conn)
      
      # write
      data.table::fwrite(c(z.meta, "."=NA, complete.row[reordered]), 
             file = csv,
             append = T)
    })
    
    # - - measure file size - -
    
    disk_size = file.info(csv)$size
    size <- utils:::format.object_size(disk_size, "Mb")
    cat(paste("... Resulting file is approximately"),size,"...")
    }
