require(classyfireR)

doClassyFire <- function(smiles = c("[H][C@](N)(CC1=CN(C)C=N1)C(O)=O"), posts.per.minute = 5, max.wait = 20){
  # credits to https://github.com/cbroeckl/RAMClustR/blob/master/R/getClassyFire.R for the SMILES
  resobj <- list()
  params <- list(label = "MetaboShiny", query_input = c(smiles), 
                 query_type = "STRUCTURE")
  submit <- httr::POST("http://classyfire.wishartlab.com/queries", 
                       body = params, encode = "json", httr::accept_json(), 
                       httr::add_headers(`Content-Type` = "application/json"))
  Sys.sleep(1)
  query_id <- jsonlite::fromJSON(httr::content(submit, 
                                               "text"))
  # TODO: implement, needs a loop
  if(any(names(query_id) == "error")) {
    if(query_id$error == "Limit exceeded") {
      Sys.sleep(60)
      params <- list(label = "ramclustR", query_input = resobj$smiles[i], 
                     query_type = "STRUCTURE")
      submit <- httr::POST("http://classyfire.wishartlab.com/queries", 
                           body = params, encode = "json", httr::accept_json(), 
                           httr::add_headers(`Content-Type` = "application/json"))
      query_id <- jsonlite::fromJSON(httr::content(submit, 
                                                   "text"))
    }
  }
  
  
  if(any(names(query_id) == "status")) {
    if(query_id$status == "500") {
      cat(" failed", '\n')
      next
    }
  }
  
  out <- list()
  out$classification_status <- "not done"
  out$number_of_elements <- 0
  Sys.sleep(1)
  skiptonext <- FALSE
  time.a <- Sys.time()
  while (out$number_of_elements == 0) {
    Sys.sleep(1)
    
    if(!RCurl::url.exists(paste0("http://classyfire.wishartlab.com/queries/",  query = query_id$id, ".json"))) {
      cat(" not done", '\n')
      break
    }
    out <- tryCatch( 
      {
        out <- jsonlite::fromJSON(paste0("http://classyfire.wishartlab.com/queries/", query = query_id$id, ".json"))
      }, 
      error = function(y) {
        out <- list()
        out$classification_status <- "not done"
        out$number_of_elements <- 0
        out
      }
    )
    
    if (round(as.numeric(difftime(Sys.time(), 
                                  time.a, units = "secs")), 3) >= max.wait) {
      cat("timed out", "\n")
      skiptonext <- TRUE
      break
    }
  }
  
  if (out$number_of_elements == 0) {
    resobj$classyfire[i, ] <- c(resobj$inchikey[i], 
                                     rep(NA, 6))
    rm(out)
  } else {
    a <- out$entities$inchikey
    if(is.null(a)) {
      a <- resobj$inchikey[i]
    } else {
      a <- gsub("InChIKey=", "", a)
    }
    b <- out$entities$kingdom$name
    if (is.null(b)) {
      b <- NA
      c <- NA
      d <- NA 
      e <- NA
      f <- NA
      g <- NA
    } else {
      c <- if(length(out$entities$superclass)>1) {out$entities$superclass$name} else {c <- NA}
      if (is.null(c)) {c <- NA}
      d <- if(length(out$entities$class)>1) {out$entities$class$name} else {d <- NA}
      if (is.null(d)) {d <- NA}
      e <- if(length(out$entities$subclass)>1) {out$entities$subclass$name} else {e <- NA}
      if (is.null(e)) {e <- NA}
      f <- if(length(out$entities$direct_parent)>1) {out$entities$direct_parent$name} else {f <- NA}
      if (is.null(f)) {f <- NA}
      g <- if(nchar(out$entities$description)>1) {out$entities$description} else {g <- NA}
      if (is.null(g)) {g <- NA}
    }
    cat(" done", '\n')
    resobj$classyfire[i, ] <- c(a, b, c, 
                                     d, e, f, g)
    rm(out)
  }
  Sys.sleep(ceiling(60/posts.per.minute))  
}
