shiny::observe({
  
  if(is.null(mSetter$do)){
    
    NULL # if not subsetting anything, nevermind
    
  }else if(!is.null(mSetter$do)){
    
    if(!is.null(mSet)){
      success = F
      try({
        updated <- doUpdate(mSet, 
                            lcl, 
                            input, 
                            do = mSetter$do)
        mSet = updated$mSet
        lcl = updated$lcl
        # ---
        success = T
      })
      
      if(success){
        if(is.ordered.mSet(updated$mSet)){
          msg = "mSet class label order still correct! :)"
          try({
            shiny::showNotification(msg) 
          })
          print(msg)
          mSet <<- mSet
          lcl$has_changed <<- TRUE
          uimanager$refresh <- c("general", "ml")
        }else{
          msg = "mSet class label order incorrect! Restoring... :("
          try({
            shiny::showNotification(msg)
          })
          print(msg)
        }
      }else{
        metshiAlert("Failed! Restoring old mSet...")
      }
      mSetter$do <- NULL
    }
  }
})