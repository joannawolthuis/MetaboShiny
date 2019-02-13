# listener for all the file pickers, trigger a window if they are clicked
observe({
  shinyFileChoose(input, 'outlist_pos', roots=global$paths$volumes, filetypes=c('csv'))
  shinyFileChoose(input, 'outlist_neg', roots=global$paths$volumes, filetypes=c('csv'))
  shinyFileChoose(input, 'custom_db', roots=global$paths$volumes, filetypes=c('csv'))
  shinyFileChoose(input, 'metadata', roots=global$paths$volumes, filetypes=c('xls', 'xlsm', 'xlsx', 'csv'))
  shinyFileChoose(input, 'database', roots=global$paths$volumes, filetypes=c('sqlite3', 'db', 'sqlite'))
  shinyFileChoose(input, 'taskbar_image_path', roots=global$paths$volumes, filetypes=c('png', 'jpg', 'jpeg', 'bmp'))
  shinyFileChoose(input, 'custom_db_img_path', roots=global$paths$volumes, filetypes=c('png', 'jpg', 'jpeg', 'bmp'))
})

# observes if a new taskbar image is chosen by user
observe({
  # - - - - 
  if(!is.list(input$custom_db_img_path)) return() # if nothing is chosen, do nothing
  img_path <- parseFilePaths(global$paths$volumes, input$custom_db_img_path)$datapath
  new_path <- file.path(getwd(), "www", basename(img_path)) # set path to copy to
  global$paths$custom.db.path <<- new_path
  # copy image to the www folder
  if(img_path != new_path) file.copy(img_path, new_path, overwrite = T)
  # - - -
  # render taskbar image preview
  output$custom_db_img <- renderImage({
    list(src = new_path, 
         width = 150,
         height = 150,
         style = "background-image:linear-gradient(0deg, transparent 50%, #aaa 50%),linear-gradient(90deg, #aaa 50%, #ccc 50%);background-size:10px 10px,10px 10px;")
  }, deleteFile = FALSE)
})

# observes if a new taskbar image is chosen by user
observe({
  # - - - - 
  if(!is.list(input$taskbar_image_path)) return() # if nothing is chosen, do nothing
  img_path <- parseFilePaths(global$paths$volumes, input$taskbar_image_path)$datapath
  new_path <- file.path(getwd(), "www", basename(img_path)) # set path to copy to
  
  # copy image to the www folder
  if(img_path != new_path) file.copy(img_path, new_path, overwrite = T)
  # - - -
  # render taskbar image preview
  output$taskbar_image <- renderImage({
    list(src = new_path, 
         width = 120,
         height = 120,
         style = "background-image:linear-gradient(0deg, transparent 50%, #aaa 50%),linear-gradient(90deg, #aaa 50%, #ccc 50%);background-size:10px 10px,10px 10px;")
  }, deleteFile = FALSE)
  # change chosen taskbar image in user option file
  setOption('user_options.txt', 
            'taskbar_image', 
            basename(new_path))
})

# observes if user is choosing a different database storage folder
observe({  
  # trigger window
  shinyDirChoose(input, "get_db_dir", 
                 roots=global$paths$volumes, 
                 session = session)
  
  if(typeof(input$get_db_dir) != "list") return() # if nothing selected or done, ignore
  
  # parse the file path given based on the possible base folders (defined in global)
  given_dir <- parseDirPath(global$paths$volumes, 
                            input$get_db_dir)
  
  if(is.null(given_dir)) return()
  # change db storage directory in user options file
  setOption("user_options.txt", "db_dir", given_dir)
  
  # render current db location in text
  output$curr_db_dir <- renderText({getOptions('user_options.txt')$db_dir})
})

# see above, but for working directory. CSV/DB files with user data are stored here.
observe({
  shinyDirChoose(input, "get_work_dir",
                 roots = global$paths$volumes,
                 session = session)
  
  if(typeof(input$get_work_dir) != "list") return()
  
  given_dir <- parseDirPath(global$paths$volumes,
                            input$get_work_dir)
  if(is.null(given_dir)) return()
  setOption("user_options.txt", "work_dir", given_dir)
  
  output$curr_exp_dir <- renderText({getOptions('user_options.txt')$work_dir})
})

# triggers if user changes their current project name
observeEvent(input$set_proj_name, {
  proj_name <<- input$proj_name
  if(proj_name == "") return(NULL) # if empty, ignore
  # change path of current db in global
  global$paths$patdb <<- file.path(getOptions("user_options.txt")$work_dir, paste0(proj_name,".db", sep=""))
  # change project name in user options file
  setOption("user_options.txt", "proj_name", proj_name)
  # print the changed name in the UI
  output$proj_name <<- renderText(proj_name)
  # change path CSV should be / is saved to in session
  global$paths$csv_loc <<- file.path(getOptions("user_options.txt")$work_dir, paste0(getOptions("user_options.txt")$proj_name,".csv"))
  
})