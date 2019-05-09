library(RSQLite)
# creates connection to SQLite db, if not exists - creates one
db <- dbConnect(SQLite(), dbname="users.db")

# create table users where logins and passwords are stored
dbSendQuery(conn = db,
            "DROP TABLE users")
dbSendQuery(conn = db,
            "CREATE TABLE users
            (username TEXT,
            password TEXT,
            role TEXT)")

# placeholder users
users = data.table(
  username = c("jobylocal", "jobydocker", paste0("advomics", c(1:30))),
  password = c("lol", "lol", rep("kirakira", 30)),
  role = c("admin", "admin", rep("noob", 30))
)

# insert some initial data to work with
DBI::dbWriteTable(db, "users", users, append=T)

DBI::dbDisconnect(db)

input <- list(username="")
local <- list()

for(i in 3:32){
  
  input$username <- users$username[i]
  print(input$username)
  
  # - - -
  
  runmode <- 'local'
  
  work_dir <- switch(runmode,
                     docker = "/userfiles/saves",
                     local = "~/MetaboShiny/saves")
  
  dbdir <- switch(runmode,
                  docker = "/userfiles/databases",
                  local = "~/MetaboShiny/databases")
  
  # check if user folder exists, otherwise make it
  userfolder = file.path(work_dir, input$username)
  
  if(!dir.exists(userfolder)){
    dir.create(userfolder)
    print("copying files...")
    cpfiles <- list.files("~/MetaboShiny/saves/jobylocal", pattern = "ADV")
    for(f in cpfiles){
      path.before <- file.path("~/MetaboShiny/saves/jobylocal", f)
      path.after <- file.path(userfolder, f)
      file.copy(path.before, path.after)
    }
  }
  
  username = input$username
  
  # check if opts file exists, otherwise make it with the proper files
  opt.loc <- file.path(userfolder, "options.txt")
  work_dir <- userfolder
  db_dir <- dbdir
  
  if(!file.exists(opt.loc)){
    contents = gsubfn::fn$paste('db_dir = /userfiles/databases
work_dir = /userfiles/saves/$username
proj_name = ADV_OMICS_19_5PPM
ppm = 2
packages_installed = Y
font1 = Open Sans
font2 = Open Sans
font3 = Open Sans
font4 = Open Sans
col1 = #1861ab
col2 = #DBDBDB
col3 = #FFFFFF
col4 = #FFFFFF
size1 = 30
size2 = 17
size3 = 13
size4 = 10
taskbar_image = gemmy_rainbow.png
gtheme = classic
gcols = #FF0004&#38A9FF&#FFC914&#2E282A&#8A00ED&#00E0C2&#95C200&#FF6BE4
gspec = RdBu')
    writeLines(contents, opt.loc)
  }
  
  
  # - - - 
}
