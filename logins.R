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
  username = c("joby", paste0("advomics", c(1:30))),
  password = c("lol", rep("kirakira", 30)),
  role = c("admin", rep("noob", 30))
)

# insert some initial data to work with
DBI::dbWriteTable(db, "users", users, append=T)

DBI::dbDisconnect(db)
