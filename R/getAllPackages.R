# # get all r files
# rfiles = list.files("~/Code/MetaDBparse/", recursive=T, pattern="\\.R", full.names = T)
# reqs = pbapply::pbsapply(rfiles, function(p){
#   requirements::req_file(p)
# })
# 
# packages = unique(unlist(reqs))
