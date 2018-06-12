args <- commandArgs(TRUE)

#args = c(1, "pos", "/hpc/cog_bioinf/ridder/users/jwolthuis/Data/Metabolomics/DSM", "BR11-20_NL1-2")

require(data.table)
require(pbapply)

i = as.numeric(args[1])
scanmode = args[2]
resdir =  gsub(x = args[3], pattern = "\\\\", replacement = "")
projname = gsub(x = args[4], pattern = "\\\\", replacement = "")
              
groups_path <- file.path(resdir, paste0(projname,"_",scanmode,"_groups.RData"))    
load(groups_path)

to_do <- 1:nrow(groups$mat)
blocksize = 2000
chunks <- split(to_do, ceiling(seq_along(to_do)/blocksize))   
mychunk = chunks[i]

outlist_path <- file.path(resdir, paste0(projname,"_",scanmode,".RData"))
load(outlist_path)

groupNames <- colnames(outlist[,-1])

outrows <- pbapply::pblapply(mychunk[[1]], function(i){
  mzidx <- groups$idx[[i]]
  # FIND PEAKS IN THIS GROUP
  members <- outlist[mzidx, ]
  # --- get intensities ---
  mzmed = as.numeric(mean(members$mzmed))
  mzmin = as.numeric(min(members$mzmed))
  mzmax = as.numeric(max(members$mzmed))
  
  ints.allsamps = rep(0, length(groupNames))
  names(ints.allsamps) = groupNames 
  # # Check for each sample if multiple peaks exists, if so take the sum!
  nrsamples <- nrow(members)
  ints.allsamps <- colSums(outlist[mzidx, -1],
                           na.rm = T)
  # --- make dt ---
  outpgrlist.a = cbind('mzmed' = mzmed, 
                       nrsamples, 
                       'mzmin' = mzmin, 
                       'mzmax' = mzmax
  )
  outpgrlist.b <- data.frame(as.list(ints.allsamps))
  outpgrlist <- cbind(outpgrlist.a, outpgrlist.b)
  
  # -----------------------
  
  as.data.table(outpgrlist)
})

outpgrlist <- rbindlist(outrows)
print(dim(outpgrlist))

dir.create(file.path(resdir, projname),showWarnings = F)
fn = file.path(resdir, projname, paste0(projname,
                              "_",
                              scanmode,
                              "_",
                              i,
                              ".csv"))
print(fn)
fwrite(outpgrlist, file=fn)