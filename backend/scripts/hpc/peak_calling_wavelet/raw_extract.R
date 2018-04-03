args = commandArgs(trailingOnly=TRUE)

outdir = args[1]
array_i = args[2]

if(!dir.exists(outdir)) dir.create(outdir)

# -----------------------

sn <- data.table::fread(file.path(outdir, "../sampleNames.txt"))
sn$batch <- as.factor(gsub(sn$File_Name,
                           pattern = "_\\d\\d\\d$",
                           replacement=""))
data.table::fwrite(sn, file.path(outdir, "sampleNames.txt"))

# --- generate breaks ---

if(array_i == 1){
  generateBreaksFwhm_jw(mzmin = 70,
                        mzmax = 600,
                        cores = cores,
                        resol=resol,
                        outdir = outdir)
}

# --- detect OS ---

get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

# --- extract data ---



