#$ -S /home/cog/jwolthuis/R-3.4.0/bin/Rscript --no-save
#$ -M J.C.Wolthuis-2@umcutrecht.nl
#$ -m beas

# --- get array num ---

args = commandArgs(trailingOnly=TRUE)
folder <- args[1]
i <- as.numeric(args[2])

# --- load libraries ---

library(xcms)

# --- set defaults ---

scale <- seq(1,64,2)

outfolder = file.path(folder, "peakcalled")

files = list.files(folder, pattern="\\.mzXML", full.names = T)

f <- files[i]

# --------------------------------------

xset_pos = NULL
attempt <- 0
startscans = c(10, 15, 20)
npeaks = 0

while( is.null(xset_pos) && attempt <= 3 ) {
    attempt <- attempt + 1
    try(
        xset_pos <- xcmsSet(method="MSW",
                        files=f,
                        scanrange = c(startscans[attempt]:150),
                        scales=scale,
                        snthresh=3,
                        nearbyPeak=T)
    )
}

if(!is.null(xset_pos)) npeaks <- nrow(xset_pos@peaks)

# --- report back ---
print(paste("Found", npeaks, "peaks in file", basename(f), "in positive mode!", sep=" "))

# --- aaaand for other mode ----

xset_neg = NULL
attempt <- 0
endscans = c(290, 285, 280)
npeaks = 0

while( is.null(xset_neg) && attempt <= 3 ) {
    attempt <- attempt + 1
    try(
    xset_neg <- xcmsSet(method="MSW",
                        files=f,
                        scanrange = c(151:endscans[attempt]),
                        scales=scale,
                        snthresh=3,
                        nearbyPeak=T)
    )
}
if(!is.null(xset_neg)) npeaks <- nrow(xset_neg@peaks)

# --- report back ---
print(paste("Found", npeaks, "peaks in file", basename(f), "in negative mode!", sep=" "))

# ---- save ----

save(xset_pos,
     xset_neg,
     file=file.path(outfolder, 
                    gsub(basename(f), 
                         pattern="\\.mzXML", 
                         replacement="_called.RData")
     ))
