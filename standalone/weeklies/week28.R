# PMIDS:
# 28130542, 27775610, 29086042, 19015140, PMC3608065, PMC5192443, 26371303, 23248816, 22978259, 19518135
# 20391594, 29225718, 28520864, 27571918, 26819880, 24947013, 28520864
# GOOD WEBSITE:
# http://casmi-contest.org/2017/index.shtml


#- - - checking the isotope scoring - - -

patdb = "/Users/jwolthuis/Analysis/SP/BrazilAndSpain_W.db"
chosen.db <- "/Users/jwolthuis/Google Drive/MetaboShiny/backend/db/internal.full.db"

conn <- RSQLite::dbConnect(RSQLite::SQLite(), patdb)
query.zero <- gsubfn::fn$paste("ATTACH '$chosen.db' AS db")

#ATTACH "/Users/jwolthuis/Google Drive/MetaboShiny/backend/db/hmdb.full.db" AS db
#cpd = 178.01787
#cpd = curr_cpd

RSQLite::dbExecute(conn, query.zero)

query.one <- gsubfn::fn$paste(strwrap(
  "SELECT base.compoundname, cpd.baseformula, cpd.fullmz, cpd.fullformula, cpd.adduct, cpd.isoprevalence, cpd.basecharge, int.* 
  FROM db.extended cpd indexed by e_idx1
  JOIN db.base base
  ON base.baseformula == cpd.baseformula AND base.charge == cpd.basecharge
  JOIN mzranges rng
  ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
  JOIN mzintensities int
  ON int.mzmed BETWEEN rng.mzmin AND rng.mzmax
  WHERE base.compoundname LIKE '%(IS)%'"
  , width=10000, simplify=TRUE))

table <- RSQLite::dbGetQuery(conn,query.one)

dim(table)

head(table)

table <- data.table::as.data.table(table[complete.cases(table),])

dim(table)

p.cpd <- split(x = table, 
              f = list(table$compoundname, table$adduct))

cpd_tab <- p.cpd$`13C6-Phenylalanine (IS).M-2H`[isoprevalence > 99.99999999]
raster::cv(cpd_tab$intensity)

cl <- parallel::makeCluster(3, "FORK")

require(data.table)

which_score = "chisq"
which_combine = "mean"
normalize = TRUE


if(nrow(table) > 0){
  res_rows <- pbapply::pblapply(p.cpd, cl=NULL, function(cpd_tab){
    
    # aggregate multimz matches
    
    cpd_tab <- setDT(cpd_tab)[, .(mzmed = mean(mzmed), intensity = sum(intensity)),
                             by=.(compoundname, baseformula, fullmz, fullformula, adduct, isoprevalence, filename)]
    
    formula = unique(cpd_tab$baseformula)
    adduct = unique(cpd_tab$adduct)
    
    if(any(cpd_tab$isoprevalence > 99.999999)){

      sorted <- as.data.table(unique(cpd_tab[order(cpd_tab$isoprevalence, 
                                                   decreasing = TRUE),]))
      
      split.by.samp <- split(sorted, 
                             sorted[,"filename"])
      # - - - - - - - - - 
      
      score <- sapply(split.by.samp, function(samp_tab){
        
        if(nrow(samp_tab) == 1){
          res1 = NA
        }else{
          
          res1 = NA   
          
          theor_mat <- samp_tab[,c("fullmz", "isoprevalence")]
          theor <- matrix(ncol = nrow(theor_mat), nrow = 2, data = c(theor_mat$fullmz, theor_mat$isoprevalence),byrow = T)
          
          obs_mat <- samp_tab[,c("mzmed", "intensity")]
          obs <- matrix(ncol = nrow(obs_mat), nrow = 2, data = c(obs_mat$mzmed, obs_mat$intensity),byrow = T)

          if(normalize){
            theor[2,] <- theor[2,]/sum(theor[2,])
            obs[2,] <- obs[2,]/sum(obs[2,])
          } 
          
          ppm = 1
          
          try({
            res1 <- switch(which_score, 
                           chisq = {
                             test <- chisq.test( obs[2,], p = theor[2,], rescale.p = T)
                             # - - -
                             as.numeric(test$p.value)
                           },
                           mape = {
                             actual = obs[2,]
                             theor = theor[2,]
                             deltaSignal = abs(theor - actual)
                             percentageDifference = deltaSignal / actual * 100# Percent by element.
                             # - - -
                             mean(percentageDifference) #Average percentage over all elements.
                           },
                           cocharm = {
                             isos <- matrix(c(samp_tab$intensity, samp_tab$isoprevalence), byrow=TRUE, nrow=2, dimnames=list(istrue=0:1, iso=0:(-1+nrow(samp_tab))))
                             # - - -
                             as.numeric(DescTools::CochranArmitageTest(x = isos,
                                                                       alternative = "two.sided")$p.value)
                           },
                           mscore = {
                             InterpretMSSpectrum::mScore(obs=obs, the=theor, dppm = ppm)
                           },
                           sirius = {
                             ## if you want the so-called 'error function'
                             erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
                             ## (see Abramowitz and Stegun 29.2.29)
                             ## and the so-called 'complementary error function'
                             erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
                             
                             # f is measured peak intensity
                             smass <- function(fs, ft, ppm){
                               ppm * 10e-6 * (fs - ft)
                             }
                             prob_dev <- erfc((abs(M - m))/(sqrt(2) * smass(f)))
                           })
          })
        }
        res1
      })
      
      #pvals <- as.numeric(unlist(score))
      scores <- as.numeric(unlist(score))
      
      #df = 2*length(pvals)
      #combined_pval <- pchisq( -2*sum(log(pvals)), df, lower.tail=FALSE)
      #combined_pval <- pgamma(-sum(log(pvals)), length(pvals), lower.tail=FALSE)

      res2 <- switch(which_combine,
                     fischer = {
                       df = 2*length(pvals)
                       combined_pval <- pchisq( -2*sum(log(scores)), df, lower.tail=FALSE)
                       #combined_pval <- pgamma(-sum(log(scores)), length(scores), lower.tail=FALSE)
                       data.table::data.table(baseformula = formula,
                                              adduct = adduct,
                                              score = combined_pval)
                       },
                     mean = {
                       combined_score <- if(all(is.na(scores))) NA else mean(scores, na.rm=T)
                       
                       data.table::data.table(baseformula = formula,
                                                      adduct = adduct,
                                                      score = combined_score)
                     },
                     all = {
                      data.table::data.table(baseformula = c(formula),
                                                      adduct = c(adduct),
                                                      score = scores)
                     })
      # - if using average .. -
     
      print(res2)
      
      res2
      
    }else{
      data.table::data.table()
    }
  })

  results <- rbindlist(res_rows)
  results <- results[complete.cases(results),]
  
  #results_cocharm_norm <- results
  #results_mape_norm <- results
  #results_mape_no_norm <- results
  #results_chisq_norm <- results
  #results_chisq_no_norm <- results
  #results_mscore_norm <- results
  #results_mscore_no_norm <- results
  
  require(ggplot2)
  
  # - for score from package - 
  ggplot(results, aes(score,
                      group = adduct,
                      fill = adduct)) +
    ggtitle("IS isotopic distribution") +
    theme(plot.title=element_text(size=30, face="bold",hjust = 0.5)) + 
    geom_density(aes(y=..scaled..), alpha=0.6) +
    #geom_rect(aes(xmin=0, xmax=50, ymin=0, ymax=1.00), color="red", alpha=0.5, fill="red") +
    #geom_vline(xintercept = 0.05, color = "black") +
    #geom_histogram(aes(y=..ncount..), alpha=0.6, bins=2) +
    scale_x_log10() +
    facet_wrap(~adduct)
  
  # - for pvals from chisq -
  ggplot(results, aes(score,
                      group = adduct,
                      fill = adduct)) +
    ggtitle("IS isotopic distribution") +
    theme(plot.title=element_text(size=30, face="bold")) +
    geom_density(aes(y=..scaled..), alpha=0.6) +
    geom_vline(xintercept = 0.001, color = "red") +
    geom_vline(xintercept = 0.05, color = "green") +
    #geom_histogram(aes(y=..ncount..), alpha=0.6, bins=2) +
    #scale_x_log10() +
    facet_wrap(~adduct)
  
}else{
  data.table::data.table() 
}


# - - - FIND INTERNAL STANDARDS THAT HAVE OTHER MATCHES - - -

int.db <- "/Users/jwolthuis/Google Drive/MetaboShiny/backend/db/internal.full.db"
other.dbs <- list.files("/Users/jwolthuis/Google Drive/MetaboShiny/backend/db/",pattern=".full.db$",full.names = T)
other.dbs <- other.dbs[-which(grepl(other.dbs,pattern="internal"))]
conn.int <- RSQLite::dbConnect(RSQLite::SQLite(), int.db)

# (attempt 2) - - search one by one in each database.

standards <- RSQLite::dbGetQuery(conn.int, "SELECT * FROM extended e
                                 JOIN base b
                                 ON b.baseformula = e.baseformula
                                 AND b.charge = e.basecharge
                                 WHERE b.compoundname LIKE '%(IS)%'")

head(standards)

ppm = 1
searcher <- data.table::data.table(mzmed = standards$fullmz,
                                   mzmin = sapply(standards$fullmz,
                                                  FUN=function(mz, ppm){
                                                    mz - mz * (ppm / 1E6)}, ppm=ppm),
                                   mzmax = sapply(standards$fullmz,
                                                  FUN=function(mz, ppm){
                                                    mz + mz * (ppm / 1E6)}, ppm=ppm))
#use hmdb as example?

conn.pat <- RSQLite::dbConnect(RSQLite::SQLite(), patdb)

search.db <- "/Users/jwolthuis/Google Drive/MetaboShiny/backend/db/hmdb.full.db"

RSQLite::dbExecute(conn.pat, gsubfn::fn$paste("ATTACH '$search.db' AS db"))

results <- pbapply::pblapply(1:nrow(searcher), function(i){
  
  row <- searcher[i,]
  mzmed = row$mzmed
  mzmin = row$mzmin
  mzmax = row$mzmax
  
  query.one <- gsubfn::fn$paste(strwrap(
    "CREATE TEMP TABLE unfiltered AS
    SELECT cpd.baseformula, cpd.adduct
    FROM mzvals mz
    JOIN mzranges rng ON rng.ID = mz.ID
    JOIN db.extended cpd indexed by e_idx2
    ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
    AND mz.foundinmode = cpd.foundinmode
    WHERE cpd.fullmz BETWEEN $mzmin AND $mzmax",width=10000, simplify=TRUE))
  # 1. Find matches in range (reasonably fast <3)
  RSQLite::dbExecute(conn.pat, query.one)
  #  2. get isotopes for these matchies (reverse search)
  
  #RSQLite::dbExecute(conn,"drop table isotopes")
  
  query.two <- gsubfn::fn$paste(strwrap(
    "SELECT cpd.baseformula,cpd.fullformula, cpd.adduct, cpd.isoprevalence, cpd.basecharge, int.* 
    FROM db.extended cpd indexed by e_idx1
    JOIN unfiltered u
    ON u.baseformula = cpd.baseformula
    AND u.adduct = cpd.adduct
    JOIN mzranges rng
    ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
    JOIN mzintensities int
    ON int.mzmed BETWEEN rng.mzmin AND rng.mzmax
    JOIN db.base b ON
    b.baseformula = cpd.baseformula
    AND b.charge = cpd.basecharge"
    , width=10000, simplify=TRUE))
  
  table <- RSQLite::dbGetQuery(conn.pat,query.two)
  
  RSQLite::dbRemoveTable(conn.pat, "unfiltered")
  
  #print(table)
  
  table
})


# (attempt 1)- - CREATE RANGE TABLE FOR INTERNAL STANDARDS - -
# 
# standards <- RSQLite::dbGetQuery(conn, "SELECT * FROM extended e
#                                  JOIN base b
#                                  ON b.baseformula = e.baseformula
#                                  AND b.charge = e.basecharge
#                                  WHERE b.compoundname LIKE '%(IS)%'")
# mzvals <- data.table::data.table(mzmed = standards$fullmz,
#                                  foundinmode = standards$foundinmode)
# 
# mzranges <- data.table::data.table(mzmin = sapply(standards$fullmz, 
#                                                   FUN=function(mz, ppm){
#                                                     mz - mz * (ppm / 1E6)}, ppm=ppm),
#                                    mzmax = sapply(standards$fullmz, 
#                                                   FUN=function(mz, ppm){
#                                                     mz + mz * (ppm / 1E6)}, ppm=ppm))
# rtree = TRUE
# 
# sql.make.meta <- strwrap("CREATE TABLE mzvals_is(
#                          ID INTEGER PRIMARY KEY AUTOINCREMENT,
#                          mzmed decimal(30,13),
#                          foundinmode text)", width=10000, simplify=TRUE)
# RSQLite::dbExecute(conn, sql.make.meta)
# RSQLite::dbExecute(conn, "create index mzfind_is on mzvals_is(mzmed, foundinmode);")
# 
# # --- write vals to table ---
# RSQLite::dbWriteTable(conn, "mzvals_is", mzvals, append=TRUE) # insert into
# # --- make range table (choose if R*tree or not) ---
# sql.make.rtree <- strwrap("CREATE VIRTUAL TABLE mzranges_is USING rtree(
#                           ID INTEGER PRIMARY KEY AUTOINCREMENT,
#                           mzmin decimal(30,13),
#                           mzmax decimal(30,13));"
#                           , width=10000, simplify=TRUE)
# RSQLite::dbExecute(conn, if(rtree) sql.make.rtree else sql.make.normal)
# 
# # --- write ranges to table ---
# RSQLite::dbWriteTable(conn, "mzranges_is", mzranges, append=TRUE) # insert into
# 
# RSQLite::dbExecute(conn, gsubfn::fn$paste("ATTACH '$patdb' AS pat"))
# 
# 
# results_in_range <- pbapply::pblapply(other.dbs[-1], function(db){
#   
#   print(db)
#   
#   RSQLite::dbExecute(conn, gsubfn::fn$paste("ATTACH '$db' AS mydb"))
#   
#   # - - join - -
#   
#   # "SELECT mz.mzmed as orig_mz, b1.compoundname as orig_name, b1.baseformula as orig_formula, e.adduct as orig_adduct, e.isoprevalence as orig_iso,
#   #           cpd.mzmed as new_mz, b2.compoundname as new_name, b2.baseformula as new_formula, cpd.adduct as new_adduct, cpd.isoprevalence as new_iso
#   # 
#   query.one <- gsubfn::fn$paste(strwrap(
#     "SELECT mz.mzmed as orig_mz, e2.fullmz as new_mz, 
#             b1.compoundname as orig_name, b2.compoundname as new_name,
#             e1.fullformula as orig_formula, e2.fullformula as new_formula, 
#             e1.adduct as orig_adduct,  e2.adduct as new_adduct,
#             e1.isoprevalence as orig_iso, e2.isoprevalence as new_iso,
#             int.*
#     FROM mzvals_is mz
#     JOIN mzranges_is rng ON rng.ID = mz.ID
#     JOIN mydb.extended e2 indexed by e_idx2
#     ON e2.fullmz BETWEEN rng.mzmin AND rng.mzmax
#     AND mz.foundinmode = e2.foundinmode
#     JOIN extended e1
#     ON e1.fullmz = mz.mzmed
#     JOIN base b1
#     ON e1.baseformula = b1.baseformula
#     AND e1.basecharge = b1.charge
#     LEFT JOIN mydb.base b2
#     ON e2.baseformula = b2.baseformula
#     AND e2.basecharge = b2.charge
#     WHERE e1.fullformula != e2.fullformula",width=10000, simplify=TRUE))
#   # 1. Find matches in range (reasonably fast <3)
#   res <- RSQLite::dbGetQuery(conn, query.one)
#   
#   res_filt <- res[complete.cases(res),]
#   # - - - - - - -
#   
#   RSQLite::dbExecute(conn, gsubfn::fn$paste("DETACH mydb"))
#   
#   res_filt
# })

RSQLite::dbListTables(conn)

head(results_in_range[[3]])

comparable.pairs <- pbapply::pblapply(results_in_range, function(dt){
  unique(dt[,c("orig_name", "orig_adduct", "new_name", "new_adduct")])
})

# - - opt 3 - - - known package... tell jeroensaskia

library(InterpretMSSpectrum)

test <- standards$fullformula[[1]]
Rdisop::getFormula(test)

# get theoretical isotopic pattern of Glucose
glc <- Rdisop::getMolecule("C6H12O6")$isotopes[[1]][,1:3]
mScore(obs=glc, the=glc)


# modify pattern by maximum allowable error (2ppm mass error, 2% int error)
glc_theoretic <- glc
glc[1,] <- glc[1,]+2*glc[1,]/10^6
glc[2,1:2] <- c(-0.02,0.02)+glc[2,1:2]
mScore(obs=glc, the=glc_theoretic)
# simulate mass and int defects
ef <- function(x, e) {runif(1,x-x*e,x+x*e)}
glc_obs <- glc
glc_obs[1,] <- sapply(glc[1,], ef, e=2*10^-6)
glc_obs[2,] <- sapply(glc[2,], ef, e=0.02)
mScore(obs=glc_obs, the=glc)
# simulate mass and int defects systematically
ef <- function(x, e) {runif(1,x-x*e,x+x*e)}
n <- 11
mz_err <- round(seq(0,5,length.out=n),3)
int_err <- round(seq(0,0.1,length.out=n),3)
mat <- matrix(NA, ncol=n, nrow=n, dimnames=list(mz_err, 100*int_err))
glc_obs <- glc
for (i in 1:n) {
  glc_obs[1,] <- sapply(glc[1,], ef, e=mz_err[i]*10^-6)
  for (j in 1:n) {
    glc_obs[2,] <- sapply(glc[2,], ef, e=int_err[j])
    mat[i,j] <- mScore(obs=glc_obs, the=glc)
  }
}
plot(x=1:n, y=1:n, type="n",axes=FALSE, xlab="mass error [ppm]", ylab="isoratio error [%]")
axis(3,at=1:n,rownames(mat),las=2); axis(4,at=1:n,colnames(mat),las=2); box()
cols <- grDevices::colorRampPalette(colors=c(2,6,3))(diff(range(mat))+1)
cols <- cols[mat-min(mat)+1]
text(x=rep(1:n,each=n), y=rep(1:n,times=n), labels=as.vector(mat), col=cols)

# - - - - CAMERA - - - - 

BiocInstaller::biocLite("CAMERA")

library(CAMERA)
file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
an   <- xsAnnotate(xs)
an   <- groupFWHM(an)
an   <- findIsotopesWithValidation(an)


# ================== SOME HPC PIPELINE EXPERIMENTS ===============

# - - - - XCMS - - - -

library(MSnbase)
library(xcms)

set.seed(123)
mzs <- seq(1, 20, 0.1)
ints1 <- abs(rnorm(length(mzs), 10))
ints1[11:20] <- c(15, 30, 90, 200, 500, 300, 100, 70, 40, 20) # add peak
ints2 <- abs(rnorm(length(mzs), 10))
ints2[11:20] <- c(15, 30, 60, 120, 300, 200, 90, 60, 30, 23)
ints3 <- abs(rnorm(length(mzs), 10))
ints3[11:20] <- c(13, 20, 50, 100, 200, 100, 80, 40, 30, 20)

## Create the spectra.
sp1 <- new("Spectrum1", mz = mzs + rnorm(length(mzs), sd = 0.01),
           intensity = ints1)
sp2 <- new("Spectrum1", mz = mzs + rnorm(length(mzs), sd = 0.01),
           intensity = ints2)
sp3 <- new("Spectrum1", mz = mzs + rnorm(length(mzs), sd = 0.009),
           intensity = ints3)

## Combine the spectra
sp_agg <- combineSpectra(list(sp1, sp2, sp3))

## Plot the spectra before and after combining
par(mfrow = c(2, 1), mar = c(4.3, 4, 1, 1))
plot(mz(sp1), intensity(sp1), xlim = range(mzs[5:25]), type = "h", col = "red")
points(mz(sp2), intensity(sp2), type = "h", col = "green")
points(mz(sp3), intensity(sp3), type = "h", col = "blue")
plot(mz(sp_agg), intensity(sp_agg), xlim = range(mzs[5:25]), type = "h",
     col = "black")


## Create a MSWParam object
mp <- MSWParam()
## Change snthresh parameter
snthresh(mp) <- 15
mp

## Loading a small subset of direct injection, single spectrum files
library(msdata)

my_spec <- readMSData("/Users/jwolthuis/Documents/umc/data/Data/Project 2018-009 DSM feed-11 (Veronesi) - Saskia v Mil_part1/RES-2018-02-07_DSM DBS Veronesi/RES_20180207_001.mzXML", msLevel. = 1, mode = "inMemory")

polarity <- sapply(as.list(my_spec@assayData), function(x) x@polarity)
times <- sapply(as.list(my_spec@assayData), function(x) x@rt)

pos_idx <- sapply(as.list(my_spec@assayData), function(x) x@polarity == 1)
pos_times <- times[pos_idx]
neg_times <- times[!pos_idx]

pos <- filterRt(my_spec, c(min(pos_times), max(pos_times)))
neg <- filterRt(my_spec, c(min(neg_times), max(neg_times)))

combined_pos <- combineSpectra(as.list(pos@assayData), mzd = 0, intensityFun = "sum")#, main = length(pos))
combined_neg <- combineSpectra(as.list(neg@assayData), mzd = 0, intensityFun = "sum")#, main = length(neg))

pos_spec <- my_spec
# change featuredata, nspectra
pos_spec@featureData <- pos_spec@featureData[1,]
name <- rownames(pos_spec@featureData[1,])

envy <- new.env()
envy[[name]] <- combined_pos
pos_spec@assayData <- envy

MSnbase::writeMSData(pos_spec, file = "/Users/jwolthuis/Analysis/test.mzXML")

fticr <- readMSData( "/Users/jwolthuis/Analysis/test.mzXML", msLevel. = 1, mode = "onDisk")

## Perform the MSW peak detection on these:
p <- MSWParam(scales = seq(1,7,0.5), 
              peakThr = NULL, 
              ampTh = 0.0000001,
              snthresh = 0.001,
              peakScaleRange = 1,
              tuneIn = TRUE,
              minNoiseLevel = 0,
              SNR.method = "data.mean", 
              winSize.noise = 5, 
              nearbyPeak = TRUE
              )

p <- MSWParam(scales = seq(1,7,0.5), 
              peakThr = NULL, 
              ampTh = 0.0000001,
              snthresh = 0.001,
              peakScaleRange = 1,
              tuneIn = TRUE,
              minNoiseLevel = 0,
              SNR.method = "data.mean", 
              winSize.noise = 5, 
              nearbyPeak = TRUE
)
fticr_id <- findChromPeaks(fticr, param = p)

rng <- c(1,100)
head(chromPeaks(fticr_id))
plot(x = mz(fticr)[[1]][rng[1]:rng[2]], y = intensity(fticr)[[1]][rng[1]:rng[2]])

?chromPeaks

# - - - -orig peak calling  - - - -

spec <- fticr_id

spec <- MSnbase::smooth(pos_spec,method = "SavitzkyGolay", halfWindowSize = 2)

df <- matrix(data = intensity(spec)[[1]],
             ncol=1)
rownames(df) <- mz(spec)[[1]]

require(MassSpecWavelet)

scales <- seq(1,7,0.5)
wCoefs <- cwt(df, 
              scales=scales, 
              wavelet='mexh')
localMax <- getLocalMaximumCWT(wCoefs,
                               amp.Th = 100,
                               minWinSize = 5)

maxInd <- which(localMax > 0)
plot(df[1:100,], type='l')
points(maxInd, df[maxInd], col='red')

ridgeList <- getRidge(localMax,gapTh = 1)
plotRidgeList(ridgeList,range = c(1,100))

majorPeakInfo <- identifyMajorPeaks(ms = df, 
                                    ridgeList = ridgeList,
                                    wCoefs = wCoefs,
                                    SNR.Th = 0)

peakIndex <- majorPeakInfo$allPeakIndex
plotPeak(df[1:100], peakIndex, main=paste('Identified peaks with SNR >',3)) 

# - - - - - - - -

mz.idx <- majorPeakInfo$peakCenterIndex

# - - - - - - - -

peaks = MALDIquant::createMassPeaks(mass = as.numeric(rownames(df))[mz.idx],
                                    intensity = majorPeakInfo$peakValue,
                                    snr = majorPeakInfo$peakSNR,
                                    metaData = list(sample = "Test"))
