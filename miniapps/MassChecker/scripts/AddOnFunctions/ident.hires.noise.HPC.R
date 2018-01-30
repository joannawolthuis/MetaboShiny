# modified identify function to also look for adducts and their isotopes
ident.hires.noise.HPC <- function(peaklist, allAdducts, scanmode="Negative", look4=c("Cl", "Ac"), identlist=NULL, resol=280000, slope=0, incpt=0, ppm.fixed=10, ppm.iso.fixed=1) {
  
  #     peaklist=outlist_Neg
  #     scanmode="Negative"
  #     identlist=noise.Neg.MZ
  #     look4=look4.addN
  #     resol=resol
  #     slope=0 
  #     incpt=0
  #     ppm.fixed=3
  #     ppm.iso.fixed=2  
  
  metlin <- assi <- iso <- rep("", nrow(peaklist))
  theormz <- nisos <- expint <- conf <- rep(0, nrow(peaklist))
  # nrH <- nrD <- nrC <- nr13C <- nrN <- nr15N <- nrO <- nrP  <- nrS <- nrCl <- rep(0, nrow(peaklist))
  
  # add adducts to identification list
  if (scanmode == "Positive") { add.mode <-  "+" } else { add.mode <- "-" }
  # identlist <- theorMZ
  identlist.orig <- identlist
  
  for (p in 1:length(look4)) { # loop over type of adduct
    #print(p)
    identlist.Adduct <- identlist.orig
    identlist.Adduct[ , "CompoundName"] <- as.character(identlist.orig[ , "CompoundName"])
    
    #add2label <- paste("[M+", look4[p], "]", add.mode, sep="")
    if (look4[p]=="H2O") {
      add2label <- paste("[M-", look4[p], "]", add.mode, sep="")
    } else {
      add2label <- paste("[M+", look4[p], "]", add.mode, sep="")
    }
    
    identlist.Adduct[ , "CompoundName"] <- paste(identlist.Adduct[ , "CompoundName"], add2label, sep=" ")
    adductInfo <- elementInfo(look4[p], allAdducts)
    if (scanmode == "Positive") {
      adductMass <- adductInfo$mass[1] + adductInfo$isotope$mass[1] - Hmass } else {
        adductMass <- adductInfo$mass[1] + adductInfo$isotope$mass[1] + Hmass
      }
    for (q in 1:nrow(identlist.Adduct)) { # loop over compounds in database
      # construct information for compound + adduct:
      if (scanmode == "Positive") {
        identlist.Adduct[q, "Mpos"] <-  as.numeric(identlist.Adduct[q, "Mpos"]) + adductMass # + Na - H
        identlist.Adduct[q, "MNeg"] <-  0
      } else {
        identlist.Adduct[q, "Mpos"] <-  0
        identlist.Adduct[q, "MNeg"] <-  as.numeric(identlist.Adduct[q, "MNeg"]) + adductMass # + Cl + H
      }
    }
    
#     # modify columns with info for mol. formula:
#     if (look4[p] == "2Na-H") {
#       identlist.Adduct[, "nrH"] <- as.numeric(identlist.Adduct[, "nrH"]) - 1
#     } else if (look4[p] == "NH4") {
#       identlist.Adduct[, "nrH"] <- as.numeric(identlist.Adduct[, "nrH"]) + 3
#       identlist.Adduct[, "nrN"] <- as.numeric(identlist.Adduct[, "nrN"]) + 1
#     } else if (look4[p] == "Cl") {
#       identlist.Adduct[, "nrCl"] <- as.numeric(identlist.Adduct[, "nrCl"]) + 1
#       identlist.Adduct[, "nrH"] <- as.numeric(identlist.Adduct[, "nrH"]) + 1
#     } else if (look4[p] == "Ac") {
#       identlist.Adduct[ , "nrC"] <- as.numeric(identlist.Adduct[ , "nrC"]) + 2
#       identlist.Adduct[ , "nrH"] <- as.numeric(identlist.Adduct[ , "nrH"]) + 3
#       identlist.Adduct[ , "nrO"] <- as.numeric(identlist.Adduct[ , "nrO"]) + 2
#     } else if (look4[p] == "CH3OH+H") {
#       identlist.Adduct[ , "nrC"] <- as.numeric(identlist.Adduct[ , "nrC"]) + 1
#       identlist.Adduct[ , "nrH"] <- as.numeric(identlist.Adduct[ , "nrH"]) + 4
#       identlist.Adduct[ , "nrO"] <- as.numeric(identlist.Adduct[ , "nrO"]) + 1
#     } else {
#       identlist.Adduct[, "nrH"] <- as.numeric(identlist.Adduct[, "nrH"]) - 1
#     }
    identlist <- rbind(identlist, identlist.Adduct) 
  } # end for p adducts in look4
  
  
  if (scanmode == "Positive") { theor.mcol <- as.numeric(identlist[ , "Mpos"]) } else {
    theor.mcol <- as.numeric(identlist[ , "MNeg"]) }
  # apply correction using regression line obtained with ISses
  theor.mcol <- (1+slope)*theor.mcol + incpt
  
  # get mz information from peaklist
  mcol <- peaklist[ , "mzmed.pgrp"]
  # if column with average intensities is missing, calculate it:
  if (!("avg.int" %in% colnames(peaklist))){
    mzmaxcol <- which(colnames(peaklist) == "mzmax.pgrp")
    endcol <- ncol(peaklist)
    peaklist[ , "avg.int"] <- apply(peaklist[ ,(mzmaxcol+1):(endcol)], 1, mean)
  }
  
#   # generate URL for Metlin:
#   for (p in 1:nrow(peaklist)) { # for each peak  # p <- 1
#     mzpeak <- as.numeric(mcol[p])
#     # resolution as function of mz:
#     resol.mz <- resol*(1/sqrt(2)^(log2(mzpeak/200)))
#     fwhm <- mzpeak/resol.mz
#     #    massmin <- mzpeak - (0.00003 + 0.000003*mzpeak) - (1.0078250321 - 0.00054858) - fwhm  # for Metlin db
#     #    massmax <- mzpeak - (0.00003 + 0.000003*mzpeak) - (1.0078250321 - 0.00054858) + fwhm  # for Metlin db
#     if (scanmode == "Positive") {
#       massmin <- mzpeak - (1.0078250321 - 0.00054858) - fwhm  # for Metlin db
#       massmax <- mzpeak - (1.0078250321 - 0.00054858) + fwhm  # for Metlin db
#     } else {
#       massmin <- mzpeak + (1.0078250321 - 0.00054858) - fwhm  # for Metlin db
#       massmax <- mzpeak + (1.0078250321 - 0.00054858) + fwhm  # for Metlin db
#     }
#     metlin[p] <- paste("http://metlin.scripps.edu/metabo_list.php?mass_min=", massmin, "&mass_max=", massmax, sep="")
#   }
  

  # do indentification using own database:
  for (t in 1:nrow(identlist)) { # theoretical mass  # t <- 45
    # for (t in 1:169) {
    #print(as.character(identlist[t,"CompoundName"]))
    theor.mz <- theor.mcol[t]
    
    theor.comp <- as.character(identlist[t, "Composition"])
    #theor.comp <- mol.formula(identlist[t, ])
    
#     # if there's Deuterium, Tritium, 13C or 15N in the composition:
#     mass.incr <- 0 + (as.numeric(identlist[t,"nrD"])*Dmass) +
#       (as.numeric(identlist[t,"nrT"])*Tmass) +
#       (as.numeric(identlist[t,"nrC13"])*C13mass) +
#       (as.numeric(identlist[t,"nrN15"])*N15mass)
#     theor.comp <- strsplit(theor.comp, "iso")[[1]][1]
    mass.incr <- 0
    
    # resolution as function of mz:
    resol.mz <- resol*(1/sqrt(2)^(log2(theor.mz/200)))
    # calculate fine-grained isotopic distribution using MIDAs
    fwhm <- round(theor.mz/resol.mz,6)
    
    #system(paste("C:/Users/awillem7/tools/MIDAs_New/MIDAs_Example ", theor.comp, " 2  C00 \"\"  \"\"  ", fwhm, " 1 0 1e-50 2 tmp", sep=""),ignore.stderr=TRUE)
    #system(paste("C:/Users/mraves/Metabolomics/MIDAs_New/MIDAs_Example ", theor.comp, " 2  C00 \"\"  \"\"  ", fwhm, " 1 0 1e-50 2 tmp", sep=""),ignore.stderr=TRUE)
    #system(paste(path2MIDAS, theor.comp, " 2  C00 \"\"  \"\"  ", fwhm, " 1 0 1e-50 2 tmp", sep=""),ignore.stderr=TRUE)

    #system(paste("/data/home/luyf/Metabolomics/MIDAs/MIDAs_Example ", theor.comp, " 2  C00 \"\"  \"\"  ", fwhm, " 1 0 1e-50 2 tmp", sep=""),ignore.stderr=TRUE)
    options(stringsAsFactors = FALSE)
    #fgid <- read.table(file="tmp_Fine_Grained_Isotopic_Distribution", header=FALSE)
    
#     res <- try(fgid <- read.table(file="tmp_Fine_Grained_Isotopic_Distribution", header=FALSE))
#     if(inherits(res, "try-error"))
#     {
#       #error handling code, maybe just skip this iteration using
#       message("Skipped")
#       next
#     }
#     
#     # correct mass for D, T, 13C and 15N
#     fgid[ ,1] <- as.numeric(fgid[ ,1]) + mass.incr
#     # calculate percentage intensities from relative intensities
#     firstone <- as.numeric(fgid[1,2])
#     fgid[ ,3] <- as.numeric(fgid[ , 2]) / firstone
#     #fgid <- as.matrix(fgid, ncol=3)
#     # the mz in the MIDAs file are of the neutral molecule
#     if (scanmode == "Positive") { mz.iso <- as.numeric(fgid[ , 1]) + Hmass - electron }
#     if (scanmode == "Negative") { mz.iso <- as.numeric(fgid[ , 1]) - Hmass + electron }
#     fgid <- cbind(fgid, mz.iso)
#     colnames(fgid) <- c("mz","rel.int", "perc.int", "mz.iso")
    
    # compensate mz for presence of adduct
#     Adduct.mass <- theor.mz - mz.iso[1]
#     fgid[ , "mz.iso"] <- as.numeric(fgid[ , "mz.iso"]) + Adduct.mass
    
    # set tolerance for mz accuracy of main peak
    mtol <- theor.mz*ppm.fixed/1000000
    # find main peak
    selp <- which(mcol > (theor.mz - mtol) & mcol < (theor.mz + mtol))
    # selp <- which(mcol > (theor.mz - 0.01) & mcol < (theor.mz + 0.01))
    # peaklist[selp, c(1:4,(endcol+1))]
    
    # set tolerance for mz accuracy of isotope peaks
    itol <- theor.mz*ppm.iso.fixed/1000000
    
#     if (length(selp) > 1) { # more than one candidate peak for main; select best one based on isotope pattern
#       #cat(as.character(identlist[t, "CompoundName"])); print(" has >1 candidate peaks")
#       conf.local <- rep(0, length(selp))
#       for (p in 1:length(selp)) { # p <- 2
#         # determine isotope pattern for each candidate peak
#         # obs.mz <- peaklist[selp[p],"mzmed.pgrp"]
#         conf.local[p] <- match.isotope.pattern(peaklist, scanmode, selp[p], fgid, ppm.iso.fixed)
#       }
#       # selp <- selp[abs(mcol[selp] - theor.mz) == min(abs(mcol[selp] - theor.mz))] }
#       selp <- na.exclude(selp[conf.local == max(conf.local)])
#     }
    if (length(selp) > 1) { # more than one candidate peak for main; select best one based on mz.diff
      selp <- selp[abs(mcol[selp] - theor.mz) == min(abs(mcol[selp] - theor.mz))] 
    }
    if (length(selp) == 1) { # match for main
      assi[selp] <- paste(assi[selp], as.character(identlist[t,"CompoundName"]), sep=";")
      theormz[selp] <- theor.mz
      # conf[selp] <- match.isotope.pattern(peaklist, scanmode, selp, fgid, ppm.iso.fixed)
      
      #       nrH[selp] <- identlist[t,"nrH"]
      #       nrD[selp] <- identlist[t,"nrD"]
      #       nrC[selp] <- identlist[t,"nrC"]
      #       nr13C[selp] <- identlist[t,"nrC13"]
      #       nrN[selp] <- identlist[t,"nrN"]
      #       nr15N[selp] <- identlist[t,"nrN15"]
      #       nrO[selp] <- identlist[t,"nrO"]
      #       nrP[selp] <- identlist[t,"nrP"]
      #       nrS[selp] <- identlist[t,"nrS"]
      #       nrCl[selp] <- identlist[t,"nrCl"]   
      
      # assign isotope peaks
#       mz.main <- peaklist[selp, "mzmed.pgrp"]  # mz of main peak
#       int.main <- peaklist[selp, "avg.int"]  # intensity of main peak (= 100%)
#       # deviation from theoretical mass:
#       diff <- theor.mz - mz.main
#       # calculate expected intensities and select isotopes with exp.int > threshold
#       fgid[ , "exp.int"] <- fgid[ , "perc.int"] * int.main      
#       fgid.subset <- fgid[(fgid[ , "exp.int"] > thresh), ]
#       nisos[selp] <- nrow(fgid.subset) - 1
#       if (nrow(fgid.subset) > 1) { # avoid error message if fgid.subset has only 1 line
#         for (f in 2:nrow(fgid.subset)) { # f <- 2
#           mz.target <- fgid.subset[f, "mz.iso"] - diff
#           int.target <- fgid.subset[f, "exp.int"]
#           # print(itarget)
#           sel.iso <- peaklist[ , "mzmed.pgrp"] > (mz.target - itol) & peaklist[ , "mzmed.pgrp"] < (mz.target + itol)
#           # sum(sel.iso)
#           # sel.iso <- peaklist[ , "mzmed.pgrp"] > (mz.target - 0.01) & peaklist[ , "mzmed.pgrp"] < (mz.target + 0.01)
#           # peaklist[sel.iso, c(1:4,23)]
#           if (sum(sel.iso) == 1) { # 2 separate if-statements because of error if sum(sel.iso) = 0
#             if (peaklist[sel.iso, "avg.int"] >  (int.target/2)) { # match
#               iso[sel.iso] <- paste(paste(iso[sel.iso], as.character(identlist[t,"CompoundName"]), "iso", f, sep=" "),";", sep="")
#               # peaklist[sel.iso, ]
#               expint[sel.iso] <- fgid.subset[f, "exp.int"]
#             }
#           } else if (sum(sel.iso) > 1) {
#             nrs.iso <- which(sel.iso)
#             nr.iso <- nrs.iso[which(abs(peaklist[sel.iso, "avg.int"] - int.target) == min(abs(peaklist[sel.iso, "avg.int"] - int.target)))]
#             if (peaklist[nr.iso, "avg.int"] >  (int.target/2)) {  # match
#               # print(peaklist[nr.iso, "avg.int"])
#               iso[nr.iso] <- paste(paste(iso[nr.iso], as.character(identlist[t,"CompoundName"]), "iso", f, sep=" "),";", sep="")
#               expint[nr.iso] <- fgid.subset[f, "exp.int"]
#             } # end if
#           } # end else if
#         } # end for f 
#       } # end if
    } # end if 
  } # end for t
  #   cbind(peaklist, nrH, nrD, nrC, nr13C, nrN, nr15N, nrO, nrP, nrS, nrCl, assi, theormz, conf, nisos, iso, expint, metlin)
  cbind(peaklist, assi, theormz, conf, nisos, iso, expint, metlin)
}
