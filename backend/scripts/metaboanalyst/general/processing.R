##################################################
## R script for MetaboAnalyst
## Description: processing raw data types
##
## Author: Jeff Xia, jeff.xia@mcgill.ca
## McGill University, Canada
##
## License: GNU GPL (>= 2)
###################################################

# basic sanity check for the content
# return 1 or 0 based on the result
SanityCheckData<-function(){

    msg = NULL;
    cls <- dataSet$orig.cls;
    dataSet$small.smpl.size <- 0;
    # check class info
    if(dataSet$cls.type == "disc"){
        if(substring(dataSet$format,4,5)=="ts"){
            if(dataSet$design.type =="time"){
                msg<-c(msg, "The data is time-series data.");
            }else if(dataSet$design.type =="time0"){
                msg<-c(msg, "The data is time-series only data.");
            }else{
                msg<-c(msg, "The data is not time-series data.");
            }
            clsA.num <- length(levels(dataSet$facA));
            clsB.num <- length(levels(dataSet$facB));
            msg<-c(msg, paste(clsA.num, "groups were detected in samples for factor", dataSet$facA.lbl));
            msg<-c(msg, paste(clsB.num, "groups were detected in samples for factor", dataSet$facB.lbl));
        }else{
            # checking if too many groups but a few samples in each group
            cls.lbl <- dataSet$orig.cls;
            min.grp.size <- min(table(cls.lbl));

            cls.num <- length(levels(cls.lbl));
            if(cls.num/min.grp.size > 3){
                dataSet$small.smpl.size <- 1;
                msg <- c(msg, "<font color='red'>Too many groups with very small number of replicates!</font>");
                msg <- c(msg, "<font color='red'>Only a subset of methods will be available for analysis!</font>");
            }

            msg<-c(msg, paste(cls.num, "groups were detected in samples."));
            dataSet$cls.num <- cls.num;
            dataSet$min.grp.size <- min.grp.size;

            if(dataSet$paired){
                msg<-c(msg,"Samples are paired.");
                # need to first set up pair information if not csv file
                if(!(dataSet$type=="conc" | dataSet$type=="specbin" | dataSet$type=="pktable" )){
                    pairs<-ReadPairFile();
                    # check if they are of the right length
                    if(length(pairs)!=nrow(dataSet$orig)){
                        AddErrMsg("Error: the total paired names are not equal to sample names.");
                        return(0);
                    }else{
                        # matching the names of the files
                        inx<-match(rownames(dataSet$orig),names(pairs));
                        #check if all matched exactly
                        if(sum(is.na(inx))>0){
                            AddErrMsg("Error: some paired names not match the sample names.");
                            return(0);
                        }else{
                            dataSet$pairs<-pairs[inx];
                        }
                    }
                }

                pairs<-dataSet$pairs;

                # check if QC samples are present
                qc.hits <- tolower(as.character(cls)) %in% "qc";
                dataSet$orig <- dataSet$orig;
                if(sum(qc.hits) > 0){
                    AddErrMsg("<font color='red'>Error: QC samples not supported in paired analysis mode.</font>");
                    AddErrMsg("You can perform QC filtering using regular two-group labels.");
                    AddErrMsg("Then re-upload your data (without QC samples) for paired analysis.");
                    return(0);
                }else{
                    pairs<-as.numeric(pairs);
                }

                label<-as.numeric(pairs);
                cls <- as.factor(ifelse(label>0,1,0));
                dataSet$pairs<-label;

                lev<-unique(pairs);
                uni.cl<-length(lev);
                uni.cl.abs<-uni.cl/2;             
                sorted.pairs<-sort(pairs,index=TRUE);

                if(!all(sorted.pairs$x==c(-uni.cl.abs:-1,1:uni.cl.abs))){
                    AddErrMsg("There are some problems in paired sample labels! ");
                    if(uni.cl.abs != round(uni.cl.abs)){
                        AddErrMsg("The total samples must be of even number!");
                    }else{
                        AddErrMsg(paste("And class labels between ",-uni.cl.abs,
                            " and 1, and between 1 and ",uni.cl.abs,".",sep=""));
                    }
                    return(0);
                }else{  
                    msg<-c(msg,"The labels of paired samples passed sanity check.");
                    msg<-c(msg, paste("A total of", uni.cl.abs, "pairs were detected."));
                    # make sure paired samples are sorted 1:n/2 and -1:-n/2

                    x<-sorted.pairs$ix[(uni.cl.abs+1):uni.cl]
                    y<-sorted.pairs$ix[uni.cl.abs:1]
                    index<-as.vector(cbind(x,y));
                    cls<-cls[index];
                    pairs <- pairs[index];
                    dataSet$pairs<-pairs;
                    dataSet$orig.cls<-cls;
                    dataSet$orig<-dataSet$orig[index,];
                }
            }else{
                msg<-c(msg,"Samples are not paired.");
            }
        }

        #samples may not be sorted properly, need to do some sorting at the beginning 
        if(substring(dataSet$format,4,5)=="ts"){
            nfacA <- dataSet$facA;
            nfacB <- dataSet$facB;
            if(dataSet$design.type =="time" | dataSet$design.type =="time0"){
                # determine time factor and should order first by subject then by each time points
                if(tolower(dataSet$facA.lbl) == "time"){ 
                    time.fac <- nfacA;
                    exp.fac <- nfacB;
                }else{
                    time.fac <- nfacB;
                    exp.fac <- nfacA;
                }
                # update with new index
                ord.inx <- order(exp.fac);
            }else{
                ord.inx <- order(nfacA);
            }
            dataSet$orig.cls <- dataSet$orig.cls[ord.inx];
            dataSet$orig <- dataSet$orig[ord.inx, ];
            dataSet$facA <- dataSet$orig.facA <- dataSet$facA[ord.inx];
            dataSet$facB <- dataSet$orig.facB <- dataSet$facB[ord.inx];
        }else{
            ord.inx <- order(dataSet$orig.cls);
            dataSet$orig.cls <- cls[ord.inx];
            dataSet$orig<-dataSet$orig[ord.inx, ];
            if(dataSet$paired){
                dataSet$pairs<-dataSet$pairs[ord.inx];
            }
        }
     }
     msg<-c(msg,"Only English letters, numbers, underscore, hyphen and forward slash (/) are allowed.");
     msg<-c(msg,"<font color=\"orange\">Other special characters or punctuations (if any) will be stripped off.</font>");

     int.mat<-dataSet$orig;

     # check numerical matrix
     rowNms <- rownames(int.mat);
     colNms <- colnames(int.mat);
     naNms <- sum(is.na(int.mat));

     num.mat<-apply(int.mat, 2, as.numeric)

     if(sum(is.na(num.mat)) > naNms){
            # try to remove "," in thousand seperator if it is the cause
            num.mat <- apply(int.mat,2,function(x) as.numeric(gsub(",", "", x)));
            if(sum(is.na(num.mat)) > naNms){
                msg<-c(msg,"<font color=\"red\">Non-numeric values were found and replaced by NA.</font>");
            }else{
                msg<-c(msg,"All data values are numeric.");
            }
      }else{
            msg<-c(msg,"All data values are numeric.");
      }

      int.mat <- num.mat;
      rownames(int.mat)<-rowNms;
      colnames(int.mat)<-colNms;

      # check for columns with all constant (var =0)
      varCol <- apply(int.mat, 2, var, na.rm=T);

      constCol <- (varCol == 0 | is.na(varCol));
      constNum <- sum(constCol, na.rm=T);
      if(constNum > 0){
           msg<-c(msg, paste("<font color=\"red\">", constNum, "columns with constant or a single value were found and deleted.</font>"));
           int.mat <- int.mat[,!constCol];
      }

      # check zero, NA values
      totalCount <-nrow(int.mat)*ncol(int.mat);
      naCount<-sum(is.na(int.mat));
      naPercent<-round(100*naCount/totalCount,1)

      msg<-c(msg, paste("A total of ", naCount, " (", naPercent, "%) missing values were detected.", sep=""));
      msg<-c(msg, "<u>By default, these values will be replaced by a small value.</u>",
                          "Click <b>Skip</b> button if you accept the default practice",
                          "Or click <b>Missing value imputation</b> to use other methods");

      # obtain original half of minimal positive value (threshold)
      minConc<-min(int.mat[int.mat>0], na.rm=T)/2;

      dataSet$minConc<-minConc;
      dataSet$preproc <- as.data.frame(int.mat);
      dataSet$proc.cls <- dataSet$cls <- dataSet$orig.cls;
      if(substring(dataSet$format,4,5)=="ts"){
          dataSet$proc.facA <- dataSet$orig.facA;
          dataSet$proc.facB <- dataSet$orig.facB;
      }

      dataSet$check.msg <- c(dataSet$read.msg, msg);
      dataSet <<- dataSet;
      return(1);
}


GetGroupNumber<-function(){
    return(length(levels(dataSet$cls)));
}

IsSmallSmplSize<-function(){
    return(dataSet$small.smpl.size);
}

GetMinGroupSize<-function(){
    return(dataSet$min.grp.size);
}

IsDataContainsNegative<-function(){
    return(dataSet$containsNegative);
}

################################################################
# Note: the following step directly modifies the dataSet$proc
#################################################################

# replace zero/missing values by half of the minimum pos values, this is the default
# also we will call this method after all missing value imputation if conducted
ReplaceMin<-function(int.mat=as.matrix(dataSet$preproc)){

    minConc<-dataSet$minConc;

    # replace zero and missing values
    # we leave nagative values unchanged! ? not sure if this is the best way
    int.mat[int.mat==0 | is.na(int.mat)] <- minConc;

    # note, this is last step of processing, also save to proc
    dataSet$proc <- as.data.frame(int.mat);
    dataSet$replace.msg <- paste("Zero or missing variables were replaced with a small value:", minConc);
    dataSet <<- dataSet;
    rm(int.mat);
    gc();
}

# remove variable with over certain percentage values are missing
RemoveMissingPercent<-function(int.mat=dataSet$preproc, percent=perct){
    minConc<-dataSet$minConc;
    good.inx<-apply(is.na(int.mat), 2, sum)/nrow(int.mat)<percent;
    dataSet$preproc <- as.data.frame(int.mat[,good.inx]);

    # need to update cls labels
    dataSet$replace.msg <- c(dataSet$replace.msg, paste(sum(!good.inx), "variables were removed for threshold", round(100*percent, 2), "percent."));
    dataSet <<- dataSet;
}

# Replace by min/mean/median/KNN/BPCA/PPCA/SVMImpute
ImputeVar<-function(int.mat=dataSet$preproc, method="min"){

    new.mat<-NULL;
    msg<-dataSet$replace.msg;

    if(method=="exclude"){
        good.inx<-apply(is.na(int.mat), 2, sum)==0
        new.mat<-int.mat[,good.inx];
        msg <- c(msg,"Variables with missing values were excluded.");
    }else if(method=="min"){
        minConc<-dataSet$minConc;
        int.mat[int.mat==0 | is.na(int.mat)] <- minConc;
        new.mat <- int.mat;
        msg <- c(msg,"Variables with missing values were replaced with a small value.");
    }else if (method=="colmin"){
            new.mat<-apply(int.mat, 2, function(x){
                if(sum(is.na(x))>0){
                    x[is.na(x)]<-min(x,na.rm=T)/2;
                }
                x;
            });
        msg <- c(msg,"Missing variables were replaced with the half of minimum values for each feature column.");
    }else if (method=="mean"){
            new.mat<-apply(int.mat, 2, function(x){
                if(sum(is.na(x))>0){
                    x[is.na(x)]<-mean(x,na.rm=T);
                }
                x;
            });
        msg <- c(msg,"Missing variables were replaced with mean.");
    }else if (method == "median"){
            new.mat<-apply(int.mat, 2, function(x){
                if(sum(is.na(x))>0){
                    x[is.na(x)]<-median(x,na.rm=T);
                }
                x;
            });
        msg <- c(msg,"Missing variables were replaced with median.");
    }else {
        if(method == "knn"){
            suppressMessages(require(impute));
            #print("loading for KNN...");
            new.mat<-t(impute.knn(t(int.mat))$data);
        }else{
            suppressMessages(require(pcaMethods));
            if(method == "bpca"){
                new.mat<-pca(int.mat, nPcs =5, method="bpca", center=T)@completeObs;
            }else if(method == "ppca"){
                new.mat<-pca(int.mat, nPcs =5, method="ppca", center=T)@completeObs;
            }else if(method == "svdImpute"){
                new.mat<-pca(int.mat, nPcs =5, method="svdImpute", center=T)@completeObs;
            }
        }
        msg <- c(msg, paste("Missing variables were imputated using", toupper(method)));
    }
    dataSet$proc <- as.data.frame(new.mat);
    dataSet$replace.msg <- msg;
    dataSet <<- dataSet;
}

# to deal with negative values, this is after dealing with negative values
# so operate on dataSet$proc
ClearNegatives <- function(int.mat=as.matrix(dataSet$proc), method="abs"){

    if(dataSet$containsNegative){
        if(method == "min"){
            int.mat[int.mat < 0] <- dataSet$minConc;
            msg <- paste("Negative variables were replaced with a small value:", dataSet$minConc);
        }else if(method =="abs"){
            int.mat <- abs(int.mat);
            msg <- paste("Negative variables were replaced with their absolute values");
        }else{ # exclude
            good.inx<-apply(int.mat<0, 2, sum)==0
            new.mat<-int.mat[,good.inx];
            msg <- paste("Columns contains negative variables were excluded");
        }
  
    dataSet$containsNegative <- 0;
    dataSet$replace.msg <- c(dataSet$replace.msg, msg);
    dataSet$proc <- as.data.frame(int.mat);
    dataSet <<- dataSet;
  }
}

# Group peak list basede on position using xcms algorithm (align peaks wrt rt and mz)
# NMR peaks change ppm -> mz and add dummy rt
# 2-col MS need to add dummy rt
# 3-col MS can be used directly
# default mzwid MS 0.25 m/z, NMR 0.03 ppm
# bw 30 for LCMS, 5 for GCMS
GroupPeakList<-function(mzwid = 0.25, bw = 30, minfrac = 0.5, minsamp = 1, max = 50) {

    peakSet<-dataSet$peakSet;
    samples <- peakSet$sampnames;
    classlabel <- peakSet$sampclass;
    classnames <- levels(classlabel)

    classlabel <- as.vector(unclass(classlabel))
    classnum <- integer(max(classlabel))
    for (i in seq(along = classnum)){
        classnum[i] <- sum(classlabel == i)
    }

    peakmat <- peakSet$peaks;
    porder <- order(peakmat[,"mz"]);
    peakmat <- peakmat[porder,,drop=F]
    rownames(peakmat) <- NULL
    retrange <- range(peakmat[,"rt"])

    minpeakmat <- min(classnum)/2

    mass <- seq(peakmat[1,"mz"], peakmat[nrow(peakmat),"mz"] + mzwid, by = mzwid/2)
    masspos <- findEqualGreaterM(peakmat[,"mz"], mass)

    groupmat <- matrix(nrow = 512, ncol = 7 + length(classnum))
    groupindex <- vector("list", 512)

    endidx <- 0
    num <- 0
    gcount <- integer(length(classnum))
    for (i in seq(length = length(mass)-2)) {
        startidx <- masspos[i]
        endidx <- masspos[i+2]-1
        if (endidx - startidx + 1 < minpeakmat)
            next
        speakmat <- peakmat[startidx:endidx,,drop=FALSE]
        den <- density(speakmat[,"rt"], bw, from = retrange[1]-3*bw, to = retrange[2]+3*bw)
        maxden <- max(den$y)
        deny <- den$y
        gmat <- matrix(nrow = 5, ncol = 2+length(classnum))
        snum <- 0
        while (deny[maxy <- which.max(deny)] > maxden/20 && snum < max) {
            grange <- descendMin(deny, maxy)
            deny[grange[1]:grange[2]] <- 0
            gidx <- which(speakmat[,"rt"] >= den$x[grange[1]] & speakmat[,"rt"] <= den$x[grange[2]])
            gnum <- classlabel[unique(speakmat[gidx,"sample"])]
            for (j in seq(along = gcount))
                gcount[j] <- sum(gnum == j)
            if (! any(gcount >= classnum*minfrac & gcount >= minsamp))
                next
            snum <- snum + 1
            num <- num + 1
            ### Double the size of the output containers if they're full
            if (num > nrow(groupmat)) {
                groupmat <- rbind(groupmat, matrix(nrow = nrow(groupmat), ncol = ncol(groupmat)))
                groupindex <- c(groupindex, vector("list", length(groupindex)))
            }
            groupmat[num, 1] <- median(speakmat[gidx, "mz"])
            groupmat[num, 2:3] <- range(speakmat[gidx, "mz"])
            groupmat[num, 4] <- median(speakmat[gidx, "rt"])
            groupmat[num, 5:6] <- range(speakmat[gidx, "rt"])
            groupmat[num, 7] <- length(gidx)
            groupmat[num, 7+seq(along = gcount)] <- gcount
            groupindex[[num]] <- sort(porder[(startidx:endidx)[gidx]])
        }
    }
    colnames(groupmat) <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax",
                            "npeaks", classnames)

    groupmat <- groupmat[seq(length = num),]
    groupindex <- groupindex[seq(length = num)]

    # Remove groups that overlap with more "well-behaved" groups
    numsamp <- rowSums(groupmat[,(match("npeaks", colnames(groupmat))+1):ncol(groupmat),drop=FALSE])
    uorder <- order(-numsamp, groupmat[,"npeaks"])
    uindex <- rectUnique(groupmat[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE],
                         uorder)

    peakSet$groups <- groupmat[uindex,];
    peakSet$groupidx<- groupindex[uindex];
    dataSet$peakSet<-peakSet;
    dataSet <<- dataSet;
}

# object is nmr.xcmsSet object
SetPeakList.GroupValues<-function() {
    peakSet <- dataSet$peakSet;
    msg<-dataSet$peakMsg;

    peakmat <- peakSet$peaks;
    groupmat <- peakSet$groups;
    groupindex <- peakSet$groupidx;

    sampnum <- seq(length = length(peakSet$sampnames))
    intcol <- match("int", colnames(peakmat))
    sampcol <- match("sample", colnames(peakmat))

    # row is peak, col is sample
    values <- matrix(nrow = length(groupindex), ncol = length(sampnum))

    for (i in seq(along = groupindex)) {
       # for each group, replace multiple peaks from the same sample by their sum
       for(m in sampnum){
            samp.inx<-which(peakmat[groupindex[[i]], sampcol]==m)
            if(length(samp.inx)>0){
                 values[i, m] <- sum(peakmat[groupindex[[i]][samp.inx], intcol]);
            }else{
                 values[i, m] <- NA;
            }
        }
    }

    msg<-c(msg, paste("A total of", length(groupindex), "peak groups were formed. "));
    msg<-c(msg, paste("Peaks of the same group were summed if they are from one sample. "));
    msg<-c(msg, paste("Peaks appear in less than half of samples in each group were ignored."));

    colnames(values) <- peakSet$sampnames;

    if(peakSet$ncol==2){
    	rownames(values) <- paste(round(groupmat[,paste("mz", "med", sep="")],5));
    }else{
    	rownames(values) <- paste(round(groupmat[,paste("mz", "med", sep="")],5), "/", round(groupmat[,paste("rt", "med", sep="")],2), sep="");
        dataSet$three.col <- T;
    }

    dataSet$orig<-t(values);
    dataSet$proc.msg<-msg;
    dataSet$orig.cls<-as.factor(peakSet$sampclass);
    dataSet <<- dataSet;
}

# retention time correction for LC/GC-MS spectra
MSspec.rtCorrection<-function(bw=30){
    xset2<-retcor(dataSet$xset.orig)
    # re-group peaks after retention time correction
    xset2<-group(xset2, bw=bw)
    dataSet$xset.rt<-xset2;
    dataSet <<- dataSet;
}

# plot rentention time corrected spectra
PlotMS.RT<-function(imgName, format="png", dpi=72, width=NA){
    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
    if(is.na(width)){
        w <- 9;
    }else if(width == 0){
        w <- 9;
        imgSet$msrt<-imgName;
        imgSet <<- imgSet;
    }else{
        w <- width;    
    }
    h <- w*7/9;

    plotrt(dataSet$xset.rt);
}

# fill in missing peaks
MSspec.fillPeaks<-function(){
    xset3<-fillPeaks(dataSet$xset.rt);
    dataSet$xset.fill<-xset3;

    msg<-paste("A total of", dim(xset3@peaks)[1],"peaks were detected from these samples");
    msg<-c(msg, paste("with an average of", round(dim(xset3@peaks)[1]/dim(xset3@phenoData)[1], 2), "peaks per spectrum."));
    dataSet$xset.msg<-msg;
    dataSet <<- dataSet;
}

# into:  integrated area of original (raw) peak
# intf:  integrated area of filtered peak
# maxo:  maximum intensity of original (raw) peak
# maxf:  maximum intensity of filtered peak

SetupMSdataMatrix<-function(intvalue = c("into","maxo","intb")){
	values <- groupval(dataSet$xset.fill, "medret", value = intvalue);
    msg<-dataSet$xset.msg;
    # transpose to make row for samples
	orig<-as.data.frame(t(values));
    msg<-dataSet$xset.msg;
    msg=c(msg, paste("These peaks were aligned", dim(orig)[2], "groups according to their mass and retention time."));
    msg=c(msg, paste("Please note, some peaks were excluded if they appear in only a few samples."));
    dataSet$xset.msg<-msg;
    dataSet$orig<-orig;
    dataSet$orig.cls<-as.factor(sampclass(dataSet$xset.fill));
    dataSet <<- dataSet;
}

IsSpectraProcessingOK<-function(){
    msg<-dataSet$xset.msg;
    res = 1;
    if(is.null(dataSet$xset.orig)){
        dataSet$xset.msg<-c(msg, "Failed to read in and process the spectra.");
        res <- 0;
    }
    if(is.null(dataSet$xset.rt)){
        dataSet$xset.msg<-c(msg, "Faiedl in retention time correction, spectra problem?");
        res <- 0;
    }
    if(is.null(dataSet$xset.fill)){
        dataSet$xset.msg<-c(msg, "Failed in filling missing peaks, spectra problem?");
        res <- 0;
    }
    dataSet <<- dataSet;
    return(res);
}

####################################################################
### ========= methods for non-specific filtering of variables====###
####################################################################


# the final variable should be less than 5000 for effective computing
FilterVariable <- function(filter, qcFilter, rsd){
    int.mat <- as.matrix(dataSet$proc);
    cls <- dataSet$proc.cls;

    # save a copy
    dataSet$filt.cls <- cls;
    if(substring(dataSet$format,4,5)=="ts"){
        dataSet$filt.facA <- dataSet$proc.facA; 
        dataSet$filt.facB <- dataSet$proc.facB; 
    }

    msg <- "";
    if(qcFilter == "T"){
        rsd <- rsd/100;
        # need to check if QC exists
        qc.hits <- tolower(as.character(cls)) %in% "qc";
        if(sum(qc.hits) > 2){ # require at least 3 QC for RSD
            qc.mat <- int.mat[qc.hits,];
            sds <- apply(qc.mat, 2, sd, na.rm=T);
            mns <- apply(qc.mat, 2, mean, na.rm=T);
            rsd.vals <- abs(sds/mns);  
            gd.inx <- rsd.vals < rsd;
            int.mat <- int.mat[,gd.inx];
            msg <- paste("Removed ", sum(!gd.inx), " features based on QC RSD values. QC samples are still kept. You can later remove them using Data Editor.");
        }else if(sum(qc.hits) > 0){
            msg <- "RSD requires at least 3 QC samples, and only non-QC based filtering can be applied.";
            return(0);
        }else{
            current.msg <<- "No QC Samples (with class label: QC) found.  Please use non-QC based filtering.";
            return(0);
        }
    }
    
    feat.num <- ncol(int.mat);
    feat.nms <- colnames(int.mat);
    nm <- NULL;
    if(filter == "none" && feat.num <= 4000){ # only allow for less than 4000
        remain <- rep(TRUE, feat.num);
        msg <- paste(msg, "No non-QC based data filtering was applied");
    }else{
        if (filter == "rsd" ){
            sds <- apply(int.mat, 2, sd, na.rm=T);
            mns <- apply(int.mat, 2, mean, na.rm=T);
            filter.val <- abs(sds/mns);
            nm <- "Relative standard deviation";
        }else if (filter == "nrsd" ){
            mads <- apply(int.mat, 2, mad, na.rm=T);
            meds <- apply(int.mat, 2, median, na.rm=T);
            filter.val <- abs(mads/meds);
            nm <- "Non-paramatric relative standard deviation";
        }else if (filter == "mean"){
            filter.val <- apply(int.mat, 2, mean, na.rm=T);
            nm <- "mean";
        }else if (filter == "sd"){
            filter.val <- apply(int.mat, 2, sd, na.rm=T);
            nm <- "standard deviation";
        }else if (filter == "mad"){
            filter.val <- apply(int.mat, 2, mad, na.rm=T);
            nm <- "Median absolute deviation";
        }else if (filter == "median"){
            filter.val <- apply(int.mat, 2, median, na.rm=T);
            nm <- "median";
        }else{ # iqr
            filter.val <- apply(int.mat, 2, IQR, na.rm=T);
            nm <- "Interquantile Range";
        }

        # get the rank of the
        rk <- rank(-filter.val, ties.method='random');

        var.num <- ncol(int.mat);
        if(var.num < 250){ # reduce 5%
            remain <- rk < var.num*0.95;
            msg <- paste(msg, "Further feature filtering based on", nm);
        }else if(ncol(int.mat) < 500){ # reduce 10%
            remain <- rk < var.num*0.9;
            msg <- paste(msg, "Further feature filtering based on", nm);
        }else if(ncol(int.mat) < 1000){ # reduce 25%
            remain <- rk < var.num*0.75;
            msg <- paste(msg, "Further feature filtering based on", nm);
        }else{ # reduce 40%, if still over 8000, then only use top 8000
            remain <- rk < var.num*0.6;
            msg <- paste(msg, "Further feature filtering based on", nm);
            if(sum(remain) > 8000){
                remain <-rk < 8000;
                msg <- paste(msg, "Reduced to 8000 features based on", nm);
            }
        }
    }
    dataSet$filt <- int.mat[, remain];
    dataSet$filter.msg <- msg
    dataSet <<- dataSet;
    current.msg <<- msg;
    return(1);
}

# mat are log normalized, diff will be ratio
CalculatePairwiseDiff <- function(mat){
    f <- function(i, mat) {
       z <- mat[, i-1] - mat[, i:ncol(mat), drop = FALSE]
       colnames(z) <- paste(colnames(mat)[i-1], colnames(z), sep = "/")
       z
    }
    res <- do.call("cbind", sapply(2:ncol(mat), f, mat));
    round(res,5);
}