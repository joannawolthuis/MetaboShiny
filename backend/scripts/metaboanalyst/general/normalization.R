##################################################
## R script for MetaboAnalyst
## Description: perform various normalization
##
## Author: Jeff Xia, jeff.xia@mcgill.ca
## McGill University, Canada
##
## License: GNU GPL (>= 2)
###################################################

# get the dropdown list for sample normalization view
GetPrenormSmplNms <-function(){
  if(is.null(dataSet$prenorm)){
        if(is.null(dataSet$filt)){
            dataSet$prenorm <- dataSet$proc;
        }else{
            dataSet$prenorm <- dataSet$filt;
        }
        dataSet <<- dataSet;
    }
    return(rownames(dataSet$prenorm));
}

GetPrenormFeatureNms <- function(){
    if(is.null(dataSet$prenorm)){
        if(is.null(dataSet$filt)){
            dataSet$prenorm <- dataSet$proc;
        }else{
            dataSet$prenorm <- dataSet$filt;
        }
        dataSet <<- dataSet;
    }
    return(colnames(dataSet$prenorm));
}

GetPrenormClsNms <- function(){
    if(is.null(dataSet$prenorm.cls)){
        if(is.null(dataSet$filt)){
            dataSet$prenorm.cls <- dataSet$proc.cls;
            if(substring(dataSet$format,4,5)=="ts"){
                dataSet$prenorm.facA <- dataSet$proc.facA;
                dataSet$prenorm.facB <- dataSet$proc.facB;
            }
        }else{
            dataSet$prenorm.cls <- dataSet$filt.cls;
            if(substring(dataSet$format,4,5)=="ts"){
                dataSet$prenorm.facA <- dataSet$filt.facA;
                dataSet$prenorm.facB <- dataSet$filt.facB;
            }
        }
        dataSet <<- dataSet;
    }
    return(levels(dataSet$prenorm.cls));
}

###############################################################
# remove the sample or feature from data 
# Note: this should happen after processing and before normalization 
# dataSet$proc dataSet$proc.cls (make a copy of this pair for restore)
########################################################

UpdateGroupItems<-function(){
    if(!exists("grp.nm.vec")){
        current.msg <<- "Cannot find the current group names!";
        return (0);
    }

    if(is.null(dataSet$filt)){
        data <- dataSet$proc;
        cls <- dataSet$proc.cls;
        if(substring(dataSet$format,4,5)=="ts"){
            facA <- dataSet$proc.facA;
            facB <- dataSet$proc.facB;
        }
    }else{
        data <- dataSet$filt;
        cls <- dataSet$filt.cls;
        if(substring(dataSet$format,4,5)=="ts"){
            facA <- dataSet$filt.facA;
            facB <- dataSet$filt.facB;
        }
    }
    hit.inx <- cls %in% grp.nm.vec;
    dataSet$prenorm <- CleanDataMatrix(data[hit.inx,]);
    dataSet$prenorm.cls <- factor(cls[hit.inx], levels=grp.nm.vec);

    if(substring(dataSet$format,4,5)=="ts"){
        dataSet$prenorm.facA <- droplevels(factor(facA[hit.inx]));
        dataSet$prenorm.facB <- droplevels(factor(facB[hit.inx]));
    }
    dataSet <<- dataSet;
    current.msg <<- "Successfully updated the group items!";
    return(length(levels(dataSet$prenorm.cls)));
}

CleanDataMatrix <- function(ndata){
    # make sure no costant columns crop up
    varCol <- apply(ndata, 2, var, na.rm=T);
    constCol <- (varCol == 0 | is.na(varCol));
    return(ndata[,!constCol]);
}

UpdateSampleItems<-function(){
    if(!exists("smpl.nm.vec")){
        current.msg <<- "Cannot find the current sample names!";
        return (0);
    }

    if(is.null(dataSet$filt)){
        data <- dataSet$proc;
        cls <- dataSet$proc.cls;
        if(substring(dataSet$format,4,5)=="ts"){
            facA <- dataSet$proc.facA;
            facB <- dataSet$proc.facB;
        }
    }else{
        data <- dataSet$filt;
        cls <- dataSet$filt.cls;
        if(substring(dataSet$format,4,5)=="ts"){
            facA <- dataSet$filt.facA;
            facB <- dataSet$filt.facB;
        }
    }

    hit.inx <- rownames(data) %in% smpl.nm.vec;
    dataSet$prenorm <- CleanDataMatrix(data[hit.inx,]);
    dataSet$prenorm.cls <- as.factor(as.character(cls[hit.inx]));
    if(substring(dataSet$format,4,5)=="ts"){
        dataSet$prenorm.facA <- as.factor(as.character(facA[hit.inx]));
        dataSet$prenorm.facB <- as.factor(as.character(facB[hit.inx]));
    }
    dataSet <<- dataSet;
    current.msg <<- "Successfully updated the sample items!";
    return(length(levels(dataSet$prenorm.cls)));
}

UpdateFeatureItems<-function(){

    if(!exists("feature.nm.vec")){
        current.msg <<- "Cannot find the selected feature names!";
        return (0);
    }

    if(is.null(dataSet$filt)){
        data <- dataSet$proc;
        cls <- dataSet$proc.cls;
        if(substring(dataSet$format,4,5)=="ts"){
            facA <- dataSet$proc.facA;
            facB <- dataSet$proc.facB;
        }
    }else{
        data <- dataSet$filt;
        cls <- dataSet$filt.cls;
        if(substring(dataSet$format,4,5)=="ts"){
            facA <- dataSet$filt.facA;
            facB <- dataSet$filt.facB;
        }
    }

    hit.inx <- colnames(data) %in% feature.nm.vec;
    dataSet$prenorm <- CleanDataMatrix(data[,hit.inx]);
    dataSet$prenorm.cls <- cls; # this is the same
    dataSet <<- dataSet;
    current.msg <<- "Successfully updated the sample items!";
    return (1);
}


Normalization<-function(rowNorm, transNorm, scaleNorm, ref=NULL, ratio=FALSE, ratioNum=20){

    data <- dataSet$proc;
    print(names(dataSet))
    
    if(is.null(dataSet$prenorm.cls)){ # can be so for regression 
        dataSet$prenorm.cls <- dataSet$proc.cls;
    }

    cls <- dataSet$prenorm.cls;

    # note, setup time factor
    if(substring(dataSet$format,4,5)=="ts"){
        nfacA <- dataSet$prenorm.facA;
        nfacB <- dataSet$prenorm.facB;
        dataSet$facA <- nfacA;
        dataSet$facB <- nfacB;
        if(dataSet$design.type =="time" | dataSet$design.type =="time0"){
            # determine time factor and should order first by subject then by each time points
            if(tolower(dataSet$facA.lbl) == "time"){ 
                time.fac <- nfacA;
                exp.fac <- nfacB;
            }else{
                time.fac <- nfacB;
                exp.fac <- nfacA;
            }
            dataSet$time.fac <- time.fac;
            dataSet$exp.fac <- exp.fac;
        }
    }

    colNames <- colnames(data);
    rowNames <- rownames(data);
    print(names(dataSet))
    # row-wise normalization
    if(rowNorm=="QuantileNorm"){
        data<-QuantileNormalize(data);
        # this can introduce constant variables if a variable is 
        # at the same rank across all samples (replaced by its average across all)

        varCol <- apply(data, 2, var, na.rm=T);
        constCol <- (varCol == 0 | is.na(varCol));
        constNum <- sum(constCol, na.rm=T);
        if(constNum > 0){
            print(paste("After quantile normalization", constNum, "columns with constant value were found and deleted."));
            data <- data[,!constCol];
            colNames <- colnames(data);
            rowNames <- rownames(data);
        }
        rownm<-"Quantile Normalization";
    }else if(rowNorm=="ProbNormT"){
        grp.inx <- cls == ref;
        ref.smpl <- apply(data[grp.inx, ], 2, mean);
        data<-t(apply(data, 1, ProbNorm, ref.smpl));
        rownm<-"Probabilistic Quotient Normalization";
    }else if(rowNorm=="ProbNormF"){
        ref.smpl <- data[ref,];
        data<-t(apply(data, 1, ProbNorm, ref.smpl));
        rownm<-"Probabilistic Quotient Normalization";
    }else if(rowNorm=="CompNorm"){
        data<-t(apply(data, 1, CompNorm, ref));
        rownm<-"Normalization by a reference feature";
    }else if(rowNorm=="SumNorm"){
        data<-t(apply(data, 1, SumNorm));
        rownm<-"Normalization to constant sum";
    }else if(rowNorm=="MedianNorm"){
        data<-t(apply(data, 1, MedianNorm));
        rownm<-"Normalization to sample median";
    }else if(rowNorm=="SpecNorm"){
        if(!exists("norm.vec")){
            norm.vec <- rep(1,nrow(data)); # default all same weight vec to prevent error
            print("No sample specific information were given, all set to 1.0");
         }
         rownm<-"Normalization by sample-specific factor";
         data<-data/norm.vec;
    }else{
        # nothing to do
        rownm<-"N/A";
    }

   # use apply will lose dimesion info (i.e. row names and colnames)
   rownames(data)<-rowNames;
   colnames(data)<-colNames;

   # note: row-normed data is based on biological knowledge, since the previous
   # replacing zero/missing values by half of the min positive (a constant) 
   # now may become different due to different norm factor, which is artificial
   # variance and should be corrected again
   #
   # stopped, this step cause troubles
   # minConc<-round(min(data)/2, 5);
   # data[dataSet$fill.inx]<-minConc;

   # if the reference by feature, the feature column should be removed, since it is all 1
    if(rowNorm=="CompNorm" && !is.null(ref)){
        inx<-match(ref, colnames(data));
        data<-data[,-inx];
        colNames <- colNames[-inx];
    }

   # record row-normed data for fold change analysis (b/c not applicable for mean-centered data)
   dataSet$row.norm<-as.data.frame(CleanData(data));

   # this is for biomarker analysis only (for compound concentraion data)
    if(ratio){
        min.val <- min(abs(data[data!=0]))/2;
        norm.data <- log2((data + sqrt(data^2 + min.val))/2);
        transnm<-"Log Normalization";
        ratio.mat <- CalculatePairwiseDiff(norm.data);

        fstats <- Get.Fstat(ratio.mat, cls);
        hit.inx <- rank(-fstats) < ratioNum;  # get top n
        ratio.mat <- ratio.mat[, hit.inx];

        data <- cbind(norm.data, ratio.mat);
        colNames <- colnames(data);
        rowNames <- rownames(data);
    }

    if(!ratio){
        # transformation
        if(transNorm=='LogNorm'){
            min.val <- min(abs(data[data!=0]))/10;
            data<-apply(data, 2, LogNorm, min.val);
            transnm<-"Log Normalization";
        }else if(transNorm=='CrNorm'){
            norm.data <- abs(data)^(1/3);
            norm.data[data<0] <- - norm.data[data<0];
            data <- norm.data;
            transnm<-"Cubic Root Transformation";
        }else{
            transnm<-"N/A";
        }
    }

    # scaling
    if(scaleNorm==''){
            data<-apply(data, 2, MeanCenter);
            scalenm<-"MeanCenter";
    }else if(scaleNorm=='AutoNorm'){
            data<-apply(data, 2, AutoNorm);
            scalenm<-"Autoscaling";
    }else if(scaleNorm=='ParetoNorm'){
            data<-apply(data, 2, ParetoNorm);
            scalenm<-"Pareto Scaling";
    }else if(scaleNorm=='RangeNorm'){
            data<-apply(data, 2, RangeNorm);
            scalenm<-"Range Scaling";
    }else{
            scalenm<-"N/A";
    }

    # note after using "apply" function, all the attribute lost, need to add back
    rownames(data)<-rowNames;
    colnames(data)<-colNames;

    # need to do some sanity check, for log there may be Inf values introduced
    data <- CleanData(data, T, F);

    dataSet$norm <- as.data.frame(data);
    dataSet$cls <- cls;

    dataSet$rownorm.method<-rownm;
    dataSet$trans.method<-transnm;
    dataSet$scale.method<-scalenm;
    dataSet$combined.method<-FALSE;
    dataSet$norm.all <- NULL; # this is only for biomarker ROC analysis
    dataSet <<- dataSet;
    return(1);
}

########################################
###row-wise norm methods, x is a row ###
########################################

# normalize by a sum of each sample, assume constant sum (1000)
# return: normalized data
SumNorm<-function(x){
    1000*x/sum(x, na.rm=T);
}

# normalize by median
MedianNorm<-function(x){
    x/median(x, na.rm=T);
}

# normalize by a reference sample (probability quotient normalization)
# ref should be the name of the reference sample
ProbNorm<-function(x, ref.smpl){
    x/median(as.numeric(x/ref.smpl), na.rm=T)
}

# normalize by a reference reference (i.e. creatinine)
# ref should be the name of the cmpd
CompNorm<-function(x, ref){
	1000*x/x[ref];
}

# perform quantile normalization on the raw data (can be log transformed later by user)
# https://stat.ethz.ch/pipermail/bioconductor/2005-April/008348.html
QuantileNormalize <- function(data){
    require('preprocessCore');
    return(t(normalize.quantiles(t(data), copy=FALSE)));
}

##############################################
###column-wise norm methods, x is a column ###
##############################################

# generalize log, tolerant to 0 and negative values
LogNorm<-function(x,min.val){
	 log2((x + sqrt(x^2 + min.val^2))/2)
}

# normalize to zero mean and unit variance
AutoNorm<-function(x){
	(x - mean(x))/sd(x, na.rm=T);
}

# normalize to zero mean but varaince/SE
ParetoNorm<-function(x){
	(x - mean(x))/sqrt(sd(x, na.rm=T));
}

# normalize to zero mean but varaince/SE
MeanCenter<-function(x){
	x - mean(x);
}

# normalize to zero mean but varaince/SE
RangeNorm<-function(x){
    if(max(x) == min(x)){
        x;
    }else{
        (x - mean(x))/(max(x)-min(x));
    }
}


##############################################
################## Summary plot ##############
##############################################

# plot two summary plot, one b4 normalization, one after
# for each plot top is box plot, bottom is a density plot
PlotNormSummary<-function(imgName=NA, format="png", dpi=72, width=NA){

    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
    if(is.na(width)){
        w <- 10.5; h <- 12;
    }else if(width == 0){
        w <- 7.2;h <- 9;
        imgSet$norm<-imgName;
        imgSet <<- imgSet;
    }else{
        w <- 7.2; h <- 9;
    }

    #
    layout(matrix(c(1,2,2,2,3,4,4,4), 4, 2, byrow = FALSE))

    # since there may be too many compounds, only plot a subsets (50) in box plot
    # but density plot will use all the data

    pre.inx<-GetRandomSubsetIndex(ncol(dataSet$proc), sub.num=50);
    namesVec <- colnames(dataSet$proc[,pre.inx]);

    # only get common ones
    nm.inx <- namesVec %in% colnames(dataSet$norm)
    namesVec <- namesVec[nm.inx];
    pre.inx <- pre.inx[nm.inx];

    norm.inx<-match(namesVec, colnames(dataSet$norm));
    namesVec <- substr(namesVec, 1, 12); # use abbreviated name

    rangex.pre <- range(dataSet$proc[, pre.inx], na.rm=T);
    rangex.norm <- range(dataSet$norm[, norm.inx], na.rm=T);

    x.label<-GetValueLabel();
    y.label<-GetVariableLabel();

    # fig 1
    op<-par(mar=c(4,7,4,0), xaxt="s");
    plot(density(apply(dataSet$proc, 2, mean, na.rm=TRUE)), col='black',lty=3,las =2, lwd=2, main="", xlab="", ylab="");
    mtext("Density", 2, 5);
    mtext("Before Normalization",3, 1)

    # fig 2
    op<-par(mar=c(7,7,0,0), xaxt="s");
    boxplot(dataSet$proc[,pre.inx], names= namesVec, ylim=rangex.pre,las = 2, col="pink", horizontal=T);
    mtext(x.label, 1, 5);

    # fig 3
    op<-par(mar=c(4,7,4,2), xaxt="s");
    plot(density(apply(dataSet$norm, 2, mean, na.rm=TRUE)), col='black',lty=3, las=2, lwd =2, main="", xlab="", ylab="");
    mtext("After Normalization",3, 1);

    # fig 4
    op<-par(mar=c(7,7,0,2), xaxt="s");
    boxplot(dataSet$norm[,norm.inx], names=namesVec, ylim=rangex.norm, las = 2, col="pink", horizontal=T);
    mtext(paste("Normalized",x.label),1, 5);

    #
}

PlotSampleNormSummary<-function(imgName=NA, format="png", dpi=72, width=NA){
    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
    if(is.na(width)){
        w <- 10.5; h <- 12;
    }else if(width == 0){
        w <- 7.2;h <- 9;
        imgSet$norm <-imgName;
        imgSet <<- imgSet;
    }else{
        w <- 7.2; h <- 9;
    }

    #
    layout(matrix(c(1,1,1,2,3,3,3,4), 4, 2, byrow = FALSE))

    # since there may be too many samples, only plot a subsets (50) in box plot
    # but density plot will use all the data

    pre.inx<-GetRandomSubsetIndex(nrow(dataSet$proc), sub.num=50);
    namesVec <- rownames(dataSet$proc[pre.inx,]);

    # only get common ones
    nm.inx <- namesVec %in% rownames(dataSet$norm)
    namesVec <- namesVec[nm.inx];
    pre.inx <- pre.inx[nm.inx];

    norm.inx<-match(namesVec, rownames(dataSet$norm));
    namesVec <- substr(namesVec, 1, 12); # use abbreviated name

    rangex.pre <- range(dataSet$proc[pre.inx,], na.rm=T);
    rangex.norm <- range(dataSet$norm[norm.inx,], na.rm=T);

    x.label<-GetValueLabel();
    y.label<-"Samples";

    # fig 1
    op<-par(mar=c(5,7,4,0), xaxt="s");
    boxplot(t(dataSet$proc[pre.inx, ]), names= namesVec, ylim=rangex.pre, las = 2, col="pink", horizontal=T);
    mtext("Before Normalization",3, 1);

    # fig 2
    op<-par(mar=c(7,7,0,0), xaxt="s");
    plot(density(apply(dataSet$proc, 1, mean, na.rm=TRUE)),col='black',lty=3, las =2, lwd=2, main="", xlab="", ylab="");
    mtext("Density", 2, 5);
    mtext(x.label, 1, 5);

    # fig 3
    op<-par(mar=c(5,7,4,2), xaxt="s");
    boxplot(t(dataSet$norm[norm.inx,]), names=namesVec, ylim=rangex.norm, las = 2, col="pink", horizontal=T);
    mtext("After Normalization",3, 1);

    # fig 4
    op<-par(mar=c(7,7,0,2), xaxt="s");
    plot(density(apply(dataSet$norm, 1, mean, na.rm=TRUE)),col='black',lty=3, las=2, lwd =2, main="", xlab="", ylab="");
    mtext(paste("Normalized",x.label),1, 5);

    #
}