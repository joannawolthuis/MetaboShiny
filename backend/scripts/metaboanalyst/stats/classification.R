#########################################################
## R script for MetaboAnalyst
## Description: perform RandomForest and SVM
##
## Author: Jeff Xia, jeff.xia@mcgill.ca
## McGill University, Canada
##
## License: GNU GPL (>= 2)
###################################################

#######################################
########### Random Forest #############
#######################################

# random forests
# -1 is contant 0: fix current random, 1 random
RF.Anal<-function(treeNum=500, tryNum=7, randomOn=1){
    suppressMessages(library(randomForest));

    # set up random numbers
    if(is.null(dataSet$random.seeds)){
        dataSet$random.seeds <- GetRandomNumbers();
        dataSet$cur.inx <- 0;
        dataSet$rn.seed <- dataSet$random.seeds[1];
    }

    if(randomOn == -1){
        rn.sd <- 123456;
    }else if(randomOn == 0){ # keep current
        rn.sd <- dataSet$rn.seed;
    }else{ # random on
        cur.inx <- dataSet$cur.inx + 1;
        rn.sd <- dataSet$random.seeds[cur.inx];        
        dataSet$cur.inx <- cur.inx;
    }
    set.seed(rn.sd);
    # save the 
    dataSet$rn.seed <- rn.sd;
    dataSet <<- dataSet;

    rf_out<-randomForest(dataSet$norm, dataSet$cls, ntree = treeNum, mtry = tryNum, importance = TRUE, proximity = TRUE);

    # set up named sig table for display
    impmat<-rf_out$importance;
    impmat<-impmat[rev(order(impmat[,"MeanDecreaseAccuracy"])),]
    sigmat<-impmat[,"MeanDecreaseAccuracy", drop=F];
    sigmat<-signif(sigmat, 5);

    write.csv(sigmat,file="randomforests_sigfeatures.csv");
    analSet$rf<-rf_out;
    analSet$rf.sigmat<-sigmat;
    analSet <<- analSet;
}

# plot variable importance ranked by MeanDecreaseAccuracy
PlotRF.Classify<-function(imgName, format="png", dpi=72, width=NA){
    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");

    if(is.na(width)){
        w <- 8;
    }else if(width == 0){
        w <- 8;
        imgSet$rf.cls<-imgName;
       imgSet <<- imgSet;
    }else{
        w <- width;
    }
    h <- w*5/8;

    #par(mfrow=c(2,1));
    par(mar=c(4,4,3,2));
    cols <- rainbow(length(levels(dataSet$cls))+1);
    plot(analSet$rf, main="Random Forest classification", col=cols);
    legend("topright", legend = c("Overall", levels(dataSet$cls)), lty=2, lwd=1, col=cols);

    #PlotConfusion(analSet$rf$confusion);
    

}

# plot variable importance ranked by MeanDecreaseAccuracy
PlotRF.VIP<-function(imgName, format="png", dpi=72, width=NA){
	vip.score <- rev(sort(analSet$rf$importance[,"MeanDecreaseAccuracy"]));
    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
    if(is.na(width)){
        w <- 8;
    }else if(width == 0){
        w <- 7;
        imgSet$rf.imp<-imgName;
        imgSet <<- imgSet;
    }else{
        w <- width;    
    }
    h <- w*7/8;

    PlotImpVar(vip.score,"MeanDecreaseAccuracy");
}

PlotRF.Outlier<-function(imgName, format="png", dpi=72, width=NA){

    cols <- GetColorSchema();
    uniq.cols <- unique(cols);
    legend.nm <- unique(as.character(dataSet$cls));
    dist.res <- outlier(analSet$rf);

    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
    if(is.na(width)){
        w <- 9;
    }else if(width == 0){
        w <- 7.2;
        imgSet$rf.outlier<-imgName;
       imgSet <<- imgSet;
    }else{
        w <- width;
    }
    h <- w*7/9;

    layout(matrix(c(1,2), 1, 2, byrow = TRUE), width=c(4,1));

    op<-par(mar=c(5,5,4,0));
    plot(dist.res, type="h", col=cols, xlab="Samples", xaxt="n", ylab="Outlying Measures", bty="n");

    # add sample names to top 5
    rankres <- rank(-abs(dist.res), ties.method="random");

    inx.x <- which(rankres < 6);
    inx.y <- dist.res[inx.x];
    nms <- names(dist.res)[inx.x];
    text(inx.x, inx.y, nms, pos=ifelse(inx.y >= 0, 3, 1), xpd=T)
 	op<-par(mar=c(5,0,4,1));
    plot.new();
    plot.window(c(0,1), c(0,1));

    legend("center", legend =legend.nm, pch=15, col=uniq.cols);

}

# get the OOB error for the last signif
GetRFOOB<-function(){
    errors = analSet$rf$err.rate;
    nrow = dim(errors)[1];
    signif(errors[nrow, 1],3);
}

GetSigTable.RF<-function(){
    GetSigTable(analSet$rf.sigmat, "Random Forest");
}

# significance measure, double[][]
GetRFSigMat<-function(){
    return(CleanNumber(analSet$rf.sigmat))
}

GetRFSigRowNames<-function(){
    rownames(analSet$rf.sigmat);
}

GetRFSigColNames<-function(){
    colnames(analSet$rf.sigmat);
}

GetRFConf.Table<-function(){
     print(xtable(analSet$rf$confusion, 
        caption="Random Forest Classification Performance"), size="\\scriptsize");
}

# return double[][] confusion matrix
GetRFConfMat<-function(){
	signif(analSet$rf$confusion,3);
}

GetRFConfRowNames<-function(){
	rownames(analSet$rf$confusion);
}

GetRFConfColNames<-function(){
	colnames(analSet$rf$confusion);
}

#######################################
########### R-SVM #####################
#######################################

# recursive SVM for feature selection and classification
RSVM.Anal<-function(cvType){

   ladder = CreateLadder(ncol(dataSet$norm));
   svm.out <- RSVM(dataSet$norm, dataSet$cls, ladder, CVtype=cvType);

   # calculate important features
   ERInd <- max( which(svm.out$Error == min(svm.out$Error)) )
   MinLevel <- svm.out$ladder[ERInd]
   FreqVec <- svm.out$SelFreq[, ERInd]
   SelInd <- which(rank(FreqVec) >= (svm.out$ladder[1]-MinLevel));
   FreqInd<-svm.out$SelFreq[SelInd, ERInd]
   names(FreqInd)<-names(dataSet$norm)[SelInd];

   #create a sig table for display
   sig.var<- rev(sort(FreqInd));
   sig.var<-as.matrix(sig.var); # 1-column matrix
   colnames(sig.var)<-"Freqency";

   write.csv(sig.var,file="svm_sigfeatures.csv");

   # add sorted features frequencies as importance indicator
   svm.out<-append(svm.out, list(sig.mat=sig.var, best.inx=ERInd));
   analSet$svm<-svm.out;
    analSet <<- analSet;
}

# Plot plsda classification performance using different components
PlotRSVM.Classification<-function(imgName, format="png", dpi=72, width=NA){
    res<-analSet$svm$Error;
    edge<-(max(res)-min(res))/100; # expand y uplimit for text

    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
    if(is.na(width)){
        w <- 8;
    }else if(width == 0){
        w <- 7;
        imgSet$svm.class<-imgName;
        imgSet <<- imgSet;
    }else{
        w <- width;
    }
    h <- w*6/8;

    plot(res,type='l',xlab='Number of variables (levels)',ylab='Error Rate',
                ylim = c(min(res)-5*edge, max(res)+18*edge), axes=F,
                main="Recursive SVM classification")
    text(res,labels =paste(100*round(res,3),'%'), adj=c(-0.3, -0.5), srt=45, xpd=T)

    points(res, col=ifelse(1:length(res)==analSet$svm$best.inx,"red","blue"));
    axis(2);
    axis(1, 1:length(res), names(res));
}

# if too many, plot top 15
PlotRSVM.Cmpd<-function(imgName, format="png", dpi=72, width=NA){
    sigs<-analSet$svm$sig.mat;
    data<-sigs[,1];
    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
    if(is.na(width)){
        w <- 8;
    }else if(width == 0){
        w <- 7;
        imgSet$svm <-imgName;
        imgSet <<- imgSet;
    }else{
        w <- width;
    }
    h <- w*7/8;

    PlotImpVar(data,"Frequency");
}

GetSigTable.SVM<-function(){
    GetSigTable(analSet$svm$sig.mat, "Recursive SVM");
}

# significance measure, double[][]
GetSVMSigMat<-function(){
    return(CleanNumber(analSet$svm$sig.mat));
}

GetSVMSigRowNames<-function(){
    rownames(analSet$svm$sig.mat);
}

GetSVMSigColNames<-function(){
    colnames(analSet$svm$sig.mat);
}

### R-code for R-SVM
### use leave-one-out / Nfold or bootstrape to permute data for external CV
### build SVM model and use mean-balanced weight to sort genes on training set
### and recursive elimination of least important genes
### author: Dr. Xin Lu, Research Scientist
### Biostatistics Department, Harvard School of Public Health

## create a decreasing ladder for recursive feature elimination
CreateLadder <- function(Ntotal, Nmin=5 ){
    x <- vector()
    x[1] <- Ntotal
    # note SVM is very computationally intensive, large step first 
    # first descend with 0.5 -> 50 var left
    # then descend with 0.6 -> 25 var left
    # then desend with 0.75 -> 5 var

    for( i in 1:100 ){
        if(x[i]>200){
            pRatio = 0.4
        }else if(x[i]>50){
            pRatio = 0.5
        }else if(x[i]>25){
            pRatio = 0.6
        }else{
            pRatio = 0.75
        }
        pp <- round(x[i] * pRatio)
        if( pp == x[i] ){
            pp <- pp-1
        }
        if( pp >= Nmin ) {
            x[i+1] <- pp
        } else{
            break
        }
    }
    x
}

## R-SVM core code
## input:
##    x: row matrix of data
##    y: class label: 1 / -1 for 2 classes
##    CVtype:
##        integer: N fold CV
##        "LOO":    leave-one-out CV
##        "bootstrape": bootstrape CV
##    CVnum:   number of CVs
##        LOO: defined as sample size
##        Nfold and bootstrape:  user defined, default as sample size
## output: a named list
##    Error: a vector of CV error on each level
##    SelFreq: a matrix for the frequency of each gene being selected in each level
##             with each column corresponds to a level of selection
##             and each row for a gene
##          The top important gene in each level are those high-freqent ones
RSVM <- function(x, y, ladder, CVtype, CVnum=0 ){
    suppressMessages(require(e1071));
    ## check if y is binary response
    Ytype <- names(table(y))
    if( length(Ytype) != 2)
    {
        print("ERROR!! RSVM can only deal with 2-class problem")
        return(0)
    }

    ## class mean
    m1 <- apply(x[ which(y==Ytype[1]), ], 2, mean)
    m2 <- apply(x[ which(y==Ytype[2]), ], 2, mean)
    md <- m1-m2

    yy <- vector( length=length(y))
    yy[which(y==Ytype[1])] <- 1
    yy[which(y==Ytype[2])] <- -1
    y <- yy

    ## check ladder
    if( min(diff(ladder)) >= 0 )
    {
        print("ERROR!! ladder must be monotonously decreasing")
        return(0);
    }

    if( ladder[1] != ncol(x) )
    {
        ladder <- c(ncol(x), ladder)
    }

    nSample <- nrow(x)
    nGene   <- ncol(x)
    SampInd <- seq(1, nSample)

    if( CVtype == "LOO" )
    {
        CVnum <- nSample
    } else
    {
        if( CVnum == 0 )
        {
            CVnum <- nSample
        }
    }

    ## vector for test error and number of tests
    ErrVec <- vector( length=length(ladder))
    names(ErrVec) <- as.character(ladder);
    nTests <- 0

    SelFreq <- matrix( 0, nrow=nGene, ncol=length(ladder))
    colnames(SelFreq) <- paste("Level", ladder);

    ## for each CV
    for( i in 1:CVnum )
    {
        ## split data
        if( CVtype == "LOO" )
        {
            TestInd <- i
            TrainInd <- SampInd[ -TestInd]
        } else {
            if( CVtype == "bootstrape" ) {
                TrainInd <- sample(SampInd, nSample, replace=T);
                TestInd <- SampInd[ which(!(SampInd %in% TrainInd ))];
            } else {
                ## Nfold
                TrainInd <- sample(SampInd, nSample*(CVtype-1)/CVtype);
                TestInd <- SampInd[ which(!(SampInd %in% TrainInd ))];
            }
        }

        nTests <- nTests + length(TestInd)

        ## in each level, train a SVM model and record test error
        xTrain <- x[TrainInd, ]
        yTrain <- y[TrainInd]

        xTest  <- x[TestInd,]
        yTest  <- y[TestInd]

        ## index of the genes used in the
        SelInd <- seq(1, nGene)
        for( gLevel in 1:length(ladder) )
        {
            ## record the genes selected in this ladder
            SelFreq[SelInd, gLevel] <- SelFreq[SelInd, gLevel] +1

            ## train SVM model and test error
            ###################################################################################
            ## note the scale is changed to T or it never returns sometime for unscaled data ###
            ## note: the classification performance is idenpendent of about scale is T/F  #####
            ## for "LOO", the test data should be as.data.frame, matrxi will trigger error #####
            ###################################################################################
             svmres <- svm(xTrain[, SelInd], yTrain, scale=T, type="C-classification", kernel="linear" )
             if( CVtype == "LOO" ){
                 svmpred <- predict(svmres, as.data.frame(xTest[SelInd], nrow=1) )
             }else{
                 svmpred <- predict(svmres, xTest[, SelInd] )
             }
             ErrVec[gLevel] <- ErrVec[gLevel] + sum(svmpred != yTest )

            ## weight vector
             W <- t(svmres$coefs*yTrain[svmres$index]) %*% svmres$SV * md[SelInd]
             rkW <- rank(W)

             if( gLevel < length(ladder) ){
                SelInd <- SelInd[which(rkW > (ladder[gLevel] - ladder[gLevel+1]))]
             }
        }
    }
    ret <- list(ladder=ladder, Error=ErrVec/nTests, SelFreq=SelFreq);
    ret;
}

PlotConfusion <- function(clsConf){
    prior(clsConf) <- 100 
    # The above rescales the confusion matrix such that columns sum to 100.
    opar <- par(mar=c(5.1, 6.1, 2, 2))
    x <- x.orig <- unclass(clsConf)
    x <- log(x + 0.5) * 2.33
    x[x < 0] <- NA
    x[x > 10] <- 10
    diag(x) <- -diag(x)
    image(1:ncol(x), 1:ncol(x),
        -(x[, nrow(x):1]), xlab='Actual', ylab='',
        col=colorRampPalette(c(hsv(h = 0, s = 0.9, v = 0.9, alpha = 1), 
                             hsv(h = 0, s = 0, v = 0.9, alpha = 1), 
                             hsv(h = 2/6, s = 0.9, v = 0.9, alpha = 1)))(41), 
        xaxt='n', yaxt='n', zlim=c(-10, 10))
    axis(1, at=1:ncol(x), labels=colnames(x), cex.axis=0.8)
    axis(2, at=ncol(x):1, labels=colnames(x), las=1, cex.axis=0.8)
    title(ylab='Predicted', line=4.5)
    abline(h = 0:ncol(x) + 0.5, col = 'gray')
    abline(v = 0:ncol(x) + 0.5, col = 'gray')
    text(1:6, rep(6:1, each=6), labels = sub('^0$', '', round(c(x.orig), 0)))
    box(lwd=2)
    par(opar) # reset par
}
