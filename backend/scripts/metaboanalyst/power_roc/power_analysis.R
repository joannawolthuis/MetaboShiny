
InitPowerAnal <- function(clsOpts){

    require(SSPA);
    if(clsOpts == "NA"){
        grp.nms <- levels(dataSet$cls)[1:2];
    }else{
        grp.nms <- strsplit(clsOpts, " vs. ")[[1]];
    }
    inx1 <- which(dataSet$cls==grp.nms[1]);
    inx2 <- which(dataSet$cls==grp.nms[2]);
    stats <- apply(as.matrix(dataSet$norm), 2, function(x) {
        tmp <- try(t.test(x[inx1], x[inx2], paired = dataSet$paired, var.equal = T));
        if(class(tmp) == "try-error") {
            return(NA);
        }else{
            return(tmp$statistic);
        }
    })

    stats <- stats[!is.na(stats)];
    n1 <- length(inx1);
    n2 <- length(inx2);

    pdD <- pilotData(statistics = stats, 
                 samplesize = sqrt(n1+n2), 
                 distribution="t",
                 df=n1+n2-2);
    analSet$power <- list(pdD = pdD);
    analSet <<- analSet;
}

PlotPowerStat<-function(imgName, format="png", dpi=72, width=NA){

    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
    if(is.na(width)){
        w <- 10;
    }else if(width == 0){
        w <- 10;
        imgSet$power<-imgName;
        imgSet <<- imgSet;
    }else{
        w <- width;
    }
    h <- w;
    Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
    require(lattice);
    print(plot(analSet$power$pdD));
    dev.off();
}

GetSampleSizeLadder <- function(maxNum){
    Jpred <- c(3, 6, 10, 16, 24, 40, 60, 100, 150, seq(200, 1000, 100));
    inx <- which(Jpred == min(Jpred[Jpred>=maxNum]))
    return(Jpred[1:inx]);
}

PerformPowerProfiling <- function(fdr.lvl, smplSize){

    res <- round(length(analSet$power$pdD@statistics)/2);
    ssD <- sampleSize(analSet$power$pdD, method="congrad", control=list(from=-6, to=6, resolution=res));
    Jpred <- GetSampleSizeLadder(smplSize);
    N <- sqrt(Jpred/2);

    pi0 <- ssD@pi0;
    if(fdr.lvl >= pi0){
        fdr.lvl <- signif(pi0-pi0/10, 3);
    }
    pwrD <- predictpower(ssD, samplesizes=N, alpha=fdr.lvl)
    analSet$power$ssD <- ssD;
    analSet$power$Jpred <- Jpred;
    analSet$power$pwrD <- pwrD;
    analSet <<- analSet;
    return(fdr.lvl);
}


PlotPowerEffectSize<-function(imgName, format="png", dpi=72, width=NA){
    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
    if(is.na(width)){
        w <- 9;
    }else if(width == 0){
        w <- 8;
        imgSet$power<-imgName;
        imgSet <<- imgSet;
    }else{
        w <- width;
    }
    h <- w*7/9;

    Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
    print(plot(analSet$power$ssD));
    dev.off();
}

PlotPowerProfile<-function(fdr.lvl, smplSize, imgName, format="png", dpi=72, width=NA){

    Jpred <- GetSampleSizeLadder(smplSize);
    N <- sqrt(Jpred/2);
    pwrD <- predictpower(analSet$power$ssD, samplesizes=N, alpha=fdr.lvl)

    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
    if(is.na(width)){
        w <- 9;
    }else if(width == 0){
        w <- 7;
        imgSet$power<-imgName;
        imgSet <<- imgSet;
    }else{
        w <- width;
    }
    h <- w*6/9;

    Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
    plot(Jpred, pwrD, type="n", ylim=c(0,1), ylab="Predicted power", xlab="Sample Size (per group)");
    grid(col = "lightgray", lty = "dotted", lwd = 1);
    lines(Jpred, pwrD, lwd=4, col="orange");
    points(Jpred, pwrD, pch=17);
    dev.off();

    analSet$power$pwrD <- pwrD;
    analSet$power$Jpred <- Jpred;
    analSet <<- analSet;
    return(pwrD);
}

GetPowerValuesX <- function(){
    return(analSet$power$Jpred);
}
