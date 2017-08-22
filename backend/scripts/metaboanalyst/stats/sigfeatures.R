##################################################
## R script for MetaboAnalyst
## Description: perform SAM and EBAM for feature selection
##
## Author: Jeff Xia, jeff.xia@mcgill.ca
## McGill University, Canada
##
## License: GNU GPL (>= 2)
###################################################

##################################
########### SAM ##################
##################################

# SAM analysis
SAM.Anal<-function(method="d.stat", paired=FALSE, varequal=TRUE){

    suppressMessages(require(siggenes));
    mat<-t(dataSet$norm); # in sam the column is sample
    cl<-as.numeric(dataSet$cls); # change to 0 and 1 for class label
    if(dataSet$cls.num==2){
        if(paired){
            cl<-as.numeric(dataSet$pairs);
        }
        if(method == "d.stat"){
            sam_out<-sam(mat, cl, method=d.stat, var.equal=varequal, R.fold=0, rand=123);
        }else{
            sam_out<-sam(mat, cl, method=wilc.stat, R.fold=0,rand=123);
        }
    }else{
        sam_out<-sam(mat, cl, rand=123);
    }
    analSet$sam<-sam_out;
    analSet <<- analSet;
}

PlotSAM.FDR<-function(delta, imgName, format="png", dpi=72, width=NA){
    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
    if(is.na(width)){
        w <- 10;
    }else if(width == 0){
        w <- 7.2;
        imgSet$sam.fdr<-imgName;
        imgSet <<- imgSet;
    }
    h <- w*3/5;
    
	par(mfrow=c(1,2), mar=c(5,6,4,1));
	mat.fdr<-analSet$sam@mat.fdr;
	plot(mat.fdr[,"Delta"],mat.fdr[,"FDR"],xlab='Delta',ylab=NA,type="b", col='blue', las=2);
        abline(v = delta, lty=3, col="magenta");
        mtext("FDR", side=2, line=5);
        par(mar=c(5,5,4,2))
	plot(mat.fdr[,"Delta"],mat.fdr[,"Called"],xlab='Delta',ylab="Significant feaure No.",type="b", col='blue', las=2);
        abline(v = delta, lty=3, col="magenta");

        hit.inx <- mat.fdr[,"Delta"] <= delta;
        my.fdr <- signif(min(mat.fdr[,"FDR"][hit.inx]), 3);
        my.sigs <- min(mat.fdr[,"Called"][hit.inx]);
        mtext(paste("Delta:", delta, " FDR:", my.fdr, " Sig. cmpds:", my.sigs), line=-2, side = 3, outer = TRUE, font=2)
    
}

SetSAMSigMat<-function(delta){
        sam.sum<-summary(analSet$sam, delta);
        summary.mat<-sam.sum@mat.sig;

        sig.mat <-as.matrix(signif(summary.mat[,-c(1,6)],5));
        write.csv(signif(sig.mat,5),file="sam_sigfeatures.csv");
        analSet$sam.cmpds<-sig.mat;
        analSet$sam.delta<-delta;
    analSet <<- analSet;
}

GetSAMSigMat<-function(){
      return(CleanNumber(analSet$sam.cmpds));
}

GetSAMSigRowNames<-function(){
    rownames(analSet$sam.cmpds);
}

GetSAMSigColNames<-function(){
    colnames(analSet$sam.cmpds);
}

GetSigTable.SAM<-function(){
    GetSigTable(analSet$sam.cmpds, "SAM");
}

PlotSAM.Cmpd<-function(imgName, format="png", dpi=72, width=NA){
    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
    if(is.na(width)){
        w <- 8;
    }else if(width == 0){
        w <- 7;
        imgSet$sam.cmpd<-imgName;
        imgSet <<- imgSet;
    }else{
        w <- width;
    }
    h <- w;

    
	plot(analSet$sam, analSet$sam.delta);
    
}

# obtain a default delta with reasonable number
# of sig features and decent FDR
GetSuggestedSAMDelta<-function(){
    	mat.fdr<-analSet$sam@mat.fdr
        deltaVec <- mat.fdr[,"Delta"];
        fdrVec <- mat.fdr[,"FDR"];
        signumVec <- mat.fdr[,"Called"];
        for(i in 1:length(deltaVec)){
            delta = deltaVec[i];
            fdr = fdrVec[i];
            called = signumVec[i];
            if(called > 0){ # at least 1 significant cmpd
                # check fdr, default threshold 0.01
                # if too many significant compounds, tight up and vice versa
                if(fdr < 0.001){
                     return (delta);
                }else if(fdr < 0.01 & called < 100){
                     return (delta);
                }else if(fdr < 0.05 & called <50){
                    return (delta);
                }else if(fdr < 0.1 & called < 20){
                     return (delta);
                }else if(called < 10){
                    return (delta);
                }
            }
        }
        return (deltaVec[1]); # if no significant found, return the first one
}

#######################################
############# EBAM ####################
#######################################

# deteriming a0, only applicable for z.ebam (default)
EBAM.A0.Init<-function(isPaired, isVarEq){
    suppressMessages(require(siggenes));
    if(isPaired){
        cl.ebam<-as.numeric(dataSet$pairs); 
    }else{
        cl.ebam<-as.numeric(dataSet$cls)-1; # change to 0 and 1 for class label
    }
    conc.ebam<-t(dataSet$norm); # in sam column is sample, row is gene
    ebam_a0<-find.a0(conc.ebam, cl.ebam, var.equal=isVarEq, gene.names = names(dataSet$norm), rand=123);
    analSet$ebam.a0<-ebam_a0;
    analSet <<- analSet;
}

# plot ebam a0 plot also return the analSet$ebam.a0 object so that the suggested a0 can be obtained
PlotEBAM.A0<-function(imgName, format="png", dpi=72, width=NA){
    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
    if(is.na(width)){
        w <- 8;
    }else if(width == 0){
        w <- 7; 
        imgSet$ebam.a0<-imgName;
        imgSet <<- imgSet;
    }
    h <- 3*w/4;
	plot(analSet$ebam.a0);
}

# note: if method is wilcoxon, the A0 and var equal will be ignored
EBAM.Cmpd.Init<-function(method="z.ebam", A0=0, isPaired=FALSE, isVarEq=TRUE){
    if(isPaired){
        cl.ebam<-as.numeric(dataSet$pairs);
    }else{
        cl.ebam<-as.numeric(dataSet$cls)-1;
    }
    conc.ebam<-t(dataSet$norm); # in sam column is sample, row is feature
    if(method=="z.ebam"){
        ebam_out<-ebam(conc.ebam, cl.ebam, method=z.ebam, a0=A0, var.equal=isVarEq, fast=TRUE, gene.names = names(dataSet$norm), rand=123);
    }else{
        ebam_out<-ebam(conc.ebam, cl.ebam, method=wilc.ebam, gene.names = names(dataSet$norm), rand=123);
    }
    analSet$ebam <-ebam_out;
    analSet <<- analSet;
}

# return double matrix with 3 columns - z.value, posterior, local.fdr
SetEBAMSigMat<-function(delta){
    ebam.sum<-summary(analSet$ebam, delta);
    summary.mat<-ebam.sum@mat.sig;
    sig.mat <-as.matrix(signif(summary.mat[,-1],5));
    write.csv(signif(sig.mat,5),file="ebam_sigfeatures.csv");
    analSet$ebam.cmpds <- sig.mat;
    analSet$ebam.delta <- delta;
    analSet <<- analSet;
}

GetEBAMSigMat<-function(){
    return(CleanNumber(analSet$ebam.cmpds));
}

GetEBAMSigRowNames<-function(){
    rownames(analSet$ebam.cmpds);
}

GetEBAMSigColNames<-function(){
    colnames(analSet$ebam.cmpds);
}

GetSigTable.EBAM<-function(){
    GetSigTable(analSet$ebam.cmpds, "EBAM");
}

# plot ebam
PlotEBAM.Cmpd<-function(imgName, format="png", dpi=72, width=NA){
    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
    if(is.na(width)){
        w <- h <- 7;
    }else if(width == 0){
        w <- h <- 7;
        imgSet$ebam.cmpd <-imgName;
        imgSet <<- imgSet;
    }else{
        w <- h <- width;
    }
	plot(analSet$ebam, analSet$ebam.delta);
}