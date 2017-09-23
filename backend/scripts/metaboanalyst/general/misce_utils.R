##################################################
## R script for MetaboAnalyst
## Description: some misceleneous tasks
##
## Author: Jeff Xia, jeff.xia@mcgill.ca
## McGill University, Canada
##
## License: GNU GPL (>= 2)
###################################################

# merge duplicated columns or rows by their mean
# dim 1 => row,  dim 2 => column
MergeDuplicates <- function(data, dim=2){

    if(is.null(dim(data))){ # a vector
        if(is.null(names(data))){
            print("Cannot detect duplicate data without names!!!");
            return();
        }
        nm.cls <- as.factor(names(data));
        uniq.len <- length(levels(nm.cls));
        if(uniq.len == length(data)){
            return(data);
        }
        new.data <- vector (mode="numeric",length=uniq.len);
        for(i in 1:uniq.len){
            dup.inx <- nm.cls == levels(nm.cls)[i];
            new.data[i] <- mean(data[dup.inx]);
        }
        names(new.data) <- levels(nm.cls);
        rem.len <- length(data) - length(new.data);
    }else{
        if(dim == 1){
            data <- t(data);
        }
        if(is.null(colnames(data))){
            print("Cannot detect duplicate data without var names!!!");
            return();
        }

        nm.cls <- as.factor(colnames(data));
        uniq.len <- length(levels(nm.cls));

        if(uniq.len == ncol(data)){
            if(dim == 1){
                data <- t(data);
            }
            return(data);
        }

        new.data <- matrix (nrow=nrow(data), ncol=uniq.len);
        for(i in 1:uniq.len){
            dup.inx <- which(nm.cls == levels(nm.cls)[i]);
            new.data[,i] <- apply(data[,dup.inx, drop=F], 1, mean);
        }
        rownames(new.data) <- rownames(data);
        colnames(new.data) <- levels(nm.cls);

        rem.len <- ncol(data) - ncol(new.data);
        if(dim == 1){
            new.data <- t(new.data);
        }
    }
    print(paste(rem.len, "duplicates are merged to their average"));
    new.data;
}

# given a data with duplicates, dups is the one with duplicates
RemoveDuplicates <- function(data, lvlOpt="mean", quiet=T){
    
    all.nms <- rownames(data);
    colnms <- colnames(data);
    dup.inx <- duplicated(all.nms);
    dim.orig  <- dim(data);
    data <- apply(data, 2, as.numeric); # force to be all numeric
    dim(data) <- dim.orig; # keep dimension (will lost when only one item) 
    rownames(data) <- all.nms;
    colnames(data) <- colnms;
    if(sum(dup.inx) > 0){
        uniq.nms <- all.nms[!dup.inx];
        uniq.data <- data[!dup.inx,,drop=F];

        dup.nms <- all.nms[dup.inx];
        uniq.dupnms <- unique(dup.nms);
        uniq.duplen <- length(uniq.dupnms);

        for(i in 1:uniq.duplen){
            nm <- uniq.dupnms[i];
            hit.inx.all <- which(all.nms == nm);
            hit.inx.uniq <- which(uniq.nms == nm);

            # average the whole sub matrix 
            if(lvlOpt == "mean"){
                uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, mean, na.rm=T);
            }else if(lvlOpt == "median"){
                uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, median, na.rm=T);
            }else if(lvlOpt == "max"){
                uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, max, na.rm=T);
            }else{ # sum
                uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, sum, na.rm=T);
            }
        }
        if(!quiet){
            current.msg <<- paste(current.msg, paste("A total of ", sum(dup.inx), " of duplicates were replaced by their ", lvlOpt, ".", sep=""), collapse="\n");
        }
        return(uniq.data);
    }else{
        if(!quiet){
            current.msg <<- paste(current.msg, "All IDs are unique.", collapse="\n");
        }
        return(data);
    }
} 

# from two column input text to data matrix (single column data frame)
getDataFromTextInput <- function(txtInput, sep.type="space"){
    lines <- unlist(strsplit(txtInput, "\r|\n|\r\n")[1]);
    if(substring(lines[1],1,1)=="#"){
        lines <- lines[-1];
    }

    # separated by tab 
    if(sep.type=="tab"){
        my.lists <- strsplit(lines, "\\t");
    }else{ # from any space
        my.lists <- strsplit(lines, "\\s+");
    }
    my.mat <- do.call(rbind, my.lists);

    if(dim(my.mat)[2] == 1){ # add 0
        my.mat <- cbind(my.mat, rep(0, nrow(my.mat)));
    }else if(dim(my.mat)[2] > 2){
        my.mat <- my.mat[,1:2];
        current.msg <- "More than two columns found in the list. Only first two columns will be used. ";
    }
    rownames(my.mat) <- data.matrix(my.mat[,1]);
    my.mat <- my.mat[,-1, drop=F];
    return(my.mat);
}

# use single core on the public server
Perform.permutation <- function(perm.num, fun){
   print(paste("performing", perm.num, "permutations ..."));
   #suppressMessages(require('multicore'));
   #core.num <- multicore:::detectCores();

   #if(core.num > 1){ # use two CPUs only, otherwise, the server will be unresponsive for other users
   #    perm.res <- mclapply(2:perm.num, fun, mc.cores =core.num-1);
   #}else{ # just regular
       perm.res <- lapply(2:perm.num,fun);
   #}
   perm.res;
}

`%fin%` <- function(x, table) {
  fmatch(x, table, nomatch = 0L) > 0L
}

# create semitransparant colors for a given class label
CreateSemiTransColors <- function(cls){

    # note, the first color (red) is for QC
    col.nms <- rainbow(length(levels(cls)));

    # convert to semi-transparent
    semi.nms <- ToSemiTransParent(col.nms);

    # now expand to the one-to-one match to cls element
    col.vec <- vector(mode="character", length=length(cls));
    for (i in 1:length(levels(cls))){
        lv <- levels(cls)[i];
        col.vec[cls==lv] <- semi.nms[i];
    }
    return(col.vec);
}

# convert rgb color i.e. "#00FF00FF" to semi transparent
ToSemiTransParent <- function (col.nms, alpha=0.5){
    rgb.mat <- t(col2rgb(col.nms));
    rgb(rgb.mat/255, alpha=alpha);
}

# col.vec should already been created
UpdateGraphSettings <- function(){
    grpnms <- GetGroupNames();
    names(colVec) <- grpnms;
    names(shapeVec) <- grpnms;
    colVec <<- colVec;
    shapeVec <<- shapeVec;
}

GetShapeSchema <- function(show.name, grey.scale){
    if(exists("shapeVec") && all(shapeVec > 0)){
        sps <- rep(0, length=length(dataSet$cls));
        clsVec <- as.character(dataSet$cls)
        grpnms <- names(shapeVec);
        for(i in 1:length(grpnms)){
            sps[clsVec == grpnms[i]] <- shapeVec[i];
        }
        shapes <- sps;
    }else{
        if(show.name | grey.scale){
            shapes <- as.numeric(dataSet$cls)+1;
        }else{
            shapes <- rep(19, length(dataSet$cls));
        }
    }
    return(shapes);
}

GetColorSchema <- function(grayscale=F){
    # test if total group number is over 9
     grp.num <- length(levels(dataSet$cls));

     if(grayscale){
        dist.cols <- colorRampPalette(c("grey90", "grey30"))(grp.num);
        lvs <- levels(dataSet$cls);
        colors <- vector(mode="character", length=length(dataSet$cls));
        for(i in 1:length(lvs)){
            colors[dataSet$cls == lvs[i]] <- dist.cols[i];
        }
     }else if(grp.num > 9){
        pal12 = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                    "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
                    "#FFFF99", "#B15928");
        dist.cols <- colorRampPalette(pal12)(grp.num);
        lvs <- levels(dataSet$cls);
        colors <- vector(mode="character", length=length(dataSet$cls));
        for(i in 1:length(lvs)){
            colors[dataSet$cls == lvs[i]] <- dist.cols[i];
        }
     }else{
        if(exists("colVec") && !any(colVec =="#NA") ){
            cols <- vector(mode="character", length=length(dataSet$cls));
            clsVec <- as.character(dataSet$cls)
            grpnms <- names(colVec);
            for(i in 1:length(grpnms)){
                cols[clsVec == grpnms[i]] <- colVec[i];
            }
            colors <- cols;
        }else{
            colors <- as.numeric(dataSet$cls)+1;
        }
     }
    return (colors);
}

# unzip the uploaded .zip files, remove the uploaded file, check for success
UnzipUploadedFile<-function(inPath, outPath, rmFile=T){
     #  a<-unzip(inPath, exdir=outPath);

       a<-try(system(paste("unzip",  "-o", inPath, "-d", outPath), intern=T));
       if(class(a) == "try-error" | !length(a)>0){
            AddErrMsg("Failed to unzip the uploaded files!");
            AddErrMsg("Possible reason: file name contains space or special characters.");
            AddErrMsg("Use only alphabets and numbers, make sure there is no space in your file name.");
            AddErrMsg("For WinZip 12.x, use \"Legacy compression (Zip 2.0 compatible)\"");
            return (0);
       }
       if(rmFile){
            RemoveFile(inPath);
       }
       return(1);
}

# clean data and remove -Inf, Inf, NA, negative and 0
CleanData <-function(bdata, 
                     removeNA=T, 
                     removeNeg=T, 
                     removeConst=T){

    if(sum(bdata==Inf, na.rm=TRUE)>0){
        inx <- bdata == Inf;
        bdata[inx] <- NA;
        bdata[inx] <- max(bdata, na.rm=T)*2
    }
    if(sum(bdata==-Inf, na.rm=TRUE)>0){
        inx <- bdata == -Inf;
        bdata[inx] <- NA;
        bdata[inx] <- min(bdata, na.rm=T)/2
    }
    if(removeNA){
        if(sum(is.na(bdata))>0){
            bdata[is.na(bdata)] <- min(bdata, na.rm=T)/2
        }
    }
    if(removeNeg){
        if(sum(bdata<=0) > 0){
            inx <- bdata <= 0;
            bdata[inx] <- NA;
            bdata[inx] <- min(bdata, na.rm=T)/2
        }
    }
    if(removeConst){
        varCol <- apply(bdata, 2, var, na.rm=T);
        constCol <- (varCol == 0 | is.na(varCol));
        constNum <- sum(constCol, na.rm=T);
        if(constNum > 0){
           bdata <- bdata[,!constCol];
        }
    }
    bdata;
}

# replace -Inf, Inf to 99999 and -99999
CleanNumber <-function(bdata){
    if(sum(bdata==Inf)>0){
        inx <- bdata == Inf;
        bdata[inx] <- NA;
        bdata[inx] <- 999999;
    }
    if(sum(bdata==-Inf)>0){
        inx <- bdata == -Inf;
        bdata[inx] <- NA;
        bdata[inx] <- -999999;
    }
    bdata;
}

# remove file
RemoveFolder<-function(folderName){
     # a<-unzip(inPath, exdir=outPath);
       a<-system(paste("rm",  "-r", folderName), intern=T);
       if(!length(a)>0){
            AddErrMsg(paste("Could not remove file -", folderName));
            return (0);
       }
       return(1);
}

# remove files
RemoveFile<-function(fileName){
    if(file.exists(fileName)){
        file.remove(fileName);
    }
}

# clear the current folder and objects in memory
ClearUserDir<-function(){
    # remove physical files
    unlink(dir(), recursive=T);
    dataSet <<- analSet <<- list();
    gc();
}

# utils to remove from
# within, leading and trailing spaces
ClearStrings<-function(query){
    # kill multiple white space
    query <- gsub(" +"," ",query);
    # remove leading and trailing space
    query<- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", query, perl=TRUE);
    return (query);
}

# remove HTML tag
# add escape for %
PrepareLatex <- function(stringVec){
    stringVec <- gsub("<(.|\n)*?>","",stringVec);
    stringVec <- gsub("%", "\\\\%", stringVec);
    stringVec;
}

# get last command from the Rhistory.R file
GetCMD<-function(regexp){
    # store all lines into a list object
    all.lines<-readLines("Rhistory.R");

    all.matches<-grep(regexp, all.lines, value=T);
    if(length(all.matches)==0){
        return(NULL);
    }else{
        # only return the last command
        return(all.matches[length(all.matches)]);
    }
}

# determine value label for plotting
GetValueLabel<-function(){
	if(dataSet$type=="conc"){
	    return("Concentration");
	}else {
	    return("Intensity");
	}
}

# determine variable label for plotting
GetVariableLabel<-function(){
	if(dataSet$type=="conc"){
	    return("Compounds");
	}else if(dataSet$type=="specbin"){
	    return("Spectra Bins");
	}else if(dataSet$type=="nmrpeak"){
	    return("Peaks (ppm)");
	}else if(dataSet$type=="mspeak"){
        if(dataSet$peakSet$ncol==2){
            return("Peaks (mass)");
        }else{
            return("Peaks (mz/rt)");
        }
    }else{
        return("Peaks(mz/rt)");
    }
}

# determine the number of rows and columns for a given total
# number of plots (used by Kmeans and SOM plots)
GetXYCluster<-function(total){
    if(total>16){
        ncol<-4;
        nrow<-5;
    }else if(total>12){
       ncol<-4;
       nrow<-4;
    }else if(total>9){
       ncol<-3;
       nrow<-4;
    }else if(total>6){
       ncol<-3;
       nrow<-3;
    }else if(total>4){
       ncol<-2;
       nrow<-3;
    }else{
       ncol<-1;
       nrow<-total;
    }
    c(nrow, ncol);
}

###################################################
### ====== utility classes for peak grouping=== ###
###################################################
rectUnique <- function(m, order = seq(length = nrow(m)), xdiff = 0, ydiff = 0) {

    nr <- nrow(m)
    nc <- ncol(m)
    if (!is.double(m))
        m <- as.double(m)
    .C("RectUnique",
       m,
       as.integer(order-1),
       nr,
       nc,
       as.double(xdiff),
       as.double(ydiff),
       logical(nrow(m)),
       DUP = FALSE, PACKAGE = "xcms")[[7]]
}

findEqualGreaterM <- function(x, values) {

    if (!is.double(x)) x <- as.double(x)
    if (!is.double(values)) values <- as.double(values)
    .C("FindEqualGreaterM",
       x,
       length(x),
       values,
       length(values),
       index = integer(length(values)),
       DUP = FALSE, PACKAGE = "xcms")$index + 1
}

descendMin <- function(y, istart = which.max(y)) {

    if (!is.double(y)) y <- as.double(y)
    unlist(.C("DescendMin",
              y,
              length(y),
              as.integer(istart-1),
              ilower = integer(1),
              iupper = integer(1),
              DUP = FALSE, PACKAGE = "xcms")[4:5]) + 1
}

# obtain a random subset of numbers from a total number
GetRandomSubsetIndex<-function(total, sub.num = 50){
    if(total < sub.num){
        1:total;
    }else{
        sample(1:total, sub.num);
    }
}

Get.Accuracy <- function(cm) {
    sum(diag(cm)) / sum(cm);
}

# Get a subsets of data ranked by their p values from t tests
GetTTSubsetIndex<-function(data = dataSet$norm, sub.num=50){
    if(ncol(data) < sub.num){
        1:ncol(data);
    }else{
        if(is.null(analSet$tt)){
            Ttests.Anal(0.75);
        }
        all.lod <- -log10(analSet$tt$p.value);
        
        sub.inx <-order(all.lod, decreasing = T)[1:sub.num];
        sel.inx <- 1:ncol(data) %in% sub.inx;
        sel.inx;
    }
}

# generate Latex table
GetSigTable<-function(mat, method){
    suppressMessages(require(xtable));
    if(!isEmptyMatrix(mat)){ # test if empty
        cap<-"Important features identified by";
        if(nrow(mat)>50){
            smat<-as.matrix(mat[1:50,]); # only print top 50 if too many
            colnames(smat)<-colnames(mat); # make sure column names are also copied
            mat<-smat;
            cap<-"Top 50 features identified by";
        }
        # change the rowname to first column
        col1<-rownames(mat);
        cname<-colnames(mat);
        cname<-c(GetVariableLabel(), cname);
        mat<-cbind(col1, mat);
        rownames(mat)<-NULL;
        colnames(mat)<-cname;
        print(xtable(mat, caption=paste(cap, method)), ,caption.placement="top", size="\\scriptsize");
    }else{
        print(paste("No significant features were found using the given threshold for", method));
    }
}

# test if a sig table matrix is empty
isEmptyMatrix<-function(mat){
    if(is.null(mat) | length(mat)==0){
        return(TRUE);
    }
    if(nrow(mat)==0 | ncol(mat)==0){
        return(TRUE);
    }
    if(is.na(mat[1,1])){
        return(TRUE);
    }
    return(FALSE);
}

# Compute BSS/WSS for each row of a matrix which may have NA
# Columns have labels
# x is a numeric vector,
# cl is consecutive integers
Get.bwss<-function(x, cl){
   K <- max(cl) - min(cl) + 1
   tvar <- var.na(x);
   tn <- sum(!is.na(x));
   wvar <- wn <- numeric(K);

   for(i in (1:K)) {
     if(sum(cl == (i + min(cl) - 1)) == 1){
        wvar[i] <- 0;
        wn[i] <- 1;
     }

     if(sum(cl == (i + min(cl) - 1)) > 1) {
        wvar[i] <- var.na(x[cl == (i + min(cl) - 1)]);
        wn[i] <- sum(!is.na(x[cl == (i + min(cl) - 1)]));
     }
   }

   WSS <- sum.na(wvar * (wn - 1));
   TSS <- tvar * (tn - 1)
   (TSS - WSS)/WSS;
}


# Compute SSQ for each row of a matrix which may have NA
# Columns have labels cl=consecutive integers
# note: this is desgined for ASCA parition data
# in which Within group (WSS) is
# zero, so, we only need TSS

Get.tss<-function(x, cl){

   K <- max(cl) - min(cl) + 1
   tvar <- apply(x, 1, var.na);
   tn <- apply(!is.na(x), 1, sum);
   wvar <- matrix(0, nrow(x), K);
   wn <- matrix(0, nrow(x), K);

   for(i in (1:K)) {
     if(sum(cl == (i + min(cl) - 1)) == 1){
        wvar[, i] <- 0;
        wn[, i] <- 1;
     }

     if(sum(cl == (i + min(cl) - 1)) > 1) {
        wvar[, i] <- apply(x[, cl == (i + min(cl) - 1)], 1, var.na);
        wn[, i] <- apply(!is.na(x[, cl == (i + min(cl) - 1)]), 1, sum);
     }
   }

   WSS <- apply(wvar * (wn - 1), 1, sum.na)
   TSS <- tvar * (tn - 1)

   return(TSS);
}

sum.na <- function(x,...){
     res <- NA
     tmp <- !(is.na(x) | is.infinite(x))
     if(sum(tmp) > 0)
           res <- sum(x[tmp])
     res
}

var.na <- function(x){
        res <- NA
        tmp <- !(is.na(x) | is.infinite(x))
        if(sum(tmp) > 1){
              res <- var(x[tmp])
         }
        res
}

#######################################################
## calculate Fisher's Least Significant Difference (LSD)
## adapted from the 'agricolae' package
##############################################
LSD.test <- function (y, trt, alpha = 0.05){
    clase<-c("aov","lm")
    name.y <- paste(deparse(substitute(y)))
    name.t <- paste(deparse(substitute(trt)))
    if("aov"%in%class(y) | "lm"%in%class(y)){
        A<-y$model
        DFerror<-df.residual(y)
        MSerror<-deviance(y)/DFerror
        y<-A[,1]
        ipch<-pmatch(trt,names(A))
        name.t <-names(A)[ipch]
        trt<-A[,ipch]
        name.y <- names(A)[1]
    }
    junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
    means <- tapply.stat(junto[, 1], junto[, 2], stat="mean") #change
    sds <- tapply.stat(junto[, 1], junto[, 2], stat="sd")     #change
    nn <- tapply.stat(junto[, 1], junto[, 2], stat="length")  #change
    std.err <- sds[, 2]/sqrt(nn[, 2])
    Tprob <- qt(1 - alpha/2, DFerror)
    LCL <- means[,2]-Tprob*std.err
    UCL <- means[,2]+Tprob*std.err
    means <- data.frame(means, std.err, replication = nn[, 2], LCL, UCL)
    names(means)[1:2] <- c(name.t, name.y)
    #row.names(means) <- means[, 1]
    ntr <- nrow(means)
    nk <- choose(ntr, 2)
    nr <- unique(nn[, 2])

    comb <- combn(ntr, 2)
    nn <- ncol(comb)
    dif <- rep(0, nn)
    LCL1<-dif
	UCL1<-dif
    sig<-NULL
    pvalue <- rep(0, nn)
    for (k in 1:nn) {
        i <- comb[1, k]
        j <- comb[2, k]
        if (means[i, 2] < means[j, 2]){
            comb[1, k]<-j
            comb[2, k]<-i
        }
        dif[k] <- abs(means[i, 2] - means[j, 2])
        sdtdif <- sqrt(MSerror * (1/means[i, 4] + 1/means[j,4]))
        pvalue[k] <- 2 * (1 - pt(dif[k]/sdtdif, DFerror));
        pvalue[k] <- round(pvalue[k],6);
        LCL1[k] <- dif[k] - Tprob*sdtdif
		UCL1[k] <- dif[k] + Tprob*sdtdif
        sig[k]<-" "
        if (pvalue[k] <= 0.001) sig[k]<-"***"
        else  if (pvalue[k] <= 0.01) sig[k]<-"**"
        else  if (pvalue[k] <= 0.05) sig[k]<-"*"
        else  if (pvalue[k] <= 0.1) sig[k]<-"."
    }
    tr.i <- means[comb[1, ],1]
    tr.j <- means[comb[2, ],1]
    output<-data.frame("Difference" = dif, pvalue = pvalue,sig,LCL=LCL1,UCL=UCL1)
    rownames(output)<-paste(tr.i,tr.j,sep=" - ");
	output;
}

tapply.stat <-function (y, x, stat = "mean"){
    cx<-deparse(substitute(x))
    cy<-deparse(substitute(y))
    x<-data.frame(c1=1,x)
    y<-data.frame(v1=1,y)
    nx<-ncol(x)
    ny<-ncol(y)
    namex <- names(x)
    namey <- names(y)
    if (nx==2) namex <- c("c1",cx)
    if (ny==2) namey <- c("v1",cy)
    namexy <- c(namex,namey)
    for(i in 1:nx) {
        x[,i]<-as.character(x[,i])
    }
    z<-NULL
    for(i in 1:nx) {
        z<-paste(z,x[,i],sep="&")
    }
    w<-NULL
    for(i in 1:ny) {
        m <-tapply(y[,i],z,stat)
        m<-as.matrix(m)
        w<-cbind(w,m)
    }
    nw<-nrow(w)
    c<-rownames(w)
    v<-rep("",nw*nx)
    dim(v)<-c(nw,nx)
    for(i in 1:nw) {
        for(j in 1:nx) {
            v[i,j]<-strsplit(c[i],"&")[[1]][j+1]
        }
    }
    rownames(w)<-NULL
    junto<-data.frame(v[,-1],w)
    junto<-junto[,-nx]
    names(junto)<-namexy[c(-1,-(nx+1))]
    return(junto)
}

ot.helmert <- function(k){

    if(missing(k)) stop("The number of time points is missing.")

    if (is.numeric(k) && length(k) == 1)
          if(k > trunc(k)) stop("The number of time points is not an integer.")

    
    levels <- 1:k

    T0 <- matrix(rep(1/sqrt(k), k), byrow=TRUE, ncol=k)

    T1 <- matrix(rep(0,(k-1)*k), ncol=k, byrow=TRUE)
    T1 <- array(1/sqrt(diag(outer(row(T1)[,1]+1, row(T1)[,1], "*"))),c(k-1,k))
    T1[col(T1) > row(T1)] <- 0
    T1[col(T1) == row(T1)+1] <- -(row(T1)[,1])/sqrt(diag(outer(row(T1)[,1]+1, row(T1)[,1], "*")))

    OT <- rbind(T0, T1)

    OT
}

###################################################
## Utilities for create pathway maps for MetPA
#################################################


# a function to deal with long string names
Wrap.Names<-function(cName, wrap.len=10, tol.len=5){

    nc <- nchar(cName);
    long.inx <- nc > (wrap.len+tol.len);
    long.nms <- cName[long.inx];

    # first get positions of the natural breaks space or hyphen
    pos.list <- gregexpr("[ -]", long.nms);

    for(i in 1:length(pos.list)){
        current.nm <- long.nms[i];
        pos <- pos.list[[i]]+1;
        start.pos<- c(0, pos);
        end.pos <- c(pos, nchar(current.nm)+1);
        splits <- sapply(1:(length(pos)+1), function(x) substring(current.nm, start.pos[x], end.pos[x]-1));
        long.nms[i]<-CheckMergeSplittedNames(splits);
    }

    cName[long.inx] <- long.nms;
    return (cName);
}

# given a vector with naturally splitted string elements
# check if a particular element is too long and need to be
# break by brute force
CheckMergeSplittedNames<-function(nms, wrap.len=10, tol.len=5){
     clean.nm <- "";
     current.nm <- "";
     for(i in 1:length(nms)){
        current.nm <- paste(current.nm, nms[i], sep="");
        current.len <- nchar(current.nm);

        # if too long, break into halves
        if(current.len > wrap.len + tol.len){
             break.pt <- round(current.len/2);
             current.nm <- paste(substr(current.nm, 0, break.pt), "-", "\n",
                                  substr(current.nm, break.pt+1, current.len),  sep="");
             clean.nm <- paste(clean.nm, "\n", current.nm, sep="");
             current.nm <- "";
        }else if(current.len > tol.len){
             clean.nm <- paste(clean.nm, "\n", current.nm, sep="");
             current.nm <- "";
        }else{
            if(i == length(nms)){
                clean.nm <- paste(clean.nm, current.nm, sep=ifelse(nchar(current.nm)<tol.len, "", "\n"));
            }
        }
    }
    return(clean.nm);
}

# break up a long name with brute force
BreakLongNames<-function(long.nm, wrap.len=10, tol.len=5){
      splits <- sapply(seq(1,nchar(long.nm),by=wrap.len), function(x) substr(long.nm, x, x+wrap.len-1));
      tot.len <- length(splits);
      pre.nms<- paste(splits[-tot.len], collapse="\n");
      last.nm <- splits[tot.len];
      if(nchar(last.nm) < tol.len){
            pre.nms<- paste(pre.nms, last.nm, sep="");
            last.nm <- "";
      }
      return(c(pre.nms,last.nm));
}

# get all the KEGG compounds from the pathway databases
getCmpdID<-function(dirName){
	require(KEGGgraph);
	folds<-dir(dirName);
	all.nms <- "";

	for(m in 1:length(folds)){
		files <- dir(paste(dirName, "/", folds[m], sep=""));
		cmpd.nms <- "";
		for(i in 1:length(files)){
			f <- paste(dirName, "/", folds[m],"/",files[i], sep="");
			print(f);
			g <- KEGGpathway2reactionGraph(parseKGML(f));
		   	nms <- nodes(g);
		   	start.pos <- unlist(gregexpr(":", nms))+1;
    	           	nms <- substr(nms, start.pos, nchar(nms));
    			cmpd.nms <- c(cmpd.nms, nms);
		}
		all.nms <- c(all.nms, unique(cmpd.nms));
	}
	write.csv(unique(all.nms), file="kegg_uniq.csv", row.names=F)
}

getPathName<-function(dirName, saveName){
	require(KEGGgraph);
	files<-dir(dirName);
    	nm.mat<-matrix("NA", nrow=length(files), ncol=2);
	for(i in 1:length(files)){
		f <- files[i];
		print(f);
    		path <- parseKGML(paste(dirName,"/",f, sep=""));
		nm.mat[i,]<-c(f, path@pathwayInfo@title);
	}
	write.csv(nm.mat, file=saveName);
}

# extend the axis range to both end
# vec is the values for that axis
# unit is the width to extend, 10 will increase by 1/10 of the range
GetExtendRange<-function(vec, unit=10){
    var.max <- max(vec, na.rm=T);
    var.min <- min(vec, na.rm=T);
    exts <- (var.max - var.min)/unit;
    c(var.min-exts, var.max+exts);
}

# to return a shorter names
# break long names at space, append "..." to indicate
# the abbrev
GetShortNames<-function(nm.vec, max.len= 45){
    new.nms <- vector(mode="character", length=length(nm.vec));
    for(i in 1:length(nm.vec)){
        nm <- nm.vec[i];
        if(nchar(nm) <= max.len){
            new.nms[i] <- nm;
        }else{
            wrds <- strsplit(nm, "[[:space:]]+")[[1]];
            new.nm <- "";
            if(length(wrds)>1){
                   for(m in 1:length(wrds)){
                        wrd <- wrds[m];
                        if(nchar(new.nm)+4+nchar(wrd) <= max.len){
                            new.nm <- paste(new.nm, wrd);
                        }else{
                            new.nms[i] <- paste (new.nm, "...", sep="");
                            break;
                        }
                   }
            }else{
                new.nms[i] <- paste (substr(nm, 0, 21), "...", sep="");
            }
        }
    }
    return (new.nms);
}

# count the number of digits in the values
getndp <- function(x, tol=2*.Machine$double.eps){
  ndp <- 0
  while(!isTRUE(all.equal(x, round(x, ndp), tol=tol))) ndp <- ndp+1
  if(ndp > -log10(tol)) {
        warning("Tolerance reached, ndp possibly underestimated.")
  }
  ndp
}

### convert usr coords (as used in current plot) to pixels in a png
## adapted from the imagemap package
usr2png <- function(xy,im){
  xy <- usr2dev(xy,dev.cur())
  cbind(
        ceiling(xy[,1]*im$Width),
        ceiling((1-xy[,2])*im$Height)
        )
}

usr2plt <- function(xy,dev=dev.cur()){
  olddev <- dev.cur()
  dev.set(dev)
  usr <- par("usr")
  dev.set(olddev)
  xytrans(xy,usr)
}

plt2fig <- function(xy,dev=dev.cur()){
  olddev <- dev.cur()
  dev.set(dev)
  plt <- par("plt")
  dev.set(olddev)
  xytrans2(xy,plt)
}

fig2dev <- function(xy,dev=dev.cur()){
  olddev <- dev.cur()
  dev.set(dev)
  fig <- par("fig")
  dev.set(olddev)
  xytrans2(xy,fig)
}

usr2dev <- function(xy,dev=dev.cur()){
  fig2dev(plt2fig(usr2plt(xy,dev),dev),dev)
}

xytrans2 <- function(xy,par){
  cbind(par[1]+((par[2]-par[1])*xy[,1]),
        par[3]+((par[4]-par[3])*xy[,2]))
}

xytrans <- function(xy,par){
  cbind((xy[,1]-par[1])/(par[2]-par[1]),
        (xy[,2]-par[3])/(par[4]-par[3]))
}

#  VENN DIAGRAM COUNTS AND PLOTS
getVennCounts <- function(x,include="both") {
	x <- as.matrix(x)
	include <- match.arg(include,c("both","up","down"))
	x <- sign(switch(include,
		both = abs(x),
		up = x > 0,
		down = x < 0
	))
	nprobes <- nrow(x)
	ncontrasts <- ncol(x)
	names <- colnames(x)
	if(is.null(names)) names <- paste("Group",1:ncontrasts)
	noutcomes <- 2^ncontrasts
	outcomes <- matrix(0,noutcomes,ncontrasts)
	colnames(outcomes) <- names
	for (j in 1:ncontrasts)
		outcomes[,j] <- rep(0:1,times=2^(j-1),each=2^(ncontrasts-j))
	xlist <- list()
	for (i in 1:ncontrasts) xlist[[i]] <- factor(x[,ncontrasts-i+1],levels=c(0,1))
	counts <- as.vector(table(xlist))
	structure(cbind(outcomes,Counts=counts),class="VennCounts")
}

# function to calculate tick mark based on Heckbert algorithm
# available in the "labeling" package implemented by Justin Talbot
#' Heckbert's labeling algorithm
#' Heckbert, P. S. (1990) Nice numbers for graph labels, Graphics Gems I, Academic Press Professional, Inc.
#' @author Justin Talbot \email{jtalbot@@stanford.edu}
heckbert <- function(dmin, dmax, m)
{
    range <- .heckbert.nicenum((dmax-dmin), FALSE)
    lstep <- .heckbert.nicenum(range/(m-1), TRUE)
    lmin <- floor(dmin/lstep)*lstep
    lmax <- ceiling(dmax/lstep)*lstep
    seq(lmin, lmax, by=lstep)
}

.heckbert.nicenum <- function(x, round)
{
	e <- floor(log10(x))
	f <- x / (10^e)
	if(round)
	{
		if(f < 1.5) nf <- 1
		else if(f < 3) nf <- 2
		else if(f < 7) nf <- 5
		else nf <- 10
	}
	else
	{
		if(f <= 1) nf <- 1
		else if(f <= 2) nf <- 2
		else if(f <= 5) nf <- 5
		else nf <- 10
	}
	nf * (10^e)
}

# borrowed from Hmisc
all.numeric <- function (x, what = c("test", "vector"), extras = c(".", "NA")){
    what <- match.arg(what)
    old <- options(warn = -1)
    on.exit(options(old));
    x <- sub("[[:space:]]+$", "", x);
    x <- sub("^[[:space:]]+", "", x);
    inx <- x %in% c("", extras);
    xs <- x[!inx];
    isnum <- !any(is.na(as.numeric(xs)))
    if (what == "test") 
        isnum
    else if (isnum) 
        as.numeric(x)
    else x
}

GetFileContentAsString <- function(file.nm){
    content <- paste(readLines(file.nm), collapse="\n");
    return(content);
}

ClearNumerics <-function(dat.mat){
     dat.mat[is.na(dat.mat)] <- -777;
     dat.mat[dat.mat == Inf] <- -999;
     dat.mat[dat.mat == -Inf] <- -111;
     dat.mat;
}

###########
# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) {
                           capture.output(format(utils::object.size(x), units = "auto")) })
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}

# shorthand
ShowMemoryUse <- function(..., n=20) {
    print(.ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n));
    print(warnings());
}

CleanMemory <- function(){
    for (i in 1:10){ 
        gc(reset = T);
    }
}

# need to obtain the full path to convert (from imagemagik) for cropping images
GetConvertFullPath<-function(){
    path <- system("which convert", intern=TRUE);
    if((length(path) == 0) && (typeof(path) == "character")){
        print("Could not find convert in the PATH!");
        return("NA");
    }
    return(path);
}

# need to obtain the full path to convert (from imagemagik) for cropping images
GetBashFullPath<-function(){
    path <- system("which bash", intern=TRUE);
    if((length(path) == 0) && (typeof(path) == "character")){
        print("Could not find bash in the PATH!");
        return("NA");
    }
    return(path);
}

########################################
#### Scatterplot3D
#### adapted for better visualization
#######################################

Plot3D <- function(x, y = NULL, z = NULL, color = par("col"), pch = NULL,
     main = NULL, sub = NULL, xlim = NULL, ylim = NULL, zlim = NULL,
     xlab = NULL, ylab = NULL, zlab = NULL, scale.y = 1, angle = 40,
     axis = TRUE, tick.marks = TRUE, label.tick.marks = TRUE,
     x.ticklabs = NULL, y.ticklabs = NULL, z.ticklabs = NULL,
     y.margin.add = 0, grid = TRUE, box = TRUE, lab = par("lab"),
     lab.z = mean(lab[1:2]), type = "p", highlight.3d = FALSE,
     mar = c(5, 3, 4, 3) + 0.1, col.axis = par("col.axis"),
     col.grid = "grey", col.lab = par("col.lab"), cex.symbols = par("cex"),
     cex.axis = 0.8 * par("cex.axis"), cex.lab = par("cex.lab"),
     font.axis = par("font.axis"), font.lab = par("font.lab"),
     lty.axis = par("lty"), lty.grid = 2, lty.hide = 1,
     lty.hplot = par("lty"), log = "", ...)
     # log not yet implemented
{

    mem.par <- par(mar = mar)
    x.scal <- y.scal <- z.scal <- 1
    xlabel <- if (!missing(x)) deparse(substitute(x))
    ylabel <- if (!missing(y)) deparse(substitute(y))
    zlabel <- if (!missing(z)) deparse(substitute(z))
    ## verification, init, ...
    if(highlight.3d && !missing(color))
        warning("color is ignored when highlight.3d = TRUE")


    ## color as part of `x' (data.frame or list):
    if(!is.null(d <- dim(x)) && (length(d) == 2) && (d[2] >= 4))
        color <- x[,4]
    else if(is.list(x) && !is.null(x$color))
        color <- x$color

    ## convert 'anything' -> vector
    xyz <- xyz.coords(x=x, y=y, z=z, xlab=xlabel, ylab=ylabel, zlab=zlabel,
                      log=log)
    if(is.null(xlab)) { xlab <- xyz$xlab; if(is.null(xlab)) xlab <- "" }
    if(is.null(ylab)) { ylab <- xyz$ylab; if(is.null(ylab)) ylab <- "" }
    if(is.null(zlab)) { zlab <- xyz$zlab; if(is.null(zlab)) zlab <- "" }

    if(length(color) == 1)
        color <- rep(color, length(xyz$x))
    else if(length(color) != length(xyz$x))
        stop("length(color) ", "must be equal length(x) or 1")

    angle <- (angle %% 360) / 90
    yz.f <- scale.y * abs(if(angle < 1) angle else if(angle > 3) angle - 4 else 2 - angle)
    yx.f <- scale.y * (if(angle < 2) 1 - angle else angle - 3)
    if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
        temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
        temp <- xlab;  xlab <- ylab;   ylab <- temp
        temp <- xlim;  xlim <- ylim;   ylim <- temp
    }
    angle.1 <- (1 < angle && angle < 2) || angle > 3
    angle.2 <- 1 <= angle && angle <= 3
    dat <- cbind(as.data.frame(xyz[c("x","y","z")]), col = color)

    ## xlim, ylim, zlim -- select the points inside the limits
    if(!is.null(xlim)) {
        xlim <- range(xlim)
        dat <- dat[ xlim[1] <= dat$x & dat$x <= xlim[2] , , drop = FALSE]
    }
    if(!is.null(ylim)) {
        ylim <- range(ylim)
        dat <- dat[ ylim[1] <= dat$y & dat$y <= ylim[2] , , drop = FALSE]
    }
    if(!is.null(zlim)) {
        zlim <- range(zlim)
        dat <- dat[ zlim[1] <= dat$z & dat$z <= zlim[2] , , drop = FALSE]
    }
    n <- nrow(dat)
    if(n < 1) stop("no data left within (x|y|z)lim")

    y.range <- range(dat$y[is.finite(dat$y)])

### 3D-highlighting / colors / sort by y
    if(type == "p" || type == "h") {
        y.ord <- rev(order(dat$y))
        dat <- dat[y.ord, ]
        if(length(pch) > 1)
            if(length(pch) != length(y.ord))
                stop("length(pch) ", "must be equal length(x) or 1")
            else pch <- pch[y.ord]
        daty <- dat$y
        daty[!is.finite(daty)] <- mean(daty[is.finite(daty)])
        if(highlight.3d && !(all(diff(daty) == 0)))
            dat$col <- rgb(seq(0, 1, length = n) * (y.range[2] - daty) / diff(y.range), g=0, b=0)
    }

### optim. axis scaling
    p.lab <- par("lab")
    ## Y
    y.range <- range(dat$y[is.finite(dat$y)], ylim)
    y.prty <- pretty(y.range, n = lab[2],
        min.n = max(1, min(.5 * lab[2], p.lab[2])))
    y.scal <- round(diff(y.prty[1:2]), digits = 12)
    y.add <- min(y.prty)
    dat$y <- (dat$y - y.add) / y.scal
    y.max <- (max(y.prty) - y.add) / y.scal
    if(!is.null(ylim)) y.max <- max(y.max, ceiling((ylim[2] - y.add) / y.scal))
#    if(angle > 2) dat$y <- y.max - dat$y  ## turn y-values around
    ## X
    x.range <- range(dat$x[is.finite(dat$x)], xlim)
    x.prty <- pretty(x.range, n = lab[1],
        min.n = max(1, min(.5 * lab[1], p.lab[1])))
    x.scal <- round(diff(x.prty[1:2]), digits = 12)
    dat$x <- dat$x / x.scal
    x.range <- range(x.prty) / x.scal
    x.max <- ceiling(x.range[2])
    x.min <-   floor(x.range[1])
    if(!is.null(xlim)) {
        x.max <- max(x.max, ceiling(xlim[2] / x.scal))
        x.min <- min(x.min,   floor(xlim[1] / x.scal))
    }
    x.range <- range(x.min, x.max)
    ## Z
    z.range <- range(dat$z[is.finite(dat$z)], zlim)
    z.prty <- pretty(z.range, n = lab.z,
        min.n = max(1, min(.5 * lab.z, p.lab[2])))
    z.scal <- round(diff(z.prty[1:2]), digits = 12)
    dat$z <- dat$z / z.scal
    z.range <- range(z.prty) / z.scal
    z.max <- ceiling(z.range[2])
    z.min <-   floor(z.range[1])
    if(!is.null(zlim)) {
        z.max <- max(z.max, ceiling(zlim[2] / z.scal))
        z.min <- min(z.min,   floor(zlim[1] / z.scal))
    }
    z.range <- range(z.min, z.max)

### init graphics
    plot.new()
    if(angle.2) {x1 <- x.min + yx.f * y.max; x2 <- x.max}
    else        {x1 <- x.min; x2 <- x.max + yx.f * y.max}
    plot.window(c(x1, x2), c(z.min, z.max + yz.f * y.max))
    temp <- strwidth(format(rev(y.prty))[1], cex = cex.axis/par("cex"))
    if(angle.2) x1 <- x1 - temp - y.margin.add
    else        x2 <- x2 + temp + y.margin.add
    plot.window(c(x1, x2), c(z.min, z.max + yz.f * y.max))
    if(angle > 2) par("usr" = par("usr")[c(2, 1, 3:4)])
    usr <- par("usr") # we have to remind it for use in closures
    title(main, sub, ...)

### draw axis, tick marks, labels, grid, ...
    xx <- if(angle.2) c(x.min, x.max) else c(x.max, x.min)
    if(grid) {
	## grids
	###################
	# XY wall
        i <- x.min:x.max;
        segments(i, z.min, i + (yx.f * y.max), yz.f * y.max + z.min,
                 col = col.grid, lty = lty.grid);

        i <- 0:y.max;
        segments(x.min + (i * yx.f), i * yz.f + z.min,
                 x.max + (i * yx.f), i * yz.f + z.min,
                 col = col.grid, lty = lty.grid);

	######################
	# XZ wall
	# verticle lines
        temp <- yx.f * y.max;
        temp1 <- yz.f * y.max;
        i <- (x.min + temp):(x.max + temp);
        segments(i, z.min + temp1, i, z.max + temp1,
                 col = col.grid, lty = lty.grid);

	# horizontal lines
        i <- (z.min + temp1):(z.max + temp1);
        segments(x.min + temp, i, x.max + temp, i,
                 col = col.grid, lty = lty.grid)


	##################
	# YZ wall
	# horizontal lines
        i <- xx[2]:x.min;
	mm <- z.min:z.max;
        segments(i, mm, i + temp, mm + temp1,
                 col = col.grid, lty = lty.grid);
	# verticle lines
        i <- 0:y.max;
        segments(x.min + (i * yx.f), i * yz.f + z.min,
                 xx[2] + (i * yx.f), i * yz.f + z.max,
                 col = col.grid, lty = lty.grid)


	# make the axis into solid line
        segments(x.min, z.min, x.min + (yx.f * y.max), yz.f * y.max + z.min,
                 col = col.grid, lty = lty.hide);
        segments(x.max, z.min, x.max + (yx.f * y.max), yz.f * y.max + z.min,
                 col = col.axis, lty = lty.hide);
        segments(x.min + (y.max * yx.f), y.max * yz.f + z.min,
                 x.max + (y.max* yx.f), y.max * yz.f + z.min,
                col = col.grid, lty = lty.hide);
        segments(x.min + temp, z.min + temp1, x.min + temp, z.max + temp1,
                col = col.grid, lty = lty.hide);
        segments(x.max + temp, z.min + temp1, x.max + temp, z.max + temp1,
                col = col.axis, lty = lty.hide);
        segments(x.min + temp, z.max + temp1, x.max + temp, z.max + temp1,
                col = col.axis, lty = lty.hide);
        segments(xx[2], z.max, xx[2] + temp, z.max + temp1,
                col = col.axis, lty = lty.hide);
    }
    if(axis) {
        if(tick.marks) { ## tick marks
            xtl <- (z.max - z.min) * (tcl <- -par("tcl")) / 50
            ztl <- (x.max - x.min) * tcl / 50
            mysegs <- function(x0,y0, x1,y1)
                segments(x0,y0, x1,y1, col=col.axis, lty=lty.axis)
            ## Y
            i.y <- 0:y.max
            mysegs(yx.f * i.y - ztl + xx[1], yz.f * i.y + z.min,
                   yx.f * i.y + ztl + xx[1], yz.f * i.y + z.min)
            ## X
            i.x <- x.min:x.max
            mysegs(i.x, -xtl + z.min, i.x, xtl + z.min)
            ## Z
            i.z <- z.min:z.max
            mysegs(-ztl + xx[2], i.z, ztl + xx[2], i.z)

            if(label.tick.marks) { ## label tick marks
                las <- par("las")
                mytext <- function(labels, side, at, ...)
                    mtext(text = labels, side = side, at = at, line = -.5,
                          col=col.lab, cex=cex.axis, font=font.lab, ...)
                ## X
                if(is.null(x.ticklabs))
                    x.ticklabs <- format(i.x * x.scal)
                mytext(x.ticklabs, side = 1, at = i.x)
                ## Z
                if(is.null(z.ticklabs))
                    z.ticklabs <- format(i.z * z.scal)
                mytext(z.ticklabs, side = if(angle.1) 4 else 2, at = i.z,
                       adj = if(0 < las && las < 3) 1 else NA)
                ## Y
                temp <- if(angle > 2) rev(i.y) else i.y ## turn y-labels around
                if(is.null(y.ticklabs))
                    y.ticklabs <- format(y.prty)
                else if (angle > 2)
                    y.ticklabs <- rev(y.ticklabs)
                text(i.y * yx.f + xx[1],
                     i.y * yz.f + z.min, y.ticklabs,
                     pos=if(angle.1) 2 else 4, offset=1,
                     col=col.lab, cex=cex.axis/par("cex"), font=font.lab)
            }
        }

        ## axis and labels

        mytext2 <- function(lab, side, line, at)
            mtext(lab, side = side, line = line, at = at, col = col.lab,
                  cex = cex.lab, font = font.axis, las = 0)
        ## X
        lines(c(x.min, x.max), c(z.min, z.min), col = col.axis, lty = lty.axis)
        mytext2(xlab, 1, line = 1.5, at = mean(x.range))
        ## Y
        lines(xx[1] + c(0, y.max * yx.f), c(z.min, y.max * yz.f + z.min),
              col = col.axis, lty = lty.axis)
        mytext2(ylab, if(angle.1) 2 else 4, line= 0.5, at = z.min + y.max * yz.f)

        ## Z
        lines(xx[c(2,2)], c(z.min, z.max), col = col.axis, lty = lty.axis)
        mytext2(zlab, if(angle.1) 4 else 2, line= 1.5, at = mean(z.range))

    }

### plot points
    x <- dat$x + (dat$y * yx.f)
    z <- dat$z + (dat$y * yz.f)
    col <- as.character(dat$col)
    if(type == "h") {
        z2 <- dat$y * yz.f + z.min
        segments(x, z, x, z2, col = col, cex = cex.symbols, lty = lty.hplot, ...)
        points(x, z, type = "p", col = col, pch = pch, cex = cex.symbols, ...)
    }
    else points(x, z, type = type, col = col, pch = pch, cex = cex.symbols, ...)

### box-lines in front of points (overlay)
    if(axis && box) {
        lines(c(x.min, x.max), c(z.max, z.max),
              col = col.axis, lty = lty.axis)
        lines(c(0, y.max * yx.f) + x.max, c(0, y.max * yz.f) + z.max,
              col = col.axis, lty = lty.axis)
        lines(xx[c(1,1)], c(z.min, z.max), col = col.axis, lty = lty.axis)
    }


    # par(mem.par) # we MUST NOT set the margins back
    ### Return Function Object
    ob <- ls() ## remove all unused objects from the result's enviroment:
    rm(list = ob[!ob %in% c("angle", "mar", "usr", "x.scal", "y.scal", "z.scal", "yx.f",
        "yz.f", "y.add", "z.min", "z.max", "x.min", "x.max", "y.max",
        "x.prty", "y.prty", "z.prty")])
    rm(ob)
    invisible(list(
        xyz.convert = function(x, y=NULL, z=NULL) {
            xyz <- xyz.coords(x, y, z)
            if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
                temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
            }
            y <- (xyz$y - y.add) / y.scal
            return(list(x = xyz$x / x.scal + yx.f * y,
                y = xyz$z / z.scal + yz.f * y))
        },
        points3d = function(x, y = NULL, z = NULL, type = "p", ...) {
            xyz <- xyz.coords(x, y, z)
            if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
                temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
            }
            y2 <- (xyz$y - y.add) / y.scal
            x <- xyz$x / x.scal + yx.f * y2
            y <- xyz$z / z.scal + yz.f * y2
            mem.par <- par(mar = mar, usr = usr)
            on.exit(par(mem.par))
            if(type == "h") {
                y2 <- z.min + yz.f * y2
                segments(x, y, x, y2, ...)
                points(x, y, type = "p", ...)
            }
            else points(x, y, type = type, ...)
        },
        plane3d = function(Intercept, x.coef = NULL, y.coef = NULL,
            lty = "dashed", lty.box = NULL, ...){
            if(!is.atomic(Intercept) && !is.null(coef(Intercept))) Intercept <- coef(Intercept)
            if(is.null(lty.box)) lty.box <- lty
            if(is.null(x.coef) && length(Intercept) == 3){
                x.coef <- Intercept[if(angle > 2) 3 else 2]
                y.coef <- Intercept[if(angle > 2) 2 else 3]
                Intercept <- Intercept[1]
            }
            mem.par <- par(mar = mar, usr = usr)
            on.exit(par(mem.par))
            x <- x.min:x.max
            ltya <- c(lty.box, rep(lty, length(x)-2), lty.box)
            x.coef <- x.coef * x.scal
            z1 <- (Intercept + x * x.coef + y.add * y.coef) / z.scal
            z2 <- (Intercept + x * x.coef +
                (y.max * y.scal + y.add) * y.coef) / z.scal
            segments(x, z1, x + y.max * yx.f, z2 + yz.f * y.max, lty = ltya, ...)
            y <- 0:y.max
            ltya <- c(lty.box, rep(lty, length(y)-2), lty.box)
            y.coef <- (y * y.scal + y.add) * y.coef
            z1 <- (Intercept + x.min * x.coef + y.coef) / z.scal
            z2 <- (Intercept + x.max * x.coef + y.coef) / z.scal
            segments(x.min + y * yx.f, z1 + y * yz.f,
                x.max + y * yx.f, z2 + y * yz.f, lty = ltya, ...)
        },

        wall3d = function(Intercept, x.coef = NULL, y.coef = NULL,
            lty = "dashed", lty.box = NULL, ...){
            if(!is.atomic(Intercept) && !is.null(coef(Intercept))) Intercept <- coef(Intercept)
            if(is.null(lty.box)) lty.box <- lty
            if(is.null(x.coef) && length(Intercept) == 3){
                x.coef <- Intercept[if(angle > 2) 3 else 2]
                y.coef <- Intercept[if(angle > 2) 2 else 3]
                Intercept <- Intercept[1]
            }
            mem.par <- par(mar = mar, usr = usr)
            on.exit(par(mem.par))
            x <- x.min:x.max
            ltya <- c(lty.box, rep(lty, length(x)-2), lty.box)
            x.coef <- x.coef * x.scal
            z1 <- (Intercept + x * x.coef + y.add * y.coef) / z.scal
            z2 <- (Intercept + x * x.coef +
                (y.max * y.scal + y.add) * y.coef) / z.scal
            segments(x, z1, x + y.max * yx.f, z2 + yz.f * y.max, lty = ltya, ...)
            y <- 0:y.max
            ltya <- c(lty.box, rep(lty, length(y)-2), lty.box)
            y.coef <- (y * y.scal + y.add) * y.coef
            z1 <- (Intercept + x.min * x.coef + y.coef) / z.scal
            z2 <- (Intercept + x.max * x.coef + y.coef) / z.scal
            segments(x.min + y * yx.f, z1 + y * yz.f,
                x.max + y * yx.f, z2 + y * yz.f, lty = ltya, ...)
        },
        box3d = function(...){
            mem.par <- par(mar = mar, usr = usr)
            on.exit(par(mem.par))
            lines(c(x.min, x.max), c(z.max, z.max), ...)
            lines(c(0, y.max * yx.f) + x.max, c(0, y.max * yz.f) + z.max, ...)
            lines(c(0, y.max * yx.f) + x.min, c(0, y.max * yz.f) + z.max, ...)
            lines(c(x.max, x.max), c(z.min, z.max), ...)
            lines(c(x.min, x.min), c(z.min, z.max), ...)
            lines(c(x.min, x.max), c(z.min, z.min), ...)
        }
    ))
}

GetPdfLatexFullPath<-function(){
    path <- system("which pdflatex", intern=TRUE);
    if((length(path) == 0) && (typeof(path) == "character")){
        print("Could not find pdflatex in the PATH!");
        return("NA");
    }
    return(path);
}

GetRandomNumbers <- function(){
    rm(.Random.seed);
    runif(1);
    return(.Random.seed[3:626]);
}

GetRScriptFullPath<-function(){
    path <- system("which RScript", intern=TRUE);
    if((length(path) == 0) && (typeof(path) == "character")){
        print("Could not find pdflatex in the PATH!");
        return("NA");
    }
    return(path);
}

# note, try to use the fread, however, it has issues with 
# some windows 10 files "Line ending is \r\r\n. .... appears to add the extra \r in text mode on Windows"
# in such as, use the slower read.table method
.readDataTable <- function(fileName){
    dat <- try(data.table::fread(fileName, header=TRUE, check.names=FALSE, data.table=FALSE));
    if(class(dat) == "try-error"){
        # try to use "tr" to remove double return characters
        trFileName <- paste("tr -d \'\\r\' <", fileName);
        dat <- try(data.table::fread(trFileName, header=TRUE, check.names=FALSE, data.table=FALSE));
        if(class(dat) == "try-error"){
            print("Using slower file reader ...");
            formatStr <- substr(fileName, nchar(fileName)-2, nchar(fileName))
            if(formatStr == "txt"){
                dat <-try(read.table(fileName,header=TRUE,comment.char = "", check.names=F, as.is=T));
            }else{ # note, read.csv is more than read.table with sep=","
                dat <-try(read.csv(fileName,header=TRUE,comment.char = "", check.names=F, as.is=T));
            }  
        }
    }
    return(dat);
}