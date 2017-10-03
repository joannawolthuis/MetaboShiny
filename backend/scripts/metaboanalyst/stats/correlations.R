#########################################################
## R script for MetaboAnalyst
## Description: perform correlation analysis
##
## Author: Jeff Xia, jeff.xia@mcgill.ca
## McGill University, Canada
##
## License: GNU GPL (>= 2)
###################################################

#######################################
## Pattern hunter
##########################################

# Run template on all the high region effect genes
template.match <- function(x, template, dist.name) {
  k<-cor.test(x,template, method=dist.name);
  c(k$estimate, k$stat, k$p.value)
}

Match.Pattern<-function(dist.name="pearson", pattern=NULL){
    if(is.null(pattern)){
        pattern <- paste(1:length(levels(dataSet$cls)), collapse="-");
    }
    templ <- as.numeric(ClearStrings(strsplit(pattern, "-")[[1]]));

    if(all(templ==templ[1])){
        AddErrMsg("Cannot calculate correlation on constant values!");
        return(0);
    }

    new.template <- vector(mode="numeric", length=length(dataSet$cls))
    # expand to match each levels in the dataSet$cls
    all.lvls <- levels(dataSet$cls);

    if(length(templ)!=length(all.lvls)){
        AddErrMsg("Wrong template - must the same length as the group number!");
        return(0);
    }

    for(i in 1:length(templ)){
        hit.inx <- dataSet$cls == all.lvls[i]
        new.template[hit.inx] = templ[i];
    }

    cbtempl.results <- apply(dataSet$norm, 2, template.match, new.template, dist.name);

    cor.res<-t(cbtempl.results);

    fdr.col <- p.adjust(cor.res[,3], "fdr");
    cor.res <- cbind(cor.res, fdr.col);
    colnames(cor.res)<-c("correlation", "t-stat", "p-value", "FDR");
    ord.inx<-order(cor.res[,3]);

    sig.mat <- signif(cor.res[ord.inx,],5);

    fileName <- "correlation_pattern.csv";
    write.csv(sig.mat,file=fileName);

    analSet$corr$sig.nm<-fileName;
    analSet$corr$cor.mat<-sig.mat;
    analSet$corr$pattern <- pattern;
    analSet <<- analSet;
    return(1);
}

GenerateTemplates <- function(){

    level.len <- length(levels(dataSet$cls));

    # only specify 4: increasing, decreasing, mid high, mid low, constant
    incs <- 1:level.len;
    desc <- level.len:1;

    if(level.len > 2){
        # use ceiling, so that the peak will be right for even length
        mid.pos <- ceiling((level.len+1)/2);
        mid.high <- c(1:mid.pos, seq(mid.pos-1,by=-1,length.out=level.len-mid.pos));
        mid.low <- c(mid.pos:1, seq(2, length.out=level.len-mid.pos));

        res <- rbind(incs, desc, mid.high, mid.low); # add the constant one
    }else{
        res <- rbind(incs, desc);
    }
    # turn into string
    res <- apply(res, 1, paste, collapse="-");

    # add the ledgends
    res <- c(paste(levels(dataSet$cls), collapse="-"), res);
    return (res);
}

# calculate correlation of all other feature to a given feature name
FeatureCorrelation<-function(dist.name, varName){

    cbtempl.results <- apply(dataSet$norm, 2, template.match, dataSet$norm[,varName], dist.name);
    cor.res<-t(cbtempl.results);

    fdr.col <- p.adjust(cor.res[,3], "fdr");
    cor.res <- cbind(cor.res, fdr.col);
    colnames(cor.res)<-c("correlation", "t-stat", "p-value", "FDR");
    ord.inx<-order(cor.res[,3])
    sig.mat <-signif(cor.res[ord.inx,],5);

    fileName <- "correlation_feature.csv";
    write.csv(sig.mat,file=fileName);

    analSet$corr$sig.nm<-fileName;
    analSet$corr$cor.mat<-sig.mat;
    analSet$corr$pattern<-varName;
    analSet <<- analSet;
    return(1);
}

PlotCorr <- function(imgName, format="png", dpi=72, width=NA){

    cor.res <- analSet$corr$cor.mat;
    pattern <- analSet$corr$pattern;
    title <- paste(GetVariableLabel(), "correlated with the", pattern);
    if(nrow(cor.res) > 25){

        # first get most signficant ones (p value)
        ord.inx<-order(cor.res[,3]);
        cor.res <- cor.res[ord.inx, ];
        cor.res <- cor.res[1:25, ];

        # then order by their direction (correlation)
        ord.inx<-order(cor.res[,1]);
        if(sum(cor.res[,1] > 0) == 0){ # all negative correlation
            ord.inx <- rev(ord.inx);
        }
        cor.res <- cor.res[ord.inx, ];
        title <- paste("Top 25", tolower(GetVariableLabel()), "correlated with the", pattern);
    }

    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
    if(is.na(width)){
        w <- h <- 7.2;
    }else if(width == 0){
        w <- 7.2;
        imgSet$corr <- imgName;
        imgSet <<- imgSet;
    }else{
        w <- h <- width;
    }
    par(mar=c(5,6,4,3))
    rownames(cor.res)<-substr(rownames(cor.res), 1, 18);
    cols <- ifelse(cor.res[,1] >0, "mistyrose","lightblue");

    dotchart(cor.res[,1], pch="", xlim=c(-1,1), xlab="Correlation coefficients", main=title);
    rownames(cor.res) <- NULL;
    barplot(cor.res[,1], space=c(0.5, rep(0, nrow(cor.res)-1)), xlim=c(-1,1), xaxt="n", col = cols, add=T,horiz=T);}


GetCorrSigFileName <- function(){
    analSet$corr$sig.nm;
}

GetCorSigMat<-function(){
    as.matrix(CleanNumber(analSet$corr$cor.mat));
}

GetCorSigRowNames<-function(){
    rownames(analSet$corr$cor.mat);
}

GetCorSigColNames<-function(){
    colnames(analSet$corr$cor.mat);
}

GetSigTable.Corr<-function(){
    GetSigTable(analSet$corr$cor.mat, "Pattern search using correlation analysis");
}

PlotCorrHeatMap<-function(target="col", cor.method="pearson", 
                colors="bwm", fix.col=F, no.clst=F, top=F, topNum=50){

    main <- xlab <- ylab <- NULL;
    data <- dataSet$norm;
    if(target == 'row'){
        data <- t(data);
    }

    if(ncol(data) > topNum){
        filter.val <- apply(data.matrix(data), 2, IQR, na.rm=T);
        rk <- rank(-filter.val, ties.method='random');
        data <- as.data.frame(data[,rk <= topNum]);
        print("Data is reduced to 1000 vars ..");
    }
    corr.mat<-cor(data, method=cor.method);
    # use total abs(correlation) to select
    if(top){
        cor.sum <- apply(abs(corr.mat), 1, sum);
        cor.rk <- rank(-cor.sum);
        var.sel <- cor.rk <= topNum;
        corr.mat <- corr.mat[var.sel, var.sel];
    }

    # set up parameter for heatmap
    suppressMessages(require(RColorBrewer));
    suppressMessages(require(gplots));
    colors="heat"
    if(colors=="gbr"){
        colors <- colorRampPalette(c("green", "black", "red"), space="rgb")(256);
    }else if(colors == "heat"){
        colors <- heat.colors(256);
    }else if(colors == "topo"){
        colors <- topo.colors(256);
    }else if(colors == "gray"){
        colors <- colorRampPalette(c("grey90", "grey10"))(256);
    }else{
        colors <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256));
    }
    if(no.clst){
        rowv=FALSE;
        colv=FALSE;
        dendro= "none";
    }else{
        rowv=TRUE;
        colv=TRUE;
        dendro= "both";
    }

    require(pheatmap);
    if(fix.col){
        breaks <- seq(from = -1, to = 1, length = 257);
        res <- pheatmap(corr.mat, 
            fontsize=8, fontsize_row=8, 
            cluster_rows = colv, 
            cluster_cols = rowv,
            color = colors,
            breaks = breaks
            );
      }else{
        res <- pheatmap(corr.mat, 
            fontsize=8, fontsize_row=8, 
            cluster_rows = colv, 
            cluster_cols = rowv,
            color = colors
            );
      }

     # need to re-order according to the clusters
     new.ord <- res$tree_row$order;
     corr.mat <- corr.mat[new.ord, new.ord];
     write.csv(signif(corr.mat,5), file=file.path(exp_dir, "correlation_table.csv"))
}