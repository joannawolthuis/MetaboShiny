####################################################################
## R script for MetaboAnalyst
## Description: perform Dendrogram, Heatmap, Kmeans & SOM analysis
##
## Author: Jeff Xia, jeff.xia@mcgill.ca
## McGill University, Canada
##
## License: GNU GPL (>= 2)
###################################################

####################################
########### Dendrogram ##############
#####################################

PlotHCTree<-function(imgName, format="png", dpi=72, width=NA, smplDist, clstDist){
    # set up data set
    hc.dat<-as.matrix(dataSet$norm);
    colnames(hc.dat)<-substr(colnames(hc.dat), 1, 18) # some names are too long
    # set up distance matrix
    if(smplDist == 'euclidean'){
	dist.mat<-dist(hc.dat, method = smplDist);
    }else{
	dist.mat<-dist(1-cor(t(hc.dat), method = smplDist));
    }

    # record the paramters
    analSet$tree<-list(dist.par=smplDist, clust.par=clstDist);
    analSet <<- analSet;
    # build the tree
    hc_tree<-hclust(dist.mat, method=clstDist);

    # plot the tree
    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
    if(is.na(width)){
        w <- minH <- 630;
        myH <- nrow(hc.dat)*10 + 150;
        if(myH < minH){
            myH <- minH;
        }   
        w <- round(w/72,2);
        h <- round(myH/72,2);
    }else if(width == 0){
        w <- h <- 7.2;
        imgSet$tree<-imgName;
        imgSet <<- imgSet;
    }else{
        w <- h <- 7.2;
    }

    Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
    par(cex=0.8, mar=c(4,2,2,8));
    if(dataSet$cls.type == "disc"){
        clusDendro<-as.dendrogram(hc_tree);
        cols <- GetColorSchema();
        names(cols) <- rownames(hc.dat);
        labelColors <- cols[hc_tree$order];
        colLab <- function(n){
        if(is.leaf(n)) {
            a <- attributes(n)
            labCol <- labelColors[a$label];
            attr(n, "nodePar") <- 
                if(is.list(a$nodePar)) c(a$nodePar, lab.col = labCol,pch=NA) else
                               list(lab.col = labCol,pch=NA)
            }
            n
        }
        clusDendro<-dendrapply(clusDendro, colLab)
        plot(clusDendro,horiz=T,axes=T);
        par(cex=1);
        legend.nm <- as.character(dataSet$cls);
        legend("topleft", legend = unique(legend.nm), pch=15, col=unique(cols), bty = "n");
     }else{
        plot(as.dendrogram(hc_tree), hang=-1, main=paste("Cluster with", clstDist, "method"), xlab=NULL, sub=NULL, horiz=TRUE);
     }
     dev.off();
}

# inx has to be 1 or 2
GetClassLabel<-function(inx){
    levels(dataSet$cls)[inx]
}
############################
########### SOM #############
#############################

# SOM analysis
SOM.Anal<-function(x.dim, y.dim, initMethod, neigb = 'gaussian'){
    require(som);
    analSet$som<-som(as.matrix(dataSet$norm), xdim=x.dim, ydim=y.dim, init=initMethod, neigh=neigb);
    analSet <<- analSet;
}

# get members for given cluster index, return a character string
GetSOMClusterMembers<-function(i, j){
	clust<-analSet$som$visual;
	xTrue<-clust$x == i;
	yTrue<-clust$y == j;
    hit.inx <- xTrue & yTrue;

    all.cols <- GetColorSchema();
    paste("<font color=\"", all.cols[hit.inx], "\">", rownames(dataSet$norm)[hit.inx], "</font>",collapse =", ");
}

GetAllSOMClusterMembers<-function(){
	clust<-analSet$som$visual;
	xdim<-analSet$som$xdim;
	ydim<-analSet$som$ydim;

	clust.df = data.frame();
	rowNameVec = c();
	i = 0;
	while(i < xdim){
	    j = 0;
	    while(j < ydim){
		xTrue<-clust$x == i;
		yTrue<-clust$y == j;
		if(i==0 & j==0){ # bug in R, the first one need to be different
			clust.df <- rbind(paste(rownames(dataSet$norm)[xTrue & yTrue], collapse = " "));
			rowNameVec <- c(paste("Cluster(", i, ",", j,")"));
		}else{
			clust.df <- rbind(clust.df, paste(rownames(dataSet$norm)[xTrue & yTrue], collapse=" "));
			rowNameVec <- c(rowNameVec, paste("Cluster(", i, ",", j,")"));
		}
		j = j+1;
	   }
	   i = i+1;
	}
	row.names(clust.df)<- rowNameVec;
	colnames(clust.df)<-"Samples in each cluster";
	print(xtable(clust.df, align="l|p{8cm}", caption="Clustering result using SOM"),caption.placement="top", size="\\scriptsize");
}


# plot SOM map for  less than 20 clusters
PlotSOM <- function(imgName, format="png", dpi=72, width=NA){
	xdim<-analSet$som$xdim;
	ydim<-analSet$som$ydim;
    total<-xdim*ydim;
    if(total>20) { return();}

    ylabel<-GetValueLabel();
    clust<-analSet$som$visual;

    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
    if(is.na(width)){
        w <- 9;
    }else if(width == 0){
        w <- 7;
        imgSet$som<-imgName;
        imgSet <<- imgSet;
    }else{
        w <- width;
    }
    h <- w*8/9;

    Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
    par(mfrow = GetXYCluster(total), mar=c(5,4,2,2));
	for (i in 0:(xdim-1)) {
            xTrue<-clust$x == i;
            for (j in 0:(ydim-1)) {
               yTrue<-clust$y == j;
               sel.inx<-xTrue & yTrue; # selected row
               if(sum(sel.inx)>0){ # some cluster may not contain any member
                    matplot(t(dataSet$norm[sel.inx, ]), type="l", col='grey', axes=F, ylab=ylabel,
                            main=paste("Cluster(", i, ",", j,")", ", n=", sum(sel.inx), sep=""))
                    lines(apply(dataSet$norm[sel.inx, ], 2, median), type="l", col='blue', lwd=1);
                }else{ # plot a dummy 
                    plot(t(dataSet$norm[1, ]), type="n", axes=F, ylab=ylabel,
                            main=paste("Cluster(", i, ",", j,")",", n=", sum(sel.inx),sep=""))
                }
                axis(2);
                axis(1, 1:ncol(dataSet$norm), substr(colnames(dataSet$norm), 1, 7), las=2);
            }
	}
    dev.off();
}

##################################
########### K-means ##############
###################################

# functions for k-means analysis
Kmeans.Anal<-function(clust.num){
    analSet$kmeans<-kmeans (dataSet$norm, clust.num, nstart=100);
    analSet <<- analSet;
}

PlotKmeans<-function(imgName, format="png", dpi=72, width=NA){
	clust.num <- max(analSet$kmeans$cluster);

    if(clust.num>20) return();
	# calculate arrangement of panel
    ylabel<-GetValueLabel();
    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
    if(is.na(width)){
        w <- 9;
    }else if(width == 0){
        w <- 7;
        imgSet$kmeans<-imgName;
        imgSet <<- imgSet;
    }else{
        w <- width;
    }
    h <- w*8/9;

    Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
    par(mfrow = GetXYCluster(clust.num), mar=c(5,4,2,2));
	for (loop in 1:clust.num) {
	    matplot(t(dataSet$norm[analSet$kmeans$cluster==loop,]), type="l", col='grey', ylab=ylabel, axes=F,
            main=paste("Cluster ",loop, ", n=", analSet$kmeans$size[loop], sep=""))
            lines(apply(dataSet$norm[analSet$kmeans$cluster==loop,], 2, median), type="l", col='blue', lwd=1);
        axis(2);
        axis(1, 1:ncol(dataSet$norm), substr(colnames(dataSet$norm), 1, 7), las=2);
	}
    dev.off();
}

# get cluster member for give index
# add HTML color to the names based on its group membership
GetKMClusterMembers<-function(i){
       all.cols <- GetColorSchema();
       hit.inx <- analSet$kmeans$cluster== i;

       paste("<font color=\"", all.cols[hit.inx], "\">", rownames(dataSet$norm)[hit.inx], "</font>",collapse =", ");
       # paste(all.cols[hit.inx], rownames(dataSet$norm)[hit.inx], collapse =", ");
}

GetAllKMClusterMembers<-function(){
	clust.df = data.frame();
	rowNameVec = c();
	i = 1;
	clust.num<-max(analSet$kmeans$cluster);
	while(i<=clust.num){
		if(i==1){
			clust.df <- rbind(paste(rownames(dataSet$norm)[analSet$kmeans$cluster== i], collapse = " "));
		}else{
			clust.df <- rbind(clust.df,paste(rownames(dataSet$norm)[analSet$kmeans$cluster== i], collapse = " "));
		}
		rowNameVec <- c(rowNameVec, paste("Cluster(", i, ")"));
		i = i+1;
	}
	row.names(clust.df)<- rowNameVec;
	colnames(clust.df)<-"Samples in each cluster";
	print(xtable(clust.df, align="l|p{8cm}", caption="Clustering result using K-means"), caption.placement="top", size="\\scriptsize");
}

# plot a sub heatmap based on results from t-tests/ANOVA, VIP or randomforest
PlotSubHeatMap <- function(imgName, format="png", dpi=72, width=NA, dataOpt, scaleOpt, smplDist, clstDist, palette, method.nm, top.num, viewOpt, rowV=T, colV=T, border=T, grp.ave=F){

    var.nms = colnames(dataSet$norm);
    if(top.num < length(var.nms)){
        if(method.nm == 'tanova'){
            if(GetGroupNumber() == 2){
                if(is.null(analSet$tt)){
                    Ttests.Anal();
                }
                var.nms <- names(sort(analSet$tt$p.value))[1:top.num];
            }else{
                if(is.null(analSet$aov)){
                    ANOVA.Anal();
                }
                var.nms <- names(sort(analSet$aov$p.value))[1:top.num];
            }
        }else if(method.nm == 'cor'){
            if(is.null(analSet$cor.res)){
                Match.Pattern();
            }

            # re-order for pretty view
            cor.res <- analSet$cor.res;

            ord.inx<-order(cor.res[,3]);
            cor.res <- cor.res[ord.inx, ];

            ord.inx<-order(cor.res[,1]);
            cor.res <- cor.res[ord.inx, ];

            var.nms <- rownames(cor.res)[1:top.num];
        }else if(method.nm == 'vip'){
            if(is.null(analSet$plsda)){
                PLSR.Anal();
                PLSDA.CV();
            }
            vip.vars <- analSet$plsda$vip.mat[,1];# use the first component
            var.nms <- names(rev(sort(vip.vars)))[1:top.num];
        }else if(method.nm == 'rf'){
            if(is.null(analSet$rf)){
                RF.Anal();
            }
            var.nms <- GetRFSigRowNames()[1:top.num];
        }
    }
    var.inx <- match(var.nms, colnames(dataSet$norm));
    PlotHeatMap(imgName, format, dpi, width, dataOpt, scaleOpt, smplDist, clstDist, palette, viewOpt, rowV, colV, var.inx, border, grp.ave);
}


PlotHeatMap<-function(imgName, format="png", dpi=72, width=NA, dataOpt, scaleOpt, smplDist, clstDist, palette, viewOpt="detail", rowV=T, colV=T, var.inx=NA, border=T, grp.ave=F){

    # record the paramters
    analSet$htmap<-list(dist.par=smplDist, clust.par=clstDist);
    analSet <<- analSet;

    # set up data set
    if(dataOpt=="norm"){
        my.data <- dataSet$norm;
    }else{
        my.data <- dataSet$proc;
    }

    if(is.na(var.inx)){
        hc.dat<-as.matrix(my.data);
    }else{
        hc.dat<-as.matrix(my.data[,var.inx]);
    }

    colnames(hc.dat)<-substr(colnames(hc.dat),1,18) # some names are too long
    hc.cls <- dataSet$cls;
    if(grp.ave){ # only use group average
        lvs <- levels(hc.cls);
        my.mns <- matrix(ncol=ncol(hc.dat),nrow=length(lvs));
        for(i in 1:length(lvs)){
            inx <-hc.cls == lvs[i];
            my.mns[i,]<- apply(hc.dat[inx, ], 2, mean);
        }
        rownames(my.mns) <- lvs;
        colnames(my.mns) <- colnames(hc.dat);
        hc.dat <- my.mns;
        hc.cls <- as.factor(lvs);
    }

    # set up colors for heatmap
    if(palette=="gbr"){
        colors <- colorRampPalette(c("green", "black", "red"), space="rgb")(256);
    }else if(palette == "heat"){
        colors <- heat.colors(256);
    }else if(palette == "topo"){
        colors <- topo.colors(256);
    }else if(palette == "gray"){
        colors <- colorRampPalette(c("grey90", "grey10"), space="rgb")(256);
    }else{
        suppressMessages(require(RColorBrewer));
        colors <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256));
    }

    imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
 
    if(is.na(width)){
        minW <- 630;
        myW <- nrow(hc.dat)*18 + 150;
        if(myW < minW){
            myW <- minW;
        }   
        w <- round(myW/72,2);
    }else if(width == 0){
        w <- 7.2;
        imgSet$heatmap<-imgName;
    }else{
        w <- 7.2;
    }

    myH <- ncol(hc.dat)*18 + 150;
    h <- round(myH/72,2);

   if(viewOpt == "overview"){
        if(is.na(width)){
            if(w > 9){
                w <- 9;
            }
        }else if(width == 0){
            if(w > 7.2){
                w <- 7.2;
            }
            imgSet$heatmap<-imgName;
        }else{
            w <- 7.2;
        }
        if(h > w){
            h <- w;
        }
    }

    # make the width smaller fro group average
    if(grp.ave){
        w <- nrow(hc.dat)*25 + 300;
        w <- round(w/72,2);
    }

    imgSet <<- imgSet;

    if(border){
        border.col<-"grey60";
    }else{
        border.col <- NA;
    }

    Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
    if(dataSet$cls.type == "disc"){
        library(pheatmap);
        annotation <- data.frame(class= hc.cls);
        rownames(annotation) <-rownames(hc.dat); 

        # set up color schema for samples
        if(palette== "gray"){
            cols <- GetColorSchema(T);
            uniq.cols <- unique(cols);
        }else{
            cols <- GetColorSchema();
            uniq.cols <- unique(cols);
        }
        names(uniq.cols) <- unique(as.character(dataSet$cls));
        ann_colors <- list(class= uniq.cols);

        pheatmap(t(hc.dat), 
            annotation=annotation, 
            fontsize=8, fontsize_row=8, 
            clustering_distance_rows = smplDist,
            clustering_distance_cols = smplDist,
            clustering_method = clstDist, 
            border_color = border.col,
            cluster_rows = colV, 
            cluster_cols = rowV,
            scale = scaleOpt, 
            color = colors,
            annotation_colors = ann_colors
            );
    }else{
        heatmap(hc.dat, Rowv = rowTree, Colv=colTree, col = colors, scale="column");
    }
    dev.off();
}

