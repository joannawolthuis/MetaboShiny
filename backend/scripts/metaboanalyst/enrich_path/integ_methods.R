##################################################
## R script for INMEX
## Description: Integrative analysis methods
##
## Author: Jeff Xia, jeff@hancocklab.com
## University of British Columbia, Canada
###################################################

PerformIntegCmpdMapping <- function(cmpdIDs, org, idType){
    current.msg <<- NULL;
    dataSet$cmpd.orig <- cmpdIDs;
    dataSet$cmpd.org <- org;
    if(idType == "name"){
        cmpd.mat <- getDataFromTextInput(cmpdIDs, "tab");
    }else{ # all other id should not contains space
        cmpd.mat <- getDataFromTextInput(cmpdIDs, "space");
    }
    dataSet$cmpd <- rownames(cmpd.mat); # this is for compatible with name_match function
    dataSet$cmpd.mat <- cmpd.mat;
    dataSet <<- dataSet;
    CrossReferencing(idType, hmdb=T, pubchem=T, chebi=F, kegg=T, metlin=F);
    return(1);
}

PerformIntegGeneMapping <- function(geneIDs, org, idType){
    current.msg <<- NULL;
    dataSet$q.type.gene <- idType;
    dataSet$gene.orig <- geneIDs;
    dataSet$gene.org <- org;
    gene.mat <- getDataFromTextInput(geneIDs);
    gene.vec <- rownames(gene.mat);
    dataSet$gene.mat <- gene.mat;
    dataSet$gene <- gene.vec;
    dataSet <<- dataSet;

    enIDs <- doGeneIDMapping(gene.vec, org, idType);
    gene.name.map <- list();
    gene.name.map$hit.values <- enIDs;
    gene.name.map$match.state <- ifelse(is.na(enIDs), 0, 1);
    gene.name.map<<- gene.name.map;

    current.msg <<- c(paste("A total of ", length(unique(enIDs)), "unique genes were uploaded."));
    na.inx <- is.na(enIDs);
    if(sum(!na.inx) > 0){
        return(1);
    }else{
        return(0);
    }
}

RemoveCmpd <- function(inx){
    dataSet$cmpd <- dataSet$cmpd[-inx];
    dataSet <<- dataSet;

    name.map$hit.inx <- name.map$hit.inx[-inx];
    name.map$hit.values <- name.map$hit.values[-inx];
    name.map$match.state <- name.map$match.state[-inx];
    name.map <<- name.map;
}

RemoveGene <- function(inx){
    dataSet$gene <- dataSet$gene[-inx];
    dataSet$gene.mat <- dataSet$gene.mat[-inx, ,drop=F];
    dataSet <<- dataSet;

    gene.name.map$hit.inx <- gene.name.map$hit.inx[-inx];
    gene.name.map$hit.values <- gene.name.map$hit.values[-inx];
    gene.name.map$match.state <- gene.name.map$match.state[-inx];
    gene.name.map <<- gene.name.map;
}

PrepareIntegData <- function(){

    done <- 0;
    # prepare gene list
    if(!is.null(dataSet$gene.mat)){
        gene.mat <- dataSet$gene.mat;
        enIDs <- gene.name.map$hit.values;
        rownames(gene.mat) <- enIDs;

        na.inx <- is.na(enIDs);
        gene.mat <- gene.mat[!na.inx, ,drop=F];
        gene.mat <- RemoveDuplicates(gene.mat);
        current.msg <<- c(paste("A total of ", nrow(gene.mat), "unique genes were uploaded."));

        if(!exists("inmex.imps")){
            inmex.imps <- list();
        }
        inmex.imps$gene.mat <- gene.mat;
        inmex.imps <<- inmex.imps;
        done <- 1;
    }

    # prepare compound list
    if(!is.null(dataSet$cmpd.mat)){
        nm.map <- GetFinalNameMap();
        valid.inx <- !(is.na(nm.map$kegg)| duplicated(nm.map$kegg));
        cmpd.vec <- nm.map$query[valid.inx];
        kegg.id <- nm.map$kegg[valid.inx];

        cmpd.mat <- dataSet$cmpd.mat;
        hit.inx <- match(cmpd.vec, rownames(cmpd.mat));

        cmpd.mat <- cmpd.mat[hit.inx, ,drop=F];
        rownames(cmpd.mat) <- kegg.id;
        cmpd.mat <- RemoveDuplicates(cmpd.mat);
        current.msg <<- c(current.msg, paste("A total of ", nrow(cmpd.mat), "unique compounds were found."));
        inmex.imps$cmpd.mat <- cmpd.mat;
        done <- 1;
    }
    inmex.imps <<- inmex.imps;
    return(done);
}

# used for integrative analysis 
# as well as general pathways analysis for meta-analysis results
PerformIntegPathwayAnalysis <- function(topo="dc", enrich="hyper", libOpt="integ"){

    LoadKEGGLib(libOpt);
    set.size<-length(inmexpa$mset.list);
    ms.list <- current.mset.list;

    # prepare for the result table
    res.mat<-matrix(0, nrow=set.size, ncol=7);
    rownames(res.mat)<-names(inmexpa$path.ids);
    colnames(res.mat)<-c("Total", "Expected", "Hits", "P.Value", "Topology", "PVal.Z",  "Topo.Z");

    inmex.method <<- libOpt;
    path.mat <<- NULL;

    if(libOpt == "genetic" && !is.null(inmex.imps$gene.mat)){
        gene.sbls <- doEntrez2SymbolMapping(rownames(inmex.imps$gene.mat));
        inmex.imps$gene.mat <- cbind(Name=gene.sbls, inmex.imps$gene.mat);
        gene.vec <- paste(inmex.org, ":", rownames(inmex.imps$gene.mat), sep="");
        rownames(inmex.imps$gene.mat) <- gene.vec;
        write.csv(inmex.imps$gene.mat, file="MetaboAnalyst_result_genes.csv");

        impMat <- inmex.imps$gene.mat;
        uniq.count <- inmexpa$uniq.gene.count;
        uniq.len <- inmexpa$gene.counts;
    }else if(libOpt == "metab" && !is.null(inmex.imps$cmpd.mat)){
        cmpd.nms <- doKEGG2NameMapping(rownames(inmex.imps$cmpd.mat));
        inmex.imps$cmpd.mat <- cbind(Name=cmpd.nms, inmex.imps$cmpd.mat);
        cmpd.vec <- paste("cpd:", rownames(inmex.imps$cmpd.mat), sep="");
        rownames(inmex.imps$cmpd.mat) <- cmpd.vec;
        write.csv(inmex.imps$cmpd.mat, file="MetaboAnalyst_result_cmpds.csv");

        impMat <- inmex.imps$cmpd.mat;
        uniq.count <- inmexpa$uniq.cmpd.count
        uniq.len <- inmexpa$cmpd.counts;
    }else{ # integ
        impMat <- NULL;
        uniq.count <- uniq.len <- 0;
        if(!is.null(inmex.imps$cmpd.mat)){
            cmpd.nms <- doKEGG2NameMapping(rownames(inmex.imps$cmpd.mat));
            inmex.imps$cmpd.mat <- cbind(Name=cmpd.nms, inmex.imps$cmpd.mat);
            cmpd.vec <- paste("cpd:", rownames(inmex.imps$cmpd.mat), sep="");
            rownames(inmex.imps$cmpd.mat) <- cmpd.vec;
            write.csv(inmex.imps$cmpd.mat, file="MetaboAnalyst_result_cmpds.csv");
            impMat <- inmex.imps$cmpd.mat;
            uniq.count <- inmexpa$uniq.cmpd.count
            uniq.len <- inmexpa$cmpd.counts;
        }
        if(!is.null(inmex.imps$gene.mat)){
            gene.sbls <- doEntrez2SymbolMapping(rownames(inmex.imps$gene.mat));
            inmex.imps$gene.mat <- cbind(Name=gene.sbls, inmex.imps$gene.mat);
            gene.vec <- paste(inmex.org, ":", rownames(inmex.imps$gene.mat), sep="");
            rownames(inmex.imps$gene.mat) <- gene.vec;
            write.csv(inmex.imps$gene.mat, file="MetaboAnalyst_result_genes.csv");
            impMat <- rbind(impMat, inmex.imps$gene.mat);
            uniq.count <- uniq.count + inmexpa$uniq.gene.count;
            uniq.len <- uniq.len + inmexpa$gene.counts;
        }
    }
    # now project to pathways
    # combine results for genes and cmpds
    ora.vec <- rownames(impMat);
    colnames(impMat) <- c("Name", "logFC");
    impMat <- data.frame(Name=impMat[,1], logFC=as.numeric(impMat[,2]));
    rownames(impMat) <- ora.vec;

    # need to cut to the universe covered by the pathways, not all genes 
    ora.vec <- ora.vec[ora.vec %in% current.universe]
    q.size<-length(ora.vec);

    # note, we need to do twice one for nodes (for plotting)
    # one for query for calculating, as one node can be associated with multiple matches
    # get the matched nodes on each pathway
    hits.path <- lapply(ms.list, function(x) {unlist(lapply(x, function(var){any(var%in%ora.vec);}),use.names=FALSE)});
    names(hits.path) <- inmexpa$path.ids;

    # get the matched query for each pathway
    hits.query <- lapply(ms.list, function(x) {ora.vec%in%unlist(x);});

    hit.num<-unlist(lapply(hits.query, function(x){sum(x)}), use.names=FALSE);

    if(sum(hit.num) == 0){
        AddErrMsg("No hits found for your input!");
        return(0);
    }

    set.num<-uniq.len;
    res.mat[,1]<-set.num;
    res.mat[,2]<-q.size*(set.num/uniq.count);
    res.mat[,3]<-hit.num;

    # use lower.tail = F for P(X>x)
    if(enrich=="hyper"){
        res.mat[,4]<-phyper(hit.num-1, set.num, uniq.count-set.num, q.size, lower.tail=F);
    }else if(enrich == "fisher"){
        res.mat[,4]<-GetFisherPvalue(hit.num, q.size, set.num, uniq.count);
    }else{
        print("Not defined enrichment method!");
        print(enrich);
    }

    # adjust for multiple testing problems
    # res.mat[,5] <- p.adjust(res.mat[,4], "fdr");

    # toplogy test
    if(topo == "bc"){
        imp.list <- inmexpa$bc;
    }else if(topo == "dc"){
        imp.list <- inmexpa$dc;
    }else if(topo == "cc"){
        imp.list <- inmexpa$cc;       
    }else{
        print("Not defined topological measure!");
        print(topo);
    }

    # now, perform topological analysis		
    # calculate the sum of importance
    res.mat[,5] <- mapply(function(x, y){sum(x[y])}, imp.list, hits.path);

    # now add two more columns for the scaled values
    res.mat[,6] <- scale(-log(res.mat[,4]));
    res.mat[,7] <- scale(res.mat[,5]);

    res.mat <- res.mat[hit.num>0,];
    ord.inx<-order(res.mat[,4], res.mat[,5]);

    res.mat <- signif(res.mat[ord.inx,],5);

    #get gene symbols
    resTable <- data.frame(Pathway=rownames(res.mat), res.mat);

    # now save to different formats
    # csv
    write.csv(resTable, file="MetaboAnalyst_result_pathway.csv", row.names=F);

    path.hits <<- hits.path;

    # store results from individual analysis
    inmex.impMat <<- impMat; 
    inmex.impTopo <<- imp.list;
    path.mat <<- resTable;

    SetBarParams();
    return(1);
}

GetFisherPvalue<-function(numSigMembers, numSigAll, numMembers, numAllMembers){
    z <- cbind(numSigMembers, numSigAll-numSigMembers, numMembers-numSigMembers, numAllMembers-numMembers-numSigAll+numSigMembers);
    z <- lapply(split(z, 1:nrow(z)), matrix, ncol=2);
    z <- lapply(z, fisher.test, alternative = 'greater');
    p.values <- as.numeric(unlist(lapply(z, "[[", "p.value")));
    return(p.values);
}

# plot an scatterplot bubble chart overview of the matched pathways
# x axis is the pathway impact factor
# y axis is the p value (from ORA) 
# return the circle information
SetBarParams <- function(){

    y <- path.mat$Topology;
    x <- path.mat$P.Value;
    x = -log(x);

    x <- scale(x);
    y <- scale(y);
    base <- abs(min(c(x,y)));

    x <- x + base;
    y <- y + base;
    # names(y) <-  rownames(path.mat);

    # set circle size according to
    # sum of p and topo (since they
    # alrealy bring to same range

    # we do twice to reduce difference for plotting
    radi.vec <- sqrt(x+y);

    resTable <- data.frame(x=x, y=y); 
    rownames(resTable) <- rownames(path.mat);

    # display only top 100 sorted by p
    if(nrow(resTable) > 20){
        resTable <- resTable[1:20,];
    }

    resTable <- resTable[nrow(resTable):1,];
    bar.data <<- resTable;
}

GetBarParams<-function(){
     return(as.matrix(bar.data));
}

GetBarNames<-function(){
    # single quote apostrophe caused trouble 
    return(rownames(bar.data));
}

GetIntegResultPathIDs<-function(){
    inmexpa$path.ids[rownames(path.mat)];
}

GetIntegResultPathNames<-function(){
     return(rownames(path.mat));
}

GetIntegResultColNames<-function(){
    return(colnames(path.mat)[-1]);
}

GetIntegResultMatrix<-function(){
   return(as.matrix(path.mat[,-1]));
}

GetGeneMappingResultTable<-function(){

    qvec <- dataSet$gene;
    enIDs <- gene.name.map$hit.values;

    # style for highlighted background for unmatched names
    pre.style<-NULL;
    post.style<-NULL;

    # style for no matches
    if(dataSet$q.type.gene == "name"){
        no.prestyle<-"<strong style=\"background-color:yellow; font-size=125%; color=\"black\">";
        no.poststyle<-"</strong>";
    }else{
        no.prestyle<-"<strong style=\"background-color:red; font-size=125%; color=\"black\">";
        no.poststyle<-"</strong>";
    }

    # contruct the result table with cells wrapped in html tags
    # the unmatched will be highlighted in different background
    html.res<-matrix("", nrow=length(qvec), ncol=5);
    csv.res<-matrix("", nrow=length(qvec), ncol=5);
    colnames(csv.res)<-c("Query", "Entrez", "Symbol", "Name", "Comment");

    db.path <- paste("../../libs/", inmex.org, "/entrez.csv", sep="");
    gene.db <- .readDataTable(db.path);
    hit.inx <- match(enIDs, gene.db[, "gene_id"]);
    hit.values<-gene.name.map$hit.values;
    match.state<-gene.name.map$match.state;
    gene.name.map$hit.inx <- hit.inx;
    gene.name.map <<- gene.name.map;

    for (i in 1:length(qvec)){
        if(match.state[i]==1){
            pre.style<-"";
            post.style="";
        }else{ # no matches
            pre.style<-no.prestyle;
            post.style<-no.poststyle;
        }
        hit <-gene.db[hit.inx[i], ,drop=F];

        html.res[i, ]<-c(paste(pre.style, qvec[i], post.style, sep=""),
                       paste(ifelse(match.state[i]==0 || is.na(hit$gene_id),"-", paste("<a href=http://www.ncbi.nlm.nih.gov/gene/", hit$gene_id, " target='_blank'>",hit$gene_id,"</a>", sep="")),  sep=""),
                       paste(ifelse(match.state[i]==0 || is.na(hit$symbol), "-", paste("<a href=http://www.ncbi.nlm.nih.gov/gene/", hit$gene_id, " target='_blank'>", hit$symbol,"</a>", sep="")), sep=""),
                       paste(ifelse(match.state[i]==0 || is.na(hit$name),"-", paste("<a href=http://www.ncbi.nlm.nih.gov/gene/", hit$gene_id, " target='_blank'>",hit$name,"</a>", sep="")), sep=""),
                       ifelse(match.state[i]!=1,"View",""));
        csv.res[i, ]<-c(qvec[i],
                        ifelse(match.state[i]==0, "NA", hit$gene_id),
                        ifelse(match.state[i]==0, "NA", hit$symbol),
                        ifelse(match.state[i]==0, "NA", hit$name),
                        match.state[i]);
     }

     # store the value for report
     dataSet$gene.map.table <- csv.res;
     dataSet <<- dataSet;
     write.csv(csv.res, file="gene_name_map.csv", row.names=F);
     return(as.vector(html.res));
}

GetGeneHitsRowNumber<-function(){
    return (length(gene.name.map$match.state));
}

