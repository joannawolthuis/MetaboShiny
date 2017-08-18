################################################################
# over-representation analysis using hypergeometric tests
# The prob is calculated as obtain the equal or higher number
# of hitss using 1-phyper. Since phyper is cumulative prob,
# to get P(X>=hit.num) => P(X>(hit.num-1))
#
# Author: Jeff Xia
# jeff.xia@mcgill.ca
# McGill University
#############################################################

# method is "fisher" or "hyperg"
CalculateOraScore<-function(nodeImp, method){

    # make a clean dataSet$cmpd data based on name mapping
    # only valid kegg id will be used
    nm.map <- GetFinalNameMap();
    valid.inx <- !(is.na(nm.map$kegg)| duplicated(nm.map$kegg));
    ora.vec <- nm.map$kegg[valid.inx];
	q.size<-length(ora.vec);
    if(is.na(ora.vec) || q.size==0) {
        AddErrMsg("No valid KEGG compounds found!");
        return(0);
    }

    require(KEGGgraph);
    current.mset <- metpa$mset.list;
    uniq.count <- metpa$uniq.count;

    # check if need to be filtered against reference metabolome
    if(use.metabo.filter && exists('metabo.filter.kegg')){
        current.mset <- lapply(current.mset, function(x){x[x %in% metabo.filter.kegg]});
        analSet$ora.filtered.mset <- current.mset;
        uniq.count <- length(unique(unlist(current.mset, use.names=FALSE)));
    }

    hits <- lapply(current.mset, function(x){x[x %in% ora.vec]});
	hit.num<-unlist(lapply(hits, function(x){length(x)}), use.names=FALSE);
	set.size<-length(current.mset);
	set.num<-unlist(lapply(current.mset, length), use.names=FALSE);

	# prepare for the result table
	res.mat<-matrix(0, nrow=set.size, ncol=8);
    rownames(res.mat)<-names(current.mset);
    colnames(res.mat)<-c("Total", "Expected", "Hits", "Raw p", "-log(p)", "Holm adjust", "FDR", "Impact");

    if(nodeImp == "rbc"){
        imp.list <- metpa$rbc;
        analSet$topo.msg <- "Your selected node importance measure for topological analysis is \\textbf{relative betweenness centrality}.";
    }else{
        imp.list <- metpa$dgr;
        analSet$topo.msg <- "Your selected node importance measure for topological analysis is \\textbf{out degree centrality}.";
    }

    res.mat[,1]<-set.num;
    res.mat[,2]<-q.size*(set.num/uniq.count);
	res.mat[,3]<-hit.num;
    if(method == "fisher"){
        res.mat[,4]<-GetFisherPvalue(hit.num, q.size, set.num, uniq.count);
        analSet$rich.msg <- "The selected over-representation analysis method is \\textbf{Fishers' exact test}.";
    }else{
        # use lower.tail = F for P(X>x)
        res.mat[,4]<-phyper(hit.num-1, set.num, uniq.count-set.num, q.size, lower.tail=F);
        analSet$rich.msg <- "The selected over-representation analysis method is \\textbf{Hypergeometric test}.";
    }

    res.mat[,5] <- -log(res.mat[,4]);

    # adjust for multiple testing problems
    res.mat[,6] <- p.adjust(res.mat[,4], "holm");
    res.mat[,7] <- p.adjust(res.mat[,4], "fdr");

    # calculate the sum of importance
    res.mat[,8] <- mapply(function(x, y){sum(x[y])}, imp.list, hits);

    res.mat <- res.mat[hit.num>0,];
	ord.inx<-order(res.mat[,4], res.mat[,8]);
    analSet$ora.mat <- signif(res.mat[ord.inx,],5);
    analSet$ora.hits <- hits;
    analSet$node.imp <- nodeImp;

    analSet <<- analSet;
    save.mat <- analSet$ora.mat;
    rownames(save.mat) <- GetORA.pathNames();
    write.csv(save.mat, file="pathway_results.csv");
    return(1);
}

GetORA.pathNames<-function(){
    hit.inx <- match(rownames(analSet$ora.mat), metpa$path.ids);
    return(names(metpa$path.ids)[hit.inx]);
}

GetORA.keggIDs<-function(){
    kegg.vec<-rownames(analSet$ora.mat);
    kegg.vec <- paste("<a href=http://www.genome.jp/kegg-bin/show_pathway?",kegg.vec," target=_new>KEGG</a>", sep="");
    return (kegg.vec);
}

# only for human pathways
GetORA.smpdbIDs<-function(){
    return (SetupSMPDBLinks(rownames(analSet$ora.mat)));
}

GetFisherPvalue<-function(numSigMembers, numSigAll, numMembers, numAllMembers){
    z <- cbind(numSigMembers, numSigAll-numSigMembers, numMembers-numSigMembers, numAllMembers-numMembers-numSigAll+numSigMembers);
    z <- lapply(split(z, 1:nrow(z)), matrix, ncol=2);
    z <- lapply(z, fisher.test, alternative = 'greater');
    p.values <- as.numeric(unlist(lapply(z, "[[", "p.value"), use.names=FALSE));
    return(p.values);
}

CalculateQeaScore<-function(nodeImp, method){

    # first, need to make a clean dataSet$norm data based on name mapping
    # only contain valid kegg id will be used
    nm.map <- GetFinalNameMap();
    valid.inx <- !(is.na(nm.map$kegg)| duplicated(nm.map$kegg));
    nm.map <- nm.map[valid.inx,];
    orig.nms <- nm.map$query;

    kegg.inx <- match(colnames(dataSet$norm),orig.nms);
    hit.inx <- !is.na(kegg.inx);
    path.data<-dataSet$norm[,hit.inx];
    colnames(path.data) <- nm.map$kegg[kegg.inx[hit.inx]];

    # now, perform the enrichment analysis
    require(KEGGgraph);
    current.mset <- metpa$mset.list;
    uniq.count <- metpa$uniq.count;

    # check if a reference metabolome is applied
    if(use.metabo.filter && exists('metabo.filter.kegg')){
        current.mset<-lapply(current.mset, function(x) {x[x %in% metabo.filter.kegg]});
        analSet$qea.filtered.mset <- current.mset;
        uniq.count <- length(unique(unlist(current.mset), use.names=FALSE));
    }

    hits <- lapply(current.mset, function(x) {x[x %in% colnames(path.data)]});
	hit.inx<-unlist(lapply(hits, function(x) {length(x)}), use.names=FALSE) > 0;
    hits <- hits[hit.inx]; # remove no hits

    if(method == "gt"){
        require('globaltest');
        qea.obj <- gt(dataSet$cls, path.data, subsets=hits);
        qea.res <- result(qea.obj);
        match.num <- qea.res[,5];
        raw.p <- qea.res[,1];
        analSet$rich.msg <- "The selected pathway enrichment analysis method is \\textbf{Globaltest}.";
    }else{
        require('GlobalAncova');
        qea.obj <- NULL;
        qea.res <- GlobalAncova(xx=t(path.data), group=dataSet$cls, test.genes=hits, method="approx");
        match.num <- qea.res[,1];
        raw.p <- qea.res[,3];
        analSet$rich.msg <- "The selected pathway enrichment analysis method is \\textbf{GlobalAncova}.";
    }

    log.p <- -log(raw.p);
    # add adjust p values
    holm.p <- p.adjust(raw.p, "holm");
    fdr.p <- p.adjust(raw.p, "fdr");

    # calculate the impact values
    if(nodeImp == "rbc"){
        imp.list <- metpa$rbc[hit.inx];
        analSet$topo.msg <- "Your selected node importance measure for topological analysis is \\textbf{relative betweenness centrality}.";
    }else{
        imp.list <- metpa$dgr[hit.inx];
        analSet$topo.msg <- "Your selected node importance measure for topological analysis is \\textbf{out degree centrality}.";
    }
    imp.vec <- mapply(function(x, y){sum(x[y])}, imp.list, hits);

    set.num<-unlist(lapply(current.mset[hit.inx], length), use.names=FALSE);
    res.mat <- cbind(set.num, match.num, raw.p, log.p, holm.p, fdr.p, imp.vec);
    rownames(res.mat)<-rownames(qea.res);
    colnames(res.mat)<-c("Total Cmpd", "Hits", "Raw p", "-log(p)", "Holm adjust", "FDR", "Impact");

    ord.inx<-order(res.mat[,3], -res.mat[,7]);
    res.mat<-signif(res.mat[ord.inx,],5);

    # also calculate univariate p values
    # use lm model for t-tests (with var.equal=T), one-way anova, and linear regression (continuous);
    univ.p <- apply(as.matrix(path.data), 2, function(x) {
        tmp <- try(lm(as.numeric(dataSet$cls)~x));
        if(class(tmp) == "try-error") {
            return(NA);
        }else{
            tmp<-anova(tmp)
            return(tmp[1,5]);
        }
    });

    names(univ.p) <- colnames(path.data);

    # now store this new data
    dataSet$norm.path <- path.data;
    dataSet <<- dataSet;

    analSet$qea.hits <- hits;
    analSet$qea.univp <- signif(univ.p,7);
    analSet$qea.method <- method;
    analSet$qea.obj <- qea.obj;
    analSet$qea.mat <- res.mat;
    analSet$node.imp <- nodeImp;

    analSet <<- analSet;
    rownames(res.mat) <- GetQEA.pathNames();
    write.csv(res.mat, file="pathway_results.csv");

    return(1);
}

GetQEA.pathNames<-function(){
    hit.inx <- match(rownames(analSet$qea.mat),metpa$path.ids);
    return(names(metpa$path.ids)[hit.inx]);
}

GetQEA.keggIDs<-function(){
    kegg.vec<-rownames(analSet$qea.mat);
    kegg.vec <- paste("<a href=http://www.genome.jp/kegg-bin/show_pathway?",kegg.vec," target=_new>KEGG</a>", sep="");
    return (kegg.vec);
}

GetQEA.smpdbIDs<-function(){
    return (SetupSMPDBLinks(rownames(analSet$qea.mat)));
}

# only works for human (hsa.rda) data
SetupSMPDBLinks <- function(kegg.ids){
    smpdb.vec <- names(metpa$path.smps)[match(kegg.ids,metpa$path.smps)]
    lk.len <- length(smpdb.vec);
    all.lks <- vector(mode="character", length=lk.len);
    for(i in 1:lk.len){
          lks <- strsplit(smpdb.vec[i], "; ")[[1]];
          if(!is.na(lks[1])){
               # all.lks[i]<-paste("<a href=http://pathman.smpdb.ca/pathways/",lks,"/pathway target=_new>SMP</a>", sep="", collapse="\n");
                all.lks[i]<-paste("<a href=http://www.smpdb.ca/view/",lks," target=_new>SMP</a>", sep="", collapse="\n");
          }
    }
    return (all.lks);
}

# given a metset inx, return hmtl highlighted pathway cmpds
GetHTMLPathSet<-function(msetNm){
    pathid <- metpa$path.ids[msetNm]; 
    mset <- metpa$mset.list[[pathid]];

    hits <- NULL;
    if(analSet$type=="pathora"){
        hits <- analSet$ora.hits;
    }else{
        hits <- analSet$qea.hits;
    }

    # highlighting with different colors
    red.inx <- which(mset %in% hits[[pathid]]);
    
    # use actual cmpd names
    nms <- names(mset);
    nms[red.inx] <- paste("<font color=\"red\">", "<b>", nms[red.inx], "</b>", "</font>",sep="");
    return(cbind(msetNm, paste(nms, collapse="; ")));
}