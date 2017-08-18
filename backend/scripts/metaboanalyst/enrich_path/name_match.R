##################################################
## Description: mapping b/w names & database identifiers
##
## Author: Jeff Xia, jeff.xia@mcgill.ca
## McGill University, Canada
## License: GNU GPL (>= 2)
###################################################

# given a list of compound names or ids, find matched name or ids from selected databases
CrossReferencing <- function(q.type, hmdb=T, pubchem=T, chebi=F, kegg=T, metlin=F){

       # record the filter for 8 major databases
       return.cols <<- c(hmdb, pubchem, chebi, kegg, metlin);

       # record all the data
       if(!exists("name.map")){
            name.map <<- list();
       }

       # distribute job
       dataSet$q.type <- q.type;
       dataSet <<- dataSet;
       MetaboliteMappingExact(q.type);

       # do some sanity check
       todo.inx <-which(is.na(name.map$hit.inx));
       if(length(todo.inx)/length(name.map$hit.inx) > 0.5){
            nmcheck.msg <<- c(0, "Over half of the compound IDs could not be matched to our database. Please make 
                            sure that correct compound IDs or common compound names are used.");
       }else if (length(todo.inx) > 15){
            nmcheck.msg <<- c(2, "There are >15 compounds without matches. You can either proceed or if necessary, update these compound IDs and upload again.");        
       }else{
            nmcheck.msg <<- c(1, "Name matching OK, please inspect (and manual correct) the results then proceed.");   
       }
}

# Mapping from different metabolite IDs
# For compound names to other id, can do exact or approximate match
# For other IDs, except HMDB ID, all other may return multiple /non-unique hits
# multiple hits or non-unique hits will all users to manually select
MetaboliteMappingExact<-function(q.type){
       qvec <- dataSet$cmpd;

       # variables to record results
       hit.inx = vector(mode='numeric', length=length(qvec)); # record hit index, initial 0
       match.values = vector(mode='character', length=length(qvec)); # the best matched values (hit names), initial ""
       match.state = vector(mode='numeric', length=length(qvec));  # match status - 0, no match; 1, exact match; initial 0 

       cmpd.db <- readRDS("../../libs/compound_db.rds");
       if(q.type == "hmdb"){
            hit.inx <- match(tolower(qvec), tolower(cmpd.db$hmdb));
            match.values <- cmpd.db$name[hit.inx];
            match.state[!is.na(hit.inx)] <- 1;
       }else if(q.type == "pubchem"){
            hit.inx <- match(tolower(qvec), tolower(cmpd.db$pubchem));
            match.values <- cmpd.db$name[hit.inx];
            match.state[!is.na(hit.inx)] <- 1;
       }else if(q.type == "chebi"){
            hit.inx <- match(tolower(qvec), tolower(cmpd.db$chebi));
            match.values <- cmpd.db$name[hit.inx];
            match.state[!is.na(hit.inx)] <- 1;
       }else if(q.type == "metlin"){
            hit.inx <- match(tolower(qvec), tolower(cmpd.db$metlin));
            match.values <- cmpd.db$name[hit.inx];
            match.state[!is.na(hit.inx)] <- 1;
       }else if(q.type == "kegg"){
            hit.inx <- match(tolower(qvec), tolower(cmpd.db$kegg));
            hit.inx2 <- match(tolower(qvec), rev(tolower(cmpd.db$kegg)));

            # unique hits
            nonuniq.hits <- hit.inx + hit.inx2 != nrow(cmpd.db) + 1;
            hit.inx[nonuniq.hits] <- NA;
            match.values <- cmpd.db$name[hit.inx];
            match.state[!is.na(hit.inx)] <- 1;

       }else if(q.type == "name"){
            # first find exact match to the common compound names
            hit.inx <- match(tolower(qvec), tolower(cmpd.db$name));
            match.values <- cmpd.db$name[hit.inx];
            match.state[!is.na(hit.inx)] <- 1;

            # then try to find exact match to synanyms for the remaining unmatched query names one by one
            syn.db <- readRDS("../../libs/syn_nms.rds")
            syns.list <-  syn.db$syns.list;
            todo.inx <-which(is.na(hit.inx));
            if(length(todo.inx) > 0){
                for(i in 1:length(syns.list)){
                    syns <-  syns.list[[i]];
                    hitInx <- match(tolower(qvec[todo.inx]), tolower(syns));

                    hitPos <- which(!is.na(hitInx));
                    if(length(hitPos)>0){
                        # record matched ones
                        orig.inx<-todo.inx[hitPos];
                        hit.inx[orig.inx] <- i;                  
                        # match.values[orig.inx] <- syns[hitInx[hitPos]];  # show matched synnames
                        match.values[orig.inx] <- cmpd.db$name[i];    # show common name
                        match.state[orig.inx] <- 1;

                        # update unmatched list
                        todo.inx<-todo.inx[is.na(hitInx)];
                    }
                    if(length(todo.inx) == 0) break;
                }
            }
      }else{
            print(paste("Unknown compound ID type:", q.type));
            # guess a mix of kegg and hmdb ids
            hit.inx <- match(tolower(qvec), tolower(cmpd.db$hmdb));
            hit.inx2 <- match(tolower(qvec), tolower(cmpd.db$kegg));
            nohmdbInx <- is.na(hit.inx);
            hit.inx[nohmdbInx]<-hit.inx2[nohmdbInx]
      }
      # empty memory
      gc();

      name.map$hit.inx <- hit.inx;
      name.map$hit.values <- match.values;
      name.map$match.state <- match.state;
      name.map <<- name.map;
}

GetMappingResultTable<-function(){

       qvec <- dataSet$cmpd;
       if(is.null(qvec)){
            return();
       }
       # style for highlighted background for unmatched names
       pre.style<-NULL;
       post.style<-NULL;

       # style for no matches
       if(dataSet$q.type == "name"){
            no.prestyle<-"<strong style=\"background-color:yellow; font-size=125%; color=\"black\">";
            no.poststyle<-"</strong>";
       }else{
            no.prestyle<-"<strong style=\"background-color:red; font-size=125%; color=\"black\">";
            no.poststyle<-"</strong>";
       }

       hit.inx<-name.map$hit.inx;
       hit.values<-name.map$hit.values;
       match.state<-name.map$match.state;

       # contruct the result table with cells wrapped in html tags
       # the unmatched will be highlighted in different background
       html.res<-matrix("", nrow=length(qvec), ncol=8);
       csv.res<-matrix("", nrow=length(qvec), ncol=8);
       colnames(csv.res)<-c("Query", "Match", "HMDB", "PubChem", "ChEBI", "KEGG", "METLIN", "Comment");
       cmpd.db <- readRDS("../../libs/compound_db.rds");
       for (i in 1:length(qvec)){
           if(match.state[i]==1){
               pre.style<-"";
               post.style="";
           }else{ # no matches
               pre.style<-no.prestyle;
               post.style<-no.poststyle;
           }
           hit <-cmpd.db[hit.inx[i], ,drop=F];
           html.res[i, ]<-c(paste(pre.style, qvec[i], post.style, sep=""),
                       paste(ifelse(match.state[i]==0, "", hit.values[i]), sep=""),
                       paste(ifelse(match.state[i]==0 || is.na(hit$hmdb_id) || hit$hmdb_id=="" || hit$hmdb_id=="NA","-", paste("<a href=http://www.hmdb.ca/metabolites/", hit$hmdb_id, " target='_blank'>",hit$hmdb_id,"</a>", sep="")),  sep=""),
                       paste(ifelse(match.state[i]==0 || is.na(hit$pubchem_id) || hit$pubchem_id=="" || hit$pubchem_id=="NA", "-", paste("<a href=http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=", hit$pubchem_id," target='_blank'>", hit$pubchem_id,"</a>", sep="")), sep=""),
                       paste(ifelse(match.state[i]==0 || is.na(hit$chebi_id) || hit$chebi_id==""|| hit$chebi_id=="NA","-", paste("<a href=http://www.ebi.ac.uk/chebi/searchId.do?chebiId=", hit$chebi_id, " target='_blank'>",hit$chebi_id,"</a>", sep="")), sep=""),
                       paste(ifelse(match.state[i]==0 || is.na(hit$kegg_id) || hit$kegg_id==""|| hit$kegg_id=="NA","-",paste("<a href=http://www.genome.jp/dbget-bin/www_bget?", hit$kegg_id, " target='_blank'>", hit$kegg_id,"</a>", sep="")), sep=""),
                       paste(ifelse(match.state[i]==0 || is.na(hit$metlin_id) || hit$metlin_id==""|| hit$metlin_id=="NA","-",paste("<a href=http://metlin.scripps.edu/metabo_info.php?molid=", hit$metlin_id," target='_blank'>",hit$metlin_id,"</a>", sep="")), sep=""),
                       ifelse(match.state[i]!=1,"View",""));
           csv.res[i, ]<-c(qvec[i],
                           ifelse(match.state[i]==0, "NA", hit.values[i]),
                           ifelse(match.state[i]==0, "NA", hit$hmdb_id),
                           ifelse(match.state[i]==0, "NA", hit$pubchem_id),
                           ifelse(match.state[i]==0, "NA", hit$chebi_id),
                           ifelse(match.state[i]==0, "NA", hit$kegg_id),
                           ifelse(match.state[i]==0, "NA", hit$metlin_id),
                           match.state[i]);
     }
     # return only columns user selected

     # add query and match columns at the the beginning, and 'Detail' at the end
     return.cols <- c(TRUE, TRUE, return.cols, TRUE);
     html.res <- html.res[,return.cols, drop=F];
     csv.res <- csv.res[,return.cols, drop=F];

     # store the value for report
     dataSet$map.table <- csv.res;
     dataSet <<- dataSet;
     write.csv(csv.res, file="name_map.csv", row.names=F);
     return(as.vector(html.res));
}

GetHitsRowNumber<-function(){
    return(length(name.map$hit.inx));
}

PerformDetailMatch<-function(inx){
    q <- dataSet$cmpd[inx];
    if(dataSet$q.type == "name"){
        PerformApproxMatch(q);
    }else{
        PerformMultiMatch(q);
    }
}

PerformMultiMatch<-function(q){
    cmpd.db <- readRDS("../../libs/compound_db.rds");
    matched.inx <- which(cmpd.db$kegg %in% q);
    if(length(matched.inx) > 0) {
        # record all the candidates,
        candidates <- cbind(matched.inx, cmpd.db$name[matched.inx]);
        candidates <<- candidates;
        return();   
    }
    candidates <<- NULL;
}

PerformApproxMatch<-function(q){
    cmpd.db <- readRDS("../../libs/compound_db.rds");
    # only for none lipids
    nonLipidInx <- cmpd.db$lipid == 0;
    com.nms <- cmpd.db$name[nonLipidInx];

    syn.db <- readRDS("../../libs/syn_nms.rds")
    syns.vec <- syn.db$syns.vec[nonLipidInx];
    syns.list <- syn.db$syns.list[nonLipidInx];

    matched.dist <- NULL;
    q.length <- nchar(q);
    s <- c(0, 0.1, 0.2);
    for (j in s) {
        new.q <- q;
        if(q.length > 32){ # note: agrep fail for exact match when length over 32 characters
            new.q<-substr(q, 1, 32);
        }
        matched <- FALSE;
        matched.inx <- agrep(new.q, syns.vec, ignore.case=T, max.distance=j, useBytes=T);
        if(length(matched.inx) > 0) {
            # record all the candidates,
            # don't use cbind, since all will be converted to character mode
            # for data.frame specify "stringsAsFactors" to prevent convert value col into factor
            candidates <- data.frame(index=vector(mode = "numeric", length=length(matched.inx)),
                                value=vector(mode = "character", length=length(matched.inx)),
                                score=vector(mode = "numeric", length=length(matched.inx)),
                                stringsAsFactors = FALSE);

            for(n in 1:length(matched.inx)){
                nm.vec<-syns.list[[matched.inx[n]]];
                # try approximate match, note: in some cases, split into element will break match using whole string
                hit3.inx <- agrep(q,nm.vec,ignore.case=T, max.distance=j, useBytes=T);
                if(length(hit3.inx)>0){
                    hit3.nm <- vector(mode = "character", length=length(hit3.inx));
                    hit3.score <- vector(mode = "numeric", length=length(hit3.inx));
                    for(k in 1:length(hit3.inx)){
                        idx <- hit3.inx[k];
                        hit3.nm[k] <- nm.vec[idx];
                        hit3.score[k] <- j + abs(nchar(nm.vec[idx])-nchar(q))/(10*nchar(q));
                    }

                    # now get the best match, the rule is that the first two character should matches
                    # first check if first two character are digits or characters, otherwise will cause error
                    matches2 <- c();
                    if(length(grep("^[1-9a-z]{2}", q, ignore.case=T))>0){
                        matches2 <- grep(paste("^", substr(q, 1, 2), sep=""), hit3.nm);
                    }else if (length(grep("^[1-9a-z]", q, ignore.case=T))>0){
                        matches2 <- grep(paste("^", substr(q, 1, 1), sep=""), hit3.nm);
                    }

                    if(length(matches2)>0){
                        hit3.score[matches2] <- hit3.score[matches2] - 0.05;
                    }

                    best.inx<-which(hit3.score==min(hit3.score))[1];
                    candidates[n,1]<-matched.inx[n];
                #    candidates[n,2]<-hit3.nm[best.inx]; # show matched syn names
                    candidates[n,2]<-com.nms[matched.inx[n]] # show common names
                    candidates[n,3]<-hit3.score[best.inx];
                }      
            }
            rm.inx <- is.na(candidates[,2]) | candidates[,2]=="NA" | candidates[,2]=="";
            candidates<-candidates[!rm.inx, ];  
            candidates<-candidates[order(candidates[,3], decreasing=F), , drop=F];    
            if(nrow(candidates) > 10){
                candidates<-candidates[1:10,];
            }
            candidates <<- candidates;
            return();   
        }
    }
    candidates <<- NULL;
}


# Get all candidate compound names for a given index - 3 col - inx, name, score
GetCandidateList<-function(){
     # contruct the result table with cells wrapped in html tags
     # the unmatched will be highlighted in different background
     can.mat<-matrix("", nrow=nrow(candidates)+1, ncol= 6);
     cmpd.db <- readRDS("../../libs/compound_db.rds");

     # need to exclude lipids, to be consistent with approx matching part so that same index can be used to fetch db entries
     nonLipidInx <- cmpd.db$lipid == 0;
     cmpd.db <-cmpd.db[nonLipidInx,];

     for (i in 1:nrow(candidates)){
           hit.inx <- candidates[i, 1];
           hit.name <- candidates[i, 2];
           hit <-cmpd.db[hit.inx, ,drop=F];
           can.mat[i, ]<-c(hit.name,
                       paste(ifelse(hit$hmdb_id=="NA","", paste("<a href=http://www.hmdb.ca/metabolites/", hit$hmdb_id, " target='_blank'>",hit$hmdb_id,"</a>", sep="")), sep=""),
                       paste(ifelse(hit$pubchem_id=="NA", "", paste("<a href=http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=", hit$pubchem_id," target='_blank'>", hit$pubchem_id,"</a>", sep="")), sep=""),
                       paste(ifelse(hit$chebi_id=="NA","", paste("<a href=http://www.ebi.ac.uk/chebi/searchId.do?chebiId=", hit$chebi_id, " target='_blank'>",hit$chebi_id,"</a>", sep="")), sep=""),
                       paste(ifelse(hit$kegg_id=="NA","",paste("<a href=http://www.genome.jp/dbget-bin/www_bget?", hit$kegg_id, " target='_blank'>", hit$kegg_id,"</a>", sep="")), sep=""),
                       paste(ifelse(hit$metlin_id=="NA","",paste("<a href=http://metlin.scripps.edu/metabo_info.php?molid=", hit$metlin_id," target='_blank'>",hit$metlin_id,"</a>", sep="")), sep=""));
      }
      # add "none" option
      can.mat[nrow(candidates)+1,] <- c("None of the above", "", "", "", "", "");

      # add the hit columns
      return.cols <- c(TRUE, return.cols);

      return(as.vector(can.mat[,return.cols, drop=F]));
}

GetCanListRowNumber<-function(inx){
    return(nrow(candidates)+1); # include the "none" row
}

GetQuery<-function(inx){
    return(dataSet$cmpd[inx]);
}

# set matched name based on user selection from all potential hits
# note: to change object in the enclosing enviroment, use "<<-"
SetCandidate<-function(query_inx, can_inx){
      can_mat<- candidates;
      cmpd.db <- readRDS("../../libs/compound_db.rds");
      if(can_inx <= nrow(can_mat)){
            name.map$hit.inx[query_inx] <- can_mat[can_inx,1];
            name.map$hit.values[query_inx] <- can_mat[can_inx,2];
            name.map$match.state[query_inx] <- 1;
            # re-generate the CSV file
            hit <-cmpd.db[name.map$hit.inx[query_inx], ,drop=F];
            csv.res <- dataSet$map.table;
            if(ncol(csv.res) > 6){ # general utilities
                csv.res[query_inx, ]<-c(csv.res[query_inx, 1],
                           name.map$hit.values[query_inx],
                           hit$hmdb_id,
                           hit$pubchem_id,
                           hit$chebi_id,
                           hit$kegg_id,
                           hit$metlin_id,
                           1);
            }else{
                csv.res[query_inx, ]<-c(csv.res[query_inx, 1],
                           name.map$hit.values[query_inx],
                           hit$hmdb_id,
                           hit$pubchem_id,
                           hit$kegg_id,
                           1);
            }
            write.csv(csv.res, file="name_map.csv", row.names=F);
            dataSet$map.table <- csv.res;
            dataSet <<- dataSet;
      }else{ #no match
            name.map$hit.inx[query_inx] <- 0;
            name.map$hit.values[query_inx] <- "";
            name.map$match.state[query_inx] <- 0;
            name.map <<- name.map;
      }
}

# return the final (after user selection) map as dataframe
# three col, original name, HMDB name and KEGG ID,
# for Enrichment and pathway analysis, respectively
GetFinalNameMap<-function(){

    hit.inx<-name.map$hit.inx;
    hit.values<-name.map$hit.values;
    match.state<-name.map$match.state;

    qvec <- dataSet$cmpd;
    nm.mat<-matrix(nrow=length(qvec), ncol=3);
    colnames(nm.mat)<-c("query", "hmdb", "kegg");
    cmpd.db <- readRDS("../../libs/compound_db.rds");
    for (i in 1:length(qvec)){
        hit <-cmpd.db[hit.inx[i], ,drop=F];
        if(match.state[i]==0){
            hmdb.hit <- NA;
            kegg.hit <- NA;
        }else{
            hmdb.hit <- ifelse(nchar(hit.values[i])==0, NA, hit.values[i]);
            kegg.hit <- ifelse(nchar(hit$kegg_id)==0, NA, hit$kegg_id);
        }
        nm.mat[i, ]<-c(qvec[i],hmdb.hit, kegg.hit);
     }
     return(as.data.frame(nm.mat));
}

GetMapTable<-function(){
    require(xtable);
    print(xtable(dataSet$map.table, caption="Result from Compound Name Mapping"),
          tabular.environment = "longtable", caption.placement="top", size="\\scriptsize");
}

PathMapping<-function(qvec){
       qvec <- strsplit(qvec, "; *")[[1]];

       cmpd.db <- readRDS("../../libs/compound_db.rds");
       # load to local environment to save memory
       syn.db <- readRDS("../../libs/syn_nms.rds");

       if(!exists('path.list')){
        	LoadSmpLib();
       }

       # first find exact match to the common compound names
       hit.inx <- match(tolower(qvec), tolower(cmpd.db$name));

       match.values <- cmpd.db$name[hit.inx];
       match.ids <- cmpd.db$hmdb[hit.inx];

       # then try to find exact match to synanyms for the remaining unmatched query names one by one
       todo.inx <-which(is.na(hit.inx));
       if(length(todo.inx) > 0){
            syns.list <-  syn.db$syns.list;
            for(i in 1:length(syns.list)){
                syns <-  syns.list[[i]];
                hitInx <- match(tolower(qvec[todo.inx]), tolower(syns));

                hitPos <- which(!is.na(hitInx));
                if(length(hitPos)>0){
                    # record matched ones
                    orig.inx<-todo.inx[hitPos];
                    match.values[orig.inx] <- cmpd.db$name[i];    # show common name
                    match.ids[orig.inx] <- cmpd.db$hmdb[i];

                    # update unmatched list
                    todo.inx<-todo.inx[is.na(hitInx)];
                }
            }
       }

       na.inx<-is.na(match.values);
       qvec.nm <- match.values[!na.inx];
       qvec.id <- match.ids[!na.inx];

       # then used the matched.values to search the pathway database
       hits<-lapply(path.list, match, qvec.nm);

       set.size<-length(path.list);
       res.mat<-matrix(NA, nrow=set.size, ncol=4);

       path.nms<- names(path.list);

       for(i in 1:set.size){
            path <- path.list[[i]];
            m.inx <- hits[[i]];
            m.inx <- m.inx[!is.na(m.inx)];
            if(length(m.inx)>0){

                # note: need to set cmpd highlights in SMPDB using the hmdbid
                # syntax: http://pathman.smpdb.ca/pathways/SMP00055/pathway?reset=true&highlight[HMDB00243]=true
                # changed URL--> http://www.smpdb.ca/view/SMP00055  
                pathUrl <- paste(path.link[i], "?reset=true", sep="");
                highLights <- paste("highlight[", qvec.id[m.inx], "]=true", sep="", collapse="&");

                res.mat[i, 1] <- paste("<a href=", paste(pathUrl, "&", highLights, sep=""), ">",path.nms[i],"</a>", sep="");
                res.mat[i, 2] <- paste("<a href=http://www.hmdb.ca/metabolites/", qvec.id[m.inx], ">",qvec.nm[m.inx],"</a>", sep="", collapse="; ");
                res.mat[i, 3] <- length(m.inx);
                res.mat[i, 4] <- 1/length(path);
            }
        }

        res.mat<-res.mat[!is.na(res.mat[,1]),];
        ord.inx <- order(res.mat[,3], res.mat[,4], decreasing = T);
        res.mat <- res.mat[ord.inx, -c(3,4)];
        path.res <<- res.mat;

}

GetPathNames<-function(){
    return(path.res[,1]);
}

GetMatchedCompounds<-function(){
    return(path.res[,2]);
}