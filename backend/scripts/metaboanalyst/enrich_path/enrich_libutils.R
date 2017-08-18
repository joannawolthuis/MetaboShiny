##################################################
## R script for MetaboAnalyst
## Description: Library manage
##
## Author: Jeff Xia, jeff.xia@mcgill.ca
## McGill University, Canada
## License: GNU GPL (>= 2)
###################################################


###############################
## Metabolite set library
###############################

LoadMsetLib<-function(libname="pathway"){
    if(is.null(current.msetlib) || current.msetlib$lib.name!=libname){
        libPath <- paste("../../libs/msets/", libname, ".rda", sep="");
        load(libPath, .GlobalEnv);
    }
}

SetCachexiaSetUsed <- function(used){
    cachexia.set.used <<- used;
}

# set current user selected metset library for search
# if enrichment analysis, also prepare lib by
# creating to a list of metabolite sets
# exlude metabolite sets with compounds < excludeNum
# or exlude the metabolite set derived from test data
SetCurrentMsetLib<-function(lib.type, excludeNum=0){

    if(lib.type!="self"){
        LoadMsetLib(lib.type);
    }
    if(lib.type=="self"){
        ms.list <- user.mset;
    }else{
        # create a named list, use the ids for list names
        ms.list<-strsplit(current.msetlib[,3],"; ");
        names(ms.list)<-current.msetlib[,2];
    }

    if(excludeNum > 0){
        cmpd.count<-lapply(ms.list, length);
        sel.inx <- cmpd.count >= excludeNum;
        ms.list <- ms.list[sel.inx];
    }

    # note exclude the cachaxia signature if test data used
    if(cachexia.set.used){
        sel.inx <- which(names(ms.list) == "CANCER CACHEXIA PATIENTS");
        ms.list[sel.inx]<-NULL;
    }

    current.mset <<- ms.list;
    # total uniq cmpds in the mset lib
    uniq.count <<- length(unique(unlist(current.mset, use.names = FALSE)));
}

# methods to return the selected metset library to java for display
GetMsetNames<-function(){
    return(current.msetlib$name);
}

GetMsetMembers<-function(){
    return(current.msetlib$member);
}

GetMsetReferences<-function(){
    return(current.msetlib$reference);
}

# read user upload metabolite set library file
# two col csv file, first name, second cmpd list
Setup.UserMsetLibData<-function(filePath){
    dat <- .readDataTable(filePath);
    libCheck.msg <- NULL;
    if(class(dat) == "try-error") {
            libCheck.msg <-c(libCheck.msg, "Data format error - fail to read in the data!");
            AddErrMsg(libCheck.msg);
            return(0);
    }

    if(is.null(dim(dat)) || dim(dat)[2]!=2){
            libCheck.msg <-c(libCheck.msg, "Data format error - must have two columns!");
            AddErrMsg(libCheck.msg);
            return(0);
    }

    # create a named list, use the ids for list names
    mset.list<-strsplit(dat[,2],"; ");
    mset.ids <- paste("USER", sprintf("%04d",1:nrow(dat)), sep="");
    names(mset.list)<-dat[,1];
    names(mset.ids)<-dat[,1];

    cmpd.db <- readRDS("../../libs/compound_db.rds");

    # now need to check all metabolites match HMDB names
    # and get the statistics
    unmatched.nms <- NULL;
    unmatched.num <- 0;
    hmdb.nms <- tolower(cmpd.db$name);
    for(i in 1:length(mset.list)){
        mset <- mset.list[[i]];
        hit.inx <- match(tolower(mset), hmdb.nms);
        unmatched.nms <- c(unmatched.nms, mset[is.na(hit.inx)]);
        unmatched.num <- unmatched.num + sum(is.na(hit.inx));
    }

    # save the lib data
    user.mset <<- mset.list;
    user.mset.ids <<- mset.ids;

    if(unmatched.num > 0) {
        user.mset.info <<- paste("A total of", unmatched.num, "compounds were no matched to HMDB common names.",
                    "They are:", paste(unmatched.nms, collapse="; "), ". Please correct these names. Otherwise,",
                    "they will be ignored during the enrichment analysis.");
    }else{
        user.mset.info <<- paste("A total of", length(mset.list), "were sucessfully added to the library.");
    }
    return(1);
}

GetMsetLibCheckMsg<-function(){
    return (user.mset.info);
}

Get.ConcRef<-function(cmpd.nm){
    if(is.null(conc.db)){
        conc.db <<-  .readDataTable("../../libs/cmpd_conc.csv");
    }
    matches <- subset(conc.db, name == cmpd.nm & biotype==dataSet$biofluid, select=c(conc, pubmed, references, notes));
    if(nrow(matches)==0){
        return(NA);
    }
    return(list(concs = matches$conc, pmid = matches$pubmed, refs = matches$references, note = matches$notes));
}


# load pathway library
LoadSmpLib<-function(){
    paths <- .readDataTable("../../libs/smp_path.csv");
    path.list<-strsplit(paths[,2],"; ");
    names(path.list)<-paths[,1];
    path.list <<- path.list;
    path.link <<- paths[,3];
}

SearchMsetLibraries<-function(query, type){
       if(!exists("lib.search")){
            lib.search <<- list();
       }
       
       query <- ClearStrings(query);

       if(nchar(query)==0){
            return();
       }

       if(type=="name"){
            SearchByName(query);
       }else{
            SearchByCompound(query);
       }
}

# search compound in all member compounds of metabolite set
SearchByCompound <- function(query){

       names.vec <- current.msetlib$member;
       matched.inx <- NULL;
       matched <- FALSE;
       exact = FALSE;

       # matching from close match to more dist match (max.dist = 0.5)
       # once a closer match found, stop trying more distant one
       matched.dist <- NULL;
       s <- seq(0, 0.2, .1)
       for (i in s) {
            matched.inx <- agrep(query,names.vec,ignore.case=T, max.distance=i);
            if(length(matched.inx) > 0) {
                matched.dist <- i;
                matched <- TRUE;
                break;
            }
        }

        if(matched){
            # now break down the set into each individual metabolites and find out which one gives the best hit
            matched.list<- vector(mode = "list", length=length(matched.inx));
            for(i in 1:length(matched.inx)){
                matched.list[i]<-strsplit(current.msetlib[matched.inx[i], "member"],"; *");
            }

            # re-do the matching, and record the matched values & sort
            matched.score <- NULL;
            hit.value <- vector(mode = "character", length=length(matched.inx)); # save the exact hit
            matched.value <- vector(mode = "character", length=length(matched.inx)); # save the whole metset
            for (i in 1:length(matched.inx)) {
                matched.nm <- matched.list[[i]];
                # test if it is exact match
                if((matched.dist == 0.0) & (!is.na(hit.inx <- match(tolower(query), tolower(matched.nm))))){
                    matched.score[i] <- -1.0;
                    exact <- TRUE;
                }else{ # try approximate match, note: we only record the first match in each set
                    hit.inx <- agrep(query,matched.nm,ignore.case=T, max.distance=matched.dist)[1];
                    # matched.dist 0.0, 0.1, 0.2, with the matches of the same distance, add fine adjustment
                    # based on the length b/w query and matched name
                    # use query length for normalization
                    matched.score[i] <- matched.dist + abs(nchar(matched.nm[hit.inx])-nchar(query))/(1000*nchar(query));
                }

                # wrap up hit metabolite sets in html tags
                html.tag <- "<p>";
                for(m in 1:length(matched.list[[i]])){
                   current.cmpd <- matched.list[[i]][m];
                   if(m == hit.inx){
                        current.cmpd <- paste("<font color=\"red\">", "<b>", current.cmpd, "</b>", "</font>",sep="");
                   }
                   if(m == 1){
                        html.tag <- paste(html.tag, current.cmpd, sep="");
                   }else {
                        html.tag <- paste(html.tag, "; ", current.cmpd, sep="");
                   }
                }
                hit.value[i] <-matched.list[[i]][hit.inx] ;
                matched.value[i] <- paste(html.tag, "</p>");
            }

            matched.table <- cbind(current.msetlib$name[matched.inx],
                                   matched.value,
                                   current.msetlib$reference[matched.inx]);
            if(exact){
                exact.inx <- matched.score == -1;
                lib.search$matched.table <- matched.table[exact.inx, ];
                lib.search$best.hit <- "NA";
            }else{
                # rank results based on the matched scores
                ord.inx <- order (matched.score, decreasing=F);
                lib.search$matched.table <- matched.table[ord.inx, ];
                lib.search$best.hit <- hit.value[ord.inx][1];
            }
        }else{
            lib.search$best.hit <- "NA";
        }
        lib.search<<-lib.search;
}

getBestHit<-function(){
    return(lib.search$best.hit);
}

# since String[][] is not supported, have to return as 1D vector, 
# matrix can be directly convert to vector, note default will be column first
GetMsetLibSearchResult<-function(){
    return (as.vector(lib.search$matched.table));
}

# given a metabolite set name, search its index
SearchByName <- function(query){
       # no need for suggestions for metabolite set name search
       lib.search$best.hit <- "NA";
       names.vec <- current.msetlib$name;
       matched <- FALSE;

       # matching from exact match (max.dist = 0) to more dist match (max.dist = 0.5)
       # once a closer match found, stop trying more distant one
       matched.inx <- match(tolower(query), tolower(names.vec));
       if(is.na(matched.inx)){ # try approximate match
            s <- seq(0, 0.2, .1)
            for (i in s) {
                matched.inx <- agrep(query,names.vec,ignore.case=T, max.distance=i);
                if(length(matched.inx) > 0) {
                    matched = TRUE;
                    break;
                }
            }
        }else{
            matched = TRUE;
        }

        if(matched){
            # wrap up in html tags
            matched.names <- paste("<p><font color=\"red\"><b>", names.vec[matched.inx], "</b></font></p>",sep="");
            lib.search$matched.table <- cbind(matched.names, current.msetlib$member[matched.inx], current.msetlib$reference[matched.inx]);
        }else{
            lib.search$matched.table <-"NA";
        }
}

GetSMPDBimg<-function(msetInx){
    msetNm <- GetMetSetName(msetInx);
    inx <- which(current.msetlib$name == msetNm);
    return(current.msetlib$image[inx]);
}

# note, this process can be long, need to return a value
# to force Java to wait
SetKEGG.PathLib<-function(kegg.rda){

       dataSet$lib.msg <- paste("Your selected pathway library code is \\textbf{", kegg.rda, "}(KEGG organisms abbreviation).");
       path.dir <-paste ("../../libs/kegg/", kegg.rda, ".rda", sep="");

       dataSet <<- dataSet;
       load(path.dir, .GlobalEnv); # this one works for 32/64 bit, the above one not work in 64 bit?
       return(1);
}


# read user upload metabolome as a list of kegg ids
Setup.KEGGReferenceMetabolome<-function(filePath){
    inFile <- file(filePath, "r");
    ref.vec<-try(scan(inFile, 'character', strip.white = T, sep="\n")); # must be single column
    close(inFile);
    libCheck.msg <- NULL;
    if(class(ref.vec) == "try-error") {
            libCheck.msg <-c(libCheck.msg, "Data format error - fail to read in the data!");
            print(libCheck.msg);
            AddErrMsg(libCheck.msg);
            return(0);
    }

    cmpd.db <- readRDS("../../libs/compound_db.rds");

    # now need to check all metabolites match HMDB names
    # and get the statistics
    hits <- tolower(ref.vec)%in%tolower(cmpd.db$kegg_id);
    metabo.filter.kegg <<- ref.vec[hits];
    unmatched.num <- sum(!hits);
    if(unmatched.num > 0) {
        unmatched.nms <- ref.vec[!hits];
        metabo.ref.info <<- paste("A total of", unmatched.num, "compounds were no matched to KEGG compound IDs.",
                    "They are:", paste(unmatched.nms, collapse="; "), ". Please correct these names. Otherwise,",
                    "they will be ignored during the enrichment analysis.");
    }else{
        metabo.ref.info <<- paste("A total of", length(ref.vec), "were sucessfully added to the library.");
    }
    return(1);
}

# read user upload metabolome as a list of hmdb cmpd names
Setup.HMDBReferenceMetabolome<-function(filePath){
    inFile <- file(filePath, "r");
    ref.vec<-try(scan(inFile, 'character', strip.white = T, sep="\n")); # must be single column
    close(inFile);
    libCheck.msg <- NULL;
    if(class(ref.vec) == "try-error") {
            libCheck.msg <-c(libCheck.msg, "Data format error - fail to read in the data!");
            AddErrMsg(libCheck.msg);
            return(0);
    }

    cmpd.db <- readRDS("../../libs/compound_db.rds");

    # now need to check all metabolites match HMDB names
    # and get the statistics
    hits <- tolower(ref.vec)%in%tolower(cmpd.db$name);
    metabo.filter.hmdb <<- ref.vec[hits];
    unmatched.num <- sum(!hits);
    if(unmatched.num > 0) {
        unmatched.nms <- ref.vec[!hits];
        metabo.ref.info <<- paste("A total of", unmatched.num, "compounds were no matched to HMDB compound names.",
                    "They are:", paste(unmatched.nms, collapse="; "), ". Please correct these names. Otherwise,",
                    "they will be ignored during the enrichment analysis.");
    }else{
        metabo.ref.info <<- paste("A total of", length(ref.vec), "were sucessfully added to the library.");
    }
    return(1);
}

GetRefLibCheckMsg<-function(){
    return (metabo.ref.info);
}

SetMetabolomeFilter<-function(TorF){
    use.metabo.filter <<- TorF;
}