# === GENERAL OPTIONS ===

options(stringsAsFactors = FALSE)
if(Sys.getenv('SHINY_PORT') == "") options(shiny.maxRequestSize=10000*1024^2)

# === LOAD LIBRARIES ===

library(ggplot2)
library(DT)
library(DBI)
library(RSQLite)
library(gsubfn)
library(data.table)
library(parallel)
library(pbapply)
library(enviPat)
library(plotly)
library(jsonlite)
library(shinyFiles)
library(stringr)
library(colourpicker)

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")           
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

# - make ggplot? -

# === GET OPTIONS ===

wd <- "/Users/jwolthuis/Google Drive/MetaboShiny"

# --- laod adduct table for general use ---

load(file.path(wd, "backend/umcfiles/adducts/AdductTableWKZ.RData"))
sourceDir(file.path(wd, "backend/scripts/joanna"))
data(isotopes)

# --- beta stuff ---
session_cl <<- NA
mode <- "time"

# --- check options ---

opt_conn <- file(".conf")
options_raw <<- readLines(opt_conn)
print(options_raw)
dbDir <<- str_match(options_raw[[1]], "(?<=')(.*)(?=')")[1,1]
exp_dir <<- str_match(options_raw[[2]], "(?<=')(.*)(?=')")[1,1]
proj_name <<- str_match(options_raw[[3]], "(?<=')(.*)(?=')")[1,1]

patdb <<- file.path(exp_dir, paste0(proj_name, ".db"))

close(opt_conn)

# === SOURCE OWN CODE ===

sourceAll <- function(where, 
                      which=c("general", "stats", "time", "enrich_path", "power_roc", "utils")){
  library(compiler)
  print("sourcing R code ... ");
  for(i in 1:length(which)){
    script.loc <- file.path(where, which[i])
    print(script.loc)
    files <- list.files(file.path(where, which[i]),full.names=TRUE, pattern=".R$");
    print(files)
    for(f in files){
      print(f)
      source(f)
    }
  }
  return("TRUE");
}

############ BUILD DATABASES ###############

# === SOURCE METABOANALST CODE ===

# === LATER FUNCTIONS ===

get_exp_vars <- function(){
  conn <- dbConnect(RSQLite::SQLite(), patdb) # change this to proper var later
  dbGetQuery(conn, "PRAGMA table_info(setup)")$name
}

browse_db <- function(chosen.db){
  conn <- dbConnect(RSQLite::SQLite(), chosen.db) # change this to proper var later
  # --- browse ---
  result <- dbGetQuery(conn, "SELECT DISTINCT compoundname as Compound, baseformula as Formula, description as Description, charge as Charge FROM base")
  # --- result ---
  result
}

get_matches <- function(mz, chosen.db){
  # --- connect to db ---
  req("patdb")
  conn <- dbConnect(RSQLite::SQLite(), patdb) # change this to proper var later
  # 0. Attach db
  query.zero <- fn$paste("ATTACH '$chosen.db' AS db")

  dbExecute(conn, query.zero)
  query.one <- fn$paste(strwrap(
    "CREATE TEMP TABLE unfiltered AS
    SELECT cpd.baseformula, cpd.adduct
    FROM mzvals mz
    JOIN mzranges rng ON rng.ID = mz.ID
    JOIN db.extended cpd indexed by e_idx2
    ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
    AND mz.foundinmode = cpd.foundinmode
    WHERE ABS(mz.[mzmed.pgrp] - $mz) < 0.000000000001",width=10000, simplify=TRUE))
  # 1. Find matches in range (reasonably fast <3)
  dbExecute(conn, query.one)
  #  2. get isotopes for these matchies (reverse search)
  query.two <- fn$paste(strwrap(
    "CREATE TEMP TABLE isotopes AS
    SELECT cpd.baseformula, cpd.adduct, cpd.isoprevalence, cpd.basecharge 
    FROM db.extended cpd indexed by e_idx1
    JOIN unfiltered u
    ON u.baseformula = cpd.baseformula
    AND u.adduct = cpd.adduct
    JOIN mzranges rng
    ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax"
    , width=10000, simplify=TRUE))
  dbExecute(conn, query.two)
  query.three <-  strwrap(
    "SELECT DISTINCT base.compoundname as Compound, base.identifier as Identifier, iso.adduct as Adduct, base.description as Description 
    FROM isotopes iso
    JOIN db.base base indexed by b_idx1
    ON base.baseformula = iso.baseformula AND
    base.charge = iso.basecharge
    GROUP BY iso.baseformula, iso.adduct
    HAVING COUNT(iso.isoprevalence > 99.99999999999) > 0"
    , width=10000, simplify=TRUE)
  # 3. get the info you want
  results <- dbGetQuery(conn,query.three)
  # --- results ---
  results
  }


get_mzs <- function(baseformula, charge, chosen.db){
  # --- connect to db ---
  req("patdb")
  conn <- dbConnect(RSQLite::SQLite(), patdb) # change this to proper var later
  print(baseformula)
  print(charge)
  query.zero <- fn$paste("ATTACH '$chosen.db' AS db")
  print(query.zero)
  dbExecute(conn, query.zero)
  # search combo of baseformula and charge matching your choice and find all possible mzvals and adducts
  query.one <-  fn$paste(strwrap(
    "CREATE TEMP TABLE possible_options AS
    SELECT DISTINCT e.fullmz, e.adduct, e.isoprevalence
    FROM db.extended e
    WHERE e.baseformula = '$baseformula' 
    AND e.basecharge = $charge"
    , width=10000, simplify=TRUE))
  print(query.one)
  
  dbExecute(conn, query.one)
  
  # join with patdb
  query.two <- fn$paste(strwrap(
    "CREATE TEMP TABLE isotopes AS
    SELECT DISTINCT mz.[mzmed.pgrp], o.*
    FROM possible_options o
    JOIN mzranges rng
    ON o.fullmz BETWEEN rng.mzmin AND rng.mzmax
    JOIN mzvals mz
    ON rng.ID = mz.ID", width=10000, simplify=TRUE))
  print(query.two)
  
  dbExecute(conn, query.two)
  
  # isofilter and only in
  
  query.three <-  strwrap(
    "SELECT DISTINCT iso.[mzmed.pgrp], adduct
    FROM isotopes iso
    GROUP BY iso.adduct
    HAVING COUNT(iso.isoprevalence > 99.99999999999) > 0"
    , width=10000, simplify=TRUE)
  print(query.three)
  
  results <- dbGetQuery(conn,query.three)
  print(results)
}

ggplotMeba <- function(mz, draw.average, cols=c("Red", "Green")){
  cols <- if(is.null(cols)) c("Red", "Green") else(cols)
  profile <- getProfile(mz)
  plot <- if(draw.average){
    ggplot(data=profile) +
      geom_line(size=0.3, aes(x=Time, y=Abundance, group=Sample, color=Group), alpha=0.4) +
      stat_summary(fun.y="mean", size=1.5, geom="line", aes(x=Time, y=Abundance, color=Group, group=Group)) +
      scale_x_discrete(expand = c(0, 0)) +
      theme_minimal(base_size = 10) +
      scale_color_manual(values=cols)
  } else{
    ggplot(data=profile) +
      geom_line(size=0.7, aes(x=Time, y=Abundance, group=Sample, color=Group)) +
      scale_x_discrete(expand = c(0, 0)) +
      theme_minimal(base_size = 10) +
      scale_color_manual(values=cols)
  }
  # ---------------
  ggplotly(plot, tooltip="Sample")
}

ggplotSummary <- function(mz, cols){

    if(substring(dataSet$format,4,5)!="ts"){
      
      par(mar=c(4,4,2,2), mfrow = c(1,2), oma=c(0,0,2,0));
      
      mns <- by(as.numeric(dataSet$proc[, cmpdNm]), dataSet$proc.cls, mean, na.rm=T);
      sds <- by(as.numeric(dataSet$proc[, cmpdNm]), dataSet$proc.cls, sd, na.rm=T);
      
      ups <- mns + sds;
      dns <- mns - sds;
      
      # all concentration need start from 0
      y <- c(0, dns, mns, ups);
      
      rg <- range(y) + 0.05 * diff(range(y)) * c(-1, 1)
      pt <- pretty(y)
      
      axp=c(min(pt), max(pt[pt <= max(rg)]),length(pt[pt <= max(rg)]) - 1);
      
      # ymk <- pretty(c(0,ymax));
      x <- barplot(mns, col= unique(GetColorSchema()), las=2, yaxp=axp, ylim=range(pt));
      arrows(x, dns, x, ups, code=3, angle=90, length=.1);
      axis(1, at=x, col="white", col.tick="black", labels=F);
      box();
      mtext("Original Conc.", line=1);
      
      boxplot(dataSet$norm[, cmpdNm]~dataSet$cls,las=2, col= unique(GetColorSchema()));
      mtext("Normalized Conc.", line=1);
      title(main=cmpdNm, out=T);
      #
    }else if(dataSet$design.type =="time0"){
      #
      plotProfile(cmpdNm);
      #
    }else{
      if(dataSet$design.type =="time"){ # time trend within phenotype
        out.fac <- dataSet$exp.fac;
        in.fac <- dataSet$time.fac;
        xlab="Time";
      }else{ # factor a split within factor b
        out.fac <- dataSet$facB;
        in.fac <- dataSet$facA;
        xlab=dataSet$facA.lbl;
      }
      
      # two images per row
      img.num <- length(levels(out.fac));
      row.num <- ceiling(img.num/2)
      
      if(row.num == 1){
        h <- w*5/9;
      }else{
        h <- w*0.5*row.num;
      }
      #
      par(mar=c(3,4,4,2), mfrow=c(row.num, 2));
      # make sure all at the same range
      ylim.ext <-  GetExtendRange (dataSet$norm[, cmpdNm], 12);
      for(lv in levels(out.fac)){
        inx <- out.fac == lv;
        dat <- dataSet$norm[inx, cmpdNm];
        cls <- in.fac[inx];
        boxplot(dat ~ cls, col="#0000ff22", ylim=ylim.ext, outline=FALSE, boxwex=c(0.5, 0.5), xlab=xlab, ylab="Abundance", main=lv);
        stripchart(dat ~ cls, method = "jitter", ylim=ylim.ext, vertical=T, add = T, pch=19, cex=0.7);
        # ggplot
      }
      #
    }
}
