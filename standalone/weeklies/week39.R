
# CONVERT PERCENTAGES
probs_gluc_mplush = c(99.91, 6.61, 1.23, 0.08)
mzs_gluc_mplush = c(181.070631454524, 182.074067290619, 183.074947067451, 184.078200327884)

probs_gluc_mplush = c(99.97, 0.13, 1.03, 0.06)
mzs_gluc_mplush = c(181.070631454524, 182.074067290619, 183.074947067451, 184.078200327884)

cpd <- mzs_gluc_mplush[2]

x <- rmultinom(n = 1, size = 30000, prob = probs_gluc_mplush)

cpd <- 181.070631454524
ppmrange <- c(181.070631454524 - (2e-6 * 181.070631454524),
              181.070631454524 + (2e-6 * 181.070631454524))
# - find a similar weight option for the main isotope - 

conn <- RSQLite::dbConnect(RSQLite::SQLite(), "/Users/jwolthuis/Analysis/SP/SWINE_SEPT18.db")
RSQLite::dbExecute(conn, "DETACH db")
chosen.db <- "/Users/jwolthuis/Google Drive/MetaboShiny/backend/db/dimedb.full.db"
RSQLite::dbExecute(conn, gsubfn::fn$paste("ATTACH '$chosen.db' AS db"))

query <- gsubfn::fn$paste(strwrap(
  "SELECT DISTINCT cpd.baseformula as baseformula, 
  cpd.fullformula, 
  cpd.adduct as adduct, 
  cpd.isoprevalence as isoprevalence, 
  cpd.basecharge as charge, 
  (ABS($cpd - cpd.fullmz) / $cpd)/1e6 AS dppm
  FROM mzvals mz
  JOIN mzranges rng ON rng.ID = mz.ID
  JOIN db.extended cpd indexed by e_idx2
  ON cpd.fullmz BETWEEN rng.mzmin AND rng.mzmax
  AND mz.foundinmode = cpd.foundinmode
  WHERE $cpd BETWEEN rng.mzmin AND rng.mzmax",width=10000, simplify=TRUE))
RSQLite::dbGetQuery(conn, query)

# - - adj - - 

min = ppmrange[1]
max = ppmrange[2]

query <- gsubfn::fn$paste(strwrap(
  "SELECT DISTINCT cpd.baseformula as baseformula, 
  cpd.fullformula, 
  cpd.adduct as adduct, 
  cpd.isoprevalence as isoprevalence, 
  cpd.basecharge as charge, 
  (ABS($cpd - cpd.fullmz) / $cpd)/1e6 AS dppm
  FROM db.extended cpd indexed by e_idx2
  WHERE cpd.fullmz BETWEEN $min AND $max",width=10000, simplify=TRUE))
RSQLite::dbGetQuery(conn, query)

gluc_isos <- RSQLite::dbGetQuery(conn, "SELECT * FROM extended WHERE baseformula = 'C6H12O6' AND adduct = 'M+H' AND isoprevalence > 99.9999999999999")
gluc_isos$fullmz[1] %between% ppmrange

# find similar cpds
all_mplush_100isos <- RSQLite::dbGetQuery(conn, "SELECT * FROM extended WHERE adduct = 'M+H' AND isoprevalence > 99.9999999999999 ORDER BY fullmz ASC LIMIT 500")
#TODO: use this to find slide example... 
tail(all_mplush_100isos, n=300)

# FIND THREE PAIRS TO DO SCORING ON
# 1: 124.991349592431 (ethylene(and others with 4.3e-13 ppmdiff), sulfacetamide)
# 2: 192.074283985401 (1e-13, 2e-13, 3e-13)

cpd = 124.991349592431

cpd = 192.074283985401
ppm = 2
ppmrange <- c(cpd - (ppm * 1e-6 * cpd),
              cpd + (ppm * 1e-6 * cpd))

matchies <- get_matches(cpd, chosen.db, inshiny=F)
min = ppmrange[1]
max = ppmrange[2]

conn <- RSQLite::dbConnect(RSQLite::SQLite(), "/Users/jwolthuis/Analysis/SP/CHICKENS_SEPT18.db")
RSQLite::dbExecute(conn, "DETACH db")
chosen.db <- "/Users/jwolthuis/Google Drive/MetaboShiny/backend/db/chebi.full.db"
RSQLite::dbExecute(conn, gsubfn::fn$paste("ATTACH '$chosen.db' AS db"))

RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS unfiltered")
query <- gsubfn::fn$paste(strwrap(
  "SELECT DISTINCT cpd.baseformula as baseformula, 
  cpd.fullformula, 
  cpd.adduct as adduct, 
  cpd.isoprevalence as isoprevalence, 
  cpd.basecharge as charge, 
  (ABS($cpd - cpd.fullmz) / $cpd)/1e6 AS dppm
  FROM db.extended cpd indexed by e_idx2
  WHERE cpd.fullmz BETWEEN $min AND $max
  AND cpd.isoprevalence > 99.99999999",width=10000, simplify=TRUE))
res <- RSQLite::dbGetQuery(conn, query)

DBI::dbWriteTable(conn, "unfiltered", res, overwrite = T)

RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS isotopes")
RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS adducts")

RSQLite::dbExecute(conn,
                   "CREATE TABLE isotopes AS
                   SELECT DISTINCT cpd.*
                   FROM db.extended cpd
                   JOIN unfiltered u
                   ON u.baseformula = cpd.baseformula
                   AND u.adduct = cpd.adduct")

DBI::dbExecute(conn,"CREATE TABLE adducts AS
               SELECT DISTINCT *
               FROM db.extended cpd
               JOIN unfiltered u
               ON u.baseformula = cpd.baseformula")

isos <- score.isos("/Users/jwolthuis/Analysis/SP/CHICKENS_SEPT18.db", method = "mscore", inshiny = F)

targets <- c("C8H9N2O3S1", "C2H6O4S1")

septable <- adducts[baseformula == targets[1]]


# get isotope percentages
cpds <- res$fullformula

res <- lapply(cpds, function(cpd){
  query <- gsubfn::fn$paste(strwrap(
    "SELECT DISTINCT fullmz, isoprevalence
    FROM db.extended cpd
    WHERE cpd.fullformula = '$cpd'
    ORDER BY isoprevalence DESC",width=10000, simplify=TRUE))
  res <- RSQLite::dbGetQuery(conn, query)
  res$cpd <- cpd
  # - - - -
  res
})

cpd_tab <- p.cpd$C7H11N3O2Na1

# do next steps
samp_tab <- split.by.samp$`BR10-8`

#convert to perc


theor_mat <- samp_tab[,c("fullmz", "isoprevalence")]
theor <- matrix(ncol = nrow(theor_mat), nrow = 3, data = c(theor_mat$fullmz, theor_mat$isoprevalence),byrow = T)

obs_mat <- samp_tab[,c("mzmed", "intensity")]
obs <- matrix(ncol = nrow(obs_mat), nrow = 3, data = c(obs_mat$mzmed, obs_mat$intensity),byrow = T)

theor[2,] <- theor[2,]/sum(theor[2,]) * 100
obs[2,] <- obs[2,]/sum(obs[2,]) * 100

theor[3,] <- c("C7H11N3O2Na1")
obs[3,] <- c("sample")

isotable = as.data.table(rbind(t(theor), t(obs)))
# - - - 

obs.dt <- as.data.table(t(obs))
colnames(obs.dt) <- c("m/z", "isotope", "label")

isotable <- rbindlist(res)
isotable <- rbindlist(list(isotable, obs.dt))
colnames(isotable) <- c("m/z", "isotope", "label")
isotable$`m/z` <- factor(isotable$`m/z`)
isotable$isotope <- as.numeric(isotable$isotope)

ggplot(data=isotable, aes(x=`m/z`, y= isotope, color = label, fill=label, width=.2)) +
  geom_bar(stat="identity") +
  geom_line(aes(x=`m/z`, y=isotope, group=formula)) +
  theme(axis.text=element_text(size=10),
        aspect.ratio = .2/1,
        axis.title=element_text(size=19,face="bold"),
        legend.title=element_text(size=15, face="bold"),
        legend.text=element_text(size=12),
        text=element_text(size=16,  family="Trebuchet MS", hjust = 0.5),
        axis.text.x = element_text(angle=90, hjust=1, size=16)) +
  facet_grid(rows = vars(label),scales = "fixed")


## - - - child rows in table - - -

resu <- get_matches(curr_cpd, chosen.db = chosen.db, inshiny=F)

child.groups <- split(resu, list(res$baseformula))
library(DT)

as.character(datatable(resu))
datatable(
  cbind(' ' = '&oplus;', res), escape = -2,
  options = list(
    columnDefs = list(
      list(visible = FALSE, targets = c(2,5,6,7,8,9)),
      list(orderable = FALSE, className = 'details-control', targets = 1)
    )
  ),
  callback = JS("
                table.column(1).nodes().to$().css({cursor: 'pointer'});
                var format = function(d) {
                return '<div style=\"background-color:#eee; padding: .5em;\"> Model: ' +
                d[0] + ', mpg: ' + d[2] + ', cyl: ' + d[3] + '</div>';
                };
                table.on('click', 'td.details-control', function() {
                var td = $(this), row = table.row(td.closest('tr'));
                if (row.child.isShown()) {
                row.child.hide();
                td.html('&oplus;');
                } else {
                row.child(format(row.data())).show();
                td.html('&CircleMinus;');
                }
                });"
))
datatable(
  cbind(' ' = '&oplus;', mtcars), escape = -2,
  options = list(
    columnDefs = list(
      list(visible = FALSE, targets = c(0, 2, 3)),
      list(orderable = FALSE, className = 'details-control', targets = 1)
    )
  ),
  callback = JS("
                table.column(1).nodes().to$().css({cursor: 'pointer'});
                var format = function(d) {
                return '<div style=\"background-color:#eee; padding: .5em;\"> Model: ' +
                d[0] + ', mpg: ' + d[2] + ', cyl: ' + d[3] + '</div>';
                };
                table.on('click', 'td.details-control', function() {
                var td = $(this), row = table.row(td.closest('tr'));
                if (row.child.isShown()) {
                row.child.hide();
                td.html('&oplus;');
                } else {
                row.child(format(row.data())).show();
                td.html('&CircleMinus;');
                }
                });"
))

# - - - subgraphs for network analysis - - - 


{
  matches <- pbapply::pblapply(1:nrow(current), function(i){
    
    mz <- current$fullmz[i]
    
    #print(paste("---", mz, "---"))
    bounds <- c(mz - mz*(ppm*1e-6), mz + mz*(ppm*1e-6))
    matches <- V(g)$name[which(V(g)$name %between% bounds)]
    name <- paste0(current$baseformula[i], "\n", "[", current$adduct[i], "]", "\n", round(as.numeric(current$isoprevalence)[i],digits = 2), "%")
    name <- rep(name, times = length(matches))
    #print(name)
    #name <- paste0(current[i, c("baseformula", "adduct", "isoprevalence")], collapse = "\n")
    # - - - -
    list(name = name, match = matches)
  })
  
  mzs <- unlist(sapply(matches, function(x) if(length(x$match) > 0) x$match))
  names <- unlist(sapply(matches, function(x) if(length(x$match) > 0) x$name))
  
  subg <- induced_subgraph(g, mzs)
  
  V(subg)$color <- c("black")
  V(subg)[V(subg)$name %in% unique(as.character(ends(subg, es=E(subg))))]$color <- c("red")
  
  V(subg)$name <- names
  V(subg)$label.color <- c("black")
  
  la <- layout.circle(subg)
  n = length(names)
  
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  
  lab.locs <- radian.rescale(x=1:n, direction=-1, start=0)
  
  plot(subg,
       edge.arrow.size = .5,
       vertex.size = 5,
       vertex.label.degree = lab.locs,
       vertex.label.cex = 1.5,
       vertex.label.dist = 2.2,
       layout = la,
       edge.color="black",
       edge.width=E(subg)$weight)
  
  # - - - - - - - - - - - -
  
  print(mean(strength(subg, mode = "all")))
  print(mean_distance(subg, directed=FALSE))
  print(edge_density(subg))
}

current <- split.by.formula[[4]]

test <- cluster_walktrap(g)

# - - netClass package? - - -

# http://pablobarbera.com/big-data-upf/html/02b-networks-descriptive-analysis.html

# - - find internal standards in internal db - -
assumed_add <- c("M+H", "M-H", "M+Na", "M-Cl")

conn <- RSQLite::dbConnect(RSQLite::SQLite(), "/Users/jwolthuis/Analysis/SP/CHICKENS_SEPT18.db")
RSQLite::dbExecute(conn, "DETACH db")
chosen.db <- "/Users/jwolthuis/Google Drive/MetaboShiny/backend/db/internal.full.db"
RSQLite::dbExecute(conn, gsubfn::fn$paste("ATTACH '$chosen.db' AS db"))

internal.standards <- RSQLite::dbGetQuery(conn,
                                          "SELECT * FROM db.base base
                                          JOIN db.extended cpd
                                          ON base.baseformula = cpd.baseformula
                                          AND base.charge = cpd.basecharge
                                          WHERE base.compoundname LIKE '%IS%'")
RSQLite::dbExecute(conn, "DROP TABLE IF EXISTS unfiltered")
query <- gsubfn::fn$paste(strwrap(
  "SELECT DISTINCT cpd.baseformula as baseformula, 
  cpd.fullformula, 
  cpd.adduct as adduct, 
  cpd.isoprevalence as isoprevalence, 
  cpd.basecharge as charge
  FROM 
  JOIN db.base base
  ON cpd.baseformula = base.baseformula
  AND cpd.basecharge = base.charge
  WHERE base.compoundname LIKE '%IS%')",width=10000, simplify=TRUE))
res <- RSQLite::dbGetQuery(conn, query)


