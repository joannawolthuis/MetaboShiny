trim <- function (x) gsub("^\\s+|\\s+$", "", x)

trim.trailing <- function (x) sub("\\s+$", "", x)

#' @export
list.to.df <- function(lst,
                       type){
  dfStorage <- pbapply::pblapply(seq_along(lst), FUN=function(y, n, i){
    matches <- as.character(y[[i]])
    if((length(matches) == 0) && (typeof(matches) == "character")){ matches = NA }
    temp.df <- data.frame(stringsAsFactors = FALSE, "mz"=rep(n[[i]], length(matches)), "CompoundName"=matches[[1]], "Adduct"=matches[[2]], "Isotope"=matches[[3]], "Source"=rep(type, length(matches)))
    temp.df
  },y=lst, n=names(lst))
  df <- rbindlist(dfStorage)
  # --- return ---
  df
}

#' @export
factorize <- function(vector){
  items <- as.factor(vector)
  new.levels <- c(0:(length(levels(items)) - 1 ))
  levels(items) <- as.numeric(new.levels)
  # --- return ---
  items
}

#' @export
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

#' @export
ppm_range <- function(mz, ppm) c((mz - (ppm/1000000 * mz)), (mz + (ppm/1000000 * mz)))

# --- isotopey testing stuff ---

# lactate <- "C3H5O3"
# glucose <- "C6H12O6"
#
# glucose.2H <- mergeform.joanna(glucose, "H2")
#
# lact.isos <- isopattern(isotopes, lactate, charge = -1)
# gluc.isos <- isopattern(isotopes, glucose.2H, charge = 2)
#
# plot.isos <- function(isopat.table){
#   plot(isopat.table[, 1], isopat.table[, 2], type = "h", xlab = "m/z",
#        ylab = "Relative abundance")
# }
#
# plot.isos(lact.isos[[1]])
# plot.isos(gluc.isos[[1]])
#
# # --- testy stuff ---
#
# plot(lact.isos[[1]][, 1], lact.isos[[1]][, 2], type = "h", xlab = "m/z",
#      ylab = "Relative abundance", col='blue')
# lines(gluc.isos[[1]][, 1], gluc.isos[[1]][, 2], type = "h", xlab = "m/z",
#       ylab = "Relative abundance", col='red')
#
# lact.isos[[1]]

#' @export
download.chebi.joanna <- function (release = "latest", woAssociations = FALSE) {
  chebi_download <- tempdir()
  utils::download.file("ftp://ftp.ebi.ac.uk/pub/databases/chebi/archive/",
                paste0(chebi_download, "releases.txt"), quiet = TRUE)
  releases <- gsub("rel", "", read.table(paste0(chebi_download,
                                                "releases.txt"), quote = "\"", comment.char = "")[,
                                                                                                  9])
  message("Validating ChEBI release number ... ", appendLF = FALSE)
  if (release == "latest") {
    release <- max(releases)
  }
  else {
    release <- releases[match(release, releases)]
  }
  message("OK")
  ftp <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/chebi/archive/rel",
                release, "/Flat_file_tab_delimited/")
  message("Downloading compounds ... ", appendLF = FALSE)
  utils::download.file(paste0(ftp, "compounds.tsv.gz"), paste0(chebi_download,
                                                        "compounds.tsv"), quiet = TRUE)
  compounds <- as.data.frame.array(read.delim2(paste0(chebi_download,
                                                      "compounds.tsv")))
  message("DONE", appendLF = TRUE)
  message("Downloading synonyms ... ", appendLF = FALSE)
  utils::download.file(paste0(ftp, "names.tsv.gz"), paste0(chebi_download,
                                                    "names.tsv"), quiet = TRUE)
  names <- suppressWarnings(as.data.frame.array(read.delim2(paste0(chebi_download,
                                                                   "names.tsv"))))
  message("DONE", appendLF = TRUE)
  message("Downloading formulas ... ", appendLF = FALSE)
  utils::download.file(paste0(ftp, "chemical_data.tsv"), paste0(chebi_download,
                                                         "formulas.tsv"), quiet = TRUE)
  formulas <- suppressWarnings(as.data.frame.array(read.delim2(paste0(chebi_download,
                                                                      "formulas.tsv"))))

  message("Downloading structures ... ", appendLF = FALSE)

  utils::download.file(paste0(ftp, "structures.csv.gz"), paste0(chebi_download,
                                                                "structures.csv"), quiet = TRUE)

  message("DONE", appendLF = TRUE)

  structures <- suppressWarnings(as.data.frame.array(read.delim2(paste0(chebi_download,
                                                                        "structures.csv"), sep=",")))
  # mol.rows <- which(structures$TYPE == "mol")
  # inchi.rows <- which(structures$TYPE == "InChI")
  # inchikey.rows <- which(structures$TYPE == "InChIKey")

  smile.rows <- which(structures$TYPE == "SMILES")

  structures <- structures[smile.rows,]

  message("DONE", appendLF = TRUE)
  message("Building ChEBI ... ", appendLF = TRUE)

  #compounds <- compounds[compounds[, "STAR"] >= 3, ]
  latest <- compounds[, c("ID", "NAME", "DEFINITION")]
  old <- compounds[, c("ID", "PARENT_ID")]
  old <- merge(x = old, y = latest, by.x = "PARENT_ID", by.y = "ID")
  compounds <- rbind(latest, old[, c("ID", "NAME", "DEFINITION")])
  compounds[compounds[, "NAME"] == "null", "NAME"] <- NA
  compounds <- compounds[complete.cases(compounds), ]
  DB <- suppressWarnings((merge(compounds[, c("ID", "NAME", "DEFINITION")],
                                names[, c("COMPOUND_ID", "SOURCE", "NAME")], by.x = "ID",
                                by.y = "COMPOUND_ID", all.x = TRUE)))
  ChEBI <- unique(DB[, c("ID", "NAME.x", "DEFINITION")])
  colnames(ChEBI) <- c("ID", "ChEBI", "DEFINITION")

  message(" KEGG Associations ... ", appendLF = FALSE)
  KEGG <- unique(DB[DB[, "SOURCE"] == "KEGG COMPOUND", c("ID",
                                                         "NAME.y")])
  KEGG <- KEGG[complete.cases(KEGG), ]
  colnames(KEGG) <- c("ID", "KEGG")
  message("DONE", appendLF = TRUE)
  message(" IUPAC Associations ... ", appendLF = FALSE)
  IUPAC <- unique(DB[DB[, "SOURCE"] == "IUPAC", c("ID", "NAME.y")])
  IUPAC <- IUPAC[complete.cases(IUPAC), ]
  colnames(IUPAC) <- c("ID", "IUPAC")
  message("DONE", appendLF = TRUE)
  message(" MetaCyc Associations ... ", appendLF = FALSE)
  MetaCyc <- unique(DB[DB[, "SOURCE"] == "MetaCyc", c("ID",
                                                      "NAME.y")])
  MetaCyc <- MetaCyc[complete.cases(MetaCyc), ]
  colnames(MetaCyc) <- c("ID", "MetaCyc")
  message("DONE", appendLF = TRUE)
  message(" ChEMBL Associations ... ", appendLF = FALSE)
  ChEMBL <- unique(DB[DB[, "SOURCE"] == "ChEMBL", c("ID",
                                                    "NAME.y")])
  ChEMBL <- ChEMBL[complete.cases(ChEMBL), ]
  colnames(ChEMBL) <- c("ID", "ChEMBL")
  message("DONE", appendLF = TRUE)
  DB <- unique(merge(DB["ID"], ChEBI, by = "ID", all.x = TRUE))
  DB <- unique(merge(DB, KEGG, by = "ID", all.x = TRUE))
  DB <- unique(merge(DB, IUPAC, by = "ID", all.x = TRUE))
  DB <- unique(merge(DB, MetaCyc, by = "ID", all.x = TRUE))
  DB <- unique(merge(DB, ChEMBL, by = "ID", all.x = TRUE))
  rm(ChEBI, ChEMBL, compounds, IUPAC, KEGG, latest, MetaCyc,
     names, old)
  if ("FORMULA" %in% unique(formulas[, "TYPE"])) {
    message(" Formula Associations ... ", appendLF = FALSE)
    formula <- formulas[formulas[, "TYPE"] == "FORMULA",
                        c("COMPOUND_ID", "CHEMICAL_DATA")]
    colnames(formula) <- c("ID", "FORMULA")
    DB <- merge(DB, formula, by = "ID", all.x = TRUE)
    DB <- merge(DB, DB[, c("ChEBI", "FORMULA")], by = "ChEBI",
                all.x = TRUE)
    DB[is.na(DB[, "FORMULA.x"]), "FORMULA.x"] <- "null"
    DB[is.na(DB[, "FORMULA.y"]), "FORMULA.y"] <- "null"
    DB[DB[, "FORMULA.x"] != "null" & DB[, "FORMULA.y"] ==
         "null", "FORMULA.y"] <- DB[DB[, "FORMULA.x"] !=
                                      "null" & DB[, "FORMULA.y"] == "null", "FORMULA.x"]
    DB[DB[, "FORMULA.y"] != "null" & DB[, "FORMULA.x"] ==
         "null", "FORMULA.x"] <- DB[DB[, "FORMULA.y"] !=
                                      "null" & DB[, "FORMULA.x"] == "null", "FORMULA.y"]
    DB <- unique(DB[DB[, "FORMULA.x"] != "null" & DB[, "FORMULA.y"] !=
                      "null", c("ID", "DEFINITION","ChEBI", "KEGG", "IUPAC", "MetaCyc",
                                "ChEMBL", "FORMULA.x")])
    rm(formula)
    message("DONE", appendLF = TRUE)
  }
  else {
    message("NOT AVAILABLE FOR THIS RELEASE")
  }
  message("Downloading molecular weights ... ", appendLF = FALSE)
  if ("MASS" %in% unique(formulas[, "TYPE"])) {
    mass <- formulas[formulas[, "TYPE"] == "MASS", c("COMPOUND_ID",
                                                     "CHEMICAL_DATA")]
    colnames(mass) <- c("ID", "MASS")
    DB <- merge(DB, mass, by = "ID", all.x = TRUE)
    DB <- merge(DB, DB[, c("ChEBI", "MASS")], by = "ChEBI",
                all.x = TRUE)
    DB[is.na(DB[, "MASS.x"]), "MASS.x"] <- "null"
    DB[is.na(DB[, "MASS.y"]), "MASS.y"] <- "null"
    DB[DB[, "MASS.x"] != "null" & DB[, "MASS.y"] == "null",
       "MASS.y"] <- DB[DB[, "MASS.x"] != "null" & DB[,
                                                     "MASS.y"] == "null", "MASS.x"]
    DB[DB[, "MASS.y"] != "null" & DB[, "MASS.x"] == "null",
       "MASS.x"] <- DB[DB[, "MASS.y"] != "null" & DB[,
                                                     "MASS.x"] == "null", "MASS.y"]
    DB <- unique(DB[, c("ID", "DEFINITION", "ChEBI", "KEGG", "IUPAC",
                        "MetaCyc", "ChEMBL", "FORMULA.x", "MASS.x")])
    rm(mass)
    message("DONE", appendLF = TRUE)
  }
  else {
    message("NOT AVAILABLE FOR THIS RELEASE")
  }
  message("Downloading monoisotopic molecular weights ... ",
          appendLF = FALSE)
  if ("MONOISOTOPIC MASS" %in% unique(formulas[, "TYPE"])) {
    mmass <- formulas[formulas[, "TYPE"] == "MONOISOTOPIC MASS",
                      c("COMPOUND_ID", "CHEMICAL_DATA")]
    colnames(mmass) <- c("ID", "MONOISOTOPIC")
    DB <- merge(DB, mmass, by = "ID", all.x = TRUE)
    DB <- merge(DB, DB[, c("ChEBI", "MONOISOTOPIC")], by = "ChEBI",
                all.x = TRUE)
    DB[is.na(DB[, "MONOISOTOPIC.x"]), "MONOISOTOPIC.x"] <- "null"
    DB[is.na(DB[, "MONOISOTOPIC.y"]), "MONOISOTOPIC.y"] <- "null"
    DB[DB[, "MONOISOTOPIC.x"] != "null" & DB[, "MONOISOTOPIC.y"] ==
         "null", "MONOISOTOPIC.y"] <- DB[DB[, "MONOISOTOPIC.x"] !=
                                           "null" & DB[, "MONOISOTOPIC.y"] == "null", "MONOISOTOPIC.x"]
    DB[DB[, "MONOISOTOPIC.y"] != "null" & DB[, "MONOISOTOPIC.x"] ==
         "null", "MONOISOTOPIC.x"] <- DB[DB[, "MONOISOTOPIC.y"] !=
                                           "null" & DB[, "MONOISOTOPIC.x"] == "null", "MONOISOTOPIC.y"]
    DB <- unique(DB[, c("ID","DEFINITION", "ChEBI", "KEGG", "IUPAC",
                        "MetaCyc", "ChEMBL", "FORMULA.x", "MASS.x", "MONOISOTOPIC.x")])
    rm(mmass)
    message("DONE", appendLF = TRUE)
  }
  else {
    message("NOT AVAILABLE FOR THIS RELEASE")
  }
  message("Downloading molecular charges ... ", appendLF = FALSE)
  if ("CHARGE" %in% unique(formulas[, "TYPE"])) {
    charge <- formulas[formulas[, "TYPE"] == "CHARGE", c("COMPOUND_ID",
                                                         "CHEMICAL_DATA")]
    colnames(charge) <- c("ID", "CHARGE")
    DB <- merge(DB, charge, by = "ID", all.x = TRUE)
    DB <- merge(DB, DB[, c("ChEBI", "CHARGE")], by = "ChEBI",
                all.x = TRUE)
    DB[is.na(DB[, "CHARGE.x"]), "CHARGE.x"] <- "null"
    DB[is.na(DB[, "CHARGE.y"]), "CHARGE.y"] <- "null"
    DB[DB[, "CHARGE.x"] != "null" & DB[, "CHARGE.y"] ==
         "null", "CHARGE.y"] <- DB[DB[, "CHARGE.x"] != "null" &
                                     DB[, "CHARGE.y"] == "null", "CHARGE.x"]
    DB[DB[, "CHARGE.y"] != "null" & DB[, "CHARGE.x"] ==
         "null", "CHARGE.x"] <- DB[DB[, "CHARGE.y"] != "null" &
                                     DB[, "CHARGE.x"] == "null", "CHARGE.y"]
    DB <- unique(DB[, c("ID", "DEFINITION", "ChEBI", "KEGG", "IUPAC",
                        "MetaCyc", "ChEMBL", "FORMULA.x", "MASS.x", "MONOISOTOPIC.x",
                        "CHARGE.x")])
    message("DONE", appendLF = TRUE)
  }
  else {
    message("NOT AVAILABLE FOR THIS RELEASE")
  }
  DB[DB == "null"] <- NA

  DB = merge(DB, structures[,-1], by.x = "ID",  by.y = "COMPOUND_ID", all.x = TRUE, incomparables = "unknown")

  DB <- unique(DB[complete.cases(DB[, c("ID", "DEFINITION","ChEBI", "FORMULA.x",
                                        "MASS.x", "MONOISOTOPIC.x", "CHARGE.x", "STRUCTURE")]), ])
  colnames(DB) <- c("ID", "DEFINITION","ChEBI", "KEGG", "IUPAC", "MetaCyc",
                    "ChEMBL", "FORMULA", "MASS", "MONOISOTOPIC", "CHARGE", "STRUCTURE")

  if (woAssociations == TRUE) {
    compounds <- unique(rbind(setNames(DB[, c("ChEBI", "DEFINITION","FORMULA",
                                              "MASS", "MONOISOTOPIC", "CHARGE")], c("NAME","DEFINITION","FORMULA",
                                                                                    "MASS", "MONOISOTOPIC", "CHARGE","STRUCTURE")), setNames(DB[,
                                                                                                                                    c("KEGG","DEFINITION", "FORMULA", "MASS", "MONOISOTOPIC", "CHARGE","STRUCTURE")],
                                                                                                                                 c("NAME", "DEFINITION","FORMULA", "MASS", "MONOISOTOPIC", "CHARGE","STRUCTURE")),
                              setNames(DB[, c("IUPAC", "DEFINITION","FORMULA", "MASS", "MONOISOTOPIC",
                                              "CHARGE")], c("NAME", "DEFINITION","FORMULA", "MASS", "MONOISOTOPIC",
                                                            "CHARGE")), setNames(DB[, c("MetaCyc", "DEFINITION","FORMULA",
                                                                                        "MASS", "MONOISOTOPIC", "CHARGE","STRUCTURE")], c("NAME", "DEFINITION",
                                                                                                                              "FORMULA", "MASS", "MONOISOTOPIC", "CHARGE","STRUCTURE")),
                              setNames(DB[, c("ChEMBL","DEFINITION", "FORMULA", "MASS", "MONOISOTOPIC",
                                              "CHARGE","STRUCTURE")], c("NAME", "DEFINITION","FORMULA", "MASS", "MONOISOTOPIC",
                                                            "CHARGE","STRUCTURE"))))
    compounds <- compounds[complete.cases(compounds), ]
    return(compounds)
  }
  else {
    return(DB)
  }
}


sdfStream.joanna <- function (input, output, append = FALSE, fct, Nlines = 10000,
                              startline = 1, restartNlines = 10000, silent = FALSE, ...)
{
  require(ChemmineR)

  stop <- FALSE
  f <- file(input, "r")
  n <- Nlines
  offset <- 0
  if (startline != 1) {
    fmap <- file(input, "r")
    shiftback <- 2
    chunkmap <- scan(fmap, skip = startline - shiftback,
                     nlines = restartNlines, what = "a", blank.lines.skip = FALSE,
                     quiet = TRUE, sep = "\n")
    startline <- startline + (which(grepl("^\\${4,4}", chunkmap,
                                          perl = TRUE))[1] + 1 - shiftback)
    if (is.na(startline))
      stop("Invalid value assigned to startline.")
    dummy <- scan(f, skip = startline - 2, nlines = 1, what = "a",
                  blank.lines.skip = FALSE, quiet = TRUE, sep = "\n")
    close(fmap)
    offset <- startline - 1
  }
  counter <- 0
  cmpid <- 1
  partial <- NULL
  while (!stop) {
    counter <- counter + 1
    chunk <- readLines(f, n = n)
    if (length(chunk) > 0) {
      if (length(partial) > 0) {
        chunk <- c(partial, chunk)
      }
      inner <- sum(grepl("^\\${4,4}", chunk, perl = TRUE)) <
        2
      while (inner) {
        chunklength <- length(chunk)
        chunk <- c(chunk, readLines(f, n = n))
        if (chunklength == length(chunk)) {
          inner <- FALSE
        }
        else {
          inner <- sum(grepl("^\\${4,4}", chunk, perl = TRUE)) <
            2
        }
      }
      y <- regexpr("^\\${4,4}", chunk, perl = TRUE)
      index <- which(y != -1)
      indexDF <- data.frame(start = c(1, index[-length(index)] +
                                        1), end = index)
      complete <- chunk[1:index[length(index)]]
      if ((index[length(index)] + 1) <= length(chunk)) {
        partial <- chunk[(index[length(index)] + 1):length(chunk)]
      }
      else {
        partial <- NULL
      }
      index <- index + offset
      indexDF <- data.frame(SDFlineStart = c(offset +
                                               1, index[-length(index)] + 1), SDFlineEnd = index)
      offset <- indexDF[length(indexDF[, 2]), 2]
      try({
        sdfset <- read.SDFset(read.SDFstr(complete))
        if (length(indexDF[, 1]) == 1) {
          suppressWarnings(sdfset <- c(sdfset, sdfset))
          resultMA <- fct(sdfset, ...)
          }
        else {
          resultMA <- fct(sdfset, ...)
        }
        if (silent == FALSE) {
          print(rownames(resultMA))
        }
        if (counter == 1 & append != TRUE) {
          unlink(output)
          write.table(resultMA, output, quote = FALSE,
                      col.names = NA, sep = "\t")
        }
        else {
          write.table(resultMA, output, quote = FALSE,
                      append = TRUE, col.names = FALSE, sep = "\t")
        }
      })
    }
    if (length(chunk) == 0) {
      stop <- TRUE
      close(f)
    }
  }
}

# plot molecules in R plot window instead of separate Java window
plot.mol = function(smi,
                    width=500,
                    height=500,
                    marg=0,
                    main='',
                    style="bow") {

  if(is.null(smi)){
    return(NULL)
  }
  #curcumin = parse.smiles("O=C(\\C=C\\c1ccc(O)c(OC)c1)CC(=O)\\C=C\\c2cc(OC)c(O)cc2")[[1]]
  molecule = rcdk::parse.smiles(smi)[[1]]
  #rcdkplot(curcumin, style="cow")
  par(mar=c(marg,marg,marg,marg)) # set margins to zero since this isn't a real plot
  dept = rcdk::get.depictor(width = width, height = height, zoom = 3, style = style,
                      annotate = "off", abbr = "on", suppressh = TRUE,
                      showTitle = FALSE, smaLimit = 100, sma = NULL)
  temp1 = rcdk::view.image.2d(molecule, dept) # get Java representation into an image matrix. set number of pixels you want horiz and vertical

  # - - return - -
  # A temp file to save the output. It will be deleted after renderImage
  # sends it, because deleteFile=TRUE.
  a <- tempfile(fileext='.png')
 
  # Generate a png
  #png(a, width=500, height=500,bg=NA)

  plot(NA,NA,xlim=c(1,10),ylim=c(1,10),xaxt='n',yaxt='n',xlab='',ylab='',main=main,bg=NA) # create an empty plot
  rasterImage(temp1,1,1,10,10) # boundaries of raster: xmin, ymin, xmax, ymax. here i set them equal to plot boundaries

  #dev.off()

  #system(gsubfn::fn$paste("convert $a -transparent white $b"))
  # Return a list
  list(src = a)
}

# @export
find.formulas <- function(mzvals, cl=FALSE, ppm=3, charge=1, element.counts = list(c("C",0,50),c("H",0,50),
                                                                                   c("N",0,50),c("O",0,50),
                                                                                   c("S",0,50),c("Na", 0, 5),
                                                                                   c("Cl", 0, 5), c("P", 0,5))){
  require(rcdk)
  require(pbapply)
  require(data.table)
  # ------------------------------------------
  found.rows <- pblapply(mzvals,cl=cl, function(mz){
    window = mz * (ppm / 1e6)
    # --- generate molecular formulae ---
    found.mfs <- generate.formula(mz, window=0.3, element.counts, validation=TRUE, charge=charge)
    rows <- if(length(found.mfs) == 0) NA else(
      rows <- lapply(found.mfs, function(formula){
        # --- check for ppm range ---
        mz.found <- formula@mass
        within.ppm <- abs(mz - mz.found) < window
        if(within.ppm){
          data.table(origMZ = mz,
                     genMZ = mz.found,
                     BaseFormula = formula@string)
        } else(data.table(origMZ=mz,
                          genMZ=NA,
                          BaseFormula=NA))
      })
    )
    # --- return ---
    unique(rbindlist(rows[!is.na(rows)]))
  })
  rbindlist(found.rows[!is.na(found.rows)])
}

get.package.table <- function(){
  status <- sapply(gbl$constants$packages, FUN=function(package){
    if(package %in% rownames(installed.packages())) "Yes" else "No"
  })
  version <- sapply(gbl$constants$packages, FUN=function(package){
    if(package %in% rownames(installed.packages())){packageDescription(package)$Version} else ""
  })
  result <-data.table(
    Package = gbl$constants$packages,
    Installed = status,
    Version = version
  )
  # --- return ---
  result
}

joanna_debugger <- function(){
  lcl <<- debug_lcl
  mSet <<- debug_mSet
  input <<- debug_input
  shown_matches <<- isolate({reactiveValuesToList(debug_matches)})
  my_selection <<- isolate({reactiveValuesToList(debug_selection)})
}
  