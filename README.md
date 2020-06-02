# MetaboShiny
Welcome to the info page on MetaboShiny! We are currently on BioRXiv and the paper itself is somewhat of a guide to how to use the software. Additionally, a visual manual to all the possible actions within the app can be viewed below.
For example input files (positive and negative peaklists + metadata) please see the `examples` folder.

http://biorxiv.org/cgi/content/short/734236v1

# Table of Contents
- [Installation](#installation)
  * [Through Docker](#through-docker)
  * [Through R](#through-r)
- [Before starting](#before-starting)
  * [Building databases](#building-databases)
  * [Updating databases](#updating-databases)
  * [Adding custom databases](#adding-custom-databases)
- [Prepare input files](#prepare-input-files)
  * [m/z peak files](#mz-peak-files)
  * [Metadata file](#metadata-file)
- [Load data files](#load-data-files)
- [Data normalization](#data-normalization)
  * [Pre-matching all m/z (optional)](#Pre-matching-all-mz-optional)
  * [Statistics panel / side bar](#statistics-panel--side-bar)
  * [Dimension reduction / plot export](#dimension-reduction--plot-export)
  * [Overview plots / search bar](#overview-plots--search-bar)

# Installation
## Through Docker
1. Install Docker. If you are on Windows and not running 10 pro or enterprise (likely) please try the following tutorial:https://docs.docker.com/toolbox/toolbox_install_windows/
2. Run the following in a terminal: `docker pull jcwolthuis/metaboshiny`
3. Make a folder in your home folder, in terminal: `mkdir -p ~/MetaboShiny/{databases,saves/admin}`
3. Run the following on command line/terminal: `docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny:latest Rscript -e "MetaboShiny::start.metshi(inBrowser=F)"`


## Through R
1. Install an R version that at least is R 3.5, install [Java](https://www.java.com/en/download/) and JDK (https://www.oracle.com/technetwork/java/javase/downloads/jdk12-downloads-5295953.html, may require creating an Oracle account)
2. Copy the commands listed in the `inst/install.metshi.R` file on this repository and run them in R (terminal or Rstudio)
3. In R(studio), run: `library(MetaboShiny); start.metshi(inBrowser=T)`;
4. Your browser should open automatically for you and you're in! Otherwise browse to the "Listening on x.x.x.x" address noted in terminal in your browser.

# Before starting
## Building databases
MetaboShiny offers over **XX** metabolite databases for m/z identification. Before any other steps are taken, it is necessary to build the databases that the user is interested in. Each database only needs to be built once. To check if a database has been built, click the "check if database exists" button below the logo (Figure [Database Tab]). The database version number and download date are listed there as well.
## Updating databases
MetaboShiny does not automatically update databases. To re-build a database of interest, click on the "build database" button below the logo in the [Database Tab](#database-tab). The database version number and download date are listed below the logo.
## Adding custom databases
**HOW TO ADD CUSTOM DB**

![Database Tab](inst/www/database_panel_a.png?raw=true "Database tab")

# Prepare input files
## m/z peak files
**DESCRIBE HOW TO GET POSITIVE/NEGATIVE MODE FILES FROM RAW MS DATA?
LINK TO EXAMPLE DATA.**

## Metadata file
**DESCRIBE METADATA FILE FORMAT, NECESSARY COLUMNS, ETC.
LINK TO EXAMPLE DATA.**

# Load data files
To load your data into MetaboShiny, follow the 6 steps listed in the [File Import](#file-import) figure.
1. Enter a unique project name.
2. Set the error margin of your mass spectrometer in parts per million (ppm).
3. Set the m/z missing values percentage threshold. This is how many samples are allowed to be missing each m/z value without it being filtered out. A common rule of thumb is to set this value to 20% (**REF TO 20/80 RULE PAPER**).
4. Select input data files
 * 4a. (optional) Input a regex string to to adjust peaklist names to metadata sample names - the match is removed from each name.
 * 4b. Upload your metadata and positive and negative mode m/z peak files.
5. Click on the arrow to merge peak data and metadata and create an SQLITE database.
6. Once step 5 is completed (green tick mark), click on the arrow to the right to generate a CSV file that will serve as the input into MetaboShiny.
7. Once step 6 is completed (green tick mark), continue to the [Data normalization](#data-normalization) step.

**THIS FIGURE NEEDS TO BE UPDATED! The % filtering is missing**
![File Import](inst/www/file_import.png?raw=true "File import")

# Data normalization
The data needs to be normalized in order to compare m/z peak values between samples and batches.

## Batches and concentration
If your metadata only contains one batch and no column that represents concentration (**EXPLAIN WHAT THIS MEANS**), then you can skip this part and continue to the [Filtering and normalization](#filtering-and-normalization) step. Otherwise, follow the steps below.
 1. Click on the "Get options" button ([Normalization](#normalization) figure).
 2. If applicable, select the variable that represents concentration in your data.
 3. If applicable, select the variable that contains your multiple batch IDs.

## Filtering and normalization
In this section, you will find multiple options and methods to filter and normalize your data. The best selection will depend on each user's data and we encourage you to look into the different methods that can be applied here. After normalization, the distribution of pre- and post-normalized peak values will be plotted for a randomly selected set of m/z values and samples, so the user can see how the data distribution has changed with the normalization and adjust their parameters if needed (see [Normalization](#normalization) figure).

Select one of each of the options for the following normalization features and then press "Go" (see [Normalization](#normalization) figure).
**ADD DESCRIPTIONS/INFO ON THE DIFFERENT OPTIONS**
- Filtering **WHAT IS BEING FILTERED OUT?**
  - Interquartile range
  - Mean
  - Median absolute deviation
  - Median
  - Non-parametric relative standard deviation (stdev)
  - Relative standard deviation (stdev)
  - Standard deviation
  - None
- Normalization type
  - By reference compound
  - By reference feature
  - By sample-specific factor
  - Median
  - Quantile normalization
  - Sum
  - None
- Data transformation
  - Cubic root transform
  - Log transform
  - None
- Scaling
  - Autoscale/Z-transform
  - Mean-center
  - Pareto scaling
  - Range scaling
  - None
- Missing values
  - Half feature minimum
  - Half sample minimum
  - Total minimum
  - Random forest. It is possible to adjust the number of trees built per variable and whether to parallelize based on forests or variables, or not.
  - KNN imputation
  - SVD imputation
  - BPCA imputation
  - PPCA imputation
  - Median
  - Mean
  - Leave them out
  - Leave them alone
- Outliers
 The user can choose whether to exclude outliers from the data analysis by toggling the "Exclude outliers?" tab.

**THIS FIGURE NEEDS TO BE UPDATED! The % filtering is now in the data upload step**
![Normalization](inst/www/normalization.png?raw=true "Normalization")



## Pre-matching all m/z (optional)
![Prematching](inst/www/prematching.png?raw=true "Prematching")
## Statistics panel / side bar
![Statistics panel start](inst/www/stats.png?raw=true "Statistics Panel Start")
### Dimension reduction / plot export
![Dimension reduction and plot export](inst/www/dimred_export.png?raw=true "Dimension reduction and plot export")
### Overview plots / search bar
![Volcano plot and search bar](inst/www/stats2.png?raw=true "Volcano plot and search bar")
### Per m/z analysis / formula prediction
![Per m/z analysis and prediction sidebar](inst/www/permz_predict.png?raw=true "Per m/z analysis and prediction sidebar")
### Machine learning
![Machine learning settings and metadata bar](inst/www/ml1.png?raw=true "Machine learning settings and metadata bar")
![Machine learning results](inst/www/ml2.png?raw=true "Machine learning results")
### Match isotope scoring and abstract searching
![Isoscoring and abstract searching](inst/www/isoscore_wordcloud.png?raw=true "Isoscoring and abstract searching")
### Plot customization
![Plot customization](inst/www/colorbar.png?raw=true "Plot customization")



Please report any issues and feedback on the Issues page here, along with suggestions! =)
