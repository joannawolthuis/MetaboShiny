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
  * [Batches and concentration](#batches-and-concentration)
  * [Filtering and normalization](#filtering-and-normalization)
- [Pre-matching all m/z](#pre-matching-all-mz)
- [Analysis](#analysis)
  *

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
 1. Click on the "Get options" button (see [Normalization](#normalization) figure).
 2. If applicable, select the variable that represents concentration in your data.
 3. If applicable, select the variable that contains your multiple batch IDs.

## Filtering and normalization
In this section, you will find multiple options and methods to filter and normalize your data. The best selection will depend on each user's data and we encourage you to look into the different methods that can be applied here. After normalization, the distribution of pre- and post-normalized peak values will be plotted for a randomly selected set of m/z values and samples, so the user can see how the data distribution has changed with the normalization and adjust their parameters if needed (see [Normalization](#normalization) figure).

Select one of each of the options for the following normalization features and then press "Go" (see [Normalization](#normalization) figure). It is advised to save your data after completing this step (button on the bottom center of the screen).
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


# Pre-matching all m/z
*This step is optional*

In the pre-matching tab (see [Prematching](#prematching) figure) the user can match all m/z values with all or a subset of the available databases. This can be a time-consuming step if the dataset is large and many databases are selected, but will make searching for possible m/z metabolite matches much faster in the data analysis step.
1. Toggle the "Do matching beforehand?" button to "Yes" (see [Prematching](#prematching) figure).
2. Select the databases you wish to find matches in or click on the shopping basket to match with all databases.
3. Click on "Find matches". This can take a few minutes.
4. Save your data (button on the bottom center of the screen).
If you wish to match with other or more databases, click the "Clear matches" button and re-do the steps above.

![Prematching](inst/www/prematching.png?raw=true "Prematching")


# Analysis
The analysis tab has two sections; the statistics panel and the side bar. The statistics panel contains tabs with different statistical analysis methods and the side bar contains multiple functions regarding variable choices, data subsetting, plot aesthetics, m/z matching, and plot export.

## Side bar
The side bar contains four tabs, whose descriptions and functions you can find below.

### Switch/subset
 - *Current experiment* Shows the variable(s) and subset(s) that are currently being analysed.
 - *Change of variable of interest* Here you can choose to inspect one variable, two variables in combindation, time-series, or time-series in combination with one variable. Press "do stats on selected" to change the current experiment for analysis.
 - *Subset data* The "Current sample count" shows the number of samples that are analysed in the current experiment. To subset data, select the variable that you want to subset based on and then select the group(s) that you want to inspect. Click "click to subset" to apply changes.
 - *Load existing meta-dataset* **EXPLANATION OF THIS PART**

### Search
 - *Settings* If the data is pre-matched, this part is skipped. Otherwise, here you can select the databases that you want to match your data to. Click the shopping basket to add all available databases.
 - *Compound info* This tab will display all database matches for a selected m/z value from the statistics panel (see **XX**). The displayed table in the "mz > molecule" tab can be sorted based on m/z value, adduct, isotope percentage, or the m/z value distance from the database range (dppm). The table can be copied or exported as a .csv or .xlsx file. By clicking on the funnel icon, you can filter the matched results based on adduct, database, and main and minor isotope. In the "molecule > mz" tab you can search for a specific metabolite name and the resulting table will list all m/z values from the data that match a corresponding metabolite in any of the databases.

### Plot aesthetics
 - *Plot style* Here you can choose whether to use plotly or ggplot for the figures in the statistics panel. Plotly provides interactive plots with the ability to zoom and hover over data points for information, whereas ggplot is faster but does not have interactive features. Additionally, you can choose whether plots are exported as .png or .svg files.
 - *Style(s)* Choose the plotting style for the scatter plots for t-test and ANOVA results. Options are beeswarm, boxplot, scatterplot, and violin plots. It is possible to select multiple styles, e.g., violin plot with a beeswarm scatter.
 - *Stats shown* Select whether to show median or mean lines in the plots.
 - *Marker shape* Select the variable that determines the marker shape.
 - *Color* Select the variable that determines the marker color.
 - *Hover text* Select the variable that is displayed when hovering over data points.
 - *Plot theme* Select one of the following themes for the plot background
   - Minimal
   - Grid, gray background
   - Grid, black background
   - Grid, white background
   - No grid, white background
   - Grid, white background, gray axes
   - Line drawing
 - *Continuous data* Select a color scheme for plotting continuous data, which is used in heatmaps and volcano plots. All options of the **XX** package are possible (**LINK**).
 - *Discrete data* Select colors to use when plotting discrete data. Click on the color field to select a color or type in the HEX code.

### Metadata
In this tab you can upload new metadata. The file should be in a .csv format and contain a column with sample IDs and any new metadata as additional columns with new unique headers.



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
