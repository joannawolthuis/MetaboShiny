# MetaboShiny
Welcome to the info page on MetaboShiny! We are currently on BioRXiv and the paper itself is somewhat of a guide to how to use the software. Additionally, a visual manual to all the possible actions within the app can be viewed below.
For example input files (positive and negative peaklists + metadata) please see the `examples` folder.

http://biorxiv.org/cgi/content/short/734236v1

# Table of Contents
- [MetaboShiny](#metaboshiny)
- [Table of Contents](#table-of-contents)
- [Installation](#installation)
  * [Through Docker](#through-docker)
  * [Through R](#through-r)
- [Before starting](#before-starting)
  * [Settings](#settings)
    + [Global](#global)
    + [Project](#project)
    + [Search](#search)
    + [Adducts](#adducts)
    + [Formula prediction and lookup](#formula-prediction-and-lookup)
    + [Isotope scoring](#isotope-scoring)
    + [Aesthetic](#aesthetic)
  * [Building databases](#building-databases)
  * [Updating databases](#updating-databases)
  * [Adding custom databases](#adding-custom-databases)
- [Prepare input files](#prepare-input-files)
  * [m/z peak files](#m-z-peak-files)
  * [Metadata file](#metadata-file)
- [Load data files](#load-data-files)
- [Data normalization](#data-normalization)
  * [Batches and concentration](#batches-and-concentration)
  * [Filtering and normalization](#filtering-and-normalization)
- [Continue working on an existing data set](#continue-working-on-an-existing-data-set)
  * [Save your data](#save-your-data)
  * [Load existing data set](#load-existing-data-set)
- [Pre-matching all m/z](#pre-matching-all-m-z)
- [Analysis](#analysis)
  * [Side bar](#side-bar)
    + [Switch/subset](#switch-subset)
    + [Search](#search-1)
      - [Database search settings](#database-search-settings)
        * [Compound info](#compound-info)
    + [Compound matches filtering](#compound-matches-filtering)
      - [PubMed search and word cloud](#pubmed-search-and-word-cloud)
    + [Plot aesthetics](#plot-aesthetics)
    + [Metadata](#metadata)
  * [Statistics](#statistics)
    + [Dimension reduction](#dimension-reduction)
      - [Principal components analysis (PCA)](#principal-components-analysis--pca-)
      - [Partial least squares-discriminant analysis (PLS-DA)](#partial-least-squares-discriminant-analysis--pls-da-)
      - [t-Distributed Stochastic Neighbor Embedding (t-SNE)](#t-distributed-stochastic-neighbor-embedding--t-sne-)
    + [Per m/z](#per-m-z)
      - [t-test](#t-test)
      - [Fold-change](#fold-change)
      - [Pattern](#pattern)
      - [Analysis of variance (ANOVA)](#analysis-of-variance--anova-)
      - [Analysis of variance - simultaneous component analysis (ASCA)](#analysis-of-variance---simultaneous-component-analysis--asca-)
      - [Multivariate empirical Bayes statistical time-series analysis (MEBA)](#multivariate-empirical-bayes-statistical-time-series-analysis--meba-)
    + [Overview analyses](#overview-analyses)
      - [Volcano plot](#volcano-plot)
      - [Heatmap](#heatmap)
      - [Venn diagram](#venn-diagram)
      - [Power calculation](#power-calculation)
      - [Enrichment](#enrichment)
    + [Machine learning](#machine-learning)
      - [Machine learning settings](#machine-learning-settings)
      - [Machine learning results](#machine-learning-results)
- [How to cite](#how-to-cite)
- [Feedback](#feedback)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>

# Installation
## Through Docker
1. Install Docker. If you are on Windows and not running 10 pro or enterprise (likely) please try the following tutorial:https://docs.docker.com/toolbox/toolbox_install_windows/
2. Run the following in a terminal: `docker pull jcwolthuis/metaboshiny`
3. Make a folder in your home folder, in terminal: `mkdir -p ~/MetaboShiny/{databases,saves/admin}`
3. Run the following on command line/terminal: `docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny:latest Rscript -e "MetaboShiny::start_metshi(inBrowser=F)"`


## Through R
1. Install an R version that at least is R 3.5, install [Java](https://www.java.com/en/download/) and JDK (https://www.oracle.com/technetwork/java/javase/downloads/jdk12-downloads-5295953.html, may require creating an Oracle account)
2. Copy the commands listed in the `inst/install.metshi.R` file on this repository and run them in R (terminal or Rstudio)
3. In R(studio), run: `library(MetaboShiny); start_metshi(inBrowser=T)`;
4. Your browser should open automatically for you and you're in! Otherwise browse to the "Listening on x.x.x.x" address noted in terminal in your browser.

# Before starting
## Settings
![Settings](inst/www/settings.png?raw=true "Settings")
### Global
1. Set the number of cores that MetaboShiny can use for analysis.
2. Define what to do with samples that are missing metadata (omit/keep).
### Project
The drop-down list shows all projects that have been saved.
1. Select project name and press "Apply" to select project.
2. Click "load" at the bottom of the screen to load the saved progress on the project.
### Search
When searching databases for compounds matching your m/z values you can click on the compounds in the match list. If selected, MetaboShiny will copy the ticked compound information to the clipboard (SMILES, name, formula).
### Adducts
The "Definitions" tab shows the current adduct table. Below the table is a field to import another adduct table. The "Rules" field shows the adduct rules that are used to calculate adducts when matching m/z values with database compounds. Below the table is a field to import a new adduct rule table.
### Formula prediction and lookup
Here you can set the parameters and rules to use when predicting chemical formulas.
### Isotope scoring
Select the method to use to score compounds that have the same weight (currently only M-score available). Set the intensity imprecision (default: 2%).
### Aesthetic
Change MetaboShiny colors, fonts, and font sizes. Restart MetaboShiny to apply changes.

## Building databases
MetaboShiny offers multiple metabolite databases for m/z identification. Before any other steps are taken, it is necessary to build the databases that the user is interested in. Each database only needs to be built once. To check if a database has been built, click the "check if database exists" button below the logo (Figure [Database Tab]). The database version number and download date are listed there as well.
## Updating databases
MetaboShiny does not automatically update databases. To re-build a database of interest, click on the "build database" button below the logo in the Database Tab. The database version number and download date are listed below the logo.

## Adding custom databases
Users often have their own in-house databases. You can add these in MetaboShiny as well. To do so, scroll down on the database overview page to the large '+' button. Here, you'll be shown examples of what files are required (a logo, a base csv file with columns displayed in the pop-up, and a database name + description). After you upload these, your database will be available for building after the next restart and you can click the "build database" button as you usually would.

![Database Tab](inst/www/database_panel_a.png?raw=true "Database tab")

# Prepare input files
## m/z peak files
MetaboShiny does not accept raw peak data. We suggest using either XCMS (with the MetaboAnalyst export option) or another method of choice such as MSnbase. You can find examples of three different accepted data formats (MetaboAnalyst-like, MetaboShiny native and Metabolights) in the inst/examples folder.

## Metadata file
MetaboShiny, unless using the MetaboAnalyst format, requires an additional metadata table. This should minimally have a 'sample' column that contains the same sample identifiers used in the peak table files, an 'individual' column (since multiple samples can come from one individual in time series data) and at least one column on experimental group or something alike. Examples of metadata formats are also present in the inst/examples folder.

# Load data files
![File Import](inst/www/file_import.png?raw=true "File import")

1. Enter a unique project name.
2. Set the error margin of your mass spectrometer in parts per million (ppm).
3. Set the m/z missing values percentage threshold. This is how many samples are allowed to be missing each m/z value without it being filtered out.
4. Select input data files
 * 4a. (optional) Input a regex string to to adjust peaklist names to metadata sample names - the match is removed from each name.
 * 4b. Upload your metadata and positive and negative mode m/z peak files.
5. Click on the arrow to merge peak data and metadata and convert them to a format that will serve as the input for the analyses.
6. Once step 5 is completed (green tick mark), continue to the [Data normalization](#data-normalization) step.

# Data normalization
The data needs to be normalized in order to compare m/z peak values between samples and batches.

![Normalization](inst/www/normalization.png?raw=true "Normalization")

## Batches and concentration
If your metadata only contains one batch and no column that represents concentration, then you can skip this part and continue to the [Filtering and normalization](#filtering-and-normalization) step. Otherwise, follow the steps below.
 1. Click on the "Get options" button.
 2. If applicable, select the variable that represents concentration in your data.
 3. If applicable, select the variable that contains your multiple batch IDs.

## Filtering and normalization
In this section, you will find multiple options and methods to filter and normalize your data. The best selection will depend on each user's data and we encourage you to look into the different methods that can be applied here. After normalization, the distribution of pre- and post-normalized peak values will be plotted for a randomly selected set of m/z values and samples, so the user can see how the data distribution has changed with the normalization and adjust their parameters if needed.

Select one of each of the options for the following normalization features and then press "Go". It is advised to save your data after completing this step (button on the bottom center of the screen).
- **Filtering options**
  - Interquartile range
  - Mean
  - Median absolute deviation
  - Median
  - Non-parametric relative standard deviation (stdev)
  - Relative standard deviation (stdev)
  - Standard deviation
  - None
- **Normalization type**
  - By reference compound
  - By reference feature
  - By sample-specific factor
  - Median
  - Quantile normalization
  - Sum
  - None
- **Data transformation**
  - Cubic root transform
  - Log transform
  - None
- **Scaling**
  - Autoscale/Z-transform
  - Mean-center
  - Pareto scaling
  - Range scaling
  - None
- **Missing values**
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
- **Outliers** The user can choose whether to exclude outliers from the data analysis by toggling the "Exclude outliers?" tab.

# Continue working on an existing data set
## Save your data
On the bottom center of the screen, you will find a "save" button. The data set will be saved under the name you chose in the [file importing step](#load-data-files).
## Load existing data set
1. Go to the "settings" panel
2. Select the "Project" tab
3. Select your project from the drop-down list and press "Apply"
4. Press the "load" button on the bottom center of the screen

# Pre-matching all m/z
*This step is optional*
![Prematching](inst/www/prematching.png?raw=true "Prematching")

In the pre-matching panel the user can match all m/z values with all or a subset of the available databases. This can be a time-consuming step if the dataset is large and many databases are selected, but will make searching for possible m/z metabolite matches much faster in the data analysis step.
1. Toggle the "Do matching beforehand?" button to "Yes".
2. Select the databases you wish to find matches in or click on the shopping basket to match with all databases.
3. Click on "Find matches". This can take a few minutes.
4. Save your data (button on the bottom center of the screen).
If you wish to match with other or more databases, click the "Clear matches" button and re-do the steps above.


# Analysis
The analysis tab has two sections; the statistics panel and the side bar. The statistics panel contains tabs with different statistical analysis methods and the side bar contains multiple functions regarding variable choices, data subsetting, plot aesthetics, m/z matching, and plot export.

## Side bar
The side bar contains four tabs, whose descriptions and functions you can find below.
  * Note that the side bar can be resized by dragging the left-hand side.

### Switch/subset
*See [Statistics panel figure](#statistics)*
 - **Current experiment** Shows the variable(s) and subset(s) that are currently being analysed.
 - **Change of variable of interest** Here you can choose to inspect one variable, two variables in combindation, time-series, or time-series in combination with one variable. Press "do stats on selected" to change the current experiment for analysis.
 - **Subset data** The "Current sample count" shows the number of samples that are analysed in the current experiment. To subset data, select the variable that you want to subset based on and then select the group(s) that you want to inspect. Click "click to subset" to apply changes.
 - **Load existing meta-dataset** Every time you use the subset/switch option, it saves the results from that subset. Use the drop-down menu under the subsetting field to go back to previously defined subsets.

### Search
![Search tab](inst/www/mz_search.png?raw=true "Search tab")

#### Database search settings
   * If the data is pre-matched, this part is skipped. Otherwise, here you can select the databases that you want to match your data to. Click the shopping basket to add all available databases.
   * Click the inspector icon to match your selected m/z value to your selected databases
##### Compound info
This tab will display all database matches for a selected m/z value from the statistics panel.
   * The displayed table in the "mz > molecule" tab can be sorted based on m/z value, adduct, isotope percentage, or the m/z value distance from the database range (dppm).
   * The table can be copied or exported as a .csv or .xlsx file.
   * By clicking on the funnel icon, you can filter the matched results based on adduct, database, and main and minor isotope.
   * In the "molecule > mz" tab you can search for a specific metabolite name and the resulting table will list all m/z values from the data that match a corresponding metabolite in any of the databases.
   * In the match menu, when selecting a compound:
     - the name and SMILES or formula (specified in settings) are copied to the clipboard.
     - in the compound description field, clicking on a database icon will copy the database id to the clipboard.
### Compound matches filtering
Click the funnel above the "compound info" section to filter the database matches based on adduct, main/minor isotope, or database.
  * **Adduct** The bar chart shows the ratio of matches having a certain adduct. Hover over the slices for number of matches. Click on the slices that you are interested to filter for those adducts.
  * **Isotope** You can filter for main (100% peak) or minor (<100% peak) isotopes. Click pie slices to filter.
  * **Databases** The pie shows the ratio of results that come from each database. Hover over the slices for the number of matches. Click on the slices of the databases you want to filter for.

#### PubMed search and word cloud
![PubMed search](inst/www/pubmed_wordcloud.png?raw=true "PubMed search and word cloud")
   1. Settings
      - For a PubMed search, enter your search term (e.g. metabolite name) and specify the publishing date range and how many abstracts to use in the search.
      - To use the results from your m/z database matches, toggle "own word" to "from matches".
   2. Press 'Plot' to start your search and render the word cloud
   3. In the 'filters' tab, you can search for a second set of abstracts that you want to *exclude* from your previous search.
   4. In the 'plot' tab you will find the word cloud. Toggle "cloud" to "barchart" to render a bar chart instead of a word cloud. Select the number of words to use for your plot.
   5. Filter out words that are commonly used in the English language or metabolomics research by selecting predefined word filters. The words from your search in *3* will be listed as a separate set here.
      * **stopwords** This set contains the top 200 most common words in the English language.
      * **metabolomics** This set contains words that are common in the metabolomics field: metabolism, metabolic, metabolomic, metabolomics, biochemical, mass, spectrometry, nmr, direct, infusion, exposome, papers, compounds, and compound.
      * **default** This set contains the words "exposome", "synonyms", and all the available database names.
   6. Click on a word in the word cloud to show PubMed abstracts mentioning that word, and their PubMed IDs.


### Plot aesthetics
**Plot aesthetic changes are applied on plot re-creation**
 ![Plot customization](inst/www/colorbar.png?raw=true "Plot customization")

 - **Plot style** Here you can choose whether to use plotly or ggplot for the figures in the statistics panel. Plotly provides interactive plots with the ability to zoom and hover over data points for information, whereas ggplot is faster but does not have interactive features. Additionally, you can choose whether plots are exported as .png or .svg files.
 - **Style(s)** Choose the plotting style for the scatter plots for t-test and ANOVA results. Options are beeswarm, boxplot, scatterplot, and violin plots. It is possible to select multiple styles, e.g., violin plot with a beeswarm scatter.
 - **Stats shown** Select whether to show median or mean lines in the plots.
 - **Marker shape** Select the variable that determines the marker shape.
 - **Color** Select the variable that determines the marker color.
 - **Hover text** Select the variable that is displayed when hovering over data points.
 - **Plot theme** Select one of the following themes for the plot background
   - Minimal
   - Grid, gray background
   - Grid, black background
   - Grid, white background
   - No grid, white background
   - Grid, white background, gray axes
   - Line drawing
 - **Continuous data** Select a color scheme for plotting continuous data, which is used in heatmaps and volcano plots.
 - **Discrete data** Select colors to use when plotting discrete data. Click on the color field to select a color or type in the HEX code.

### Metadata
In this tab you can upload new metadata. The file should be in a .csv format and contain a column with sample IDs and any new metadata as additional columns with new unique headers.

**Uploading new metadata replaces the old metadata**, so make sure to include all relevant columns in the new file.

## Statistics
The statistics panel has four tabs whose contents change based on whether the current experiment is a one-factor, two-factor, or a time series analysis. The four statistics categories are dimension reduction methods), per m/z value analyses, overview analyses, and machine learning.
![Statistics panel start](inst/www/stats.png?raw=true "Statistics Panel Start")

**If you click any m/z in the result table or in the plotly scatterplots, this will be recorded in the side bar. If the [prematching step](#pre-matching-all-mz) was performed then the matches will appear in the [side bar](#compound-info). Otherwise you can manually [search for matches](#search).**

### Dimension reduction

#### Principal components analysis (PCA)
1. Go to "Settings" and press "click to start PCA" to start the analysis.
2. Once done, a "Plots" section and a "Tables" sections will appear.
3. In "Plots" you can view the scores and the loadings plots. Below the plots you can select which principal component (PC) is plotted on which axis, and whether to render a 2D or 3D version of the plot.
4. The "Tables" section has three tabs showing the %explained variance of each principal component ("table"), the scree plot, and loadings.

#### Partial least squares-discriminant analysis (PLS-DA)
1. Go to "Settings" and press "click to start PLS-DA" to start the analysis.
2. Once done, a "Plots" section and a "Tables" sections will appear.
3. In the "Plots" section, you can switch between the samples, loadings, and performance (cross-validation and permutation). Below the plots you can select which principal component (PC) is plotted on which axis, and whether to render a 2D or 3D version of the plot.
4. The "Tables" section shows the %explained variance of each principal component and the loadings.

#### t-Distributed Stochastic Neighbor Embedding (t-SNE)
1. In "Settings", select the initial number of dimensions (default 30), perplexity (default 30), and maximum number of iterations (default 1000), and press "click to start t-SNE" to start the analysis.
2. Once done, a "Plots" section will appear showing the t-SNE plot. Below the plot is a toggle button to make a 2D or 3D plot.

### Per m/z
#### t-test
1. In "Settings", select whether to make a non-parametric t-test (default "no") and equal variance (default "yes"). Press "click to start t-test" to start the analysis.
2. Once finished, a "Plots" and "Tables" sections will appear.
3. In "Plots", you will find a plot showing your m/z values ordered on the x-axis and the negative log10 of the p-value from the t-test.
4. In "Tables", you will find a table showing the t-stat, p-value, negative log10 of the p-value, and the false discovery rate (FDR).

#### Fold-change
1. In "Settings", press "click to start fold-change analysis".
2. Once finished, a "Plots" and "Tables" sections will appear.
3. In "Plots", you will find a plot showing your m/z values ordered on the x-axis and the log2 of the fold-change value on the y-axis.
4. In "Tables", you will find a table showing the fold-change value and the log2 of the fold-change.

#### Pattern
1. Drag the panels to select the pattern for correlation (low-high). The pattern corresponds to the expected pattern used to search for m/z values whose peak intensities increase linearly between the groups.
2. Choose the correlation metric (pearson, spearman, kendall)
3. Click to start pattern finding
4. Once finished, a "Plots" and "Tables" section will appear.

#### Analysis of variance (ANOVA)
1. In "Settings", press "click to start ANOVA".
2. Once finished, a "Plots" and "Tables" sections will appear.
3. In "Plots", you will find a plot showing your m/z values ordered on the x-axis and the negative log10 of the p-value on the y-axis.
4. In "Tables", you will find a table showing the p-value and false discovery rate (FDR) from a Fisher's Least Significant Difference (LSD) test.

#### Analysis of variance - simultaneous component analysis (ASCA)
1. In "settings", click to start ASCA analysis
2. A "Tables" section will appear that shows the leverage and SPE values of the m/z values.

#### Multivariate empirical Bayes statistical time-series analysis (MEBA)
1. In "settings", click to start MEBA analysis
2. A "Tables" section will appear showing the m/z value Hotelling/T2 scores.

### Overview analyses
#### Volcano plot
1. In "Settings", press "click to make volcano plot".
2. Once finished, a "Plots" and "Tables" sections will appear.
3. In "Plots", you will find a plot with the log2 of the fold-change value on x-axis and the negative log10 of the t-test p-value on the y-axis.
4. In "Tables", you will find a table showing the fold-change value, the log2 of the fold-change, the raw p-value, and the negative log10 of the p-value.

#### Heatmap
1. In "Settings", select the analysis that you want to use for your heatmap. The drop-down list will give a list of all analyses that have been performed and can be used for a heatmap. Choose wether to only use significant hits (default "no"), and whether to color based on all metabolites (default "no") or only those that will appear in the heatmap. Press "Click to make heatmap".
2. The heatmap will appear in the "Plots" section. Use the sliding scale to select how many m/z values to feature on your heatmap.

#### Venn diagram
In the venn diagram you can see which m/z values overlap between the different analyses.
1. Select the analyses you would like to compare by clicking on them and pressing the "down" arrow button. Select an analysis from the lower list and press the "up" arrow to remove the analysis from the comparison.
2. Select how many hits to include from each analysis.
3. Press "click to make venn diagram". The diagram will appear on the right.
4. Below the diagram, you can select any of the included analysis to view the m/z values that overlap between them.

#### Power calculation

#### Enrichment

### Machine learning
MetaboShiny supports most of the models included in the *caret* package and adds some functionality based on training/testing selection and up/downsampling.

#### Machine learning settings
1. Select the *caret* machine learning algorithm to use for your model building. Default is *rf* or Random Forest. Furthermore, we use *glmnet* or LASSO regression in our paper as well. 
2. Select the method you wish to use for cross validation. By default, regular cross validation is suggested (users can choose the fold amount, or leave-one-out cross validation).
3. Each algorithm has custom parameters to set - oftentimes if left empty, caret will auto-tune to find the optimal settings, but it is suggested to select them yourself. 
4. Next, select the percentage of data you want to include in the training set. 
5. *Optional* If you want to include metadata as predictive columns in your model, you can do so.
6. *Optional* If you want to test and train on specific subsets of the data, select a subset as if you were going to completely subset via the sidebar, but do not submit through the button. Instead, use the buttons on the ML settings page to register which subsets you want to train and test on.
7. *Optional* Select if you want to balance classes through up/downsampling, and if so, which method to use.
8. Select the amount of models built, and if you want to re-split in train/test on every repeat, or build multiple models on the same data split.
9. Press "click to start machine learning" 
10. Once finished, an additional "results" tab will become available.
![Machine learning settings](inst/www/machinelearning_settings.png?raw=true "Machine learning settings")

#### Machine learning results
The machine learning results tab features two main plot + table sections to work with. First is the Receiver Operating Characteristic (ROC) curve often used to visualise the effectiveness of a predictive model. Here, a non-predictive model is shown as a diagonal straight line, and a 'perfect' model is a 90 degree curve hitting the top left of the graph. True positive rate and false positive rate are the axis descriptors here. Each of the models built has its own curve - in multivariate data this becomes more complex and curves are shown for each pair of possible class comparisons (e.g. with A B C one would see curves for A vs B, A vs C, A vs B, B vs C). When clicking a curve, the table below displays the features (m/z + metadata) used in this model and their importance within the model. Higher importance means that the feature is important to do a correct prediction. Models can also be selected in the left-hand table to explore the involved features without using the plot.
The second plot available is a bar plot displaying a top amount of features based on feature importance, with more important features (averaged over all repeats) displayed first in the bar plot. Users can click these to register the m/z to use in further searching.
![Machine learning results](inst/www/machinelearning_results.png?raw=true "Machine learning results")

# Feedback
Please report any issues and feedback on the Issues page here, along with suggestions! =)
