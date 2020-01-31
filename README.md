# MetaboShiny
Welcome to the info page on MetaboShiny! I am still working on a proper readme, please have patience with me :-)
We are currently on BioRXiv and the paper itself is somewhat of a guide to how to use the software.
For example input files (positive and negative peaklists + metadata) please see the `examples` folder.

http://biorxiv.org/cgi/content/short/734236v1

# Manual #

![Database Tab](inst/www/database_panel_a.png?raw=true "Database tab")
![File Import](inst/www/file_import.png?raw=true "File import")
![Normalization](inst/www/normalization.png?raw=true "Normalization")
![Prematching](inst/www/prematching.png?raw=true "Prematching")
![Statistics panel start](inst/www/stats.png?raw=true "Statistics Panel Start")
![Dimension reduction and plot export](inst/www/dimred_export.png?raw=true "Dimension reduction and plot export")
![Volcano plot and search bar](inst/www/stats2.png?raw=true "Volcano plot and search bar")
![Per m/z analysis and prediction sidebar](inst/www/permz_predict.png?raw=true "Per m/z analysis and prediction sidebar")
![Machine learning settings and metadata bar](inst/www/ml1.png?raw=true "Machine learning settings and metadata bar")


Please report any issues and feedback on the Issues page here, along with suggestions! =)
## R PATH ##
# STEP 1: requirements (these may require admin rights)
install an R version that at least is R 3.5, install [Java](https://www.java.com/en/download/) and JDK (https://www.oracle.com/technetwork/java/javase/downloads/jdk12-downloads-5295953.html, may require creating an Oracle account)
# STEP 2
copy the commands listed in the `inst/install.metshi.R` file on this repository and run them in R (terminal or Rstudio)
# STEP 3
in R: `library(MetaboShiny); start.metshi(inBrowser=T)`;
# STEP 3
browser should open automatically for you and you're in! Otherwise browse to the "Listening on x.x.x.x" address noted in terminal in your browser.

## DOCKER PATH ##
# STEP 1
Install docker
If you are on Windows and not running 10 pro or enterprise (likely) please try the following tutorial:
https://docs.docker.com/toolbox/toolbox_install_windows/
# STEP 2:
pull from docker
`docker pull jcwolthuis/metaboshiny`
# STEP 3:
run the following on command line/terminal
```
docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny R -e 'library(MetaboShiny); start.metshi(inBrowser=F);'
```
# STEP 4
navigate to http://localhost:8080 OR if on docker toolbox http://192.168.99.100:8080/
# STEP 5
you're in!
