# MetaboShiny
Welcome to the info page on MetaboShiny! I am still working on a proper readme, please have patience with me :-)
We are currently on BioRXiv and the paper itself is somewhat of a guide to how to use the software.
For example input files (positive and negative peaklists + metadata) please see the 'examples' folder.

http://biorxiv.org/cgi/content/short/734236v1

Please report any issues and feedback on the Issues page here, along with suggestions! =)
## R PATH ##
# STEP 1
pull the repository through the terminal or through download
# STEP 2: requirements (these may require admin rights)
install an R version that at least is R 3.5
install Java 
install JDK (may require creating an Oracle account)
# STEP 3
run either start.sh or startShiny.R (the latter with "Rscript startShiny.R" from the command line)
# STEP 4
wait for all the packages to install 
# STEP 5
browser should open automatically for you and you're in!

## DOCKER PATH ##
# STEP 1
Install docker
If you are on Windows and not running 10 pro or enterprise (likely) please try the following tutorial:
https://docs.docker.com/toolbox/toolbox_install_windows/
# STEP 2:
pull from docker
docker pull jcwolthuis/metaboshiny
# STEP 3:
run the following on command line/terminal
docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny Rscript startShiny.R
# STEP 4
navigate to http://localhost:8080 OR if on docker toolbox http://192.168.99.100:8080/
# STEP 5
you're in!
