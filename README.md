# MetaboShiny
Welcome to the info page on MetaboShiny! I am still working on a proper readme, please have patience with me :-)
We are currently on BioRXiv and the paper itself is somewhat of a guide to how to use the software.

http://biorxiv.org/cgi/content/short/734236v1

Please report any issues and feedback on the Issues page here, along with suggestions! =)

# STEP 0: install docker

# STEP 1: pull from docker
docker pull jcwolthuis/metaboshiny

# STEP 3: run the following on command line/terminal
docker run -p 8080:8080 -v ~/MetaboShiny/:/root/MetaboShiny/:cached --rm -it jcwolthuis/metaboshiny Rscript startShiny.R

# STEP 4: navigate to http://localhost:8080

# STEP 5: you're in!
