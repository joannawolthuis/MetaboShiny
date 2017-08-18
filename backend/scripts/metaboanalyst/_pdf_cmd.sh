# This script is called from Java directly without using Rserve
# in order to get around the PDF issue
#! /bin/bash 

PATH="/usr/texbin:/Library/TeX/texbin:/usr/bin:/usr/local/bin:${PATH}"
export PATH
Rscript $@