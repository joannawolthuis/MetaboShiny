library(XML)
library(httr)


test_hmdb <- 'http://www.hmdb.ca/metabolites/HMDB31795.xml'

test_chebi <- 'https://www.ebi.ac.uk/webservices/chebi/2.0/test/getCompleteEntity?chebiId=132282'
# -- chebi needs additional step to get xml ---
r <- GET(test_chebi)
print(content(r, "text"))
xmlParse(content(r, "parsed"), asText = T, useInternalNodes = T)
# ---------------------------------------------
data <- xmlParse(test_hmdb)

desc_nodes <- getNodeSet(data, 
                            "/*/description|definition")

desc_values = sapply(description, xmlValue)
