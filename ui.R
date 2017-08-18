# Rely on the 'WorldPhones' dataset in the datasets
# package (which generally comes preloaded).
library(datasets)

asca.table <- read.csv("backend/data/asca_sigAB.csv")

# Use a fluid Bootstrap layout
fluidPage(
  title = "Examples of DataTables",
  fluidRow(
    column(6,DT::dataTableOutput('x1')),
    column(6, plotOutput('plot1'))
  )
)