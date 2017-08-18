# Rely on the 'WorldPhones' dataset in the datasets
# package (which generally comes preloaded).
library(datasets)

asca.table <- read.csv("backend/data/asca_sigAB.csv")

# Use a fluid Bootstrap layout
fluidPage(    
  
  # Give the page a title
  titlePanel("ASCA results"),
  
  # Generate a row with a sidebar
  sidebarLayout(      
    
    # Define the sidebar with one input
    sidebarPanel(
      selectInput("mz", "Mass/charge value:", 
                  choices=asca.table$X),
      hr(),
      helpText("...")
    ),
    
    # Create a spot for the barplot
    mainPanel(
      plotOutput("mzPlot")  
    )
    
  )
)