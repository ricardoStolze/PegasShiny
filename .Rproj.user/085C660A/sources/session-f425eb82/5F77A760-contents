source('UISidePanel.R')
source('UIMainPanel.R')

# Define UI for application that draws a histogram
ui <- fluidPage(
  #using spinners, till server respondes
  waiter::use_waiter(),
  # Application title
  titlePanel("MasterThesis"),
  
  sidebarLayout(
    mySidePanel,
    myMainPanel
    
  )
)
