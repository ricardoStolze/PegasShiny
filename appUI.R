source('UISidePanel.R')
source('UIMainPanel.R')

# Define UI for application that draws a histogram
ui <- fluidPage(
  #using spinners, till server respondes
  waiter::use_waiter(),
  
  shinyFeedback::useShinyFeedback(),
  tags$style(HTML(".custom-h5 {
                  color: gray;}
                  ")),
  
  # Application title
  titlePanel("ShinyPegas"),
  
  sidebarLayout(
    mySidePanel,
    myMainPanel
    
  )
)
