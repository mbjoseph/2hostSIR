require(deSolve)
require(shiny)
require(ggplot2)

shinyUI(pageWithSidebar(
  headerPanel("Two host SIR model"),
  
  sidebarPanel(
    selectInput("tmode", "Transmission function",
                list("Density-dependent" = "dens",
                     "Frequency-dependent" = "freq")),
    sliderInput("b1", "Birth rate, species 1",
                min=0.01, max=10, value=1, step=.1),
    sliderInput("b2", "Birth rate, species 2",
                min=0.01, max=10, value=1, step=.1),
    sliderInput("d1", "Death rate, species 1",
                min=0.01, max=10, value=.5, step=.1),
    sliderInput("d2", "Death rate, species 2",
                min=0.01, max=10, value=.5, step=.1),
    sliderInput("tmax", "Number of timesteps", 
                min = .5, max = 30, value = 3, step= .5),
    sliderInput("tint", "Timestep intervals", 
                min = 0.01, max = 1, value = .01, step= .01),
    sliderInput("beta11", "Intraspecific transmission rate: species 1",
                min=0, max=100, value=.05, step=.01),
    sliderInput("beta22", "Intraspecific transmission rate: species 2",
                min=0, max=100, value=.05, step=.01),
    sliderInput("cij21", "Strength of transmission from species 2 to 1",
                min=0, max=1, value=.05, step=.01),
    sliderInput("cij12", "Strength of transmission from species 1 to 2",
                min=0, max=1, value=.05, step=.01),
    sliderInput("m", "Virulence",
                min=0, max=50, value=1.5, step=.5)
  ),
  
  # Two panels for different types of plots
  mainPanel(
    tabsetPanel(
       tabPanel("Solutions through time", plotOutput("p1")),
       tabPanel("Phase portraits", plotOutput("p2"))
      )
  )
))