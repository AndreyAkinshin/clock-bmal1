library(shiny)

slider <- function(title, min, max, value = min) {
  sliderInput(title, title, min, max, value, 0.01)
}
sliderS <- function(title, value = 3) slider(title, 1, 6, value)
sliderA <- function(title, value = 1) slider(title, 0, 20, value)
sliderK <- function(title, value = 1) slider(title, 0, 20, value)
sliderX <- function(title, value = 1) slider(title, 0, 20, value)

shinyUI(fluidPage(

  titlePanel("CLOCK:BMAL1 Circadian Oscillator Model 6D"),

  sidebarLayout(
    sidebarPanel(
      div(style = "overflow-y:scroll; max-height: 600px; position:relative;",
      h2("Parameters"),

      sliderA("A1", 1),
      sliderA("B1"),
      sliderS("s1"),
      sliderA("A2", 5),
      sliderA("B2"),
      sliderS("s2"),
      sliderA("A3", 10),
      sliderA("B3"),
      sliderS("s3"),
      sliderA("A4", 10),
      sliderA("B4"),
      sliderS("s4"),
      sliderA("A5", 10),
      sliderA("B5"),
      sliderS("s5"),
      
      sliderA("a1"),
      sliderA("b1"),
      sliderS("sigma1"),
      sliderA("a2"),
      sliderA("b2"),
      sliderS("sigma2"),
      sliderA("a3"),
      sliderA("b3"),
      sliderS("sigma3"),
      sliderA("a4"),
      sliderA("b4"),
      sliderS("sigma4"),
      sliderA("a6"),
      sliderA("b6"),
      sliderS("sigma6"),
      
      sliderK("k1"),
      sliderK("k2"),
      sliderK("k3"),
      sliderK("k4"),
      sliderK("k5"),
      sliderK("k6"),
      sliderK("C"),
      hr(),
      
      h2("Initial condition"),
      sliderX("p0"),
      sliderX("u0"),
      sliderX("w0"),
      sliderX("z0"),
      sliderX("x0"),
      sliderX("b0"),
      sliderInput("T", "T (total simulation time)", 1, 200, 20, 0.01)),
      
      br(),
      br(),
      downloadButton("save_state", "Save parameters to file"),
      br(),
      br(),
      fileInput("restore_state", "Load parameters from file",
                placeholder = ".yaml file")
    ),

    mainPanel(
      plotOutput("timelinePlot"),
      plotOutput("projPlot"),
      htmlOutput("htmlOutput")
    )
  )
))

