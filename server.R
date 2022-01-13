library(shiny)
library(ggplot2)
library(deSolve)
library(rlist)

source("core.R")

shinyServer(function(input, output, session) {

  output$timelinePlot <- renderPlot({
    res <- simulate(reactiveValuesToList(input))
    ggplot(res$timelineDf, aes(x, y, col = type)) +
      geom_line() +
      theme_bw() +
      labs(title = "Timeline plot", x = "time", y = "value", col = "variable")
  })
  
  output$projPlot <- renderPlot({
    res <- simulate(reactiveValuesToList(input))
    ggplot(res$traj, aes(e2, e3)) +
      geom_path() +
      geom_point(x = 0, y = 0) +
      theme_bw() +
      labs(title = "Projection plot", x = "Re(e2)", y = "Im(e2)")
  })
  
  output$htmlOutput <- renderUI({
    context <- reactiveValuesToList(input)
    res <- simulate(context)
    star <- res$star
    dstar <- as.list(res$model(0, res$star, NULL)[[1]])
    M <- res$M
    Ml <- res$Ml
    Me <- t(res$Me)
    rnd <- function(x) { if (is.complex(x) & abs(Im(x)) < 1e-9) round(Re(x), 6) else round(x, 6) }
    fmt <- function(x) { if (is.na(x)) "         ?" else formatC(x, digits = 6, width = 10, format = "f") }
    hr <- "-----------------------------------------------------------------<br>"
    text <- with(context, {
      paste0(
        "<pre>",
        "Stationary point:<br>",
        "p* = ", fmt(star$p), "   (p*' = ", fmt(dstar$dp), ")<br>",
        "u* = ", fmt(star$u), "   (u*' = ", fmt(dstar$du), ")<br>",
        "w* = ", fmt(star$w), "   (w*' = ", fmt(dstar$dw), ")<br>",
        "z* = ", fmt(star$z), "   (z*' = ", fmt(dstar$dz), ")<br>",
        "x* = ", fmt(star$x), "   (x*' = ", fmt(dstar$dx), ")<br>",
        "b* = ", fmt(star$b), "   (b*' = ", fmt(dstar$db), ")<br>",
        hr,
        "Linearization matrix:<br>",
        fmt(M[1,1]), " ", fmt(M[1,2]), " ", fmt(M[1,3]), " ", fmt(M[1,4]), " ", fmt(M[1,5]), " ", fmt(M[1,6]), "<br>",
        fmt(M[2,1]), " ", fmt(M[2,2]), " ", fmt(M[2,3]), " ", fmt(M[2,4]), " ", fmt(M[2,5]), " ", fmt(M[2,6]), "<br>",
        fmt(M[3,1]), " ", fmt(M[3,2]), " ", fmt(M[3,3]), " ", fmt(M[3,4]), " ", fmt(M[3,5]), " ", fmt(M[3,6]), "<br>",
        fmt(M[4,1]), " ", fmt(M[4,2]), " ", fmt(M[4,3]), " ", fmt(M[4,4]), " ", fmt(M[4,5]), " ", fmt(M[4,6]), "<br>",
        fmt(M[5,1]), " ", fmt(M[5,2]), " ", fmt(M[5,3]), " ", fmt(M[5,4]), " ", fmt(M[5,5]), " ", fmt(M[5,6]), "<br>",
        fmt(M[6,1]), " ", fmt(M[6,2]), " ", fmt(M[6,3]), " ", fmt(M[6,4]), " ", fmt(M[6,5]), " ", fmt(M[6,6]), "<br>",
        hr,
        "Eigenvalues:<br>",
        "λ1 = ", rnd(Ml[1]), "<br>",
        "λ2 = ", rnd(Ml[2]), "<br>",
        "λ3 = ", rnd(Ml[3]), "<br>",
        "λ4 = ", rnd(Ml[4]), "<br>",
        "λ5 = ", rnd(Ml[5]), "<br>",
        "λ6 = ", rnd(Ml[6]), "<br>",
        hr,
        "Eigenvectors:<br>",
        "e1 = (", rnd(Me[1,1]), ", ", rnd(Me[1,2]), ", ", rnd(Me[1,3]), ", ", rnd(Me[1,4]), ", ", rnd(Me[1,5]), ", ", rnd(Me[1,6]), ")", "<br>",
        "e2 = (", rnd(Me[2,1]), ", ", rnd(Me[2,2]), ", ", rnd(Me[2,3]), ", ", rnd(Me[2,4]), ", ", rnd(Me[2,5]), ", ", rnd(Me[2,6]), ")", "<br>",
        "e3 = (", rnd(Me[3,1]), ", ", rnd(Me[3,2]), ", ", rnd(Me[3,3]), ", ", rnd(Me[3,4]), ", ", rnd(Me[3,5]), ", ", rnd(Me[3,6]), ")", "<br>",
        "e4 = (", rnd(Me[4,1]), ", ", rnd(Me[4,2]), ", ", rnd(Me[4,3]), ", ", rnd(Me[4,4]), ", ", rnd(Me[4,5]), ", ", rnd(Me[4,6]), ")", "<br>",
        "e5 = (", rnd(Me[5,1]), ", ", rnd(Me[5,2]), ", ", rnd(Me[5,3]), ", ", rnd(Me[5,4]), ", ", rnd(Me[5,5]), ", ", rnd(Me[5,6]), ")", "<br>",
        "e6 = (", rnd(Me[6,1]), ", ", rnd(Me[6,2]), ", ", rnd(Me[6,3]), ", ", rnd(Me[6,4]), ", ", rnd(Me[6,5]), ", ", rnd(Me[6,6]), ")", "<br>",
        hr,
        hr,
        "Differential equations:<br>",
        "p' = k1 * (Г1(u) * γ1(w) - p)<br>",
        "u' = k2 * (Г2(x) * L2(p) - u)<br>",
        "w' = k3 * (Г3(x) * L3(p) - w)<br>",
        "z' = k4 * (Г4(x) * L4(p) - z)<br>",
        "x' = k5 * (Г5(b) - x)<br>",
        "b' = k6 * (C * L6(z) - b)<br>",
        hr,
        "Aux functions:<br>",
        "Г1(u) = (A1 * u<sup>s1</sup>) / (B1 + u<sup>s1</sup>)<br>",
        "Г2(x) = (A2 * x<sup>s2</sup>) / (B2 + x<sup>s2</sup>)<br>",
        "Г3(x) = (A3 * x<sup>s3</sup>) / (B3 + x<sup>s3</sup>)<br>",
        "Г4(x) = (A4 * x<sup>s4</sup>) / (B4 + x<sup>s4</sup>)<br>",
        "Г5(b) = (A5 * b<sup>s5</sup>) / (B5 + b<sup>s5</sup>)<br>",
        "γ1(w) = (a1 * w<sup>σ1</sup>) / (b1 + w<sup>σ1</sup>)<br>",
        "L2(p) = a2 / (b2 + p<sup>σ2</sup>)<br>",
        "L3(p) = a3 / (b3 + p<sup>σ3</sup>)<br>",
        "L4(p) = a4 / (b4 + p<sup>σ4</sup>)<br>",
        "L6(z) = a6 / (b6 + z<sup>σ6</sup>)<br>",
        hr,
        "Parameters:<br>",
        "A1 = ", A1, "; B1 = ", B1, "; s1 = ", s1, "; (Г1(u))<br>",
        "A2 = ", A2, "; B2 = ", B2, "; s2 = ", s2, "; (Г2(x))<br>",
        "A3 = ", A3, "; B3 = ", B3, "; s3 = ", s3, "; (Г3(x))<br>",
        "A4 = ", A4, "; B4 = ", B4, "; s4 = ", s4, "; (Г4(x))<br>",
        "A5 = ", A5, "; B5 = ", B5, "; s5 = ", s5, "; (Г5(b))<br>",
        "a1 = ", a1, "; b1 = ", b1, "; σ1 = ", sigma1, "; (γ1(w))<br>",
        "a2 = ", a2, "; b2 = ", b2, "; σ2 = ", sigma2, "; (L2(p))<br>",
        "a3 = ", a3, "; b3 = ", b3, "; σ3 = ", sigma3, "; (L3(p))<br>",
        "a4 = ", a4, "; b4 = ", b4, "; σ4 = ", sigma4, "; (L4(p))<br>",
        "a6 = ", a6, "; b6 = ", b6, "; σ6 = ", sigma6, "; (L6(z))<br>",
        "k1 = ", k1, "; k2 = ", k2, "; k3 = ", k3, "; k4 = ", k4, "; k5 = ", k5, "; k6 = ", k6, ";<br>",
        "C = ", C, ";<br>",
        "p0 = ", p0, "; u0 = ", u0, "; w0 = ", w0, "; z0 = ", z0, "; x0 = ", x0, "; b0 = ", b0, ";<br>",
        "T = ", T, "<br>",
        "</pre>"
        )
    })
    HTML(text)
  })
  
  output$save_state <- downloadHandler(
    filename = function() {
      paste0("data-", Sys.Date(), "-", randomCodeName(), ".yaml")
    },
    content = function(file) {
      data <- reactiveValuesToList(input)
      data <- data[names(data) != "restore_state"] 
      list.save(data, file)
    } 
  )
  
  loadedData <- reactive({
    list.load(input$restore_state$datapath)
  })
  observe({
    data <- loadedData()
    updateSliderInput(session, "A1", value = data$A1)
    updateSliderInput(session, "B1", value = data$B1)
    updateSliderInput(session, "s1", value = data$s1)
    updateSliderInput(session, "A2", value = data$A2)
    updateSliderInput(session, "B2", value = data$B2)
    updateSliderInput(session, "s2", value = data$s2)
    updateSliderInput(session, "A3", value = data$A3)
    updateSliderInput(session, "B3", value = data$B3)
    updateSliderInput(session, "s3", value = data$s3)
    updateSliderInput(session, "A4", value = data$A4)
    updateSliderInput(session, "B4", value = data$B4)
    updateSliderInput(session, "s4", value = data$s4)
    updateSliderInput(session, "A5", value = data$A5)
    updateSliderInput(session, "B5", value = data$B5)
    updateSliderInput(session, "s5", value = data$s5)
    
    updateSliderInput(session, "a1", value = data$a1)
    updateSliderInput(session, "b1", value = data$b1)
    updateSliderInput(session, "sigma1", value = data$sigma1)
    updateSliderInput(session, "a2", value = data$a2)
    updateSliderInput(session, "b2", value = data$b2)
    updateSliderInput(session, "sigma2", value = data$sigma2)
    updateSliderInput(session, "a3", value = data$a3)
    updateSliderInput(session, "b3", value = data$b3)
    updateSliderInput(session, "sigma3", value = data$sigma3)
    updateSliderInput(session, "a4", value = data$a4)
    updateSliderInput(session, "b4", value = data$b4)
    updateSliderInput(session, "sigma4", value = data$sigma4)
    updateSliderInput(session, "a6", value = data$a6)
    updateSliderInput(session, "b6", value = data$b6)
    updateSliderInput(session, "sigma6", value = data$sigma6)
    
    updateSliderInput(session, "k1", value = data$k1)
    updateSliderInput(session, "k2", value = data$k2)
    updateSliderInput(session, "k3", value = data$k3)
    updateSliderInput(session, "k4", value = data$k4)
    updateSliderInput(session, "k5", value = data$k5)
    updateSliderInput(session, "k6", value = data$k6)
    updateSliderInput(session, "C", value = data$C)
    
    updateSliderInput(session, "p0", value = data$p0)
    updateSliderInput(session, "u0", value = data$u0)
    updateSliderInput(session, "w0", value = data$w0)
    updateSliderInput(session, "z0", value = data$z0)
    updateSliderInput(session, "x0", value = data$x0)
    updateSliderInput(session, "b0", value = data$b0)
    updateSliderInput(session, "T", value = data$T)
    
  })
})

