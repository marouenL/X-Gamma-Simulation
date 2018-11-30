library(shiny)
library(shinydashboard)
library(MASS)


gaussianVectorSimulation <- function(sensi, sigma1, sigma2, rho, simulationNb=1000){
  
  sigma <- diag(c(sigma1, sigma2))
  cor <- matrix(c(1,rho/100,rho/100,1), nrow=2, ncol=2, byrow = T)
  cov <- data.frame(sigma %*% cor %*% t(sigma))
  gaussianVector <- data.frame(mvrnorm(n = simulationNb, mu=rep(0,2), as.matrix(cov), tol = 1e-6, 
                                       empirical = FALSE, EISPACK = FALSE))
  
  plot(exp(gaussianVector[,1])-1, exp(gaussianVector[,2])-1, main = "Shocks Simulations",
       xlab = "RF1", ylab = "RF2")
  
}

simulateCrossXGamma<- function(sensi, sigma1, sigma2, rho, simulationNb=1000){

  sigma <- diag(c(sigma1, sigma2))
  cor <- matrix(c(1,rho/100,rho/100,1), nrow=2, ncol=2, byrow = T)
  cov <- data.frame(sigma %*% cor %*% t(sigma))
  gaussianVector <- data.frame(mvrnorm(n = simulationNb, mu=rep(0,2), as.matrix(cov), tol = 1e-6, 
                                       empirical = FALSE, EISPACK = FALSE))
  simulations <- cbind(rep(0,simulationNb),
                       sensi*(exp(gaussianVector[,1])-1)*(exp(gaussianVector[,2])-1)*100)
  
  maxValue <- max(abs(simulations))

  
  plot(c(0,1), simulations[1,], col=rgb(0,0,0,alpha = abs(simulations[1,2])/maxValue),
       ylim = c(-5000000,5000000), lty=3, main = "1-day Loss Simulations", xlab = "Time",
       ylab = "Loss")
  
  
  for (i in 2:simulationNb){
    
    lines(c(0,1), simulations[i,], col=rgb(0,0,0,alpha = abs(simulations[i,2])/maxValue), lty=3)
  }
  
  loss <- quantile(simulations[,2], 0.01)
  
  lines(c(0,1), c(loss, loss), col="red")
  
  euro <- "\u20AC"
  text(0, loss-200000, paste("99% VaR = ",toString(round(loss/1000000, digits = 2))
                             ," M", euro), col = "red", adj = 0)

  
}


ui <- fluidPage(

  titlePanel("X Gamma Simulation"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("sensi", "Sensitivity:", 6588479, 0, 100000000),
      numericInput("sigma1", "Sigma RF1:", 0.026, 0, 1, 0.01),
      numericInput("sigma2", "Sigma RF2:", 0.019, 0, 1, 0.01),
      sliderInput("correl", "Correlation:", -100, 100, -0.42) 
    ),
    mainPanel(
      plotOutput("plot1"),
      plotOutput("plot2")
    ) 
  )
    
 
)


server <- function(input, output) {
  output$plot1 <- renderPlot({
    simulateCrossXGamma(input$sensi, input$sigma1, input$sigma2, input$correl)
  })
  
  output$plot2 <- renderPlot({
    gaussianVectorSimulation(input$sensi, input$sigma1, input$sigma2, input$correl)
  })
}



shinyApp(ui, server)