require(shiny)
require(deSolve) 
require(ggplot2)


shinyServer(function(input, output) {
  
  SIRsys <- function(parms, times){
    derivs <- function(times, state, parms) {
      with(as.list(c(state, parms)), {
        N1 <- S1 + I1 + R1
        N2 <- S2 + I2 + R2
        
        dS1 <- N1 * (b[1] - dd[1] * N1) -
              d[1] * S1 -
              S1 * (beta[1, 1] * I1 + beta[1, 2] * I2)* ifelse(tmode=="freq", 1/(N1+N2), 1)
        dI1 <- S1 * (beta[1, 1] * I1 + beta[1, 2] * I2)* ifelse(tmode=="freq", 1/(N1+N2), 1) -
              I1 * (d[1] + alpha[1] + sigma[1])
        dR1 <- sigma[1] * I1 - d[1] * R1
        
        dS2 <- N2 * (b[2] - dd[2] * N2) -
              d[2] * S2 -
              S2 * (beta[2, 2] * I2 + beta[2, 1] * I1)* ifelse(tmode=="freq", 1/(N1+N2), 1)
        dI2 <- S2 * (beta[2, 2] * I2 + beta[2, 1] * I1)* ifelse(tmode=="freq", 1/(N1+N2), 1) -
              I2 * (d[2] + alpha[2] + sigma[2])
        dR2 <- sigma[2]*I2 - d[2]*R2
        
        return(list(c(dS1, dI1, dR1, dS2, dI2, dR2)))
      }
      )
    }
    
    b <- c(input$b1, input$b2)
    dd <- c(.1, .1)
    d <- c(input$d1, input$d2)
    cij21 <- input$cij21
    cij12 <- input$cij12
    beta <- array(dim=c(2, 2))
    beta[1, 1] <- input$beta11
    beta[2, 2] <- input$beta22
    beta[1, 2] <- cij21 * mean(diag(beta))
    beta[2, 1] <- cij12 * mean(diag(beta))
    alpha <- c((input$m - 1)/d[1], (input$m - 1) / d[2])
    sigma <- c(5, 1)
    tmode <- input$tmode
    parms <- list(b, dd, d, beta, alpha, sigma, tmode)
    state <- with(parms, 
                  c(S1 = (b[1]-d[1])/dd[1], I1 = 1, R1 = 0, 
                    S2 = (b[2] - d[2]) / dd[2], I2 = 0, R2 = 0))
    return(ode(state, times, derivs, parms))
  }
  
output$p1 <- reactivePlot(function() {
    times <- seq(0, input$tmax, by = input$tint)
    res <- SIRsys(parms, times)
    d <- data.frame(timesteps = rep(times, 6), 
                    classes=gl(6, length(times), 
                               labels=c("Susceptible: species 1", "Susceptible: species 2", 
                               "Infectious: species 1", "Infectious: species 2", 
                               "Recovered: species 1", "Recovered: species 2")),
                    nums=c(res[,2], res[,5], res[,3], res[,6], res[,4], res[,7]))
    p <- ggplot(d) + geom_line(aes(x=timesteps, y=nums)) + 
      facet_wrap(~classes, scales="free", as.table="T", ncol=2) +
      xlab("Time") + ylab("Number of individuals") + theme_bw()
    print(p)
}
)
  
output$p2 <- reactivePlot(function(){
    times <- seq(0, input$tmax, by = input$tint)
    res <- SIRsys(parms, times)
    par(mfrow=c(1, 3))
    plot(res[,2], res[,5], type="l", 
         xlab="Susceptible individuals: species 1", 
         ylab="Susceptible individuals: species 2")
    s <- seq(length(res[,1])-1)# one shorter than data
    options(warn=-1)
    arrows(res[s, 2], res[s, 5], res[s+1, 2], res[s+1, 5], length=.05)
    plot(res[,3], res[,6], type="l",
         xlab="Infectious individuals: species 1", 
         ylab="Infectious individuals: species 2")
    arrows(res[s, 3], res[s, 6], res[s+1, 3], res[s+1, 6], length=.05)
    plot(res[,4], res[,7], type="l", 
         xlab="Recovered individuals: species 1", 
         ylab="Recovered individuals: species 2")
    arrows(res[s, 4], res[s, 7], res[s+1, 4], res[s+1, 7], length=.05)
}
)
})
