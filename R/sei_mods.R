# models
SEI.rabies <- function(nu = 0.52, mu = 0.44, R0 = 1.2, sigma = (1/22.5*365), alpha = 1, 
                       gamma = (1/3.1*365), K = 20, START.S.prop = 0.99, START.I.prop = 0.005, START.E.prop = 0.005,
                       START.R.prop = 0, START.S.density = 15, START.I.density = 0.15, START.E.density = 0.00,
                       START.R.density = 0, years = 20, steps = 1/52, model = "frequency",...){
  
  ### Deterministic skeleton of models for rabies
  
  require(deSolve)
  
  times <- seq(from = 0, to = years, by = steps)
  
  if (model == "frequency"){
    
    START.N.prop <- START.S.prop + START.I.prop + START.E.prop
    
    beta = (R0*gamma*(sigma + mu))/sigma
      
      # This function models a time step for the SIR:
      dx.dt.SIR <- function(t, y, parms) {
        N <- y["S"] + y["I"] + y["E"] 
        S <- y["S"]
        I <- y["I"]
        E <- y["E"]
        R <- y["R"]
        
        # Calculate the change in Susceptible
        dS <-  parms["nu"]*(S+E) - parms["mu"]*S - parms["beta"]*S*(I^parms["alpha"])/N
        
        #Calculate the change in Exposed
        dE <- parms["beta"]*S*(I^parms["alpha"])/N - parms["mu"]*E - parms["sigma"]*E
        
        # Calculate the change in Infected
        dI <- parms["sigma"]*E - parms["gamma"]*I
        
        # keeping track of removed
        dR <- parms["gamma"]*I
        
        # Return a list with the changes in S, E, I, R at the current time step
        return(list(c(dS, dI, dE, dR)))
      }
      
      # Create the parameter vector
      parms <- c(nu = nu, mu = mu, beta = beta, alpha = alpha, sigma = sigma, gamma = gamma)
      inits <- c(S=START.S.prop, I=START.I.prop, E=START.E.prop, R=START.R.prop)
  }
  
  if (model == "density"){
    
    START.N.density <- START.S.density + START.I.density + START.E.density
    
    # the old school rabies model
    # This function models a time step for the SIR:
    dx.dt.SIR <- function(t, y, parms) {
      N <- y["S"] + y["I"] + y["E"]
      S <- y["S"]
      I <- y["I"]
      E <- y["E"]
      R <- y["R"]
      
      ## Calculate dmort and beta here from parms and plug them in
      dmort <- (parms["nu"] - parms["mu"])/parms["K"]
      beta <- (parms["R0"]*(parms["sigma"] + parms["nu"])*(parms["gamma"] + parms ["nu"]))/(parms["sigma"]*parms["K"])
      
      # Calculate the change in Susceptible
      dS <-  parms["nu"]*N - parms["mu"]*S - dmort*S - beta*S*I
      
      #Calculate the change in Exposed
      dE <- beta*S*I - parms["mu"]*E - dmort*N*E - parms["sigma"]*E
      
      # Calculate the change in Infected
      dI <- parms["sigma"]*E - parms["mu"]*I - dmort*N*I - parms["gamma"]*I
      
      # Calculate the removed pop
      dR <- parms["gamma"]*I
      
      # Return a list with the changes in S, E, I, R at the current time step
      return(list(c(dS, dI, dE, dR)))
    }
    
    # Create the parameter vector
    parms <- c(nu = nu, mu = mu, R0 = R0, K = K, sigma = sigma, gamma = gamma)
    inits <- c(S=START.S.density, I=START.I.density, E=START.E.density, R=START.R.density)
  }
  
  # Run the ODE solver
  SIR.output <- lsoda(y = inits, 
                      times = times, 
                      func = dx.dt.SIR, 
                      parms = parms)
  SIR.output <- as.data.frame(SIR.output)
  SIR.output$N <- SIR.output$S + SIR.output$E + SIR.output$I
  return (SIR.output)
}



## Frequency dependent
ddt_1.01 <- SEI.rabies(R0 = 1.01, years = 20, model = "density")
ddt_1.01 <- data.frame(infected = unname(tapply(ddt_1.01$I, (seq_along(ddt_1.01$I)-1) %/% 4, sum)[1:260]),
                       pop = ddt_1.01$N[seq(1, 20*52, 4)], R0 = 1.01, trans = "Density")
ddt_1.1 <- SEI.rabies(R0 = 1.1, years = 20, model = "density")
ddt_1.1 <- data.frame(infected = unname(tapply(ddt_1.1$I, (seq_along(ddt_1.1$I)-1) %/% 4, sum)[1:260]),
                       pop = ddt_1.1$N[seq(1, 20*52, 4)], R0 = 1.1, trans = "Density")
ddt_1.05 <- SEI.rabies(R0 = 1.05, years = 20, model = "density")
ddt_1.05 <- data.frame(infected = unname(tapply(ddt_1.05$I, (seq_along(ddt_1.05$I)-1) %/% 4, sum)[1:260]),
                      pop = ddt_1.05$N[seq(1, 20*52, 4)], R0 = 1.05, trans = "Density")
ddt <- bind_rows(ddt_1.01, ddt_1.1, ddt_1.05)
ddt$time <- 1:260
ggplot(data = ddt, aes(x = time, y = pop, group = R0, color = R0)) + geom_line()
ggplot(data = ddt, aes(x = time, y = infected/pop, group = R0, color = as.factor(R0))) + 
  geom_line(size = 1) +
  scale_color_manual(values = c("lightcoral", "red", "darkred"), name = "R0") +
  xlab("Months") +
  ylab("Incidence (monthly proportion \n of population infected)") +
  newtheme +
  labs(tag = "B", subtitle = "Density")


fdt_1.01 <- SEI.rabies(R0 = 1.01, years = 20, model = "frequency")
fdt_1.01 <- data.frame(infected = unname(tapply(fdt_1.01$I, (seq_along(fdt_1.01$I)-1) %/% 4, sum)[1:260]),
                       pop = fdt_1.01$N[seq(1, 20*52, 4)], R0 = 1.01, trans = "Frequency")
fdt_1.1 <- SEI.rabies(R0 = 1.1, years = 20, model = "frequency")
fdt_1.1 <- data.frame(infected = unname(tapply(fdt_1.1$I, (seq_along(fdt_1.1$I)-1) %/% 4, sum)[1:260]),
                      pop = fdt_1.1$N[seq(1, 20*52, 4)], R0 = 1.1, trans = "Frequency")
fdt_1.05 <- SEI.rabies(R0 = 1.05, years = 20, model = "frequency")
fdt_1.05 <- data.frame(infected = unname(tapply(fdt_1.05$I, (seq_along(fdt_1.05$I)-1) %/% 4, sum)[1:260]),
                       pop = fdt_1.05$N[seq(1, 20*52, 4)], R0 = 1.05, trans = "Frequency")
fdt <- bind_rows(fdt_1.01, fdt_1.1, fdt_1.05)
fdt$time <- 1:260
ggplot(data = fdt, aes(x = time, y = pop, group = R0,color = as.factor(R0))) + geom_line() +
  scale_color_manual(values = c("lightcoral", "red", "darkred"), name = "R0") +
  xlab("Months") +
  ylab("Population") +
  newtheme

ggplot(data = fdt, aes(x = time, y = infected/pop, group = R0, color = as.factor(R0))) + 
  geom_line(size = 1) +
  scale_color_manual(values = c("lightcoral", "red", "darkred"), name = "R0") +
  xlab("Months") +
  ylab("Incidence (monthly proportion \n of population infected)") +
  newtheme +
  labs(tag = "A", subtitle = "Frequency")
