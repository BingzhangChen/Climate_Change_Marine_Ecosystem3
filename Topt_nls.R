load('Merged_PHY.Rdata')
library(tidyverse)
library(tidymodels)

PHY_ag <- newdat %>%
  group_by(ID) %>%
  summarize(N_obs = n()) #Number of data points in each expt

#Pick up the experiment with the greatest N_obs
IDm <- PHY_ag[which.max(PHY_ag$N_obs), 'ID']

df <- newdat %>%
  filter(ID == as.character(IDm[[1]]))

#Function converting temperature in C to Boltzmann temperature
T.K <- function(x, Tref = 15,kb=8.62E-5) {
  1/((Tref+273)*kb) - 1/((x + 273)*kb)
}

#Nonlinear function of TPC
Johnson <-  function(temp, Topt, Lnmumax, LnEa, LnEh){
      T0 <- 273.15
      kb <- 8.62E-5
      Ea <- exp(LnEa)
      Eh <- exp(LnEh)
      Ed <- Eh - Ea
      b  <- T.K(temp) - T.K(Topt)   #Standardized temperature (ivT)
      h  <- (Ea/Ed + 1) * exp(Ea * b)/(1+ Ea/Ed * exp(Eh * b) )
      h  <- exp(Lnmumax)* h 
      return(h)
}

#Set nls setting
nls.control(maxiter = 2000, tol = 1e-03, 
    minFactor = 1/204800,
    printEval = FALSE, warnOnly = TRUE)

#Set initial parameters
mumin  <- 1e-2
muMAX  <- 10
Eamin  <- 1e-2
Eamax  <- 4
Ehmin  <- .1
Ehmax  <- 20
Tomin  <- -6
Tomax  <- 50
optT   <- df[which.max(df$Growth), 'Temperature']

#Run nls using port algorithm
z <- nls(
  Growth ~ Johnson(Temperature, Topt, Lnum, LnEa, LnEh),
  data = df,
  start = list(
    Lnum = max(log(df$Growth), na.rm = T),
    LnEa = log(0.6),
    LnEh = log(3),
    Topt = optT
  ),
  lower = list(
    Lnum = log(mumin),
    LnEa = log(Eamin),
    LnEh = log(Ehmin),
    Topt = Tomin
  ),
  upper = list(
    Lnum = log(muMAX),
    LnEa = log(Eamax),
    LnEh = log(Ehmax),
    Topt = Tomax
  ),
  algorithm = 'port'
)

#Extract the model coefficients
z %>% tidy()

#Plot
df %>%
  ggplot(aes(x = Temperature, y = Growth)) +
  geom_point(alpha = 1/2) +
  geom_smooth(method = 'nls',
              formula = y  ~ Johnson(x, Topt, Lnum, LnEa, LnEh),
              se = FALSE,
              method.args = list(
                                  start = list(
                                               Lnum = log(4),
                                               LnEa = log(0.6),
                                               LnEh = log(3),
                                               Topt = optT
                                             )
                                 )
  ) +
  xlab('Temperature') +
  ylab(expression(paste('Growth rate (' * d^-1 * ')'))) +
  theme_light()


