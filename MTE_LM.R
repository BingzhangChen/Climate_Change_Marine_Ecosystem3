rm(list = ls() ) # Remove everything from memory
library(tidyverse)
library(tidymodels)
tidymodels_prefer()
load('Corrected_Metabolism.Rdata')

#Use base R
z.lm <- lm(log(Metabolism_W) ~ log(Mass_g), W)
summary(z.lm)

#Check residuals
par(mfrow = c(2,2))
plot(z.lm, page = 1)

#Run linear model using tidymodels
lm_model <- linear_reg() %>%
  set_engine("lm")

lm_form_fit <- lm_model %>%
  fit(log(Metabolism_W) ~ log(Mass_g), data = W)

lm_form_fit %>% extract_fit_engine()


#Fit to both Mass and size together
w <- read_csv('Metabolism.csv') %>%
  select(Major_taxa, Temperature_C, Mass_g, Metabolism_W)

w[w$Major_taxa == 'Mammal', 'Temperature_C'] <- 37
w[w$Major_taxa == 'Bird',   'Temperature_C'] <- 37

#Function converting temperature in C to Boltzmann temperature
T.K <- function(x, Tref = 15,kb=8.62E-5) {
  1/((Tref+273)*kb) - 1/((x + 273)*kb)
}

z.lm2 <- lm(log(Metabolism_W) ~ log(Mass_g) + T.K(Temperature_C), w)
summary(z.lm2)

