#Code NPZ model: modeling predator and prey interaction under eutrophication
library(deSolve)                                  

#autonomous system
NPZmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    function_response <- as.integer(function_response)
    if (function_response == 2){
        G = Im * P*Z/(P+Kp)
    }else if(function_response == 3){
        G = Im * P**2*Z/(P**2+Kp**2)
    }else{
       stop('Functional response incorrect!')
    }

    uptake= mumax * N/(Kn+N)*P

    zoo_mortality <- as.integer(zoo_mortality)
    if (zoo_mortality == 1){
        Zmort = mz*Z
    }else if (zoo_mortality == 2){
        Zmort = mz*Z*Z
    }else{
       stop('Zooplankton mortality incorrect!')
    }
    dN <- (1-eps)*G - uptake +Zmort
    dP <- uptake - G
    dZ <- eps*G - Zmort
    return(list(c(dN, dP, dZ)))
  })
}

pars  <- c(Im  = 2,    # /day, maximal rate of ingestion
           Kp  = 1.0,    # Half saturation constant of zoo
           eps = 0.3,    # dimensionless, growth efficiency of zoo
           mumax= 1,    # /day, maximal growth rate of phytoplankton
           mz = 0.05,   # /day, mortality rate of predator
           Kn = 0.5,    # Half saturation constant of phytoplankton
           zoo_mortality = 1, #linear
           function_response = 2)  #Holling type 2

#Equilibrium point (TODO):
#Solve the nonlinear system using Package ‘nleqslv’
#library(nleqslv)
#equ <- function(pars){
#  with(as.list(c(pars)), {
#    Preystar = rMort/(rIng*assEff)
#    Predatorstar = rGrow/rIng
#    return(c(Prey=Preystar, Predator=Predatorstar))
#  })
#}

#When to generate output (5 years, per day)
times  <- seq(0, 360*5, by = 1)

#Generate bifurication diagram
N0     <- seq(0.1, 2, 0.01)
Ntrial <- length(N0)  # different initial nutrient conditions
OUT    <- array(NA, c(Ntrial, length(times), 4))

for (k in 1:Ntrial){
    #Initial conditions
    yini  <- c(N = N0[k] - 0.011, P = 0.01, Z = 1e-3)
    
    #Main work: integration using ode
    out   <- ode(yini, times, NPZmod, pars)
    
    #Transform to a dataframe
    OUT[k,,]   <- out
}

#Plot out
pdf('NPZ_timeseries_Holling2_zmort_linear.pdf', width = 6, height = 4*2)

op <- par(font.lab = 1,
           family  = "serif",
           cex.lab = 1.2,
           cex.axis= 1.2,
           mgp     = c(2.2,1,0),
           mar     = c(4,3.5,3,1),
           mfrow   = c(Ntrial,1))

#Plot time series of nutrient, phyto, and zoo
for (i in c(1, Ntrial)){
    out <- as.data.frame(OUT[i,,])

    #Add column name
    colnames(out) <- c('time','N','P','Z')

    #Adjust Y scale
    plot(out$time, out$N,ylim=c(0, N0[i]),
        type='l', xlab='Days', ylab='Biomass') 
    points(out$time, out$P, type='l', col=2)
    points(out$time, out$Z, type='l', col=3)
    
    #Plot legend only once
    if(i==1) legend('topright', c('Nutrient','Phyto', 'Zoo'), col=1:3, lty=1)
}
dev.off()

pdf('NPZPhase.pdf', width=4, height=4)
op <- par(font.lab = 1,
           family  = "serif",
           cex.lab = 1.2,
           cex.axis= 1.2,
           mgp     = c(2.2,1,0),
           mar     = c(4,3.5,3,1),
           mfrow   = c(1,1))

for (i in c(1,Ntrial)){
    out <- as.data.frame(OUT[i,,])
    colnames(out) <- c('time','N','P','Z')
    if (i == 1){
      plot(out$P,out$Z,
        xlim=c(0, max(N0)), 
        ylim=c(0, max(N0)), type='l', xlab='Phyto', ylab='Zoo') 
    }else{
      points(out$P,out$Z,  type='l', col=i) 
    }
}
legend('topright', c('Low nutrient','High nutrient'), col=1:Ntrial, lty=1)
dev.off()

#Generate bifurication diagram
pdf('NPZ_bifurcation.pdf', width=5, height=5)
op <- par(font.lab = 1,
           family  = "serif",
           cex.lab = 1.2,
           cex.axis= 1.2,
           mgp     = c(2.2,1,0),
           mar     = c(4,3.5,3,1),
           mfrow   = c(1,1))


plot(N0, N0, 
    ylim = c(0, max(N0)),
    main = 'Bifurcation diagram for NPZ model',
    xlab = 'Total nitrogen', ylab = 'Phytoplankton',  type = 'n')
#Take the results from the last year
for (i in 1:Ntrial){
    out <- as.data.frame(OUT[i,,])
    colnames(out) <- c('time','N','P','Z')
    out <- out[361:nrow(out),]
    points(rep(N0[i],length(out$P)), out$P, pch=16, cex=.5) 
}
dev.off()

#Test different formulations of zooplankton functional response and mortality closure term
#Initial conditions
N0     <- seq(0.1, 2, 0.01)
Ntrial <- length(N0)  # different initial nutrient conditions
OUT    <- array(NA, c(Ntrial, 4, length(times),4))

for (k in 1:Ntrial){
    #Initial conditions
    yini  <- c(N = N0[k] - 0.011, P = 0.01, Z = 1e-3)
    
    #Main work: integration using ode
    for (m in 1:2){
        for (n in 1:2){
            #If change to Holling type 3
            pars$function_response <- n+1
            pars$zoo_mortality     <- m

            out   <- ode(yini, times, NPZmod, pars)
            
            #Transform to a dataframe
            OUT[k,n+(m-1)*2,,]   <- out
        }
    }
}

#Plot out as an example
for (m in 1:2){
    for (n in 1:2){
      filename <- paste0('NPZ_timeseries_Holling',n+1, 'zmort',m,'.pdf')
      pdf(filename, width = 4, height = 4)
      op <- par(font.lab = 1,
                 family  = "serif",
                 cex.lab = 1.2,
                 cex.axis= 1.2,
                 mgp     = c(2.2,1,0),
                 mar     = c(4,3.5,3,1),
                 mfrow   = c(1,1))
      
      #Plot time series of nutrient, phyto, and zoo
      
      #Add column name
      out <- as.data.frame(OUT[Ntrial,n+(m-1)*2,,])
      colnames(out) <- c('time','N','P','Z')
      
      #Adjust Y scale
      plot(out$time, out$N,ylim=c(0, sum(yini)+.5),
          type='l', xlab='Days', ylab='Biomass') 
      points(out$time, out$P, type='l', col=2)
      points(out$time, out$Z, type='l', col=3)
      
      #Plot legend only once
      legend('topright', c('Nutrient','Phyto', 'Zoo'), col=1:3, lty=1)
      
      dev.off()
    }
}

#Generate bifurication diagram
for (m in 1:2){
    for (n in 1:2){
      filename <- paste0('NPZ_bifurcation_Holling',n+1, 'zmort',m,'.pdf')

      pdf(filename, width = 5, height = 5)
      plot(N0, N0, 
          ylim = c(0, max(N0)),
          xlab = 'Total nitrogen', ylab = 'Phytoplankton',  type = 'n')
      #Take the results from the last year
      for (i in 1:Ntrial){
          out <- as.data.frame(OUT[i,n+(m-1)*2,,])
          colnames(out) <- c('time','N','P','Z')
          out <- out[361:nrow(out),]
          points(rep(N0[i],length(out$P)), out$P, pch=16, cex=.2) 
      }
      dev.off()
    }  
}

#Check the relationships between N, P, Z and total N
n <- 2
for (m in 1:2){
    NPZstar  <- OUT[,n+(m-1)*2, length(times),]
    filename <- paste0('NPZ_nutrient_Holling3_zooMort_',m,'.pdf')
    pdf(filename, width=5, height=5)
   
    plot(N0, N0, type = 'n', ylim=c(0, 1.4),
        xlab='Total nitrogen', ylab='N, P, Z')
    for (i in 1:3){
        points(N0, NPZstar[, i+1], col=i, type='l')
    }
    legend('topleft', c('NO3', 'PHY', 'ZOO'), lty=1,  col=1:3)
    dev.off()
}