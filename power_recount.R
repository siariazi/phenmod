library(deSolve)
library(optimx)
library(ggplot2)
#save.image("linearfits.RData")
#load("linearfits.RData")
setwd("C:/Siavash/Codes/phenmod/plots_data")

growth_viability <- read.csv("formaldehyde_growth_experiments_Dec2015_CFUs_final_recalc160923.csv")
growth_nash <- read.csv("formaldehyde_growth_experiments_Dec2015_NashAssay.csv")
death_viability <- read.csv("killingcurve_8615_results_recalc_160923_CFUs_final_recalc160923.csv")
death_nash <- read.csv("RECOUNT_killingcurve_8615_alldata_final_nash_final_data2.csv")

tmp1 <- merge(growth_viability,growth_nash)
tmp1 <-tmp1[tmp1$substrate !="succinate",]
tmp1 <-tmp1[tmp1$strain !="CM3745",]
growth_total <- subset(tmp1, select = -c(substrate,strain))
mean_growth <- mean(growth_total[growth_total$time==0,]$cells)
#growth_total[growth_total$time == 0,]$cells <- mean_growth
growth_total$cells[growth_total$time == 0] <- mean_growth

tmp2 <- merge(death_viability,death_nash, all.x = TRUE, all.y = TRUE)
tmp2 <-tmp2[tmp2$strain !="CM3745",]
death_total <- subset(tmp2, select = -c(strain,NashA412,standard_curve))
mean_death <- mean(death_total[death_total$time==0,]$cells)
death_total$cells[death_total$time==0] <- mean_death
colnames(death_total)[5] <- "formaldehyde_mM"

combined <- rbind(growth_total,death_total)
form.dat3 <- combined
form.dat3 <- form.dat3[order(form.dat3$time),]
names(form.dat3) <- c("flask", "time", "conc", "cells", "nash")
form.dat3$nash <- as.numeric(form.dat3$nash)
#form.dat3$cells[form.dat3$cells == 0.0] <- 1
form.dat3$cells <- form.dat3$cells + 1
form.dat3[form.dat3$conc>4,]$time <- form.dat3[form.dat3$conc>4,]$time/60 

deltaefgA<-function(t,state,parms) {
  with(as.list(c(state,parms)), {
    #rate of change
    dfex = -theta*vv*(fex-fin)*n
    dm = - vv*n*(vmax*m/(km+m))
    dfin =  theta*(fex-fin) + (vmax*m/(km+m)) - (vmax*gamma)*(fin/(kf+fin)) 
    dn = r*n*(fin/(kf+fin)) - n*alpha*(fin^2) 
    #return the rate of change
    list(c(dfex,dm,dfin,dn))
  }) 
}

fixedparms <- c(
  r = 0.23,          # growth rate from Jannell data (1/h)
  theta = 9e+04,     # diffusion coefficient (1/h)
  km = 2e-02,        # km for MDH 
  kf = 0.04,          # km for Fin, I started with 0.5  
  vv = 2e-13,        # volume of the cell to the volume of the media 
  vmax = 9e+04
)

IHS <- function(x, theta){  # IHS transformation
  (1/theta)*asinh(theta * x)
}



################ make two different functions for the residual, one for WT and one for dEfgA
## first for delta EfgA
all.resids <- NULL
resid.fun.dEfgA <- function(parms){
  #### Need to solve the DEs for each strain and concentration
  for(theconc in unique(form.dat3$conc)){
    dat <- subset(form.dat3, conc==theconc & time < 1000 & cells > 0)
    times <- unique(dat$time)
    state <- c(fex = theconc, m = 15, fin = 0, n = dat[dat$time==0,]$cells)
    working.parms <- c(fixedparms, 
                       alpha = parms[1],
                       gamma = parms[2]           # half saturation for negative hill function for growth 
    )
    out <- as.data.frame(ode(y = state, times = times, func = deltaefgA, 
                             parms= working.parms, method="lsoda", rtol=10^-10, atol=10^-10))
    out <- zapsmall(out)
    #resid_cells <- sum(((out$n - dat$n)/dat$n)^2, na.rm=T)
    #resid_cells <- sum((out$n - dat$n)^2, na.rm=T)
    #resid_cells <- sum((log(out$n)-log(dat$n))^2, na.rm = T)
    resid_cells <- sum((IHS(out$n, 1) - IHS(dat$cells, 1))^2)
    #resid_nash <- sum(((out$fex - dat$nash)/dat$nash)^2, na.rm=T)
    resid <- resid_cells #+ resid_nash
    
    all.resids <- rbind(all.resids, resid) 
  }
  sum(all.resids) + 10^5*sum(parms[working.parms<0])^2
}

#### Need to find reasonable starting parameters..
guess.parms.dEfgA <- c(
  0.002213925,    
  1.47 
)
bestparms.dEfgA <- optimx(guess.parms.dEfgA, resid.fun.dEfgA, method="Nelder-Mead")

###########################################
#### check out the fits.
# solving ODEs now with the bestparms.dEfgA
parms.dEfgA <- c(
  r = 0.23,          # growth rate from Jannell data (1/h)
  vmax = 9e+04,
  theta = 9e+04,     # diffusion coefficient (1/h)
  km = 2e-02,        # km for MDH 
  kf = 0.04,          # km for Fin, I started with 0.5  
  vv = 2e-13,        # volume of the cell to the volume of the media 
  alpha = bestparms.dEfgA$p1,
  gamma = bestparms.dEfgA$p2
)

form.dat3$n.model <- NA
outies <- NULL
for(theconc in unique(form.dat3$conc)){
  dat <- subset(form.dat3, conc==theconc & time < 1000)
  times <- dat$time
  state <- c(fex = theconc, m = 15, fin = 0, n = dat[dat$time==0,]$cells)
  out <- as.data.frame(ode(y = state, times = times, func = deltaefgA, 
                           parms= parms.dEfgA, method="lsoda"))
  form.dat3[form.dat3$conc==theconc &  form.dat3$time < 1000,]$n.model <- out$n
  out$conc <- theconc
  outies <- rbind(outies, out)
}

# ploting 

par(mfrow=c(2,3))
for(theconc in unique(form.dat3$conc)){
  tmp <- subset(form.dat3, conc==theconc & time < 1000)
  plot(cells+1~time, form.dat3, subset= time < 250, xlab="time", ylab="WT", log="y", col="white")
  title(main = c("conc:",theconc))
  lines(cells+1~time, tmp, type="b")
  lines(n.model+1~time, tmp, type="b", col="red")
  par(ask=TRUE)
}

# fex 

form.dat3$fex.model <- NA
outies <- NULL
for(theconc in unique(form.dat3$conc)){
  dat <- subset(form.dat3, conc==theconc & time < 1000)
  times <- dat$time
  state <- c(fex = theconc, m = 15, fin = 0, n = dat[dat$time==0,]$cells)
  out <- as.data.frame(ode(y = state, times = times, func = deltaefgA, 
                           parms= parms.dEfgA, method="lsoda"))
  form.dat3[form.dat3$conc==theconc &  form.dat3$time < 1000,]$fex.model <- out$fex
  out$conc <- theconc
  outies <- rbind(outies, out)
}

# ploting 

par(mfrow=c(2,3))
for(theconc in unique(form.dat3$conc)){
  tmp <- subset(form.dat3, conc==theconc & time < 1000)
  plot(nash~time, form.dat3, subset= time < 250, xlab="time", ylab="WT", log="y", col="white")
  title(main = c("conc:",theconc))
  lines(nash~time, tmp, type="b")
  lines(fex.model~time, tmp, type="b", col="red")
  par(ask=TRUE)
}
form.dat3$n.model <- NA

outies <- NULL
for(theconc in unique(form.dat3$conc)){
  dat <- subset(form.dat3, conc==theconc & time < 1000)
  times <- seq(0, 180, 0.5)
  state <- c(fex = theconc, m = 15, fin = 0, n = 10^7)
  out <- as.data.frame(ode(y = state, times = times, func = deltaefgA, 
                           parms= parms.dEfgA, method="lsoda"))
  out$conc <- theconc
  outies <- rbind(outies, out)
}

