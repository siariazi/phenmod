library(ggplot2)
setwd("C:/Siavash/EfgA/phenmod/plots_data")
library(deSolve)

vmax <- 9e+04        # kcat for methanol (1/h), was 7.2e+03, if I change it to 70000 I get regrowth 
r <- 0.23            # growth rate from Jannell data (1/h)
theta <- 90000       # diffusion coefficient (1/h) try 16000
km <- 2e-02          # km for MDH 
kf <- 0.04           # km for Fin, I started with 0.02 
vv <- 2e-13          # volume of the cell to the volume of the media 
gamma <- 2.30077     #1.85# 0.02,0.032 combined metabolic and maitenance use of formaldehyde try 1.4
alpha <- 0.03159032  #0.021 # 2,1.82 maxium killing rate of fin, metabolic independent (1/h) # starting with 0.32
sigma <- 2

parms <- c(     vmax = vmax,
                km = km,         
                kf = kf,         
                theta = theta,  
                alpha = alpha,       
                r = r, 
                vv = vv, 
                gamma = gamma, 
                sigma = sigma
)

h <- function(fin){alpha*(fin^sigma)}
times=seq(0,200,by=0.05)
efgAmodel<-function(t,state,parms) {
  with(as.list(c(state,parms)), {
    #rate of 5change
    dfex = -theta*vv*(fex-fin)*n
    dm = - vv*(vmax*m/(km+m))*n
    dfin =  theta*(fex-fin) + (vmax*m/(km+m)) - (vmax*gamma)*(fin/(kf+fin)) 
    dn = r*(fin/(kf+fin))*n - h(fin)*n
    #return the rate of change
    list(c(dfex,dm,dfin,dn))
  }) 
}
# here's a function to run the model on a range of different levels of initial external formaldehyde (fex),
# and output all the results into several summary tables
model_multi_levels<-function(model, time_range, baseline_state, parameters, fex_init_range){
  out<-data.frame() # create an empty dataframe for putting all of the outputs into
  for (i in fex_init_range) { # for-loop to run through all initial external formaldehyde levels
    state_i<-baseline_state # start off with your baseline state variables, of which all except fex will remain the same
    state_i[1]<-i # replace the fex from the baseline set of state variables with the fex value you want to model
    out_i <- as.data.frame(ode(y=state_i, times=time_range, func=model, parms=parameters, method="lsoda")) # run the model
    treatment<-rep(i, length(time_range)) # create a column called "fex_init" to append to the data tables, so you know which initial fex value yielded these outputs
    out<-rbind(out, cbind(out_i, treatment)) # append the results of this model run, along with a column indicating the initial formaldehyde concentration, to the master output dataframe
  }
  return(out)
}

times_growth=seq(0,70,by=0.05)
growth_fex_range<-c(0,2,3,4) # use one range of low fex levels for growth
#growth_fex_range<-seq(0,10,1) # use one range of low fex levels for growth
state_growth<-c(fex = 20, m = 15, fin = 0, n = 388500)
growth_out<-model_multi_levels(model=efgAmodel, time_range=times_growth, baseline_state=state_growth, parameters=parms, fex_init_range=growth_fex_range)

times_death=seq(0,25,by=0.05)
death_fex_range<-c(5, 7.5, 10, 12.5, 15, 20) # and a range of higher fex levels for death
state_death<-c(fex = 20, m = 15, fin = 0, n = 77166668)
death_out<-model_multi_levels(model=efgAmodel, time_range=times_death, baseline_state=state_death, parameters=parms, fex_init_range=death_fex_range)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

# alternatively, use the following color schemes for yellow-red or blue-green palettes:
YlOrRd<-c('#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#b10026')
PuBuGrn<-c('#bdc9e1','#67a9cf','#1c9099','#016c59')


# death counts
death_count<-read.csv('killingcurve_8615_results_recalc_160923_CFUs_final_recalc160923.csv', header=T) #from death curve 8
death_count<-subset(death_count, strain=='CM2730')
death_count<-count_data[complete.cases(death_count),] # this is to get rid of lines with "NA," because I Use "NA" to denote cases that I think are plating errors
death_count_data <- death_count[,-4]
colnames(death_count_data)[2:3] <- c("treatment","n")
death_count_data$"type" <- "data"
death_count_data$time <- death_count_data$time/60

death_count_model <- death_out[,c(1,5,6)]
death_count_model$"type" <- "model"
death_count_total <- rbind(death_count_data,death_count_model)

# death_nash
death_nash <- read.csv("killingcurve_8615_results_recalc_160923_nash_final_data.csv",header = T)
death_nash <- subset(death_nash,strain=="CM2730")
death_nash_data <- death_nash[,c(3,4,6)]
death_nash_data <- death_nash_data[complete.cases(death_nash_data),]
colnames(death_nash_data)[3] <- "fex"
death_nash_data$"type" <- "data"
death_nash_data$time <- death_nash_data$time/60

death_nash_model <- death_out[,c(1,2,6)]
death_nash_model$"type" <- "model"
death_nash_total <- rbind(death_nash_data,death_nash_model)

# growth counts
growth_count<-read.csv('formaldehyde_growth_experiments_Dec2015_CFUs_final_recalc160923.csv', header=T) # from my December formaldehyde growth curves
growth_count<-subset(growth_count, strain=='CM2730')
growth_count<-subset(growth_count, substrate=='methanol')
growth_count_data <- growth_count[,c(2,5,6)]
growth_count_data$"type" <- "data"
colnames(growth_count_data)[2:3] <- c("treatment","n")

growth_count_model <- growth_out[,c(1,5,6)]
growth_count_model$"type" <- "model"
growth_count_total <- rbind(growth_count_model,growth_count_data)
 
# growth_nash
growth_nash <- read.csv('formaldehyde_growth_experiments_Dec2015_NashAssay.csv',header = T)
growth_nash <- subset(growth_nash,strain=='CM2730')
growth_nash <- subset(growth_nash,substrate=='succinate')
growth_nash_data <- growth_nash[,c(3,5,6)]
colnames(growth_nash_data) <- c("treatment","time","fex")
growth_nash_data$"type" <- "data"
growth_nash_model <- growth_out[,c(1,2,6)]
growth_nash_model$"type" <- "model"

growth_nash_total <- rbind(growth_nash_data,growth_nash_model)

growth_count_total$treatment <- as.factor(growth_count_total$treatment)
plot_growth_count <- ggplot() + geom_line(data=growth_count_total, size=1.2,
                                       aes(x=time, y=n, color=treatment, group=treatment))
plot_growth_count <- plot_growth_count + theme_bw() + scale_color_manual(values=gg_color_hue(6), name='initial \n[formaldehyde] (mM)')
plot_growth_count <- plot_growth_count + scale_y_log10(limits=c(10, 1e+09), breaks=c(1e+02, 1e+04, 1e+06, 1e+08)) #+ scale_x_continuous(limits=c(0,3))
plot_growth_count <- plot_growth_count + xlab('time (hours)') + ylab('viable cells (CFU/mL)') + facet_grid(~type)
plot_growth_count

growth_nash_total$treatment <- as.factor(growth_nash_total$treatment)
plot_growth_nash <- ggplot() + geom_line(data = growth_nash_total,size=1.2,aes(x=time,y=fex,color=treatment,group=treatment))
plot_growth_nash <- plot_growth_nash + theme_bw() + scale_color_manual(values=gg_color_hue(6), name='initial \n[formaldehyde] (mM)')
plot_growth_nash <- plot_growth_nash + facet_grid(~type)
plot_growth_nash

death_count_total$treatment <- as.factor(death_count_total$treatment)
plot_death_count <- ggplot() + geom_line(data=death_count_total, size=1.2,
                                          aes(x=time, y=n, color=treatment, group=treatment))
plot_death_count <- plot_death_count + theme_bw() + scale_color_manual(values=gg_color_hue(6), name='initial \n[formaldehyde] (mM)')
plot_death_count <- plot_death_count + scale_y_log10(limits=c(10, 1e+09), breaks=c(1e+02, 1e+04, 1e+06, 1e+08)) + scale_x_continuous(limits=c(0,3))
plot_death_count <- plot_death_count + xlab('time (hours)') + ylab('viable cells (CFU/mL)') + facet_grid(~type)
plot_death_count

death_nash_total$treatment <- as.factor(death_nash_total$treatment)
plot_death_nash <- ggplot() + geom_line(data = death_nash_total,size=1.2,aes(x=time,y=fex,color=treatment,group=treatment))
plot_death_nash <- plot_death_nash + theme_bw() + scale_color_manual(values=gg_color_hue(6), name='initial \n[formaldehyde] (mM)')
plot_death_nash <- plot_death_nash + facet_grid(~type)
plot_death_nash

