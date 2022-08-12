library(deSolve)
vmax <- 9e+04        # kcat for methanol (1/h), was 7.2e+03
r <- 0.23            # growth rate from Jannell data (1/h)
theta <- 90000       # diffusion coefficient (1/h) try 16000
km <- 2e-02          # km for MDH 
kf <- 0.04           # km for Fin, I started with 0.02 
vv <- 2e-13          # volume of the cell to the volume of the media 
gamma <- 2.3         #1.68      # combined metabolic and maitenance use of formaldehyde try 1.4
alpha <- 0.03159032  #0.02    # maxium killing rate of fin, metabolic independent (1/h) # starting with 0.32
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
    dfex = -theta*vv*(fex-(1/2)*(-vmax*gamma+fex*theta-kf*theta+vmax+sqrt((fex*theta-gamma*vmax-kf*theta+vmax)^2+4*theta*(fex*kf*theta+kf*vmax)))/theta)*n
    dn = (1/2)*r*n*(fex*theta-vmax*gamma-kf*theta+vmax+sqrt(fex^2*theta^2-2*fex*gamma*theta*vmax+2*fex*kf*theta^2+gamma^2*vmax^2+2*gamma*kf*theta*vmax+kf^2*theta^2+2*fex*theta*vmax-2*gamma*vmax^2+2*kf*theta*vmax+vmax^2))/(theta*(kf+(1/2)*(fex*theta-vmax*gamma-kf*theta+vmax+sqrt(fex^2*theta^2-2*fex*gamma*theta*vmax+2*fex*kf*theta^2+gamma^2*vmax^2+2*gamma*kf*theta*vmax+kf^2*theta^2+2*fex*theta*vmax-2*gamma*vmax^2+2*kf*theta*vmax+vmax^2))/theta))-alpha*((1/2)*(fex*theta-vmax*gamma-kf*theta+vmax+sqrt(fex^2*theta^2-2*fex*gamma*theta*vmax+2*fex*kf*theta^2+gamma^2*vmax^2+2*gamma*kf*theta*vmax+kf^2*theta^2+2*fex*theta*vmax-2*gamma*vmax^2+2*kf*theta*vmax+vmax^2))/theta)^sigma*n
    
    #return the rate of change
    list(c(dfex,dn))
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
    fex_init<-rep(i, length(time_range)) # create a column called "fex_init" to append to the data tables, so you know which initial fex value yielded these outputs
    fex_endpoint<-as.numeric(c(out_i$fex[-1], 'NA')) # links each timepoint to the fex value of the following timepoint, for purposes of plotting arrows wtih geom_segment 
    n_endpoint<-as.numeric(c(out_i$n[-1], 'NA')) # links each timepoint to the n value of the following timepoint, for purposes of plotting arrows wtih geom_segment 
    out<-rbind(out, cbind(out_i, fex_init, fex_endpoint, n_endpoint)) # append the results of this model run, along with a column indicating the initial formaldehyde concentration, to the master output dataframe
  }
  return(out)
}

times_growth=seq(0,200,by=10)
#growth_fex_range<-c(0,2,3,4,5,7.5,12.5,15,20) # use one range of low fex levels for growth
growth_fex_range<-c(1,3,4,4.2,4.4,4.6,5) # use one range of low fex levels for growth
state_growth<-c(fex = 20, n = 1e+06)
growth_out<-model_multi_levels(model=efgAmodel, time_range=times_growth, baseline_state=state_growth, parameters=parms, fex_init_range=growth_fex_range)
growth_out$fex_init<-factor(growth_out$fex_init) # need to make the initial formaldehyde level into a factor (categorical variable rather than continuous), for easy plotting

library(ggplot2)

gg_color_hue <- function(n) { # this is a cute function I found online, for making evenly-spaced color codes
  hues = seq(15, 375, length=n+1) # n = the number of factors you're plotting
  hcl(h=hues, l=65, c=100)[1:n] 
}

theme_set(theme_bw(base_size = 20)) 

n_nullcline <- (1/2)*(alpha^2*gamma*kf^2*vmax+2*alpha*gamma*r*vmax-alpha*kf*r*theta-2*alpha*r*vmax+sqrt(alpha^4*gamma^2*kf^4*vmax^2+4*alpha^3*gamma^2*kf^2*r*vmax^2-2*alpha^3*gamma*kf^3*r*theta*vmax-8*alpha^2*gamma*kf*r^2*theta*vmax+alpha^2*kf^2*r^2*theta^2+4*alpha*r^3*theta^2))/(alpha*r*theta)
fex_nullcline <- kf/(gamma-1)

plot_growth_arrows<-ggplot() + geom_segment(data=growth_out, size=1.2, aes(x=fex, y=n, xend=fex_endpoint, yend=n_endpoint, color=fex_init), 
                                                      arrow = arrow(length=unit(0.25, "cm")))
plot_growth_arrows<-plot_growth_arrows + scale_color_manual(values=gg_color_hue(length(growth_fex_range)), name='formaldehyde \ntreatment (mM)') + ylab("viable cells (CFU/mL)") + xlab("external formaldehyde")
plot_growth_arrows<-plot_growth_arrows + geom_vline(xintercept = n_nullcline,colour='Red') + geom_vline(xintercept = fex_nullcline)
plot_growth_arrows<-plot_growth_arrows + scale_y_log10(limits=c(1e+3,1e+10)) + scale_x_continuous(limits=c(0,5))
plot_growth_arrows

