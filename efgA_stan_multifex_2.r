library(deSolve)
library(rstan)
library(ggplot2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#### Genrating a data 

vmax <- 9e+04        # kcat for methanol (1/h)
r <- 0.23            # growth rate from Jannell data (1/h)
theta <- 90000       # diffusion coefficient (1/h)
km <- 2e-02          # km for MDH 
kf <- 0.04           # km for Fin
vv <- 2e-13          # volume of the cell to the volume of the media 
gamma <- 3.6    #combined metabolic and maitenance use of formaldehyde
alpha <- 8.7    # maxium killing rate of fin, metabolic independent (1/h) 
kd <- 8        # half saturation 

parms <- c(   
  theta = theta,  
  vv = vv, 
  vmax = vmax,
  km = km, 
  kf = kf, 
  r=r,
  gamma=gamma,
  alpha=alpha,
  kd=kd
)

# the model simulates the growth of a bacterium on formaldehyde and methanol 
efgAmodel<-function(t,state,parms) {
  with(as.list(c(state,parms)), {
    #rate of change
    dfex = -theta*vv*(fex-fin)*n   # change in external formaldehyde 
    dm = - vv*n*(vmax*m/(km+m))    # change in methanol 
    dfin =  theta*(fex-fin) + (vmax*m/(km+m)) - (vmax*gamma)*(fin/(kf+fin)) # change in internal formaldehyde
    dn = r*n*(fin/(kf+fin)) - n*alpha*(fin^2)/(fin^2+kd^2) # change in number of cells 
    #return the rate of change
    list(c(dfex,dm,dfin,dn))
  }) 
}

model_multi_levels<-function(model, time_range, baseline_state, parameters, fex_init_range){
  out<-data.frame() # create an empty dataframe for putting all of the outputs into
  for (i in fex_init_range) { # for-loop to run through all initial external formaldehyde levels
    state_i<-baseline_state # start off with your baseline state variables, of which all except fex will remain the same
    state_i[1]<-i # replace the fex from the baseline set of state variables with the fex value you want to model
    out_i <- as.data.frame(ode(y=state_i, times=time_range, func=model, parms=parameters, method="lsoda")) # run the model
    fex_init<-rep(i, length(time_range)) # create a column called "fex_init" to append to the data tables, so you know which initial fex value yielded these outputs
    out<-rbind(out, cbind(out_i, fex_init)) # append the results of this model run, along with a column indicating the initial formaldehyde concentration, to the master output dataframe
  }
  return(out)
}

model_multi_state<-function(baseline_state, fex_init_range){
  state_multi<-data.frame()
  for (i in fex_init_range) { # for-loop to run through all initial external formaldehyde levels
    state_i<-baseline_state # start off with your baseline state variables, of which all except fex will remain the same
    state_i[1]<-i # replace the fex from the baseline set of state variables with the fex value you want to model
    state_multi <- rbind(state_multi,state_i)
  }
  names(state_multi) <- c("fex","m","fin","n")
  return(state_multi)
}


times_growth=seq(0,70,by=1)
growth_fex_range<-c(0,2,3) # use one range of low fex levels for growth
#growth_fex_range<-seq(0,10,1) # use one range of low fex levels for growth
state_growth<-c(fex = 20, m = 15, fin = 0, n = 388500)
multi_growth_state <- model_multi_state(baseline_state = state_growth, fex_init_range = growth_fex_range)
growth_out<-model_multi_levels(model=efgAmodel, time_range=times_growth, baseline_state=state_growth, parameters=parms, fex_init_range=growth_fex_range)

'for (f in c(0,2,3,4)){
  growth_out$n[growth_out$fex_init==f] <- rpois(71,growth_out$n[growth_out$fex_init==f])
}'


growth_out$fex_init<-factor(growth_out$fex_init) # need to make the initial formaldehyde level into a factor (categorical variable rather than continuous), for easy plotting

gg_color_hue <- function(n) { # this is a cute function I found online, for making evenly-spaced color codes
  hues = seq(15, 375, length=n+1) # n = the number of factors you're plotting
  hcl(h=hues, l=65, c=100)[1:n] 
}

plot_growth_model<-ggplot() + geom_line(data=growth_out, size=1.2, 
                                        aes(x=time, y=n, col=fex_init, group=fex_init)) # I'm multiplying n by 10^5 so that the y axis looks right
plot_growth_model<-plot_growth_model + scale_color_manual(values=gg_color_hue(length(growth_fex_range)), name='initial \n[formaldehyde] (mM)')
plot_growth_model<-plot_growth_model + theme_bw() + xlab('time (hrs)') + ylab('viable cells (CFU/mL)')
plot_growth_model<-plot_growth_model + scale_y_log10(limits=c(1e+05, 5e+8), breaks=c(1e+04, 1e+06, 1e+08, 1e+10)) + scale_x_continuous(limits=c(0,70))
plot_growth_model

efgA_stan <- '
functions {real[] efgA(real t, real[] y, real[] eta,real[] x_r,int[] x_i) {
//intialize parameters

real gamma;
real alpha;
real kd;

real dydt[4];    

gamma = eta[1];
alpha = eta[2];
kd = eta[3];

//calculate the ode values, the equations are same as model in the r code 
dydt[1] = -90000*(2e-13)*(y[1]-y[3])*y[4];
dydt[2] = - (2e-13)*y[4]*(90000*y[2]/(0.02+y[2]));
dydt[3] = 90000*(y[1]-y[3]) + (90000*y[2]/(0.02+y[2])) - (90000*10*gamma)*(y[3]/(0.04+y[3])); 
dydt[4] = 0.23*y[4]*(y[3]/(0.04+y[3])) - y[4]*10*alpha*(y[3]^2)/(y[3]^2+(10*kd)^2);
return dydt;
}
}

data {
int<lower=1> T; // time 
int<lower=1> C; //number of fex initial conditions
real f0[C]; // pass a vector of length C initial fex values
real y0[3];     // initial values of state variables; changed to 3, fex_0 in f0 now
real t0;        // initial time
real ts[T];     // real array that has n time steps
int y[T,C];       // getting the data, here I pass cell count 
real rel_tol;
real abs_tol;
int max_steps;
}

transformed data {
real x_r[0];      
int x_i[0];
}

parameters { 
real<lower=0,upper=1> gamma; 
real<lower=0,upper=1> alpha;
real<lower=0,upper=1> kd;
} 

model {
real y_hat[T,4];
real eta[3];
real s0[4]; //this will now have the initial conditions

// priors, I want to fit gamma, alpha and kd. so they have a high sd, the rest supposed to be constant 
// gamma ~ lognormal(10,5); 
// alpha ~ lognormal(10,5); 
//  kd ~ lognormal(10,5); 

eta[1] = gamma;
eta[2] = alpha;
eta[3] = kd;

//initial conditions for m, fin, and n are not changing
//so we define them outside of the loop
s0[2] = y0[1];
s0[3] = y0[2];
s0[4] = y0[3];



for (c in 1:C){
//at the beginning of each c we create the correct fex initial condition
s0[1] = f0[c];
y_hat = integrate_ode_bdf(efgA, s0, t0, ts, eta, x_r, x_i, rel_tol, abs_tol, max_steps); //note we pass s0 now, I can enter the values of rel_tol... directly here 
y[,c] ~ poisson(y_hat[,4]); //note that I vectorized sampling and made sure to get column c of data
}
}
'
# compiling the model 
efgA_stan.comp<-stan_model(model_code=efgA_stan)

############ now Stan

y_data <- matrix(nrow = 70, ncol = 3) #make an empty matrix for the data that has the right dimensions
with(growth_out,
     for(i in 1:70) y_data[i,] <<- round(n[time == i])
)

# data: y says a matrix of 2:71 not counting 0, the last column which is N data 
efgA_data <- list(T = 70, C=3, f0 = growth_fex_range, y0 = c(15,0,388500), t0 = 0, ts = 1:70, y = y_data, rel_tol = 1e-6, abs_tol = 1e-8,max_steps = 1e6)

# sampling 
fit <- sampling(efgA_stan.comp, efgA_data, chains=4, iter=100, verbose = T, control = list(adapt_delta = 0.99))
pairs(fit)

