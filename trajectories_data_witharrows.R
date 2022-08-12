library(ggplot2)
setwd("C:/Siavash/Codes/phenmod/plots_data")

growth_viability <- read.csv("formaldehyde_growth_experiments_Dec2015_CFUcounts_final.csv")
growth_nash <- read.csv("formaldehyde_growth_experiments_Dec2015_NashAssay.csv")
death_viability <- read.csv("RECOUNT_killingcurve_8615_alldata_final_nash_final_data.csv")
death_nash <- read.csv("RECOUNT_killingcurve_8615_alldata_final_viability.csv")

tmp1 <- merge(growth_viability,growth_nash)
tmp1 <-tmp1[tmp1$substrate !="succinate",]
tmp1 <-tmp1[tmp1$strain !="CM3745",]
growth_total <- subset(tmp1, select = -c(substrate,strain))

tmp2 <- merge(death_viability,death_nash)
tmp2 <-tmp2[tmp2$strain !="CM3745",]
#tmp2 <-tmp2[tmp2$formaldehyde_nash !="NA",]
death_total <- subset(tmp2, select = -c(strain,NashA412,standard_curve))
colnames(growth_total)[c(4,6)] <- c("CFU_mL","formaldehyde_nash")
growth_total = growth_total[,-3]
combined <- rbind(growth_total,death_total)

names(combined) <- c("flask", "time", "conc", "cells","nash")
# making data in time order 
combined <- combined[order(combined$time),]
# saying conc as factor for plotting with different color 
combined$conc<-factor(combined$conc)


#combined$nash <- as.numeric(combined$nash) making nash data numeric 

plot_growth_model<-ggplot() + geom_path(data=combined, size=1.2, 
                                        aes(x=nash, y=cells, col=conc, group=conc),arrow=arrow(ends = "last", type = "closed", length = unit(0.1, "inches"))) # I'm multiplying n by 10^5 so that the y axis looks right
#plot_growth_model<-plot_growth_model + theme_bw() + xlab('Fex') + ylab('viable cells (CFU/mL)')
plot_growth_model<-plot_growth_model + scale_y_log10(limits=c(1, 5e+8), breaks=c(1e+02,1e+04, 1e+06, 1e+08, 1e+10)) + scale_x_continuous(limits=c(0,20)) + theme_bw()
plot_growth_model

#####################################################

### making a plot with arrows

# we're going to make a new dataframe that's just like combined,
# but also has columns to hold nash and cells values for the endpoints of the segments in the plot
# We'll use a for-loop to break "combined" apart by concentration
# and for each concentration, add the cells_endpoint and nash_endpoint columns
# then recombine them with rbind.

# first, make an empty dataframe:
combined2<-data.frame(flask=numeric(), time=numeric(), conc=factor(), cells=numeric(), 
                      nash=numeric(), cells_endpoint=numeric(), nash_endpoint=numeric())
# then, loop through the concentrations
for(concentration in levels(combined$conc)){
  subset<-combined[combined$conc==concentration,]
  subset<-subset[order(subset$time),] # subset by concentration, and also make sure the entries are in order by timepoint
  subset$cells_endpoint<-as.numeric(c(subset$cells[-1],'NA')) # cells_endpoint is simply the cells column, offset by 1
  subset$nash_endpoint<-as.numeric(c(subset$nash[-1],'NA')) # same with nash_endpoint
  combined2<-rbind(combined2, subset) # rbind them all
}
theme_set(theme_bw(base_size = 20)) 
# I made just a few changes to the plot format, to make it more similar to the plot of model data
plot_data_arrows<-ggplot() + geom_segment(data=combined2, size=1.2, 
                                        aes(x=nash, y=cells, xend=nash_endpoint, yend=cells_endpoint, col=conc, group=conc),
                                        arrow=arrow(length = unit(0.25, "cm")))
plot_data_arrows<-plot_data_arrows + xlab('Fex') + ylab('cells')
plot_data_arrows<-plot_data_arrows + scale_y_log10(limits=c(1, 5e+8), breaks=c(1e+02,1e+04, 1e+06, 1e+08, 1e+10)) + scale_x_continuous(limits=c(0,20)) 
plot_data_arrows

