## indicators-- portfolio effects

library(TMB)
library(VAST)
library(FishStatsUtils)
library(dplyr)
library(ggplot2)
library(DHARMa)

# directory
rootdir = "~/Desktop/indicators"
Date = Sys.Date()
datedir = file.path( rootdir, Date)
dir.create(datedir)

## Model 1: portfolio effects (JDSDM) ##

#dir
mod<-"jdsdm_port_guilds_NEUS"
rundir = file.path( datedir, mod )
dir.create(rundir)

#data
dat<-read.csv(file.path(rootdir,"/nefsc_bts.csv"))
species<-read.csv(file.path(rootdir,"/species_groups.csv"))
#year_set = 2000:2019
wanted<-c("piscivore","planktivore","benthivore","benthos")
epu<-"All"
dat<-dat[dat$group %in% wanted,]

for(i in 1:length(dat$id)){
dat$guild[i]<-if(dat$group[i] == "piscivore" & dat$season[i]=="FALL"){1
          }else if(dat$group[i] == "piscivore" & dat$season[i]=="SPRING"){2
          }else if(dat$group[i] == "planktivore" & dat$season[i]=="FALL"){3
          }else if(dat$group[i] == "planktivore" & dat$season[i]=="SPRING"){4
          }else if(dat$group[i] == "benthivore" & dat$season[i]=="FALL"){5
          }else if(dat$group[i] == "benthivore" & dat$season[i]=="SPRING"){6
          }else if(dat$group[i] == "benthos" & dat$season[i]=="FALL"){7
          }else if(dat$group[i] == "benthos" & dat$season[i]=="SPRING"){8}
}
                                                          


dat$bottom_depth<-(dat$bottom_depth - mean(dat$bottom_depth))/sd(dat$bottom_depth)

covariate_data <- data.frame(Year = NA,
                             Lat = dat$lat,
                             Lon = dat$lon,
                             bottom_depth = dat$bottom_depth)

#settings
strata_limits = "EPU"
strata_limits_all <- list('All_areas' = 1:1e5)


formula <- ~ bottom_depth + I(bottom_depth^2)


Options <- c("SD_site_density"          = 0,
             "SD_site_logdensity"       = 0,
             "Calculate_Range"          = 1, # Center of gravity
             "Calculate_evenness"       = 0,
             "Calculate_effective_area" = 1, # Area occupied
             "Calculate_Cov_SE"         = 0,
             "Calculate_Synchrony"      = 1,
             "Calculate_proportion"     = 1)

settings<-make_settings( n_x=100,
                         Region="northwest_atlantic",
                         fine_scale=FALSE,
                         use_anisotropy = FALSE,
                         n_categories = 3,
                         Options=Options,
                         strata.limits = strata_limits_all,
                         purpose= "ordination",
                         ObsModel = c(1,1))
settings$RhoConfig[c('Beta1','Beta2')] = 2

#settings$RhoConfig[c('Epsilon1','Epsilon2')] = 4 #

settings$epu_to_use <- switch(epu,
                              "MAB" = "Mid_Atlantic_Bight",
                              "GOM" = "Gulf_of_Maine",
                              "GB" = "Georges_Bank",
                              "SS" = "Scotian_Shelf",
                              "All" = "All")
#run
fit = fit_model( "settings"=settings,
                 epu_to_use = settings$epu_to_use,
                 "Lat_i"=dat[,'lat'],
                 "Lon_i"=dat[,'lon'],
                 "t_i"=dat[,'year'],
                 "c_i"=as.numeric(dat[,'guild']-1),
                 "b_i"=dat[,'catch_kg'],
                 "a_i"=dat[,'areaswept_km2'],
                 "working_dir"=rundir,
                 test_fit=FALSE,
                 knot_method ="grid",
                 newtonsteps=0,
                 getsd=TRUE,
                 Use_REML = TRUE,
                 optimize_args = list("lower" = -Inf,
                                      "upper" = Inf),
                 covariate_data=covariate_data,
                 X1_formula=formula,
                 "bias.correct.control"=list("nsplit"=100) )


saveRDS(fit, file = paste0(rundir, "/fit.rds"))

category_names<- unique(dat$guild)

#plot
results = plot_results( settings=settings, fit=fit,working_dir=paste0(moddir,"/"),check_residuals=TRUE )

plot_factors( Report=fit$Report, ParHat=fit$ParHat, Data=fit$data_list, SD=fit$parameter_estimates$SD,  RotationMethod = "Varimax",
              mapdetails_list=results$map_list, Year_Set=levels(dat$year), plotdir=paste0(rundir,"/") )
saveRDS(results, file = paste0(rundir, "/results.rds"))
save.image(file = paste0(rundir, "/output.RData"))
