## indicators --condition/density

library(TMB)
library(VAST)
library(FishStatsUtils)
library(dplyr)




# directory
rootdir = "~/Desktop/indicators"
Date = Sys.Date()
datedir = file.path( rootdir, Date)
dir.create(datedir)

## Model 3: groundfish condition (condition-density)
mod<-"cond"
rundir = file.path(datedir, mod)
dir.create(rundir)
     
# epu <- "neus"
spp_list<-c("Atlantic cod","haddock","pollock", "spiny dogfish","silver hake","winter flounder","yellowtail flounder")

species<-read.csv(file.path(rootdir,"/species_groups.csv"))

species<-species[species$common_name %in% spp_list,]
species$species_number<-1:length(species$common_name)


season <- "fall"
region <- "northwest_atlantic"
epu = "All"

n_x = 200   # Specify number of knots for predictive process
#strata_limits = data.frame('Georges_Bank'=c(1130, 1140, 1150, 1160, 1170, 1180, 1190, 1200, 1210, 1220, 1230, 1240, 1250, 1290, 1300))
strata_limits_all <- list('All_areas' = 1:1e5)
strata_limits = "EPU"


nefsc_dat <- readRDS(url("https://github.com/Laurels1/Condition/raw/master/data/NEFSC_survey_data_02-13-20.rds")) %>%
  mutate(id = as.character(paste0(CRUISE6, "_", STATION)))

nefsc_dat$bottom_temp<-(nefsc_dat$BOTTEMP - mean(nefsc_dat$BOTTEMP,na.rm=T))/sd(nefsc_dat$BOTTEMP,na.rm=T)



bot<-read.csv(file.path(rootdir,"/nefsc_bts.csv"))
bot<-bot[,2:5]
colnames(bot)<-c("bottom_depth","BEGLAT","BEGLON","YEAR")
bot<-distinct(bot)

nefsc_dat<-merge(nefsc_dat, bot)

station_dat <- nefsc_dat %>%
  filter(SEASON == toupper(season)) %>%
  select(id,
         year = YEAR,
         latitude = BEGLAT,
         longitude = BEGLON,
         bottom_temp = bottom_temp) %>%
  distinct(.keep_all = TRUE)

cpue_dat <- nefsc_dat %>%
  filter(SEASON == toupper(season),
         SVSPP %in% species$SVSPP) %>%
  select(id,
         latitude = BEGLAT,
         bottom_temp = bottom_temp,
         bottom_depth = bottom_depth,
         spp = SVSPP,
         longitude = BEGLON,
         year = YEAR,
         cpue_kg = EXPCATCHWT) %>%
  distinct(.keep_all = TRUE) %>%
  right_join(station_dat) %>%
  mutate(cpue_kg = ifelse(is.na(cpue_kg),
                          0,
                          cpue_kg))

lw_dat <- nefsc_dat %>%
  filter(SEASON == toupper(season),
         SVSPP %in% species$SVSPP) %>%
  select(id,
         latitude = BEGLAT,
         longitude = BEGLON,
         spp = SVSPP,
         bottom_temp = bottom_temp,
         bottom_depth = bottom_depth,
         year = YEAR,
         weight_g = INDWT,
         length_mm = LENGTH) %>%
  mutate(weight_g = weight_g*1000,
         length_mm = length_mm*10) %>%
  right_join(station_dat) %>%
  na.omit()

nspecies<-length(spp_list) 

plotyears<-list(23:49)
plotyears[[2]]<-23:49
plotyears[[3]]<-23:49
plotyears[[4]]<-1:27
plotyears[[5]]<-c(13:38, 40)
plotyears[[6]]<-20:46
plotyears[[7]]<-c(23:33,35:49)

nefsc_dat$bottom_depth<-(nefsc_dat$bottom_depth-mean(nefsc_dat$bottom_depth, na.rm=T))/sd(nefsc_dat$bottom_depth, na.rm=T)


covariate_data <- data.frame(Year = nefsc_dat$YEAR,
                             Lat = nefsc_dat$BEGLAT,
                             Lon = nefsc_dat$BEGLON,
                             bottom_depth = nefsc_dat$bottom_depth)
covariate_data<-covariate_data[complete.cases(covariate_data),]

i<-7
for(i in 1:nspecies){
  
  species.i<- species$SVSPP[species$common_name==spp_list[i]]
  name =  paste(spp_list[i], mod, "loop","NEUS",sep="-")
  moddir = file.path( rundir, name)
  dir.create(moddir)
  
  lw_dat.spp<-lw_dat[lw_dat$spp==species.i,]
  cpue_dat.spp<-cpue_dat[cpue_dat$spp==species.i,]
  
  lw_dat.spp<-lw_dat.spp[complete.cases(lw_dat.spp[,1:2]),]
  cpue_dat.spp<-cpue_dat.spp[complete.cases(cpue_dat.spp[,1:2]),]
  
  
  n<-lw_dat.spp %>% 
    group_by(year) %>%
    summarise(no_rows = length(year))
  
  
  all_dat <- bind_rows(lw_dat.spp,cpue_dat.spp) %>%
    filter(year >= 1970,
           year <= 2019,
           !is.na(latitude),
           !is.na(longitude)) %>%
    mutate(b_i = ifelse(!is.na(cpue_kg),
                        cpue_kg,
                        weight_g),
           c_i = ifelse(!is.na(cpue_kg), 
                        0,
                        1),
           Q_i = ifelse(!is.na(cpue_kg),
                        0,
                        log(length_mm/10)))
  
  
  
  catchability_data = data.frame( "length_cm" = ifelse(!is.na(all_dat[,'cpue_kg']),
                                                    1, all_dat[,'length_mm']/10 ))
  
  
  Options <- c("SD_site_density"          = 0,
               "SD_site_logdensity"       = 0,
               "Calculate_Range"          = 1, # Center of gravity
               "Calculate_evenness"       = 0,
               "Calculate_effective_area" = 1, # Area occupied
               "Calculate_Cov_SE"         = 1,
               "Calculate_Synchrony"      = 0,
               "Calculate_proportion"     = 0)
  
  #settings
  settings <- make_settings(
                            n_x = n_x,
                            Region = region,
                            strata.limits = strata_limits_all,
                            # strata.limits = "EPU",
                            Options = Options,
                            purpose = "condition_and_density",
                            bias.correct = FALSE,
                            knot_method = "grid")
  
  # settings$epu_to_use <- "All"
 settings$FieldConfig[c("Omega","Epsilon"),"Component_1"] <-"IID"
  
  Expansion_cz <- matrix(c(0, 0, 2, 0), ncol = 2, byrow = TRUE)
  settings$ObsModel <- matrix(c(2, 4, 1, 4), ncol = 2, byrow = TRUE)
  
  covariate_data$bottom_depth<-(covariate_data$bottom_depth-mean(covariate_data$bottom_depth, na.rm=T))/sd(covariate_data$bottom_depth, na.rm=T)
  formula <- ~ bottom_depth + I(bottom_depth^2)
  
  
  #run
  fit <- fit_model(settings = settings,
                   Lat_i = all_dat[,'latitude'],
                   Lon_i = all_dat[,'longitude'],
                   # epu_to_use = settings$epu_to_use,
                   t_i = all_dat[,'year'],
                   c_i = all_dat[, 'c_i'],
                   b_i = all_dat[, 'b_i'],
                   a_i = rep(1, nrow(all_dat)),
                   Q_ik = matrix(all_dat[, 'Q_i'], ncol = 1),
                   Expansion_cz = Expansion_cz,
                   working_dir = moddir,
                   catchability_data = catchability_data,
                   Q2_formula = ~log(length_cm),
                   X1_formula = formula,
                   covariate_data = covariate_data,
                   getsd = TRUE,
                   knot_method = "grid",
                   getJointPrecision=TRUE,
                   # newtonsteps = 1,
                   # test_fit = TRUE,
                   # optimize_args = list("lower" = -Inf,
                   #                      "upper" = Inf),
                   # bias.correct.control = list("nsplit" = 100),
                   build_model = FALSE)
  
  Map = fit$tmb_list$Map
  Map$lambda2_k = factor(NA)
  
  
  fit <- fit_model(settings = settings,
                   Lat_i = all_dat[,'latitude'],
                   Lon_i = all_dat[,'longitude'],
                   # epu_to_use = settings$epu_to_use,
                   t_i = all_dat[,'year'],
                   c_i = all_dat[, 'c_i'],
                   b_i = all_dat[, 'b_i'],
                   a_i = rep(1, nrow(all_dat)),
                   Q_ik = matrix(all_dat[, 'Q_i'], ncol = 1),
                   Expansion_cz = Expansion_cz,
                   catchability_data = catchability_data,
                   Q2_formula = ~log(length_cm),
                   X1_formula = formula,
                   covariate_data = covariate_data,
                   working_dir = moddir,
                   Map = Map,
                   Use_REML=TRUE,
                   newtonsteps=1,
                   getsd = TRUE,
                   getJointPrecision=TRUE,
                   check_fit=FALSE, 
                   optimize_args = list("lower" = -Inf,
                   "upper" = Inf))

  
  saveRDS(fit, file = paste0(moddir,"/",name,".rda",sep=""))
  save.image(file = paste0(moddir, "/output.RData"))
  
  results = plot_results( settings=settings, fit=fit, 
                          years_to_plot = plotyears[[i]],
                          #years_to_plot = c(23,37,49),
                          plot_set = 3,
                          working_dir=paste0(rootdir,"/"),
                          check_residuals=FALSE)
  plot_factors( Report=fit$Report, ParHat=fit$ParHat, Data=fit$data_list, SD=fit$parameter_estimates$SD,
                mapdetails_list=results$map_list, Year_Set=fit$years_to_plot, plotdir=paste0(moddir,"/") )
  
  
}
