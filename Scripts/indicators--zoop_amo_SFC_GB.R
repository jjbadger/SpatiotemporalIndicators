## indicators--zoop index

library(TMB)
library(VAST)
library(FishStatsUtils)
library(dplyr)

# directory
rootdir = "~/Desktop/indicators/zoop EOF"
Date = Sys.Date()
datedir = file.path( rootdir, Date)
dir.create(datedir)


## Model 2: zooplankton index (EOF) ##

#data
ecomon_dat <- readRDS(file.path(rootdir,"ecomon.rds"))
epu<- "MAB"
season="spring"


spp_list <- ecomon_dat %>%
  filter(!spp == "volume",
         EPU == epu,
         season == !!season) %>%
  group_by(spp) %>%
  summarize(count = mean(abundance, na.rm = TRUE)) %>%
  arrange(desc(count)) %>%
  head(5) %>%
  arrange(spp) %>%
  pull(spp)

# spp_list <- c(spp_list, "cirr", "hyper", "volume")
#spp_list <- c("calfin", "ctyp", "cham", "tlong")


# bottom salinity-- NOT WORKING
load(url(sprintf("https://github.com/NOAA-EDAB/ECSA/blob/master/data/gridded/sal_bottom_%s_spdf.rdata?raw=true", season)))
names(ecsa_dat) <- stringr::str_extract(names(ecsa_dat), "[0-9]{4}")
sal_bottom <- ecsa_dat
sal_years <- as.numeric(gsub("X", "", names(sal_bottom))) #1992-2017
extract_bs <- function(xyear, longitude, latitude, ...) raster::extract(sal_bottom[[xyear]], y = cbind(longitude, latitude))
possibly_extract_bs <- purrr::possibly(extract_bs, otherwise = NA_real_)

yrs<-2000:2017


all_dat <- ecomon_dat %>%
  dplyr::filter(spp %in% spp_list,
                EPU == epu,
                season == !!season,
                year >= 2000) %>%
  dplyr::mutate(species_number = group_indices(., spp)-1,
                areaswept_km2 = 0.01,
                catch_ab = log(abundance + 1),
                bottom_depth = bottom_depth/100, 
                bottom_sal = as.numeric(purrr::pmap(list(year, lon, lat), possibly_extract_bs)),
                bottom_sal = as.numeric(abs(scale(bottom_sal)))) %>%
  dplyr::select(species_number,
                year,
                day,
                sfc_temp,
                abundance,
                catch_ab,
                areaswept_km2,
                bottom_depth,
                bottom_sal,
                EPU,
                lat,
                lon) %>%
  arrange(species_number, year, lat, lon) %>%
  data.frame()

######## AMO
library(reshape2)
amo<-read.table(file.path(rootdir,"amo.txt"), header=T)
colnames(amo)<-c("year", 1:12)
amo<-melt(amo, id.vars = "year", measure.vars =colnames(amo)[2:13])
colnames(amo) =  c("year","month","index")
amo$month<-as.integer(as.character(amo$month))

all_dat$month <- cut(all_dat$day,
                     breaks=c(1,31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365),
                     labels = 1:12)
all_dat$month = as.integer(as.character(all_dat$month))

all_dat$amo<-NA

for(i in 1:length(all_dat$year)){
  year.i<-all_dat$year[i]
  month.i<-all_dat$month[i]
  all_dat$amo[i]<-amo$index[amo$year==year.i & amo$month==month.i]
}

########NAO

nao<-read.csv(file.path(rootdir,"nao.csv"), header=F)
nao<-nao[,c(4:6)]
colnames(nao) =  c("year","month","index")

all_dat$nao<-NA

for(i in 1:length(all_dat$year)){
  year.i<-all_dat$year[i]
  month.i<-all_dat$month[i]
  all_dat$nao[i]<-nao$index[nao$year==year.i & nao$month==month.i]
}

##############################
all_dat$sfc_temp<-(all_dat$sfc_temp-mean(all_dat$sfc_temp, na.rm=T))/sd(all_dat$sfc_temp, na.rm=T)




covariate_data <- data.frame(Year = all_dat$year,
                             Lat = all_dat$lat,
                             Lon = all_dat$lon,
                             amo = all_dat$amo,
                             nao = all_dat$nao,
                             sfc = all_dat$sfc_temp,
                             bs  = all_dat$bottom_sal)

#covariate_data<-covariate_data[complete.cases(covariate_data),]


#dir
mod<-"eof_MAB_zoop_nocov"
rundir = file.path( datedir, mod )
dir.create(rundir)

#all_dat$catch_ab[all_dat$species_number==3 & all_dat$year %in% c(2000,2001,2003,2007,2009,2010,2014:2017)]<-NA
#all_dat$catch_ab[all_dat$species_number==2 & all_dat$year==2003]<-NA
strata.list = "EPU"
#settings
settings = make_settings( n_x = 1000,
                          Region = "northwest_atlantic",
                          fine_scale = TRUE,
                          strata.limits = strata.list,
                          purpose="EOF2",
                          n_categories = 2,
                          use_anisotropy = FALSE, # corresponds to ln_H_input params
                          ObsModel = c(1,1),
                          RhoConfig = c("Beta1" = 2, "Beta2" = 2,
                                        "Epsilon1" = 0, "Epsilon2" = 0))
# settings$FieldConfig["Omega","Component_1"] <- 0
settings$Options[c("Calculate_Range","Calculate_effective_area")] = TRUE
settings$Options[c("zerosum_penalty")]=100


settings$epu_to_use <- switch(epu,
                              "MAB" = "Mid_Atlantic_Bight",
                              "GOM" = "Gulf_of_Maine",
                              "GB" = "Georges_Bank",
                              "SS" = "Scotian_Shelf",
                              "All" = "All")

#run
fit = fit_model( settings=settings,
                 epu_to_use = settings$epu_to_use,
                 Lat_i=all_dat[,'lat'],
                 Lon_i=all_dat[,'lon'],
                 t_i=all_dat[,"year"],
                 c_i=all_dat[,"species_number"],
                 b_i= all_dat[, "abundance"],
                 a_i=all_dat[,'areaswept_km2'],
                 working_dir = rundir,
                 anisotropy = FALSE, # corresponds to ln_H_input params
                 test_fit = FALSE,
                 newtonsteps = 0,
                 getsd = TRUE,
                 run_model = TRUE,
                 Use_REML = TRUE,
                 knot_method = "grid",
                 # Method = Method,
                 optimize_args = list("lower" = -Inf,
                                      "upper" = Inf),
                 bias.correct.control = list("nsplit" = 100))

category_names = spp_list

#plot

results = plot( fit, category_names=category_names,
                check_residuals=TRUE,
                plot_set=c(3,11,14,16),
                working_dir=paste0(rundir,"/"))

plot_factors( Report=fit$Report,
              ParHat=fit$ParHat,
              Data=fit$data_list,
              category_names = category_names,
              mapdetails_list=results$map_list,
              Year_Set=fit$year_labels, plotdir=paste0(rundir,"/"))

saveRDS(fit, file = paste0(rundir, "/fit.rds"))
saveRDS(results, file = paste0(rundir, "/results.rds"))
save.image(file = paste0(rundir, "/output.RData"))


#####################################################
# sfc 
#####################################################

formula <- ~ sfc + I(sfc^2) 
rm(fit)

#dir
mod<-"eof_MAB_zoop_sfc"
rundir = file.path( datedir, mod )
dir.create(rundir)


#run
fit = fit_model( settings=settings,
                 epu_to_use = settings$epu_to_use,
                 Lat_i=all_dat[,'lat'],
                 Lon_i=all_dat[,'lon'],
                 t_i=all_dat[,"year"],
                 c_i=all_dat[,"species_number"],
                 b_i= all_dat[, "abundance"],
                 a_i=all_dat[,'areaswept_km2'],
                 working_dir = rundir,
                 covariate_data = covariate_data,
                 X1_formula=formula,
                 anisotropy = FALSE, # corresponds to ln_H_input params
                 test_fit = FALSE,
                 newtonsteps = 0,
                 getsd = TRUE,
                 run_model = TRUE,
                 Use_REML = TRUE,
                 knot_method = "grid",
                 # Method = Method,
                 optimize_args = list("lower" = -Inf,
                                      "upper" = Inf),
                 bias.correct.control = list("nsplit" = 100))

category_names = spp_list

#plot

results = plot( fit, category_names=category_names,
                check_residuals=TRUE,
                plot_set=c(3,11,14,16),
                working_dir=paste0(rundir,"/"))

plot_factors( Report=fit$Report,
              ParHat=fit$ParHat,
              Data=fit$data_list,
              category_names = category_names,
              mapdetails_list=results$map_list,
              Year_Set=fit$year_labels, plotdir=paste0(rundir,"/"))

saveRDS(fit, file = paste0(rundir, "/fit.rds"))
saveRDS(results, file = paste0(rundir, "/results.rds"))
save.image(file = paste0(rundir, "/output.RData"))


#####################################################
# sfc +amo
#####################################################

formula <- ~ sfc + I(sfc^2) + amo
rm(fit)
#dir
mod<-"eof_MAB_zoop_sfc_amo"
rundir = file.path( datedir, mod )
dir.create(rundir)

#run
fit = fit_model( settings=settings,
                 epu_to_use = settings$epu_to_use,
                 Lat_i=all_dat[,'lat'],
                 Lon_i=all_dat[,'lon'],
                 t_i=all_dat[,"year"],
                 c_i=all_dat[,"species_number"],
                 b_i= all_dat[, "abundance"],
                 a_i=all_dat[,'areaswept_km2'],
                 working_dir = rundir,
                 covariate_data = covariate_data,
                 X1_formula=formula,
                 anisotropy = FALSE, # corresponds to ln_H_input params
                 test_fit = FALSE,
                 newtonsteps = 0,
                 getsd = TRUE,
                 run_model = TRUE,
                 Use_REML = TRUE,
                 knot_method = "grid",
                 # Method = Method,
                 optimize_args = list("lower" = -Inf,
                                      "upper" = Inf),
                 bias.correct.control = list("nsplit" = 100))

category_names = spp_list

#plot

results = plot( fit, category_names=category_names,
                check_residuals=TRUE,
                plot_set=c(3,11,14,16),
                working_dir=paste0(rundir,"/"))

plot_factors( Report=fit$Report,
              ParHat=fit$ParHat,
              Data=fit$data_list,
              category_names = category_names,
              mapdetails_list=results$map_list,
              Year_Set=fit$year_labels, plotdir=paste0(rundir,"/"))

saveRDS(fit, file = paste0(rundir, "/fit.rds"))
saveRDS(results, file = paste0(rundir, "/results.rds"))
save.image(file = paste0(rundir, "/output.RData"))

#####################################################
# sfc + nao
#####################################################

formula <- ~ sfc + I(sfc^2) + amo +nao
rm(fit)
#dir
mod<-"eof_MAB_zoop_sfc_nao"
rundir = file.path( datedir, mod )
dir.create(rundir)

#run
fit = fit_model( settings=settings,
                 epu_to_use = settings$epu_to_use,
                 Lat_i=all_dat[,'lat'],
                 Lon_i=all_dat[,'lon'],
                 t_i=all_dat[,"year"],
                 c_i=all_dat[,"species_number"],
                 b_i= all_dat[, "abundance"],
                 a_i=all_dat[,'areaswept_km2'],
                 working_dir = rundir,
                 covariate_data = covariate_data,
                 X1_formula=formula,
                 anisotropy = FALSE, # corresponds to ln_H_input params
                 test_fit = FALSE,
                 newtonsteps = 0,
                 getsd = TRUE,
                 run_model = TRUE,
                 Use_REML = TRUE,
                 knot_method = "grid",
                 # Method = Method,
                 optimize_args = list("lower" = -Inf,
                                      "upper" = Inf),
                 bias.correct.control = list("nsplit" = 100))

category_names = spp_list

#plot

results = plot( fit, category_names=category_names,
                check_residuals=TRUE,
                plot_set=c(3,11,14,16),
                working_dir=paste0(rundir,"/"))

plot_factors( Report=fit$Report,
              ParHat=fit$ParHat,
              Data=fit$data_list,
              category_names = category_names,
              mapdetails_list=results$map_list,
              Year_Set=fit$year_labels, plotdir=paste0(rundir,"/"))

saveRDS(fit, file = paste0(rundir, "/fit.rds"))
saveRDS(results, file = paste0(rundir, "/results.rds"))
save.image(file = paste0(rundir, "/output.RData"))

#####################################################
# sfc +nao +amo
#####################################################

formula <- ~ sfc + I(sfc^2) + amo + nao
rm(fit)
#dir
mod<-"eof_MAB_zoop_sfc_amo_nao"
rundir = file.path( datedir, mod )
dir.create(rundir)

#run
fit = fit_model( settings=settings,
                 epu_to_use = settings$epu_to_use,
                 Lat_i=all_dat[,'lat'],
                 Lon_i=all_dat[,'lon'],
                 t_i=all_dat[,"year"],
                 c_i=all_dat[,"species_number"],
                 b_i= all_dat[, "abundance"],
                 a_i=all_dat[,'areaswept_km2'],
                 working_dir = rundir,
                 covariate_data = covariate_data,
                 X1_formula=formula,
                 anisotropy = FALSE, # corresponds to ln_H_input params
                 test_fit = FALSE,
                 newtonsteps = 0,
                 getsd = TRUE,
                 run_model = TRUE,
                 Use_REML = TRUE,
                 knot_method = "grid",
                 # Method = Method,
                 optimize_args = list("lower" = -Inf,
                                      "upper" = Inf),
                 bias.correct.control = list("nsplit" = 100))

category_names = spp_list

#plot

results = plot( fit, category_names=category_names,
                check_residuals=TRUE,
                plot_set=c(3,11,14,16),
                working_dir=paste0(rundir,"/"))

plot_factors( Report=fit$Report,
              ParHat=fit$ParHat,
              Data=fit$data_list,
              category_names = category_names,
              mapdetails_list=results$map_list,
              Year_Set=fit$year_labels, plotdir=paste0(rundir,"/"))

saveRDS(fit, file = paste0(rundir, "/fit.rds"))
saveRDS(results, file = paste0(rundir, "/results.rds"))
save.image(file = paste0(rundir, "/output.RData"))
