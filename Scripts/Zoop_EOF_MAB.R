## indicators--zoop index

library(TMB)
library(VAST)
library(FishStatsUtils)
library(dplyr)

# directory
rootdir = "~/Desktop/indicators"
Date = Sys.Date()
datedir = file.path( rootdir, Date)
dir.create(datedir)


## Model 2: zooplankton index (EOF) ##

#dir
mod<-"eof_MAB_zoop_5_amo_SFC"
rundir = file.path( datedir, mod )
dir.create(rundir)

#data
ecomon_dat <- readRDS(file.path(rootdir,"ecomon.rds"))
epu_list <- "MAB"
season="spring"

spp_list <- ecomon_dat %>%
  filter(!spp == "volume",
         EPU == epu_list,
         season == season) %>%
  group_by(spp) %>%
  summarize(count = mean(abundance, na.rm = TRUE)) %>%
  arrange(desc(count)) %>%
  head(5) %>%
  arrange(spp) %>%
  pull(spp)


# spp_list <- c(spp_list, "cirr", "hyper", "volume")
#spp_list <- c("calfin", "cham", "ctyp", "tlong")



all_dat <- ecomon_dat %>%
  dplyr::filter(spp %in% spp_list,
                EPU %in% epu_list,
                season == !!season,
                year >= 2000) %>%
  dplyr::mutate(species_number = group_indices(., spp)-1,
                areaswept_km2 = 0.01,
                catch_ab = log(abundance + 1),
                bottom_depth = bottom_depth/100) %>%
  dplyr::select(species_number,
                year,
                day,
                sfc_temp,
                abundance,
                catch_ab,
                areaswept_km2,
                bottom_depth,
                EPU,
                lat,
                lon) %>%
  arrange(species_number, year, lat, lon) %>%
  data.frame()


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

epu<- "MAB"
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
                 b_i= all_dat[, "catch_ab"],
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

