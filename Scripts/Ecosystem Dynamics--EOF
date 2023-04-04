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


epus<-c("GOM","GB","MAB")

## Model 1: zooplankton dynamics  (EOF) ##

mod<-"zoop"
rundir = file.path( datedir, mod )
dir.create(rundir)


ecomon_dat <- readRDS(file.path(rootdir,"ecomon.rds"))


for(i in 1:length(epus)){
  
  #dir
  epu.i<- epus[i]
  name =  paste(epu.i,mod,"no cov",sep="-")
  moddir = file.path(rundir, name)
  dir.create(moddir)
  
  #data
  epu_list <- epu.i
  season="spring"
  
  
  # spp_list <- ecomon_dat %>%                    #pull out top 5 species
  #   filter(!spp == "volume",
  #          EPU == epu_list,
  #          season == !!season) %>%
  #   group_by(spp) %>%
  #   summarize(count = mean(abundance, na.rm = TRUE)) %>%
  #   arrange(desc(count)) %>%
  #   head(5) %>%
  #   arrange(spp) %>%
  #   pull(spp)
  

  spp_list <- c("calfin", "ctyp", "cham", "tlong", "oithspp")     #chosen by experts
  
  
  
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
  
  
  strata.list = "EPU"
  #settings
  settings = make_settings( n_x = 100,
                            Region = "northwest_atlantic",
                            fine_scale = TRUE,
                            strata.limits = strata.list,
                            purpose="eof3",
                            n_categories = 2,
                            use_anisotropy = FALSE, # corresponds to ln_H_input params
                            ObsModel = c(1,1),
                            RhoConfig = c("Beta1" = 2, "Beta2" = 2,
                                          "Epsilon1" = 0, "Epsilon2" = 0))
  # settings$FieldConfig["Omega","Component_1"] <- 0
  settings$Options[c("Calculate_Range","Calculate_effective_area")] = TRUE
  settings$Options[c("zerosum_penalty")]=100
  
  
  settings$epu_to_use <- switch(epu.i,
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
                   working_dir = paste0(moddir,"/",sep=""),
                   getJointPrecision=TRUE,
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
  
  #plot
  
  map_list = make_map_info( Region = settings$Region,
                            spatial_list = fit$spatial_list,
                            Extrapolation_List = fit$extrapolation_list )
  
  results = plot( fit, what = "results",
                  map_list=map_list,
                  category_names = spp_list,
                  check_residuals=TRUE,
                  plot_set=c(3,16),
                  working_dir=paste0(moddir,"/"))
  
  # plot_factors( Report=fit$Report,
  #               ParHat=fit$ParHat,
  #               Data=fit$data_list,
  #               category_names = spp_list,
  #               mapdetails_list=results$map_list,
  #               Year_Set=fit$year_labels, plotdir=paste0(moddir,"/"))
  # 
 
  saveRDS(fit, file = paste0(moddir, "/fit.rds"))
  saveRDS(results, file = paste0(moddir, "/results.rds"))
  save.image(file = paste0(moddir, "/output.RData"))
  
}
