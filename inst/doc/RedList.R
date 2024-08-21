## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE,
  fig.width=8, 
  fig.height=5,
  warning = FALSE,
  message = FALSE)

## ----load packages------------------------------------------------------------
#  
#  library(intSDM)
#  library(INLA)
#  

## ----initialize workflow------------------------------------------------------
#  
#  workflow <- startWorkflow(
#          Projection = '+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=km +no_defs',
#          Species = c("Fraxinus_excelsior", "Ulmus_glabra", "Arnica_montana"),
#          saveOptions = list(projectName =  'Vascular'), Save = FALSE
#          )
#  

## ----addArea------------------------------------------------------------------
#  
#  yes <- FALSE
#  if (yes) {
#  Norway <- giscoR::gisco_get_countries(country = 'Norway', resolution = 60)
#  Norway <- st_cast(st_as_sf(Norway), 'POLYGON')
#  Norway <- Norway[which.max(st_area(Norway)),]
#  Norway <- rmapshaper::ms_simplify(Norway, keep = 0.8)
#  
#  workflow$addArea(Object = Norway)
#  workflow$plot()
#  }
#  
#  Norway <- readRDS('IntegratedLakefish/Norway.rds')
#  workflow$addArea(Object = Norway)

## ----addGBIF------------------------------------------------------------------
#  
#  workflow$addGBIF(datasetName = 'NTNU',
#                   datasetType = 'PA',
#                   limit = 10000,
#                   coordinateUncertaintyInMeters = '0,50',
#                   generateAbsences = TRUE,
#                   datasetKey = 'd29d79fd-2dc4-4ef5-89b8-cdf66994de0d')
#  
#  workflow$addGBIF(datasetName = 'UiO',
#                   datasetType = 'PA',
#                   limit = 10000,
#                   coordinateUncertaintyInMeters = '0,50',
#                   generateAbsences = TRUE,
#                   datasetKey = 'e45c7d91-81c6-4455-86e3-2965a5739b1f')
#  
#  workflow$addGBIF(datasetName = 'CZ',
#                   datasetType = 'PO',
#                   coordinateUncertaintyInMeters = '0,50',
#                   limit = 10000,
#                   datasetKey = 'b124e1e0-4755-430f-9eab-894f25a9b59c')
#  
#  workflow$plot(Species = TRUE)
#  

## ----addCovariates, eval = FALSE----------------------------------------------
#  
#  workflow$addCovariates(worldClim = 'tavg', res = 5, Function = scale)
#  workflow$plot(Covariates = TRUE)
#  

## ----metadata-----------------------------------------------------------------
#  
#  workflow$obtainMeta()
#  

## ----INLA---------------------------------------------------------------------
#  
#  workflow$addMesh(cutoff = 20 * 0.25,
#                   max.edge = c(60, 80)*0.5, #0.25
#                   offset= c(30, 40))
#  
#  workflow$plot(Mesh = TRUE)
#  

## ----Priors-------------------------------------------------------------------
#  
#  workflow$specifySpatial(prior.range = c(100, 0.1),
#                          prior.sigma = c(1, 0.1),
#                          constr = FALSE)
#  

## ----Fixed priors-------------------------------------------------------------
#  
#  workflow$specifyPriors(effectNames = 'Intercept',
#                         Mean = 0, Precision = 1)
#  
#  workflow$specifyPriors('tavg', Mean = 0, Precision = 1)
#  

## ----Bias---------------------------------------------------------------------
#  
#  workflow$biasFields('CZ',
#                      prior.range = c(100, 0.1), #1 #0.1
#                      prior.sigma = c(1, 0.1), #1 #0.1
#                      constr = FALSE)
#  

## ----options------------------------------------------------------------------
#  
#  #workflow$crossValidation(Method = 'Loo')
#  workflow$workflowOutput(c('Maps', 'Model', 'Bias'))
#  

## ----Maps---------------------------------------------------------------------
#  
#  Maps <- sdmWorkflow(workflow,inlaOptions = list(control.inla=list(int.strategy = 'ccd',
#                                                      strategy = 'gaussian',
#                                                      cmin = 0,
#                                                      diagonal = 0.1,
#                                                      control.vb=list(enable = FALSE)),
#                                    safe = TRUE,
#                                    verbose = TRUE,
#                                    inla.mode = 'experimental'),
#                predictionDim = c(400, 400))

## ----MapsOut------------------------------------------------------------------
#  
#  Maps$Fraxinus_excelsior$Maps
#  Maps$Ulmus_glabra$Maps
#  Maps$Arnica_montana$Maps
#  
#  saveRDS(Maps, 'Maps.rds')
#  

