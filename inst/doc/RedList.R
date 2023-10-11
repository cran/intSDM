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
#          Projection = '+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs',
#          Species = c("Fraxinus_excelsior", "Ulmus_glabra", "Arnica_montana"),
#          saveOptions = list(projectName =  'Vascular'), Save = FALSE
#          )
#  

## ----addArea------------------------------------------------------------------
#  
#  workflow$addArea(countryName = 'Norway', resolution = '60')
#  workflow$plot()
#  

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
#  workflow$addMesh(cutoff = 20000,
#                   max.edge = c(60000, 80000),
#                   offset= 100000)
#  
#  workflow$plot(Mesh = TRUE)
#  

## ----Priors-------------------------------------------------------------------
#  
#  workflow$specifySpatial(prior.range = c(300000, 0.05),
#                          prior.sigma = c(50, 0.2))
#  

## ----Bias---------------------------------------------------------------------
#  
#  workflow$biasFields('CZ')
#  

## ----options------------------------------------------------------------------
#  
#  workflow$workflowOutput('Maps')
#  
#  workflow$modelOptions(INLA = list(control.inla=list(int.strategy = 'eb',
#                                                      cmin = 0),
#                                    safe = TRUE,
#                                    inla.mode = 'experimental'))
#  

## ----Maps---------------------------------------------------------------------
#  
#  Maps <- sdmWorkflow(workflow)
#  Maps$Fraxinus_excelsior$Maps
#  Maps$Ulmus_glabra$Maps
#  Maps$Arnica_montana$Maps
#  

