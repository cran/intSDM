## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE,
  fig.width=8, 
  fig.height=5,
  warning = FALSE,
  message = FALSE)

## ----setup--------------------------------------------------------------------
#  
#  library(intSDM)
#  library(USAboundaries)
#  

## ----get data-----------------------------------------------------------------
#  
#  BBA <- PointedSDMs::SetophagaData$BBA
#  BBA$Species_name <- paste0('Setophaga_', BBA$Species_name)
#  BBS <- PointedSDMs::SetophagaData$BBS
#  BBS$Species_name <- paste0('Setophaga_', BBS$Species_name)
#  

## ----startWorkflow------------------------------------------------------------
#  
#  workflow <- startWorkflow(
#    Projection = "WGS84",
#    Species = c("Setophaga_caerulescens", "Setophaga_fusca", "Setophaga_magnolia"),
#    saveOptions = list(projectName =  'Setophaga'), Save = FALSE
#  )
#  

## ----addArea------------------------------------------------------------------
#  
#  workflow$addArea(Object = USAboundaries::us_states(states = "Pennsylvania"))
#  

## ----download Data------------------------------------------------------------
#  
#  workflow$addGBIF(datasetName = 'eBird', datasetType = 'PO', limit = 5000,
#                   datasetKey = '4fa7b334-ce0d-4e88-aaae-2e0c138d049e')
#  
#  workflow$addStructured(dataStructured = BBA, datasetType = 'PA',
#                         responseName = 'NPres',
#                         speciesName = 'Species_name')
#  
#  workflow$addStructured(dataStructured = BBS, datasetType = 'Counts',
#                         responseName = 'Counts',
#                         speciesName = 'Species_name')
#  
#  workflow$plot(Species = TRUE)
#  

## ----addCovariates------------------------------------------------------------
#  
#  covariates <- scale(terra::rast(system.file('extdata/SetophagaCovariates.tif',
#                                        package = "PointedSDMs")))
#  names(covariates) <- c('elevation', 'canopy')
#  
#  workflow$addCovariates(Object = covariates)
#  
#  workflow$plot(Covariates = TRUE)
#  

## ----biasFields---------------------------------------------------------------
#  
#  workflow$biasFields(datasetName  = 'eBird')
#  
#  workflow$addMesh(cutoff = 0.2,
#                   max.edge = c(0.1, 0.24),
#                   offset = c(0.1, 0.4))
#  

## ----outcomes-----------------------------------------------------------------
#  
#  workflow$workflowOutput('Model')
#  
#  workflow$modelOptions(INLA = list(control.inla=list(int.strategy = 'eb',
#                                                      cmin = 0),
#                                    safe = TRUE,
#                                    inla.mode = 'experimental'))
#  

## ----sdmWorkflow--------------------------------------------------------------
#  
#  Models <- sdmWorkflow(Workflow = workflow)
#  
#  lapply(unlist(Models, recursive = FALSE), summary)
#  

