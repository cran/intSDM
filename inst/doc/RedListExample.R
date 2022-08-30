## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  eval = FALSE,
  collapse = TRUE,
  comment = "#>",
  fig.width=8, fig.height=5
)

## ----Install package, eval = FALSE--------------------------------------------
#  
#  if (requireNamespace(intSDM)) devtools::install_github('PhilipMostert/intSDM')
#  

## ----Load packages, warning = FALSE, message = FALSE--------------------------
#  
#  library(ggplot2)
#  library(sf)
#  library(sp)
#  library(dplyr)
#  library(ggmap)
#  library(maps)
#  library(PointedSDMs)
#  library(spatstat)
#  library(maptools)
#  library(INLA)
#  library(intSDM)
#  library(rgeos)
#  library(fields)
#  library(viridis)
#  library(ggpolypath)
#  library(RColorBrewer)
#  

## ----Norway map, warning = FALSE, message = FALSE-----------------------------
#  
#  Projection <- CRS('+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')
#  
#  norwayfill <- map("world", "norway", fill=TRUE, plot=FALSE,
#                    ylim=c(58,72), xlim=c(4,32))
#  IDs <- sapply(strsplit(norwayfill$names, ":"), function(x) x[1])
#  norway.poly <- map2SpatialPolygons(norwayfill, IDs = IDs,
#                                     proj4string = Projection)
#  

## ----Read PA data, warning = FALSE, message = FALSE---------------------------
#  
#  data("PA_redlist")
#  

## ----Plot of PA, warning = FALSE, message = FALSE-----------------------------
#  
#  ggplot() +
#    gg(norway.poly) +
#    gg(PA_redlist, aes(col = factor(individualCount))) +
#    facet_grid(~species) +
#    labs(x = 'Longitude', y = 'Latitude', col = 'Grid Observation') +
#    scale_color_manual(labels = c('Absent', "Present"), values = c("#d11141", "#00aedb")) +
#    ggtitle('Vascular Plant Field Notes') +
#    theme_classic() +
#    theme(legend.position="bottom",
#          plot.title = element_text(hjust = 0.5))
#  
#  

## ----Structured data, warning = FALSE, message = FALSE------------------------
#  
#  structured <- structured_data(PA_redlist, datasetType = 'PA',
#                                speciesName = 'species',
#                                responsePA = 'individualCount',
#                                coordinateNames = colnames(PA_redlist@coords))
#  

## ----Mesh construction, warning = FALSE, message = FALSE----------------------
#  
#  mesh <- species_model(boundary = norway.poly,
#                        return = 'mesh', limit = 5000, meshParameters = list(cutoff=0.08, max.edge=c(1, 3), offset=c(1,1)))
#  
#  ggplot() +
#    gg(mesh) +
#        ggtitle('inla.mesh object') +
#    theme_classic() +
#    theme(legend.position="bottom",
#          plot.title = element_text(hjust = 0.5))
#  

## ----All species plot, warning = FALSE, message = FALSE-----------------------
#  
#  species_plot <- species_model(speciesNames = unique(structured@dataPA$PA_redlist$species),
#                                 structuredData = structured,
#                                 boundary = norway.poly,
#                                 return = 'species plot', limit = 5000, mesh = mesh)
#  species_plot +
#        ggtitle('Plot of the species data') +
#    theme_classic() +
#    theme(legend.position="bottom",
#          plot.title = element_text(hjust = 0.5))
#  

## ----Prediction maps, warning = FALSE, message = FALSE, eval = FALSE, fig.width=8, fig.height=5----
#  
#  prediction_maps <- species_model(speciesNames = unique(structured@dataPA$PA_redlist$species),
#                                 scale = TRUE, structuredData = structured,
#                                 worldclimCovariates = 'Annual Mean Temperature',
#                                 boundary = norway.poly,
#                                 return = 'predictions map', limit = 5000,
#                                 mesh = mesh)
#  
#  prediction_maps
#  

