rm(list=ls())

# # To install the latest stable version of the R-INLA package:
# install.packages("INLA",repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

library(dplyr)
library(fastDummies)
library(INLA)
library(sf)
library(spatialreg)
library(spdep)


########################################
## 1) Load data and cartography files ##
########################################
load("CancerData_SpainPROV.Rdata")

head(carto)
str(data,2)

J <- length(data)
S <- length(unique(data[[1]]$ID))
T <- length(unique(data[[1]]$Year))
t.from <- min(data[[1]]$Year)
t.to <- max(data[[1]]$Year)

## Combine data sets by row ##
data <- do.call(rbind,data)
data$Disease <- rep(1:J, each=S*T)
str(data)


##############################################
## 2) Fit spatio-temporal M-models (Leroux) ##
##############################################
# Model <- "M1"
Model <- "M2"

if(Model=="M1") source("MCAR_INLA_ST_Model1.R")
if(Model=="M2") source("MCAR_INLA_ST_Model2.R")

## Set the "compact" mode of INLA ##
inla.setOption(inla.mode="compact")
inla.getOption()$inla.mode

## INLA models ##
LCAR.t1 <- MCAR_INLA_ST(carto=carto, data=data, ID.area="ID", ID.year="Year", ID.disease="Disease",
                        O="O", E="E", spatial="Leroux", temporal="rw1", interaction="TypeI")

LCAR.t2 <- MCAR_INLA_ST(carto=carto, data=data, ID.area="ID", ID.year="Year", ID.disease="Disease",
                        O="O", E="E", spatial="Leroux", temporal="rw1", interaction="TypeII")

LCAR.t3 <- MCAR_INLA_ST(carto=carto, data=data, ID.area="ID", ID.year="Year", ID.disease="Disease",
                        O="O", E="E", spatial="Leroux", temporal="rw1", interaction="TypeIII")

LCAR.t4 <- MCAR_INLA_ST(carto=carto, data=data, ID.area="ID", ID.year="Year", ID.disease="Disease",
                        O="O", E="E", spatial="Leroux", temporal="rw1", interaction="TypeIV")

## Save the models ##
save(list=c("LCAR.t1","LCAR.t2","LCAR.t3","LCAR.t4"),
     file=paste("INLAmodels_MCAR_Leroux_",eval(Model),".Rdata",sep=""))
