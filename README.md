# Multivariate models to uncover hidden relationships between different cancer sites
This repository contains the R code to fit with INLA the M-models for multivariate spatio-temporal areal data described in the book chapter entitled _"Multivariate disease mapping models to uncover hidden relationships between different cancer sites"_ (Adin et al., 2023). 

We illustrate the methodology with a joint analysis of male mortality data for lung, colorectal, stomach and LOCP ((lip,
oral cavity and pharynx) cancer in continental Spain for the period 2006-2020. 

This repository also contains the functions to reproduce all the figures and tables of the chapter.


## Table of contents

- [Data](#Data)
- [R code](#R-code)
- [References](#References)


# Data
Lung cancer (ICD-10 codes C33 and C34), colorectal cancer (ICD-10 codes C17-C21), stomach cancer (ICD-10 codes C16) and LOCP cancer (ICD-10 codes C00-C14) deaths for male population in the 47 provinces of continental Spain during the period 2006-2020. The data are publicly available online without any form of restriction or copyright.

The [**CancerData_SpainPROV.Rdata**](https://github.com/spatialstatisticsupna/BookChapter_STMmodels/blob/master/R/CancerData_SpainPROV.Rdata) file contains the following objects:
  - **W**: spatial adjacency binary matrix of the provinces of continental Spain
  - **carto**: `sf` object containing the polygons of the provinces of continental Spain and 3 variables
    - **_ID_**: character vector of geographic identifiers
    - **_NAME_**: character vector of province names
    - **_geometry_**: sfc_GEOMETRY
  - **data**: list of `data.frames` objects corresponding to four different cancer types. Each data.frame contains the following variables:
    - **_ID_**: character vector of geographic identifiers  
    - **_Year_**: numeric vector of year’s identifiers
    - **_O_**: observed number of cancer deaths
    - **_E_**: expected number of cancer deaths
    - **_SMR_**: standardized mortality ratio


Use the following commands to load the data
```r 
> load("R/CancerData_SpainPROV.Rdata")

> head(carto)
Simple feature collection with 6 features and 2 fields
Geometry type: GEOMETRY
Dimension:     XY
Bounding box:  xmin: 161384 ymin: 4059652 xmax: 781118 ymax: 4838774
CRS:           +proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs
  ID     NAME                       geometry
1 01    Álava POLYGON ((539704 4705573, 5...
2 02 Albacete POLYGON ((539462 4215448, 5...
3 03 Alicante MULTIPOLYGON (((688545 4195...
4 04  Almería POLYGON ((523167 4060037, 5...
5 33 Asturias MULTIPOLYGON (((262649 4764...
6 05    Ávila POLYGON ((283218 4452970, 2...

> str(data,2)
List of 4
 $ Lung      :'data.frame':	705 obs. of  5 variables:
  ..$ ID  : chr [1:705] "01" "01" "01" "01" ...
  ..$ Year: int [1:705] 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 ...
  ..$ O   : num [1:705] 92 107 120 119 91 130 105 128 121 115 ...
  ..$ E   : num [1:705] 108 111 114 117 120 ...
  ..$ SMR : num [1:705] 0.851 0.963 1.053 1.017 0.758 ...
 $ Colorectal:'data.frame':	705 obs. of  5 variables:
  ..$ ID  : chr [1:705] "01" "01" "01" "01" ...
  ..$ Year: int [1:705] 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 ...
  ..$ O   : num [1:705] 63 59 65 65 72 68 70 69 57 47 ...
  ..$ E   : num [1:705] 53.2 55.1 56.9 58.8 60.7 ...
  ..$ SMR : num [1:705] 1.18 1.07 1.14 1.1 1.19 ...
 $ Stomach   :'data.frame':	705 obs. of  5 variables:
  ..$ ID  : chr [1:705] "01" "01" "01" "01" ...
  ..$ Year: int [1:705] 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 ...
  ..$ O   : num [1:705] 27 26 33 32 24 44 33 37 40 36 ...
  ..$ E   : num [1:705] 20.5 21.2 21.9 22.5 23.2 ...
  ..$ SMR : num [1:705] 1.31 1.23 1.51 1.42 1.03 ...
 $ LOCP      :'data.frame':	705 obs. of  5 variables:
  ..$ ID  : chr [1:705] "01" "01" "01" "01" ...
  ..$ Year: int [1:705] 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 ...
  ..$ O   : num [1:705] 15 14 12 9 17 17 9 13 10 10 ...
  ..$ E   : num [1:705] 10.8 11.1 11.4 11.6 11.9 ...
  ..$ SMR : num [1:705] 1.384 1.262 1.056 0.774 1.428 ...
```

# R code
R code to fit the M-models for multivariate spatio-temporal areal data with R-INLA (http://www.r-inla.org/) considered in the present chapter, and code to reproduce all the figures and tables. All the R files are written by the authors of the paper.

- [**Fit_models.R**](https://github.com/spatialstatisticsupna/BookChapter_STMmodels/blob/master/R/Fit_models.R)

  Main function to fit with R-INLA the different multivariate spatio-temporal models described in the chapter (see Section 2 for further details about _Model 1_ and _Model 2_).
  
- [**MCAR_INLA_ST_Model1.R**](https://github.com/spatialstatisticsupna/BookChapter_STMmodels/blob/master/R/MCAR_INLA_ST_Model1.R) and [**MCAR_INLA_ST_Model2.R**](https://github.com/spatialstatisticsupna/BookChapter_STMmodels/blob/master/R/MCAR_INLA_ST_Model2.R)

  R functions to fit our M-model proposals using different prior distributions for the main spatial and temporal random effects and four different types of space-time interactions. A similar implementation was originally proposed by [Vicente et al. (2020)](https://doi.org/10.1007/s00477-020-01808-x).
  
  It defines the `MCAR_INLA_ST()` function with the following arguments:
  
  * `carto`: object of class `sf` that must contain at least the variable with the identifiers of the spatial areal units specified in the `ID.area` argument.
  * `data`: object of class `data.frame` that must contain the target variables of interest specified in the arguments `ID.area`, `ID.year`, `O` and `E`.
  * `ID.area`: character; name of the variable that contains the IDs of spatial areal units. The values of this variable must match those given in the `carto` and `data` variable.
  * `ID.year`: character; name of the variable that contains the IDs of time points.
  * `ID.disease`: character; name of the variable that contains the IDs of the diseases.
  * `O`: character; name of the variable that contains the observed number of cases for each areal unit, time point and disease
  * `E`: character; name of the variable that contains the expected number of cases for each areal unit, time point and disease
  * `W`: optional argument with the binary adjacency matrix of the spatial areal units. If NULL (default), this object is computed from the carto argument (two areas are considered as neighbours if they share a common border).
  * `spatial`: one of either "Leroux" (default), "intrinsic" or "proper", which specifies the prior distribution considered for the spatial random effect.
  * `temporal`: one of either "rw1" (default) or "rw2", which specifies the prior distribution considered for the temporal random effect.
  * `interaction`: one of either "none", "TypeI", "TypeII", "TypeIII" or "TypeIV" (default), which specifies the prior distribution for the space-time interaction random effect.
  * `strategy`: one of either "gaussian", "simplified.laplace" (default), "laplace" or "adaptive", which specifies the approximation strategy considered in the inla
function. Only valid when using the inla.mode="classic" approximation technique.
  
- [**Figures_and_Tables.R**](https://github.com/spatialstatisticsupna/BookChapter_STMmodels/blob/master/R/Figures_and_Tables.R)

  R code that contains the necessary functions to replicate the figures and tables of the present chapter.
  
Additional auxiliary functions to define the M-model for the spatial/temporal random effects using the `rgeneric`construction of R-INLA are also given in the [**R/functions/**](https://github.com/spatialstatisticsupna/BookChapter_STMmodels/blob/master/R/functions/) folder. See [Vicente et al. (2022)](https://arxiv.org/abs/2210.14849) for further details about the internal parameterization of the corresponding between-disease covariance matrices.
  

# Acknowledgements
The authors would like to thank the Spanish Statistical Office (INE) and the Spanish National Epidemiology Center (area of Environmental Epidemiology and Cancer) for providing the data. This work has been supported by the project PID2020-113125RBI00/MCIN/AEI/10.13039/501100011033.

![plot](https://github.com/spatialstatisticsupna/bigDM/blob/master/micin-aei.jpg)


# References
Adin, A., Goicoa, T., and Ugarte, M.D. (2023). Multivariate disease mapping models to uncover hidden relationships between different cancer sites. _Statistical Methods at the Forefront of Biomedical Advances._

[Vicente, G., Goicoa, T., and Ugarte, M.D. (2020). Bayesian inference in multivariate spatio-temporal areal models using INLA: analysis of gender-based violence in small areas. _Stochastic Environmental Research and Risk Assessment_, vol. 34, 1421-1440, 2020.](https://doi.org/10.1007/s00477-020-01808-x).

[Vicente, G., Adin, A., Goicoa, T., and Ugarte, M.D. (2022). High-dimensional order-free multivariate spatial disease mapping. _arXiv preprint_.](https://arxiv.org/abs/2210.14849)
