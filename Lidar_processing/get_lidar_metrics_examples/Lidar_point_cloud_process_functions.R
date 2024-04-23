
library(lidR)
library(raster)
library(sf)
library(foreach)
library(doParallel)
library(dplyr)
library(Lmoments)


#' derive many lidar metrics for the input las file (using area-based approach (ABA), i.e. the input las file is for a given plot or area or AOI)
#' the lidar metrics include: height density of all returns, height percentile of first returns, height dispersion, canopy height, 
#' percentage of canopy cover (above 2m, 6m, 10m, 14m), Lmoments, leaf area density, and rumple index. In addition to lidar metrics, it also
#' return information like year and day of year the file is created and the density of point could 
#'  
#' @param las, las for the AOI
#' @param idnm, string, id of the AOI
#' @return a dataframe with idnm as row names, and all lidar metrics as columns, use colnames(output) to get names of all lidar metrics, where output is the output dataframe, notice Year, Day_Year, and Density columns are not lidar metrics


aba_metrics = function(las, idnm){
  
  hds = header(las)
  year = hds$`File Creation Year`
  day_year = hds$`File Creation Day of Year`
  density = hds$`Number of point records`/(lidR::area(las))
  las_f = filter_first(las) 
  
  ## basic1 is a list of lidar metrics based on all returns, including lidar density metric  and relative height metric  (refer publication here: https://www.sciencedirect.com/science/article/pii/S0034425713002125)
  basic1= unlist(cloud_metrics(las, ~ld_metrics(Z)))
  
  ## basic2 is a lit of lidar metrics based on first returns only (i.e. canopy top), including percentile of returned heights
  basic2 = unlist(cloud_metrics(las_f, ~ld_metrics_cp(Z)))
  
  ## refer here https://github.com/wanwanliang/lidRmetrics/blob/main/R/metrics_dispersion.R for explanation on dispersion metrics
  disp_metrics = unlist(cloud_metrics(las, ~metrics_dispersion(Z)))
  
  
  chm = mean(las_f$Z, na.rm=T)
  ccs = unlist(cloud_metrics(las_f, ~canopy_cover(Z)))

  ## refer here https://github.com/wanwanliang/lidRmetrics/blob/main/R/metrics_Lmoments.R for explanation on Lmoments metrics
  lmoments = unlist(cloud_metrics(las, ~metrics_Lmoments(Z)))
  
  ## refer here https://github.com/wanwanliang/lidRmetrics/blob/main/R/metrics_lad.R for explanation on LAD metrics
  lads = unlist(cloud_metrics(las, ~metrics_lad(Z)))
  #rump = unlist(cloud_metrics(las, ~metrics_rumple(X,Y,Z,0.5))) 
  #names(rump) = 'rumple'


  metrics_all = c(chm, basic1, basic2, disp_metrics,  lmoments, lads,  ccs, year, day_year, density)
  names(metrics_all)[1] = 'CHM'
  nl = length(metrics_all)
  names(metrics_all)[c((nl-2):nl)] = c('Year','Day_Year','Density')
  nms = names(metrics_all)
  ld_raA = as.data.frame(metrics_all)
  ld_raA = t(ld_raA)
  colnames(ld_raA) = nms
  rownames(ld_raA) = idnm
  
  return(ld_raA)
}


#'similar with aba_metrics, but the function below derives lidar metrics as rasterstack with given resolution
#'  
#' @param las, las for the AOI
#' @param res, integer, resolution of output lidar metrics 
#' @return a raster stack with all 48 lidar metrics


get_metricsAll_raster = function(las, res){
  #las <- readLAS(chunk)

  ## remove noise
  las = las[!las$Classification%in%c(6, 7,9, 10, 11, 12,13,14,15,16,17,18),]
  las =  classify_noise(las, ivf(3,2))
  las <- filter_poi(las, (Classification != LASNOISE))

  ## get normalized height after removing DTM
  las <- normalize_height(las, tin(extrapolate = knnidw(10, 2, 5)))
  ## remove some potential noise or outliers
  las =  filter_poi(las, Z >= 0 & Z < min(100, quantile(Z, .9999)))
  ## get the first return (i.e. canopy top returns)
  las_f = filter_first(las)

  ## refer to function aba_metrics to get more details on meaning of lidar metrics below; grid_metrics is a lidR function to return lidar metrics as grid
  basic1 <- grid_metrics(las, res= res, ~ld_metrics(Z))
  basic2 <- grid_metrics(las_f, res= res, ~ld_metrics_cp(Z))
  disp_metrics = grid_metrics(las, res=res, ~metrics_dispersion(Z))
  
  
  chm = raster(rasterize_canopy(las, res=res, algorithm = p2r()))
  names(chm) = 'CHM'
  ccs = grid_metrics(las_f,res= res, ~canopy_cover(Z))
  
  
  lmoments = grid_metrics(las, res= res,~metrics_Lmoments(Z))
  lads = grid_metrics(las, res=res, ~metrics_lad(Z))
  #rump = grid_metrics(las, res= res,~metrics_rumple(X,Y,Z,0.5))  # this metric takes too much time
  #names(rump) = 'rumple'
 
  
  var_nms = c(c('CHM'),names(basic1), names(basic2), names(disp_metrics), names(lmoments),
              names(lads),  names(ccs))

  tryCatch(
    
    {
      ld_raA = stack(chm, basic1, basic2, disp_metrics,  lmoments, lads,   ccs)
      return(ld_raA)
    },
    error = function(error_message) {
      basic1 = projectRaster(basic1, chm,method = 'bilinear')
      basic2 = projectRaster(basic2, chm,method = 'bilinear')
      disp_metrics = projectRaster(disp_metrics, chm,method = 'bilinear')
      lmoments = projectRaster(lmoments, chm,method = 'bilinear')
      lads = projectRaster(lads, chm,method = 'bilinear')
      #rump = projectRaster(rump, chm,method = 'bilinear')
      ccs = projectRaster(ccs, chm,method = 'bilinear')
      
      ld_raA = stack(chm, basic1, basic2, disp_metrics,  lmoments, lads,   ccs)
      return(ld_raA)
    }
    
  )
  
}


#' filter and normalize las file
#' remove points that are identifed as noise, buildings, etc. based on classification here: https://desktop.arcgis.com/en/arcmap/latest/manage-data/las-dataset/lidar-point-classification.htm
#' remove noise identified by the classify_noise function
#' @param las a las object read by readLAS or las_catalog or  clip_roi output
#' @return L-moments and L-moment ratios

normalize_las = function(las){
  las <- filter_duplicates(las)
  las = las[!las$Classification%in%c(6, 7,9, 10, 11, 12,13,14,15,16,17,18),]
  las =  classify_noise(las, ivf(3,2))
  las <- filter_poi(las, (Classification != LASNOISE))
  las <- normalize_height(las, tin(extrapolate = knnidw(10, 2, 5)))
  las =  filter_poi(las, Z >= 0 & Z < min(100, quantile(Z, .99999)))
  return(las)
}



##### Below are functions that are used by aba_metrics and get_metricsAll_raster functions to calculate different types of lidar metrics

canopy_cover = function(rasterLayer){
  pixelsGreaterThan2 <- rasterLayer > 2
  pixelsGreaterThan6 <- rasterLayer > 6
  pixelsGreaterThan14 <- rasterLayer > 14
  pixelsGreaterThan10 <- rasterLayer > 10
  pixelsGreaterThan20 <- rasterLayer > 20
  
  
  # Calculate the number of pixels greater than 2
  numPixelsGreaterThan2 <- sum(pixelsGreaterThan2[], na.rm = TRUE)
  numPixelsGreaterThan6 <- sum(pixelsGreaterThan6[], na.rm = TRUE)
  numPixelsGreaterThan10 <- sum(pixelsGreaterThan10[], na.rm = TRUE)
  numPixelsGreaterThan14 <- sum(pixelsGreaterThan14[], na.rm = TRUE)
  #numPixelsGreaterThan20 <- sum(pixelsGreaterThan20[], na.rm = TRUE)
  
  # Calculate the total number of pixels (excluding NA values)
  totalPixels <- sum(!is.na(rasterLayer[]))
  
  # Calculate the percentage
  cc2 = (numPixelsGreaterThan2 / totalPixels) 
  cc6 = (numPixelsGreaterThan6 / totalPixels) 
  cc10 = (numPixelsGreaterThan10 / totalPixels) 
  cc14 = (numPixelsGreaterThan14 / totalPixels) 
  #cc20 = (numPixelsGreaterThan20 / totalPixels) 
  outs = list(canopy_cover2=cc2, canopy_cover6= cc6, canopy_cover10=cc10,
              canopy_cover14= cc14)
  return(outs)
}

ld_metrics = function(Z){
  ah = mean(Z, na.rm=T)
  ahstd = sd(Z, na.rm=T)
  ach = mean(Z[Z>1.3], na.rm=T)
  achstd = sd(Z[Z>1.3], na.rm=T)
  qach = ach^2
  if (is.na(ach)){
    ach = 0
    achstd = 0
    qach = 0
  }
  
  bins = seq(1.3, max(Z), length.out=11)
  
  d30 = length(Z[Z>bins[4]])/length(Z)
  d40 = length(Z[Z>bins[5]])/length(Z)
  d60 = length(Z[Z>bins[7]])/length(Z)
  d80 = length(Z[Z>bins[9]])/length(Z)
  
  probs= c(0.1, 0.3, 0.5, 0.7, 0.8,0.98)
  qts = quantile(Z, probs= probs)
  rhs = list()
  for (i in 1:length(qts)){
    nm = paste('RH_', probs[i]*100, sep="")
    rhs[[nm]] = qts[i]
  }
  
  outs = list(meanZ= ah, sd_Z= ahstd, average_canopy= ach, sd_canopy= achstd, q_ach= qach, d30= d30, d40= d40,  d60= d60, d80= d80,
              rh10= rhs[['RH_10']], rh30= rhs[['RH_30']], rh50= rhs[['RH_50']], rh70= rhs[['RH_70']],
              rh80= rhs[['RH_80']], rh98= rhs[['RH_98']])
  return(outs)
}

ld_metrics_cp = function(Z){
  probs= c(0.1,0.2,0.3, 0.5,0.6,0.7,0.9,0.95)
  qts = quantile(Z, probs= probs)
  outs = list()
  for (i in 1:length(qts)){
    nm = paste('h_', probs[i]*100, sep="")
    outs[[nm]] = qts[i]
  }
  
  return(outs)
}

metrics_rumple <- function(x, y, z, pixel_size, zmin=NA) {
  
  if (!is.na(zmin)) z <- z[z>zmin]
  
  r <- NA_real_
  
  if (length(z) > 2) {
    
    D <-  data.table::data.table(X=x, Y=y, Z=z)
    
    D <- LAS(D, header = rlas::header_create(D), check=F)
    
    D <- lidR::decimate_points(D, lidR::highest(pixel_size))
    
    r <- lidR::rumple_index(x = D$X, y = D$Y, z = D$Z)
  }
  
  return(list(rumple=r))
  
}


### from https://github.com/ptompalski/lidRmetrics
metrics_dispersion <- function(z, dz=1, zmin=0) {
  
  
  m = list(
    
    ziqr = IQR(z), # Interquartile distance 
    
    
    
    #MAD around the mean
    zMADmean = mean(abs(z - mean(z))),
    
    # MAD around the median - Median absolute deviation (median of the absolute deviations from the data's median)
    zMADmedian = median(abs(z - median(z))),
    
    # Canopy relief ratio ((mean - min) / (max – min))
    CRR = ((mean(z) - min(z)) / (max(z) - min(z))),
    
    zentropy = entropy(z, dz), #forcing z to be always above 0, otherwise entropy does not work
    
    VCI = VCI(z[z>0], zmax = max(z), by = dz) #forcing z to be always above 0, otherwise entropy does not work
  )
  return(m)
}

### from https://github.com/ptompalski/lidRmetrics
metrics_lad <- function(z, zmin=NA, dz = 1, k = 0.5, z0 = 2) {
  
  if (!is.na(zmin)) z <- z[z>zmin]
  
  lad_max <- lad_mean <- lad_cv <- lad_min <- lai <- NA_real_
  
  if(length(z) > 2) {
    
    ladprofile <- lidR::LAD(z, dz = dz, k = k, z0 = z0)
    
    lad_max <- with(ladprofile, max(lad, na.rm = TRUE))
    lad_mean <- with(ladprofile, mean(lad, na.rm = TRUE))
    lad_cv <- with(ladprofile, sd(lad, na.rm=TRUE)/mean(lad, na.rm = TRUE))
    lad_min <- with(ladprofile, min(lad, na.rm = TRUE))
    lai <- with(ladprofile, sum(lad, na.rm = TRUE))
    
  }
  
  lad_metrics <- list(lad_max = lad_max,
                      lad_mean = lad_mean,
                      lad_cv = lad_cv,
                      lad_min = lad_min,
                      lai = lai)
  
  return(lad_metrics)
}


### from https://github.com/ptompalski/lidRmetrics
#' L-moments and L-moment ratios
#' 
#' Calculates L-moments and L-moment ratios of point cloud heihts.
#' 
#' @param z Z coordinate of the point cloud
#' @param zmin Minimum height. If set, heights below are ignored in calculations.
#' @return L-moments and L-moment ratios
#' @export

metrics_Lmoments <- function(z, zmin=NA) {
  # Lmoments - code from Murray Woods. Modified.
  
  if (!requireNamespace("Lmoments", quietly = TRUE)) {
    stop("Package \"Lmoments\" needed for this function to work.",
         call. = FALSE)
  }
  
  
  if (!is.na(zmin)) z <- z[z>zmin]
  
  lmom_temp <- list(L1 = as.numeric(NA), L2 = as.numeric(NA), L3 = as.numeric(NA), L4 = as.numeric(NA),
                    Lskew = as.numeric(NA), Lkurt = as.numeric(NA), Lcoefvar = as.numeric(NA))
  
  if (length(z) >= 1) {
    
    lmom <- as.numeric(Lmoments::Lmoments(z))
    
    # an error occurs in some cases and Lmoments are not calcualted or only some are calculated
    # Code below is to check the output of the Lmoments and store it only if 4 values are present
    
    if(length(lmom)==4) {
      names(lmom) <- paste0("L",1:4)
      
      #calculate Lmoments ratios:
      Lskew = lmom[3] / lmom[2]
      Lkurt = lmom[4] / lmom[2]
      Lcoefvar = lmom[2] / lmom[1]
      
      lmom <- as.list(lmom)
      
      lmom[["Lskew"]] <- as.numeric(Lskew)
      lmom[["Lkurt"]] <- as.numeric(Lkurt)
      lmom[["Lcoefvar"]] <- as.numeric(Lcoefvar)
      
    } else {
      lmom <- lmom_temp
    }
  } else {
    lmom <- lmom_temp
  }
  return(lmom)
}


### from https://github.com/ptompalski/lidRmetrics
#' Basic metrics
#' 
#' Most common descriptive statistics used to characterize the vertical distribution of a point cloud.
#' 
#' @param z Z coordinate of the point cloud
#' @param zmin Minimum height. If set, heights below are ignored in calculations.
#' @return A set of descriptive statistics including: total number of points, maximum height, minimum height, mean height, 
#' standard deviation of height, coefficient of variation of height, skewness and kurtosis of height
#' @export
#' 
#' @examples
#' library(lidR)
#' library(lidRmetrics)
#' LASfile <- system.file("extdata", "Megaplot.laz", package="lidR")
#' las <- readLAS(LASfile, select = "*", filter = "-keep_random_fraction 0.5")
#' 
#' m1 <- cloud_metrics(las, ~metrics_basic(z = Z))
#' 
#' m2 <- grid_metrics(las, ~metrics_basic(z = Z), res = 20)

metrics_basic <- function(z, zmin=NA) {
  
  if (!is.na(zmin)) z <- z[z>zmin]
  
  n <- length(z)
  zmax <- max(z)
  zminimum <- min(z) #this is to confirm if any threshold was applied
  zmean <- mean(z)
  zsd <- stats::sd(z)
  zcv <- zsd / zmean * 100
  zskew <- (sum((z - zmean)^3)/n)/(sum((z - zmean)^2)/n)^(3/2)
  zkurt <- n * sum((z - zmean)^4)/(sum((z - zmean)^2)^2)
  
  return(list(
    n=n,
    zmax=zmax,
    zmin = zminimum,
    zmean = zmean, 
    zsd = zsd, 
    zcv = zcv,
    zskew = zskew,
    zkurt = zkurt
  )
  )
  
}


### code below all from https://github.com/ptompalski/lidRmetrics
#' Dispersion metrics
#' 
#' Metrics characterizing variation of point cloud heights. 
#' 
#' @param z Z coordinate of the point cloud
#' @param dz layer thickness to use when calculating entropy and VCI.
#' @param zmin Minimum height. If set, heights below are ignored in calculations.
#' @return Interquartile distance, mean absolute deviation (MAD) around the mean and median, canopy relief ratio
#' @export
#' 
#' @examples
#' library(lidR)
#' library(lidRmetrics)
#' LASfile <- system.file("extdata", "Megaplot.laz", package="lidR")
#' las <- readLAS(LASfile, select = "*", filter = "-keep_random_fraction 0.5")
#' 
#' m1 <- cloud_metrics(las, ~metrics_dispersion(z = Z))
#' 
#' m2 <- grid_metrics(las, ~metrics_dispersion(z = Z), res = 40)


metrics_dispersion <- function(z, dz=1, zmin=NA) {
  
  if (!is.na(zmin)) z <- z[z>zmin]
  
  m = list(
    
    ziqr = IQR(z), # Interquartile distance 
    
    
    # AAD (Average Absolute Deviation):
    # https://en.wikipedia.org/wiki/Average_absolute_deviation
    
    #MAD around the mean
    zMADmean = mean(abs(z - mean(z))),
    
    # MAD around the median - Median absolute deviation (median of the absolute deviations from the data's median)
    zMADmedian = median(abs(z - median(z))),
    
    # Canopy relief ratio ((mean - min) / (max – min))
    CRR = ((mean(z) - min(z)) / (max(z) - min(z))),
    
    zentropy = entropy(z[z>0], dz), #forcing z to be always above 0, otherwise entropy does not work
    
    VCI = VCI(z[z>0], zmax = max(z), by = dz) #forcing z to be always above 0, otherwise entropy does not work
  )
  return(m)
}




#' LAD metrics
#' 
#' Metrics based on the leaf area density. \code{lidR::LAD()} used to calculate the leaf area density. 
#' 
#' @param z Z coordinate of the point cloud
#' @param zmin Minimum height. If set, heights below are ignored in calculations.
#' @param dz numeric. The thickness of the layers used (height bin)
#' @param k numeric. is the extinction coefficient
#' @param z0 numeric. The bottom limit of the profile
#' @return min, max, mean, cv and sum (LAI), of the leaf area density profile.
#' @export
#' 
#' @examples
#' library(lidR)
#' library(lidRmetrics)
#' LASfile <- system.file("extdata", "Megaplot.laz", package="lidR")
#' las <- readLAS(LASfile, select = "*", filter = "-keep_random_fraction 0.5")
#' 
#' m1 <- cloud_metrics(las, ~metrics_lad(z = Z))
#' 
#' m2 <- grid_metrics(las, ~metrics_lad(z = Z), res = 40)


metrics_lad <- function(z, zmin=NA, dz = 1, k = 0.5, z0 = 2) {
  
  if (!is.na(zmin)) z <- z[z>zmin]
  
  lad_max <- lad_mean <- lad_cv <- lad_min <- lai <- NA_real_
  
  if(length(z) > 2) {
    
    ladprofile <- lidR::LAD(z, dz = dz, k = k, z0 = z0)
    
    lad_max <- with(ladprofile, max(lad, na.rm = TRUE))
    lad_mean <- with(ladprofile, mean(lad, na.rm = TRUE))
    lad_cv <- with(ladprofile, sd(lad, na.rm=TRUE)/mean(lad, na.rm = TRUE))
    lad_min <- with(ladprofile, min(lad, na.rm = TRUE))
    lai <- with(ladprofile, sum(lad, na.rm = TRUE))
    
  }
  
  lad_metrics <- list(lad_max = lad_max,
                      lad_mean = lad_mean,
                      lad_cv = lad_cv,
                      lad_min = lad_min,
                      lai = lai)
  
  return(lad_metrics)
}



#' Canopy volume classes
#' 
#' Canopy volume classes based on Lefsky et al 1999 (see references), modified. A voxel rerprenetation of a forest stand
#' is divided into four classes including: open gap space, closed gap space, euphotic zone, and oligophotic zone. 
#' This function is meant to be used within metrics_voxels.
#' 
#' @param x,y,z  X, Y, Z coordinate of the voxels
#' @param n Point count inside each voxel. Used to distinguish filled and empty voxels.
#' @return Percentage of voxels in each class
#' @references 
#' Lefsky, M. A., Cohen, W. B., Acker, S. A., Parker, G. G., Spies, T. A., & Harding, D. (1999). Lidar Remote Sensing of the Canopy Structure and Biophysical Properties of Douglas-Fir Western Hemlock Forests. Remote Sensing of Environment, 70(3), 339-361. doi:10.1016/S0034-4257(99)00052-8
#' 


metrics_lefsky <- function(x, y, z, n) {
  
  OpenGapSpace <- ClosedGapSpace <- Euphotic <- Oligophotic <- NA_real_
  
  if (length(z) > 2) {
    
    dvox <-  data.table::data.table(X=x, Y=y, Z=z, n=n)
    
    dvox_top <- dvox %>% 
      dplyr::filter(!is.na(n)) %>%
      dplyr::group_by(X, Y) %>%
      dplyr::summarise(Zmax = max(Z), .groups = "keep")
    
    
    vox_stats <- dvox %>% 
      dplyr::left_join(dvox_top, by = c("X", "Y")) %>%
      dplyr::mutate(class = 
                      dplyr::case_when(
                        Z > Zmax ~ "OpenGapSpace",                  #empty voxels above canopy
                        Z < Zmax & is.na(n) ~ "ClosedGapSpace",     #empty voxels below canopy
                        Z >= 0.65 * Zmax ~ "Euphotic",              #filled voxels in the top 0.65 portion of the canopy
                        TRUE ~ "Oligophotic"                        #filled voxels in the lower section of the canopy
                      )) %>% 
      dplyr::group_by(class) %>%
      dplyr::count() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(perc = n/sum(n) * 100) %>%
      dplyr::select(-n)
    
    
    #ensure all classes are reported, NA if empty:
    all_classes <- c("OpenGapSpace", "ClosedGapSpace", "Euphotic", "Oligophotic")
    
    vox_stats <- dplyr::tibble(class = all_classes) %>%
      dplyr::left_join(vox_stats, by = "class")
    
    out <- vox_stats %>% tidyr::pivot_wider(names_from = class, values_from = perc) %>% #, names_glue = "{.value}_{class}") %>%
      as.list()
    
  } else {
    out <- list(OpenGapSpace=OpenGapSpace, ClosedGapSpace=ClosedGapSpace, Euphotic=Euphotic, Oligophotic=Oligophotic)
  }
  
  return(out)
}



#' Calculate rumple index
#' 
#' A wrapper of the \code{lidR::rumple_index} function that allows to calculate rumple index without the need for CHM, and 
#' can be used directly in the e.g. \code{grid_metrics} function. The function combines the two required steps, i.e. creating a surface model, and calculating rumple index, into one.
#' Top surface is created using highest points within each pixel.
#' 
#' 
#' @param x,y,z  X, Y, Z coordinates of a point cloud
#' @param pixel_size pixel size
#' @param zmin Minimum height. If set, heights below are ignored in calculations.
#' @return Same as in \code{lidR::rumple_index} - the calculated rumple index
#' @export
#' @examples
#' library(lidR)
#' library(lidRmetrics)
#' LASfile <- system.file("extdata", "Megaplot.laz", package="lidR")
#' las <- readLAS(LASfile, select = "*", filter = "-keep_random_fraction 0.5")
#' 
#' m1 <- cloud_metrics(las, ~metrics_rumple(x = X, y = Y, z = Z, pixel_size = 1))
#' 
#' m2 <- grid_metrics(las, ~metrics_rumple(x = X, y = Y, z = Z, pixel_size = 1), res = 40)


metrics_rumple <- function(x, y, z, pixel_size, zmin=NA) {
  
  if (!is.na(zmin)) z <- z[z>zmin]
  
  r <- NA_real_
  
  if (length(z) > 2) {
    
    D <-  data.table::data.table(X=x, Y=y, Z=z)
    
    D <- LAS(D, header = rlas::header_create(D), check=F)
    
    D <- lidR::decimate_points(D, lidR::highest(pixel_size))
    
    r <- lidR::rumple_index(x = D$X, y = D$Y, z = D$Z)
  }
  
  return(list(rumple=r))
  
}
