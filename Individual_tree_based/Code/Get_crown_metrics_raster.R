

library(lidR)
library(raster)
library(sf)
library(dplyr)
library(exactextractr)
library(terra)

#' Wrapper function to get rasterized crown metrics
#' 
#' @param las The input LAS file (un-normalized or filtered).
#' @param res An integer to indicate raster pixel size (default is 5).
#' 
#' @return A list with the first item is crowns in geospatial polygons and 2nd item is the raster brick of canopy height model and crown metrics.
#' @export 
#' 
#' @examples
#' # Example usage:
#' # las <- readLAS("path/to/your/las/file.las")
#' # crown_metrics <- get_crown_metrics(las, res = 10)
get_crown_metrics = function(las, res=5) {
  
  las = normalize_las(las)
  
  ## derive canopy height model (CHM)
  chm0 =rasterize_canopy(las, res=res, algorithm = p2r())
  chm = raster(chm0)
  
  ## get crown from normalized las, min_z is used to define the minimum height of the tree
  crowns = get_individual_tree_crown(las, chm_res= 0.5,  min_z = 1.3) ## crowns is a spatial polygon dataframe
  crown_metrics = rasterize_crown_metrics(crowns, chm)
  
  crown_metrics = raster::brick(crown_metrics)
  chm_crown_metrics = stack(chm, crown_metrics)
  names(chm_crown_metrics)[1] = 'CHM'
  return(list(crowns = crowns, chm_crown_metrics = chm_crown_metrics))
}




### below are helper functions 

## used by get_crown_metrics to get normalized las and filter out potential noise
## ** you can edit this function to fit your way of normalizing/cleaning las file *****##

#' Helper function to get normalized las and filter out potential noise
#' @param las The input las file (un-normalized or filtered).
#' @return A normalized las without potential noise.
normalize_las = function(las){
  ## remove potential noise
  las <- filter_duplicates(las)
  las = las[!las$Classification%in%c(6, 7,9, 10, 11, 12,13,14,15,16,17,18),]
  #las = las[!las$Classification%in%c(6,18),]
  las =  classify_noise(las, ivf(3,2))
  las <- filter_poi(las, (Classification != LASNOISE))
  ## get normalized las and remove potential noise
  las <- normalize_height(las, tin(extrapolate = knnidw(10, 2, 5)))
  las =  filter_poi(las, Z >= 0 & Z < min(100, quantile(Z, .9999)))
  
  return(las)
}


## used by get_crown_metrics function to return tree crowns as geospatial polygons
## ** you can edit this function to add more individual tree based metrics *****##

#' Helper function to return tree crowns as geospatial polygons using normalized las
#' @param las_normalized normalized las file.
#' @param chm_res the spatial resolution of canopy height map for crown segmentation.
#' @param min_z minimize height for a tree to be included
#' @return all tree crowns in as geospatial polygons with individual tree based features
get_individual_tree_crown = function(las_normalized, chm_res= 0.5,  min_z = 1.3){
  
  # can run below function to make the point density to a given density, np, and this can potential if lidar tiles have very different density
  # las_normalized = decimate_points(las_normalized, random(10))
  
  #### below parameters are from Yahampath ###
  thr <- c(0,2,5,10,15)
  edg <- c(0, 1.5)
  
  # get a 0.5 meter resolution canopy height for crown delineation, you can change to 
  chm <- rasterize_canopy(las_normalized, chm_res, pitfree(subcircle = 0.2, thr, edg))
  
  # smooth with median filter for chm, below parameters are from Yahampath
  kernel <- matrix(1,3,3)
  custom_ws <- function(x) {
    y <- (1/8) * x 
    y[x < 32] <- 4
    y[x > 80] <- 10
    return(y)
  }
  
  # uncomment this to enable canopy smooth
  #chm_smoothed <- terra::focal(chm, w = kernel, fun = median, na.rm = TRUE)
  
  # comment this use smoothed chm
  chm_smoothed =chm
  
  # get tree tops, hmin: the minimum tree height
  ttops <- locate_trees(las = chm_smoothed, algorithm = lmf(ws = custom_ws, hmin=min_z)) 
  
  # get several individual tree based metric if there are trees
  ## ** feel free to add any potential metrics below ** ##
  if (nrow(ttops) > 0 ) {
    
    ## segment individual tree 
    tree_segments <- segment_trees(las = las_normalized, silva2016(chm_smoothed , ttops, max_cr_factor = 0.7))
    
    ## count number of tree
    num_of_trees <- length(unique(tree_segments$treeID) |> na.omit())
    
    ## get the standard crown metrics 
    crowns <- crown_metrics(tree_segments, func = .stdtreemetrics, concaveman = c(3, 0), geom = "concave")
    
    ## calculate crown area using sf pacakge, 
    # this area is very different from the crown_hull area calculated from the stdtreemetrics above
    crowns$area = st_area(crowns)
    
    ## get a rought estimate of crown diameter 
    crowns$D1 = 2*sqrt(crowns$area/pi) # assuming crown shape is circular
    ## try to calculate crown volume by multiple tree height and crown area
    crowns$V1 = crowns$Z*crowns$area
    ## get tree_height * crown diameter; Inspired by Jucker et al (2016), https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.13388
    crowns$H_CD1 = crowns$Z*crowns$D1
    
  }
  
  return(crowns)
  
}


## ** you can edit this function to rasterize additioinal crown metrics *****##
#' Used by get_crown_metrics to rasterize all crown metrics.
#' @param crowns the geospatial polygons for crowns with individual tree based metrics from get_individual_tree_crown funtion.
#' @param chm  the chm with desired output pixel size to be used as a raster template for rasterization.
#' @return rasterized crown metrics.
rasterize_crown_metrics = function(crowns, chm){
  
  # you may want to filter out trees with small crowns
  #crowns = crowns[crowns$area>3,]
  
  # use chm as template for rasterizing crown metrics
  template = chm*0 + 1
  
  ## for each crown, extract the coordinates that overlap with each crown, and the correponsding percentage of pixel overlapping by each crown
  ## crowns_coverageP is a list, each list element is for a tree; each list element is a dataframe including pixel coordinates, and percentage of each pixel covered by this tree
  crowns_coverageP <- exact_extract(template, crowns, include_xy = T, default_value=0)
  
  ## get percentage of each crown overlapping with each pixel
  crowns_coverageP =  Map(function(df) {df$coverage_fraction2 <- df$coverage_fraction/sum(df$coverage_fraction); return(df)}, crowns_coverageP)
  ## to remove some edge effect, for each pixel remove trees with less than 5% overlapping with this pixel
  crowns_coverageP =  Map(function(df) {df <- df[df$coverage_fraction2>=0.05,]; return(df)}, crowns_coverageP)
  
  ## get crown metric for each pixel by each tree
  ## here is the logic to get pixel-level crown metric:  
  ## for each pixel, multiple the crown metric (e.g. crown area) by the percetage of this crown overlapping with this pixel.
  ## example, if the crown area is 50, and 50% of this crown is overlapping with this pixel, the crown area of this tree for this pixel is 50
  
  crowns_coverageP_Area =  Map(function(df, factor) {df$value <- df$value * df$coverage_fraction2 * factor/sum(df$coverage_fraction2); return(df)}, crowns_coverageP, crowns$area)
  crowns_coverageP_V1 =  Map(function(df, factor) {df$value <- df$value * df$coverage_fraction2 * factor/sum(df$coverage_fraction2); return(df)}, crowns_coverageP, crowns$V1)
  crowns_coverageP_HCD1 =  Map(function(df, factor) {df$value <- df$value * df$coverage_fraction2 * factor/sum(df$coverage_fraction2); return(df)}, crowns_coverageP, crowns$H_CD1)
  
  
  ## now combine list of trees (each list element is a dataframe) into 1 dataframe; 
  crowns_coverageP_Area = do.call(rbind, crowns_coverageP_Area)
  crowns_coverageP_V1 = do.call(rbind, crowns_coverageP_V1)
  crowns_coverageP_HCD1 = do.call(rbind, crowns_coverageP_HCD1)
  
  
  ## now aggregate crown metrics by coordinates rather than its orginal format by  each crown, you can define your own aggregation function
  fun = mean
  crowns_coverageP_Area_agg1 = aggregate(value ~ x + y, FUN = fun, data= crowns_coverageP_Area, na.rm=T)
  crowns_coverageP_V1_agg1 = aggregate(value ~ x + y, FUN = fun, data= crowns_coverageP_V1, na.rm=T)
  crowns_coverageP_HCD1_agg1 = aggregate(value ~ x + y, FUN = fun, data= crowns_coverageP_HCD1, na.rm=T)
  
  fun = max
  crowns_coverageP_Area_agg2 = aggregate(value ~ x + y, FUN = fun, data= crowns_coverageP_Area, na.rm=T)
  crowns_coverageP_V1_agg2 = aggregate(value ~ x + y, FUN = fun, data= crowns_coverageP_V1, na.rm=T)
  crowns_coverageP_HCD1_agg2 = aggregate(value ~ x + y, FUN = fun, data= crowns_coverageP_HCD1, na.rm=T)
  
  
  fun = sum
  crowns_coverageP_Area_agg3 = aggregate(value ~ x + y, FUN = fun, data= crowns_coverageP_Area, na.rm=T)
  crowns_coverageP_V1_agg3 = aggregate(value ~ x + y, FUN = fun, data= crowns_coverageP_V1, na.rm=T)
  crowns_coverageP_HCD1_agg3 = aggregate(value ~ x + y, FUN = fun, data= crowns_coverageP_HCD1, na.rm=T)
  
  # count number of trees for each pixel
  tree_count = aggregate(value ~ x + y, FUN = length, data= crowns_coverageP_Area) 
  
  area_mean = rasterize_each_metric(crowns_coverageP_Area_agg1, chm )
  V1_mean = rasterize_each_metric(crowns_coverageP_V1_agg1, chm )
  HCD1_mean = rasterize_each_metric(crowns_coverageP_HCD1_agg1, chm )
  
  area_max = rasterize_each_metric(crowns_coverageP_Area_agg2,  chm)
  V1_max = rasterize_each_metric(crowns_coverageP_V1_agg2,  chm)
  HCD1_max = rasterize_each_metric(crowns_coverageP_HCD1_agg2, chm)
  
  area_sum = rasterize_each_metric(crowns_coverageP_Area_agg3,  chm)
  V1_sum = rasterize_each_metric(crowns_coverageP_V1_agg3, chm)
  HCD1_sum = rasterize_each_metric(crowns_coverageP_HCD1_agg3,  chm)
  
  tree_number = rasterize_each_metric(tree_count, chm)
  
  
  stacks = c(area_mean, V1_mean, HCD1_mean, area_max, V1_max, HCD1_max, area_sum, V1_sum, HCD1_sum, tree_number)
  names(stacks) = c('Crown_area_Avg', 'Tree_volume_Avg', 'Tree_height_diameter_Avg','Crown_area_Max', 'Tree_volume_Max', 'Tree_height_diameter_Max',
                    'Crown_area_Sum', 'Tree_volume_Sum', 'Tree_height_diameter_Sum','N_tree')
  
  return(stacks)
}

#'Used by rasterize_crown_metrics function to rasterize each crown metric.
#' @param data dataframe with coordinates and value for interested metric.
#' @param cords_columns list of strings for columns names for longitude and latitude.
#' @param field character to indicate the column name for the interested metric.
#' @param template_raster template raster for rasterization.
#' @return one raster for the interested metric
#
rasterize_each_metric = function(data,  chm, cords_columns = c("x", "y"), field='value'){
  vect_df <- vect(data, geom=c("x", "y"))
  metric_raster =rasterize(vect_df, rast(chm), field, na.rm=T, background=0)
}
