

#########################################################################
# This code cleans the point cloud and generates Individual Tree Crown  #
# segments using parallel                                                #
# input: a point cloud.                                                 #
# outputs: cleaned point cloud for outliers: normalized to topography:  #
#         canopy height model, tree tops, tree segmentation, save to    #
#         shape files                                                   #
#########################################################################

# clean up memory in global environment
rm(list = ls(globalenv()))

# import libraries
library(lidR)
library(sf)
library(terra)
library(raster)
library(foreach)
library(doParallel)

# keep commented for bulk operations
#library(rgl) 

# 3D window open: keep comment for HPC
#rgl.open()

# function for local maxima filter variable custom window size based on tree height
custom_ws <- function(x) {
  y <- (1/8) * x 
  y[x < 32] <- 4
  y[x > 80] <- 10
  return(y)
}

# setup project directory
prjdir <- 'your work dir with las or laz files'
outdir <- 'work dir for the output'



# listing data dir
laslist <- list.files(prjdir,'.las$|.laz$')
print(laslist)


# define function to process lidar data for tree segmentation and crown extraction
tree_seg_crow = function(item){
  # create file path
  #print(item)
  file_path <- paste(prjdir, item, sep = '/')
  # open Lidar data
  las <- readLAS(files = file_path, filter = "-drop_z_below 0")
  # mean and SD for threshold based filtering
  meanZ <- mean(las$Z, na.rm= TRUE)
  sdZ <- sd(las$Z, na.rm = TRUE)
  # thresholds filtering for bad Z (points)
  threshold1 <- meanZ + 5 * sdZ
  threshold2 <- meanZ - 5 * sdZ
  las_filter1 <- las[las$Z <= threshold1]
  las_filter2 <- las_filter1[las_filter1$Z >= threshold2]
  # generate DTM
  dtm_idw <- rasterize_terrain(las_filter2, res = 1, algorithm = knnidw())
  # height normalization(lidr documentation)
  las_normalized <- normalize_height(las_filter2, dtm_idw, na.rm = TRUE)
  # drop z more than 120
  las_normalized <- las_normalized[las_normalized$Z <= 120]
  # CHM pit free Khosravipour et al.2014
  thr <- c(0,2,5,10,15)
  edg <- c(0, 1.5)
  chm <- rasterize_canopy(las_normalized, 1, pitfree(subcircle = 0.2, thr, edg))
  # smooth with median filter for chm
  kernel <- matrix(1,3,3)
  chm_smoothed <- terra::focal(chm, w = kernel, fun = median, na.rm = TRUE)
  # tree tops
  ttops <- locate_trees(las = chm_smoothed, algorithm = lmf(ws = custom_ws, hmin=5)) # lms 2.5
  
  # condition to filter no tree polygons(there are some files no trees higher than 5m)
  if (nrow(ttops) > 0 ) {
    # crown delineation / segmentation (works well from seed=0.5, cr =0.6)
    tree_segments <- segment_trees(las = las_normalized,
                                   algorithm = dalponte2016(chm = chm_smoothed,
                                                            treetops = ttops,
                                                            th_tree =  5,
                                                            th_seed = 0.5,
                                                            th_cr = 0.6,
                                                            max_cr = 10))
    num_of_trees <- length(unique(tree_segments$treeID) |> na.omit())
    # check for partial trees
    if (num_of_trees > 0){
      # draw crowns
      crowns <- crown_metrics(tree_segments, func = .stdtreemetrics, geom = "concave")
      #plot(crowns["convhull_area"], main = "Crown area") # uncomment to plot
      
      # convert to shp files
      #cvx_hulls_sf <- st_as_sf(crowns["convhull_area"], coords = c("X", "Y"), crs = st_crs(las))
      cvx_hulls_sf <- st_as_sf(crowns, coords = c("X", "Y"), crs = st_crs(las))
      # create an empty vector to store extracted values
      max_values <- numeric(length = nrow(cvx_hulls_sf))
      
      # Loop through each canopy polygon
      for (i in 1:nrow(cvx_hulls_sf)) {
        # Extract raster values within the current polygon
        values_within_polygon <- extract(chm_smoothed, cvx_hulls_sf[i, , drop = FALSE])
        num_vect <- values_within_polygon$focal_median
        num_vect <- as.numeric(num_vect)
        # Find maximum value for the current polygon
        max_values[i] <- max(num_vect, na.rm = TRUE)
      }
      # rename the column
      cvx_hulls_sf$TreeZ <- max_values
      
      
      # modified name for saving files
      mod_name <- gsub(".laz", "", item)
      out_shp_path <- file.path(outdir, paste0(mod_name, ".shp"))
      # save as shape files
      st_write(cvx_hulls_sf, out_shp_path, append=FALSE )
    }else{
      print("No full trees to deliniate **********************************")
    }
  }else{
    print("No tress h > 5m       *****************************************")
  }
  

  
}

## parallel the code over all the las/laz tiles
registerDoParallel(cores = 6)

results <- foreach(i = laslist,  .packages = c('lidR', 'raster', 'sf','terra') , .errorhandling='pass') %dopar% {
  tree_seg_crow(i)
}

stopImplicitCluster()
