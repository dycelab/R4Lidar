
## set to your working directory with las data, there is an example laz file in Data folder
setwd("...")

las = readLAS('las4test.laz')
las

## run the function, make sure you loaded all the functions in the R code (in Code folder)
outs= get_crown_metrics(las, res=5) 

crowns= outs[[1]]
crown_metrics = outs[[2]]

## the code below will give you the name of each metric
names(crown_metrics)

par(mfrow=c(1,3))
plot(crown_metrics[[1]])
plot(st_geometry(crowns), add=T)
plot(crown_metrics[[2]])
plot(crown_metrics[[3]])
