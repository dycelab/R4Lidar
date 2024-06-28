

setwd("/uufs/chpc.utah.edu/common/home/dycelab/data/Lidar/Canada/British_columbia/")

las = readLAS('bc_082f028_3_3_4_xyes_6_utm11_2017.laz')
las
outs= get_crown_metrics(las, res=10)
crowns= outs[[1]]
crown_metrics = outs[[2]]
names(crown_metrics)

par(mfrow=c(1,3))
plot(crown_metrics[[1]])
plot(st_geometry(crowns), add=T)
plot(crown_metrics[[2]])
plot(crown_metrics[[3]])
