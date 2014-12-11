plot.lagr = function(obj, target, type=c("raw", "coef", "is.zero"), id="id") {
    
    
    #Follow here for a one-dimensional effect-modifying parameter:
    if (obj$dim == 1) {
        #Put the target in a plotting slot:
        if (type == "coef") {
            val = obj$coefs[[target]]
        } else if (type == "is.zero") {
            val = obj$is.zero[[target]]
        } else if (type == "raw") {
            val = obj$data[[target]]
        }
        
        plot(x=obj$coords, y=val, xlab="location", ylab=target, type='l', bty='n')
    } 
    
    #Follow here for a two-dimensional effect-modifying parameter:
    if (obj$dim == 2) {
    
        
        if (is(obj$data, "Spatial")) {
            polygons = obj$data
        } else {
            #If the data was not specified as a spatial data frame, make a Voronoi diagram:
            crds = obj$coords
            z = deldir(crds[,1], crds[,2])
            w = tile.list(z)
            polys = vector(mode='list', length=length(w))
            
            for (i in seq(along=polys)) {
                pcrds = cbind(w[[i]]$x, w[[i]]$y)
                pcrds = rbind(pcrds, pcrds[1,])
                polys[[i]] = Polygons(list(Polygon(pcrds)), ID=as.character(i))
            }
            SP = SpatialPolygons(polys)
            polygons = SpatialPolygonsDataFrame(SP, data=data.frame(x=crds[,1], 
                y=crds[,2], row.names=sapply(slot(SP, 'polygons'), 
                function(x) slot(x, 'ID'))))
            
            polygons@data$id = rownames(polygons@data)
        } 
        
        #Put the target in a plotting slot:
        if (type == "coef") {
            polygons@data[[paste(".", target, sep="")]] = obj$coefs[[target]]
        } else if (type == "is.zero") {
            polygons@data[[paste(".", target, sep="")]] = obj$is.zero[[target]]
        } else if (type == "raw") {
            polygons@data[[paste(".", target, sep="")]] = obj$data[[target]]
        }
        
        #Repair any holes in the shapefile:
        slot(polygons, "polygons") <- lapply(slot(polygons, "polygons"), checkPolygonsHoles)
        polygons2 <- unionSpatialPolygons(polygons, as.character(polygons[[id]]))
        
        points = fortify(polygons2, region=id)
        names(points)[which(names(points)=='id')] = id
        df = join(points, polygons@data, by=id)      
        
        #Plot the shapes:
        ggplot(df, aes_string("long", "lat", group=id, fill=paste(".", target, sep=""))) + geom_polygon() +
            scale_fill_gradient2(low=muted("blue"), mid="white", high="orange", limits=range(polygons@data[[paste(".", target, sep="")]], na.rm=TRUE), name="")

    }


}
