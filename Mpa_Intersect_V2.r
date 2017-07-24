# Algorithm to perform analyses of contribution of geomorphic features
# for an Exclusive Economic Zone (EEZ), and a set of MPAs intersecting this
# EEZ.
#
# Authors:
#   Emmanuel Blondel (FAO),
#   Levi Westerveld (UNEP GRID-ARENDAL)
#

#=============================================================================
# PACKAGE DEPENDENCIES
#=============================================================================
cat("Starting...\n")
require(RFigisGeo)
require(parallel)
require(jsonlite)

#=============================================================================
# INPUTS DEFINITION
#=============================================================================

#type of area (EEZ or ECOREGION)
areaType <- "EEZ"
#areaType <- "ECOREGION"

#area identifier (EEZ mrgid or ECOREGION provinc)
areaId <- "5677"
#areaId <- "Eastern Caribbean"


#=============================================================================
# INTERNAL VARIABLES
#=============================================================================

#check sys locales
cat("Checking system locales \n")
sys <- Sys.getlocale()
cat(sprintf("Default System locale = '%s' \n", sys))
if(sys == "C"){
  syslocale <- "en_US.UTF-8"
  cat(sprintf("Setting 'LC_ALL' category to '%s' \n", syslocale))
  Sys.setlocale(category = "LC_ALL", locale = syslocale)
  sys <- Sys.getlocale()
  cat(sprintf("New System locale = '%s' \n", sys))
}

#options
options(stringsAsFactors = FALSE)

#debugging
debugMode <- TRUE

#to run in parallel
runParallel <- TRUE
totalCores <- 8
securityCoreNb <- 1
#totalCores <- detectCores()
#if(totalCores > 8) securityCoreNb <- 2
runCores <- totalCores - securityCoreNb
cat(sprintf("Running with %s cores on %s cores availables \n", runCores, totalCores))

#spatial data
selected_area <- NULL
selected_intersect <- NULL

#bbox used to fetch geomorphic features
target_bbox <- NULL

#geomorphic types
baselayers <- c("shelf", "slope", "abyss", "hadal")
featurelayers <- c("seamounts", "guyots", "canyons", "ridges",
                   "plateaus", "spreading_ridges", "rift_valleys",
                   "glacial_troughs", "shelf_valleys", "trenches",
                   "troughs", "terraces", "fans", "rises",
                   "bridges","escarpments")
gtypes <- c(baselayers, featurelayers)

#list of geomorphic types for which bbox should be ignored
gtypesIgnoringBbox <- c("abyss")

#use Eckert IV equal area projection for area calculation
areaCRS <- CRS("+proj=eck4 +lon_0=Central Meridian +x_0=False Easting +y_0=False Northing")


#=============================================================================
# BUSINESS FUNCTIONS
#=============================================================================

#method to fetch a set of geomorphic features
# @param baseUrl object of class "character" the OWS baseUrl
# @param typeName object of class "character" layer name
# @param bbox object of class "matrix" representing a spatial bbox
# @param filter object of class "character" giving a CQL filter string
# @param verbose object of class "logical" Default is TRUE to print logs
# @returns an object of class "SpatialPolygonsDataFrame"
fetchFeatures <- function(baseUrl, typeName, bbox = NULL, filter, verbose = TRUE){
  #base config
  config <- list(
    baseUrl = baseUrl,
    service = "WFS",
    version = "1.0.0",
    request = "GetFeature",
    typeName = typeName
  )
  
  #optional bbox parameter
  if(!missing(bbox) && !is.null(bbox)){
    config$bbox = paste(bbox, collapse=",")
  }
  
  #optional filter
  if(!missing(filter)){
    config$cql_filter = filter;
  }
  
  request <- paste0(config$baseUrl,paste(lapply(2:length(config), function(i){paste(names(config)[i],config[[i]],sep="=")}),collapse="&"))
  if(verbose) cat("Reading WFS ", request,"\n")
  return(readWFS(request, target.dir = getwd(), verbose = verbose))
  
}

#method to fetch a set of geomorphic features
# @param typeName the layer name
# @param bbox a spatial bbox
# @param verbose TRUE to print logs
# @returns an object of class "SpatialPolygonsDataFrame"
fetchGeomorphicFeatures <- function(typeName, bbox, verbose = TRUE){
  return(
    fetchFeatures(
      baseUrl = "http://geoserver-dev.d4science-ii.research-infrastructures.eu/geoserver/ows?",
      typeName = typeName, bbox = bbox,
      verbose = verbose
    )
  )
}

#method to fetch a complete dataset of geomorphic features
# @param gtypes object of class "character" the list of geomorphic types to retrieves
# @param gtypesIgnoringBbox object of class "character" the list of geomorphic types
#        for which bbox should be ignored for data fetching
# @param bbox object of class "matrix" giving a spatial bbox
# @param runParallel object of class "logical", if code has to be parallelized. Default is FALSE
# @param runCores object of class "integer". Number of cores to use. Default is 16
# @param verbose object of class "logical". Default is TRUE to print logs
# @returns an object of class "SpatialPolygonsDataFrame"
fetchGeomorphicFeaturesAll <- function(gtypes, gtypesIgnoringBbox, bbox,
                                       runParallel = FALSE, runCores = NULL,
                                       verbose = TRUE){
  
  #doFetch
  doFetch <- function(chunk_of_gtypes, gtypesIgnoringBbox, bbox, verbose){
    cat(sprintf("Fetching data for chunk [%s] \n", paste(chunk_of_gtypes,collapse=',')))
    cat("----------------------------------------------------------------------------\n");
    out <- lapply(
      chunk_of_gtypes,
      function(gtype){
        cat(sprintf("Fetching data for geomorphic type %s \n",gtype))
        
        #typeName
        typeName <- paste0("W_mpa:geo_fea_", gtype)
        
        #bbox
        layerBbox <- bbox
        if(gtype %in% gtypesIgnoringBbox){
          cat(sprintf("Ignoring bbox parameter for %s \n",gtype))
          layerBbox <- NULL
        }
        
        #fetching
        sp <- fetchGeomorphicFeatures(typeName, bbox, verbose = TRUE)
        if(!is.null(sp)) sp@data$gtype <- gtype
        return(sp)
      }
    )
    out <- out[!sapply(out,is.null)]
    out <- do.call("rbind", out)
    return(out)
  }
  
  #chunks of gtypes
  chunks <- list(gtypes)
  if(runParallel) chunks <- suppressWarnings(split(gtypes, 1:runCores))
  
  #run fetching
  sp <- switch(as.character(runParallel),
               "FALSE" = lapply(chunks,
                                doFetch, gtypesIgnoringBbox, bbox, verbose),
               "TRUE"  = mclapply(chunks,
                                  doFetch, gtypesIgnoringBbox, bbox, verbose,
                                  mc.preschedule = TRUE, mc.cores = runCores)
  )
  sp <- sp[!sapply(sp,is.null)]
  sp <- do.call("rbind", sp)
  return(sp)
}

#method to intersect selected areas with geomorphic features
# @param areas1 object of class "SpatialPolygonsDataFrame"
# @param areas2 object of class "SpatialPolygonsDataFrame"
# @param crs object of class "CRS" giving the CRS to use fo surface area computation
# @param runParallel object of class "logical", if code has to be parallelized. Default is FALSE
# @param runCores object of class "integer". Number of cores to use. Default is 16
# @param verbose object of class "logical". Default is TRUE to print logs
# @returns an object of class "SpatialPolygonsDataFrame"
intersectAll <- function(areas1, areas2, crs,
                         runParallel = FALSE, runCores = NULL,
                         verbose = TRUE){
  
  #doIntersect
  doIntersect <- function(chunk_of_areas, areas1, areas2, crs){
    areas1_to_intersect = areas1[chunk_of_areas,]
    out <- intersection(areas1_to_intersect, areas2, areaCRS = crs)
    return(out)
  }
  
  #chunks of areas
  if(runParallel){
    chunks <- suppressWarnings(split(row.names(areas1), 1:runCores))
    chunks <- chunks[sapply(chunks, function(x){length(x)>0})]
    nbChunks <- length(chunks)
    if(nbChunks < runCores){
      runCores <- nbChunks
    }
  }
  
  #run intersect
  sp <- switch(as.character(runParallel),
               "FALSE" = intersection(areas1, areas2, areaCRS = crs),
               "TRUE"  = mclapply(chunks,
                                  doIntersect, areas1, areas2, crs,
                                  mc.preschedule = TRUE, mc.cores = runCores)
  )
  if(class(sp) == "list"){
    sp <- sp[!sapply(sp,is.null)]
    sp <- do.call("rbind", sp)
  }
  
  return(sp)                  
}

#=============================================================================
# BUSINESS ALGORITHM
#=============================================================================

#Area selection
#--------------
cat(sprintf("Select %s '%s' \n", areaType, areaId))
system.time(
  selected_area <- switch(areaType,
                          "EEZ" = fetchFeatures(
                            baseUrl = "http://geoserver-dev.d4science-ii.research-infrastructures.eu/geoserver/ows?",
                            typeName = "W_mpa:eez",
                            filter = sprintf("mrgid_eez = '%s'", areaId, verbose = debugMode)
                          ),
                          "ECOREGION" = fetchFeatures(
                            baseUrl = "http://geoserver-dev.d4science-ii.research-infrastructures.eu/geoserver/ows?",
                            typeName = "W_mpa:meow_ppow",
                            filter = sprintf("ecoregion = '%s'", areaId, verbose = debugMode)
                          )
  )
)

if(class(selected_area) != "SpatialPolygonsDataFrame" || nrow(slot(selected_area,"data")) == 0){
  stop(sprintf("Error with selection of %s '%s'", areaType, areaId))
}

#compute area (use same area CRS for everything)
selected_area@data$surface = gArea(spTransform(selected_area, areaCRS), byid=TRUE)

#target bbox for the analysis
target_bbox <- slot(selected_area, "bbox")

#Intersect selection
#-------------------
cat(sprintf("Select MPAs intersecting  %s with identifier = '%s' \n", areaType, areaId))
system.time(
  selected_intersect <- switch(areaType,
                               "EEZ" = fetchFeatures(
                                 baseUrl = "http://geoserver-dev.d4science-ii.research-infrastructures.eu/geoserver/ows?",
                                 typeName = "W_mpa:intersect_mpa_eez_v1",
                                 filter = sprintf("mrgid_eez = '%s'", areaId, verbose = debugMode)
                               ),
                               "ECOREGION" = fetchFeatures(
                                 baseUrl = "http://geoserver-dev.d4science-ii.research-infrastructures.eu/geoserver/ows?",
                                 typeName = "W_mpa:intersect_mpa_ecoregions_v1",
                                 filter = sprintf("ecoregion = '%s'", areaId, verbose = debugMode)
                               )
  )
)

if(class(selected_intersect) != "SpatialPolygonsDataFrame" || nrow(slot(selected_intersect, "data")) == 0){
  stop(sprintf("Error with selection of %s - MPAs intersects. No intersects found!",areaType))
}

#compute area (use same area CRS for everything)
selected_intersect@data$surface = gArea(spTransform(selected_intersect, areaCRS), byid=TRUE)

#Merge areas (EEZ, and MPAs as homogenized geospatial areas)
source.area.names <- switch(areaType,
                            "EEZ" = c("gml_id", "mrgid_eez", "geoname", "surface"),
                            "ECOREGION" = c("gml_id", "ecoregion", "ecoregion", "surface")
)
target.area.names <- c("gml_id", "id", "name", "surface", "type")
slot(selected_area, "data") <- cbind(selected_area@data[,source.area.names], type = areaType)
colnames(selected_area@data) <- target.area.names
slot(selected_intersect, "data") <- cbind(selected_intersect@data[,c("gml_id", "wdpaid", "name", "surface")], type = "MPA")
colnames(selected_intersect@data) <- target.area.names

#merge
selected_areas <- rbind(selected_area, selected_intersect)

#Geomorphic features selection
#-----------------------------
#note: likely to be parallelized for optimization

#load shapefiles of global scale geomorphic features to be intersected with EEZ and selected MPA
cat("Fetching geomorphic features\n")
system.time(
  geomorphicFeatures <- fetchGeomorphicFeaturesAll(
    gtypes, gtypesIgnoringBbox, target_bbox,
    runParallel = runParallel, runCores = runCores,
    verbose = debugMode
  )
)
cat("Geomorphic features fetched\n")


#Compute intersections
#----------------------
#note: likely to be parallelized for optimization

#intersect between selected_mpa and geomorphic features, use areaCRS argument in RFigisGeo intersection to compute area percentage
if(length(geomorphicFeatures) > 0 && length(selected_areas) > 0){
  cat("Intersecting areas\n")
  system.time(
    intersects <- intersectAll(selected_areas, geomorphicFeatures, areaCRS,
                               runParallel = runParallel, runCores = runCores,
                               verbose = debugMode)
  )
}

cat("Intersection done\n")

#Compile report
#--------------
#report to give for each geomorphic feature the percentage of area covering
#the map on the area covering the entire EEZ
cat("Reporting\n")

df <- slot(intersects, "data")
agg <- aggregate(
  df$geo_area,
  by=list(gtype=df$gtype, id = df$id, name = df$name, type = df$type),
  FUN="sum"
)
colnames(agg) <- c("gtype", "id", "name", "type", "area")
agg1 <- agg
agg2 <- agg[agg$type == areaType,]
report <- merge(agg1, agg2, by = "gtype")
colnames(report) <- c("gtype", "id", "name", "type", "surface", "area_id", "area_name", "area_type", "surface_in_area")
report <- report[order(as.factor(report$type),as.factor(report$name)),]
report <- report[,c("gtype","id","name","type","surface")]
report <- reshape(report, v.names = "surface", idvar = colnames(report)[2:4], direction = "wide", timevar="gtype")
report <- merge(selected_areas@data, report, by = c("id","name","type"))
report <- report[,-which(colnames(report)=="gml_id")]

#replace NAs
invisible(sapply(colnames(report), function(x){if(class(report[,x])=="numeric"){ report[is.na(report[,x]),x] <<- 0 }}))

#add allMPAs
#20170426 Deactivate this part (now added before intersecting for managing MPA overlaps)
##allMPAs <- aggregate(report[,4:ncol(report)], by=list(report$type), FUN="sum")
##allMPAs <- cbind(data.frame(id = 0, name = "All MPAs", type = "MPA"), allMPAs[allMPAs[,1] == "MPA",2:ncol(allMPAs)])
##report <- rbind(report, allMPAs)

#add all geomorphic features as columns (even if no values)
for(geomorphic.type in gtypes){
  report.gtype <- paste("surface",geomorphic.type,sep=".")
  if(!(report.gtype %in% colnames(report))){
    report[report.gtype] <- 0
  }
}

#order geomorphic features by UI preferences
report.fields <- colnames(report)
report <- report[,c(report.fields[regexpr("surface.",report.fields) == -1], paste("surface",gtypes,sep="."))]

#rename surface. columns
columns <- colnames(report)
surf.columns <- columns[regexpr("surface.", columns) > 0]
colnames(report)[regexpr("surface.", colnames(report)) > 0] <- as.character(sapply(surf.columns, function(x){as.character(unlist(strsplit(x,"surface."))[2])}))


#all MPAS
system.time(
  allMPAs <- data.frame(
    id = 0,
    name = "All MPAs",
    type = "MPA",
    surface = gArea(spTransform(gUnaryUnion(selected_intersect), areaCRS)),
    do.call("cbind", lapply(colnames(report)[5:length(colnames(report))], function(x){
      sp <- intersects[intersects@data$gtype == x & intersects@data$type == "MPA",]
      out <- data.frame(var = 0)
      colnames(out) <- x
      if(length(sp)>0){
        sp <- gUnaryUnion(sp)
        out[1L,] <- raster::area(spTransform(sp, areaCRS))
      }
      return(out)
    }))
  )
)

#merge to report
report <- rbind(report, allMPAs)


#=============================================================================
# OUTPUT DEFINITION
#=============================================================================

cat("Producing output\n")
#output
#------
outputFormat <- "json"
outputFile<- paste0("output_report_", format(Sys.time(),"%Y%m%d%H%M%S"),".", outputFormat)
switch(outputFormat,
       "csv"   = {
         write.table(report, outputFile, row.names = FALSE, sep=";")
       },
       "json"  = {
         file.create(outputFile)       
         conn <- file(outputFile)
         json <- toJSON(report, pretty = TRUE)
         writeLines(as.character(json), con = conn, sep = "\n", useBytes = FALSE)
       }
)
