###############################################################################
# Script: 05_elevation_graph.R
#' Purpose: Colour rawgraph.40k nodes by mean land elevation
#' Author: Dominik Jud
#'
#' Input:  rawgraph.40k
#' Output: elevGraph (world graph with mean elevation)
###############################################################################

library(geoGraph)
library(elevatr)
library(dplyr)

########################################
# 1. Elevation raster over the world
########################################
world_sfc <- st_as_sfc(st_bbox(
  c(xmin = -179.99, xmax = 179.99, ymin = -90, ymax = 90),
  crs = 4326
))
world <- st_sf(geometry = world_sfc)

elevation_data <- elevatr::get_elev_raster(
  locations = world, z = 1, clip = "locations"
)

names(elevation_data)[1] <- "elevation"

spatr <- terra::rast(elevation_data)

########################################
# 2. Assign mean elevation onto the rawgraph.40k grid and remove sea connections
########################################
landGraph <- worldgraph.40k

costs <- data.frame(
  habitat = c("sea", "land"),
  cost    = c(1000,   1)
)

#landGraph <- setCosts(landGraph, attr.name = "habitat", method = "mean", cost.rules = costs)

landGraph <- dropDeadEdges(landGraph, thres = 10)

plotEdges(landGraph)

elevGraph <- assignByRaster(
  landGraph,
  raster    = spatr,
  fun       = mean,
  attr.name = "elevation"
)

########################################
# 3. Rename the raster output, mask sea to 0
########################################

node_attr <- getNodesAttr(elevGraph)
node_attr$elevation      <- node_attr$raster_points
node_attr$raster_points  <- NULL

is_sea <- as.character(node_attr$habitat) == "sea"
node_attr$elevation[is_sea] <- 0

node_attr$elevation[is.na(node_attr$elevation)] <- 0

elevGraph@nodes.attr <- node_attr

########################################
# 4. Plot
########################################

elev  <- node_attr$elevation
n_col <- 256                                       
ramp  <- colorRampPalette(
  c("#f7fcb9", "#addd8e", "#78c679",
    "#41ab5d", "#238443", "#005a32")
)(n_col)

## normalise land elevations to [1, n_col] indices
land_max      <- max(elev, na.rm = TRUE)
elev_idx      <- pmax(1, ceiling(elev / land_max * n_col))
node_colours  <- ramp[elev_idx]

## sea nodes get their own colour, outside the ramp
node_colours[is_sea] <- "lightblue"

plot(elevGraph, col = node_colours, reset = TRUE, edges = FALSE)

########################################
# 6. Set costs based on elevation difference
########################################

## Cost of moving between two adjacent nodes increases with their elevation
## difference: exp(-|dz| * cost.coeff) gives a weight in (0, 1] that decays

diff.cost <- function(x1, x2, cost.coeff) {
  1 + abs(x1 - x2) * cost.coeff        
}

elevGraph <- setCosts(
  elevGraph,
  node.values = getNodesAttr(elevGraph)$elevation,
  method      = "function",
  FUN         = diff.cost,
  cost.coeff  = 0.01
)

plot(elevGraph, col = node_colours, edges = TRUE)

########################################
# 7. assign graph to the hgdp object for later use
########################################

hgdpElevGraph <- setGraph(hgdp, elevGraph)
isConnected(hgdpElevGraph)
plot(hgdpElevGraph, col.gGraph = node_colours, reset = TRUE)



addis <- cbind(38, 9)
ori <- closestNode(elevGraph, addis)

paths <- dijkstraFrom(hgdpElevGraph, ori)

addis <- as.vector(addis)
plot(hgdpElevGraph, reset = TRUE, col.gGraph = node_colours)
plot(paths)
points(addis[1], addis[2], pch = "x", cex = 2)
text(addis[1] + 35, addis[2], "Addis Ababa", cex = .8, font = 2)
points(hgdp, col.nodes = "black")
