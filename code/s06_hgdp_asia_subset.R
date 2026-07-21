##############################################################################
# Custom gGraph for Central/South Asia + East Asia HGDP populations
##############################################################################


# here we make a custom gGraph for or hgdp subset to show that for these populations adding mountains
# is worth doing 

library(geoGraph)
library(terra)
library(sf)
library(elevatr)
library(readr)

##################################
# creat a hexGraph for the area 
##################################

area.of.interest <- st_as_sfc(
  st_bbox(c(
    xmin = 55,
    xmax = 145,
    ymin = 5,
    ymax = 70
  ), crs = 4326)
)

aoi.bound <- st_sf(geometry = area.of.interest)

hexGraph <- makeHexGrid(geo.box = aoi.bound, spacing = 80)

hexGraph <- findLand(hexGraph)

colors <- data.frame(
  habitat = c("sea", "land"),
  color = c("#b0d4e8", "#90c090")
)

hexGraph <- setColors(hexGraph, colors)

plot(hexGraph, reset = TRUE)

# set the costs for sea cells to 100 and land cells to 1
cost.rules <- data.frame(
  habitat = c("sea", "land"),
  cost = c(100, 1)
)

hexGraph <- setCosts(hexGraph, method = "mean", attr.name = "habitat", cost.rules = cost.rules)
hexGraph <- dropDeadEdges(hexGraph, thres = 100)
plotEdges(hexGraph)

##################################
# add elevation sd for each cell
##################################

elevation.data <- get_elev_raster(
  locations = aoi.bound,
  z = 4,
  clip = "locations"
)

spatr <- rast(elevation.data)

elevGraph <- assignByRaster(graph = hexGraph, raster = spatr, layer.name = "elevation", fun = "sd")

node_attr <- getNodesAttr(elevGraph)

is_sea <- as.character(node_attr$habitat) == "sea"
node_attr$elevation[is_sea] <- 0

node_attr$elevation[is.na(node_attr$elevation)] <- 0

elevGraph@nodes.attr <- node_attr

########################################
# Plot
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


#########################################
# create a gData with the hgdp subset
#########################################

# fst
Fst_pairwise <- read_csv("data/raw/original/hgdp/Fst_pairwise.csv")
fst_wide <- as.data.frame(Fst_pairwise)
rownames(fst_wide) <- fst_wide[[1]]     
fst_wide <- fst_wide[, -1]
fst_mat  <- as.matrix(fst_wide)


#subset the hgdp data to only include East Asia and Central/South Asia populations
asia <- hgdp[hgdp@data$Region %in% c("EAST_ASIA", "CENTRAL_SOUTH_ASIA")]
asia <- asia[!asia@data$Population %in% c("Japanese", "TundraNentsi"), ]
asia@data$Population[asia@data$Population == "Han (pooled)"] <- "Han"

data <- asia@data
pops <- data$Population

fst_shared  <- fst_mat[pops, pops]  

lat = data$Latitude
lon = data$Longitude
coords <- cbind.data.frame(lon, lat)
colnames(coords) <- c("lon", "lat")

asiaGD <- new("gData", coords = coords, data = data, gGraph.name = "elevGraph")
plot(asiaGD, col.gGraph = node_colours)

#########################################
# search for the best cost surface 
#########################################

attr <- getNodesAttr(elevGraph)
sd_ele <- attr$elevation
base_habitat <- getNodesAttr(elevGraph)$habitat   # snapshot before any changes

run_threshold <- function(threshold) {
  # Reclassify cells based on this threshold
  habitat_new <- dplyr::case_when(
    base_habitat == "sea"         ~ "sea",
    is.na(sd_ele)                 ~ "land",
    sd_ele > threshold            ~ "rugged",
    TRUE                          ~ "land"
  )
  
  g_local <- setNodesAttr(elevGraph, attr.name = "habitat", values = habitat_new)
  
  cost.rules <- data.frame(
    habitat = c("sea", "land", "rugged"),
    cost    = c(100,   1,      100)
  )
  
  g_local <- setCosts(g_local, method = "mean",
                      attr.name = "habitat",
                      cost.rules = cost.rules)
  
  # Assign the modified graph so the gData can find it by name
  base::assign(".gg_current_graph", g_local, envir = globalenv())
  on.exit(rm(".gg_current_graph", envir = globalenv()), add = TRUE)
  
  gd <- asiaGD
  gd@gGraph.name <- ".gg_current_graph"
  
  tryCatch({
    m <- dijkstraBetween(gd)
    d <- as.matrix(gPath2dist(m))
    
    pops_order <- unique(getData(gd)$Population)
    pops_order[pops_order == "Han (pooled)"] <- "Han"
    rownames(d) <- colnames(d) <- pops_order
    
    d <- d[rownames(fst_shared), colnames(fst_shared)]
    
    cor(as.numeric(as.dist(fst_shared)),
        as.numeric(as.dist(d)),
        method = "pearson")
  }, error = function(e) {
    message("Error at threshold = ", threshold, ": ", conditionMessage(e))
    NA_real_
  })
}

run_threshold(max(sd_ele, na.rm = TRUE))


# Thresholds chosen from the SD distribution using quantiles
thresholds <- quantile(sd_ele,
                       rep(seq(0.7, 0.999, by = 0.001)),
                       na.rm = TRUE)


profile <- data.frame(
  quantile   = names(thresholds),
  threshold  = as.numeric(thresholds),
  mantel_r   = vapply(as.numeric(thresholds), run_threshold, numeric(1))
)

print(profile)

library(ggplot2)

ggplot(profile, aes(threshold, mantel_r)) +
  geom_line() +
  geom_point(size = 2) +
  geom_text(aes(label = quantile), vjust = -1, size = 3) +
  labs(x = "SD threshold (m)",
       y = "Mantel r vs pairwise FST",
       title = "Grid search: SD threshold defining 'rugged' cells (cost = 100)",
       subtitle = "Labels show quantile of SD distribution") +
  theme_minimal()

best_idx <- which.max(profile$mantel_r)
cat("Best threshold:", round(profile$threshold[best_idx], 1), "m",
    "(", profile$quantile[best_idx], "of SD distribution)\n",
    "Mantel r =", round(profile$mantel_r[best_idx], 4), "\n")

#########################################
# compare the elevGraph with best threshold to the original graph
#########################################

best_threshold <- round(profile$threshold[best_idx], 1)

# Reclassify cells at the winning threshold
habitat_final <- dplyr::case_when(
  base_habitat == "sea"       ~ "sea",
  is.na(sd_ele)               ~ "land",
  sd_ele > best_threshold     ~ "rugged",
  TRUE                        ~ "land"
)

# Build the final graph
optiGraph <- setNodesAttr(hexGraph,
                                 attr.name = "habitat",
                                 values    = habitat_final)

# Set colours: sea light blue, land pale green, rugged brown
colors_final <- data.frame(
  habitat = c("sea",       "land",     "rugged"),
  color   = c("#b0d4e8",   "#90c090",  "#8b6340")
)
optiGraph <- setColors(optiGraph, colors_final)

# Apply the same costs used in the grid search (for consistency)
cost.rules <- data.frame(
  habitat = c("sea", "land", "rugged"),
  cost    = c(100,   1,      100)
)
optiGraph <- setCosts(optiGraph,
                             method     = "mean",
                             attr.name  = "habitat",
                             cost.rules = cost.rules)

plot(optimal_hexGraph, reset = TRUE)

myGraph  <- dropCosts(optiGraph)
hgdp.sub.raw <- setGraph(asiaGD, "myGraph")
hgdp.sub.ele <- setGraph(asiaGD, "optiGraph")

# get geodesic distances
coords  <- getCoords(asiaGD)
gc_dist <- sp::spDists(coords, longlat = TRUE)
rownames(gc_dist) <- colnames(gc_dist) <- asiaGD@data$Population

# get raw geographic distances
m <- dijkstraBetween(hgdp.sub.raw)
raw_geog_dist <- gPath2dist(m)

# get geographic distances with elevation
m_elev <- dijkstraBetween(hgdp.sub.ele)
geog_dist <- gPath2dist(m_elev)

########################################
# 5. compare the geographic and geodetic distances with fst
########################################

mantel_r_gc <- vegan::mantel(as.dist(fst_shared),
                             as.dist(gc_dist),
                             method = "pearson")

mantel_r_geog <- vegan::mantel(as.dist(fst_shared),
                               as.dist(geog_dist),
                               method = "pearson")

mantel_r_raw <- vegan::mantel(as.dist(fst_shared),
                               as.dist(raw_geog_dist),
                               method = "pearson")

mantel_r_gc$statistic
mantel_r_raw$statistic
mantel_r_geog$statistic

