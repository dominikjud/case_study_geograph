##############################################################################
# Custom gGraph for Central/South Asia + East Asia HGDP populations
##############################################################################


# here we make a custom gGraph for or hgdp subset to show that for these populations adding mountains
# is worth doing 

library(geoGraph)
library(terra)
library(sf)
library(elevatr)

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

# here we use a hinge function (alternatively we could also set just a threshhold)

hinge.cost <- function(x1, x2, threshold, cost.coeff) {
  mean_sd <- (x1 + x2) / 2
  1 + pmax(0, mean_sd - threshold) * cost.coeff
}

run_hinge <- function(threshold, cost.coeff) {
  g_local <- setCosts(
    g,
    node.values = getNodesAttr(g)$elevation,
    method      = "function",
    FUN         = hinge.cost,
    threshold   = threshold,
    cost.coeff  = cost.coeff
  )
  
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
    message("Error at threshold=", round(threshold, 1),
            ", coeff=", signif(cost.coeff, 3),
            ": ", conditionMessage(e))
    NA_real_
  })
}


run_hinge(threshold = 0, cost.coeff = 0)


# Thresholds chosen from the actual SD distribution
thresholds <- quantile(sd_vec, c(0.50, 0.70, 0.80, 0.85, 0.90, 0.95, 0.98, 0.99),
                       na.rm = TRUE)

# Coefficients log-spaced
coeffs <- 10^seq(-4, 1, length.out = 12)

# Full grid
grid <- expand.grid(threshold = as.numeric(thresholds),
                    cost.coeff = coeffs)

grid$mantel_r <- mapply(run_hinge,
                        threshold  = grid$threshold,
                        cost.coeff = grid$cost.coeff)


head(grid[order(-grid$mantel_r), ], 10)

library(ggplot2)

ggplot(grid, aes(x = threshold, y = cost.coeff, fill = mantel_r)) +
  geom_tile() +
  scale_y_log10() +
  scale_fill_viridis_c(option = "magma", name = "Mantel r") +
  labs(x = "SD threshold (m)",
       y = "Cost coefficient (log scale)",
       title = "Hinge cost function: (threshold, coefficient) landscape") +
  theme_minimal()


best_idx  <- which.max(grid$mantel_r)
best_pair <- grid[best_idx, ]

cat("Best grid point:\n",
    "  threshold  =", round(best_pair$threshold, 1), "m\n",
    "  cost.coeff =", signif(best_pair$cost.coeff, 3), "\n",
    "  Mantel r   =", round(best_pair$mantel_r, 4), "\n")

# Local refinement, parameterising in log space for the coefficient
nm <- optim(
  par     = c(threshold = best_pair$threshold,
              log_coeff = log(best_pair$cost.coeff)),
  fn      = function(p) -run_hinge(threshold = p[1], cost.coeff = exp(p[2])),
  method  = "Nelder-Mead",
  control = list(maxit = 60, reltol = 1e-4)
)

cat("\nRefined optimum:\n",
    "  threshold  =", round(nm$par[1], 1), "m\n",
    "  cost.coeff =", signif(exp(nm$par[2]), 3), "\n",
    "  Mantel r   =", round(-nm$value, 4), "\n")




