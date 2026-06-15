# ──────────────────────────────────────────────────────────────────
# geoGraph: comparing square vs hexagonal grid distances from Cairo
# across Eurasia, replicating the Di Santo et al. (2026) setup
# ──────────────────────────────────────────────────────────────────

library(geoGraph)
library(terra)
library(sf)

# ──────────────────────────────────────────────────────────────────
# 1. Define the study area (to match Figure 4B of Di Santo et al.)
# ──────────────────────────────────────────────────────────────────

geo.box <- c(xmin = -15, xmax = 170, ymin = 0, ymax = 75)

# ──────────────────────────────────────────────────────────────────
# 2. Build a SQUARE-grid gGraph 
# ──────────────────────────────────────────────────────────────────

g_sq <- makeSquareGrid(30000, lon.range = c(-15, 170),
                       lat.range = c(0, 75))

# Assign land/sea status via rnaturalearth coastline overlay
g_sq <- findLand(g_sq)

# Habitat: land = 1, sea = impassable
cost_rules <- data.frame(
  habitat = c("land", "sea"),
  cost    = c(1,      1e6),
  stringsAsFactors = FALSE
)

g_sq <- setCosts(g_sq, attr.name = "habitat",
                 cost.rules = cost_rules, method = "mean")
g_sq <- dropDeadNodes(g_sq)
#g_sq <- setDistCosts(g_sq)

color.rules = data.frame(
  habitat = c("land", "sea"),
  color   = c("lightgreen",  "lightblue"))


g_sq <- setColors(g_sq, color.rules)
plot(g_sq)


# ──────────────────────────────────────────────────────────────────
# 3. Build a HEXAGONAL gGraph at ~110 km spacing
# ──────────────────────────────────────────────────────────────────

g_hex <- makeHexGrid(geo.box = geo.box, spacing = 110)

g_hex <- findLand(g_hex)

g_hex <- setCosts(g_hex, attr.name = "habitat",
                  cost.rules = cost_rules, method = "mean")
g_hex <- dropDeadNodes(g_hex)
#g_hex <- setDistCosts(g_hex)

g_hex <- setColors(g_hex, color.rules)
plot(g_hex)
plotEdges(g_hex)


# ──────────────────────────────────────────────────────────────────
# 4. Define target locations spanning Eurasia
# ──────────────────────────────────────────────────────────────────
# Cairo as source just like in Di Santo et al. (2026)

targets <- data.frame(
  name = c("Cairo (source)",
           "Istanbul", "Athens", "Rome", "Madrid", "Paris",
           "London", "Berlin", "Stockholm", "Helsinki",
           "Murmansk",         # high latitude, ~69°N
           "Tehran", "Tbilisi", "Baku",
           "Tashkent", "Almaty",
           "Novosibirsk",      # mid-Siberia, ~55°N
           "Yakutsk",          # NE Siberia, ~62°N, very far east
           "Vladivostok",      # Pacific coast, ~43°N
           "Magadan",          # NE Russia, ~60°N, far east
           "Beijing",
           "Ulaanbaatar"),
  lon  = c(31.24,
           28.98, 23.73, 12.50, -3.70, 2.35,
           -0.13, 13.40, 18.07, 24.94,
           33.08,
           51.39, 44.79, 49.87,
           69.27, 76.89,
           82.93,
           129.74,
           131.89,
           150.81,
           116.41,
           106.92),
  lat  = c(30.04,
           41.01, 37.98, 41.90, 40.42, 48.86,
           51.51, 52.52, 59.33, 60.17,
           68.97,
           35.69, 41.72, 40.41,
           41.31, 43.24,
           55.01,
           62.04,
           43.12,
           59.56,
           39.90,
           47.89),
  stringsAsFactors = FALSE
)

# ──────────────────────────────────────────────────────────────────
# 5. Build a gData for each gGraph
# ──────────────────────────────────────────────────────────────────

gd_sq  <- new("gData", coords = targets[, c("lon", "lat")],
              data = targets["name"], gGraph.name = "g_sq")
gd_sq <- closestNode(gd_sq, attr.name = "habitat", attr.value = "land")


gd_hex <- new("gData", coords = targets[, c("lon", "lat")],
              data = targets["name"], gGraph.name = "g_hex")
gd_hex <- closestNode(gd_hex, attr.name = "habitat", attr.value = "land")

# Sanity check: all targets should map to land nodes
isConnected(g_sq,  gd_sq)
isConnected(g_hex, gd_hex)

# ──────────────────────────────────────────────────────────────────
# 6. Compute least-cost paths FROM Cairo to all other targets
# ──────────────────────────────────────────────────────────────────

cairo_idx <- 1   # first row in `targets`

paths_sq  <- dijkstraFrom(gd_sq,  getNodes(gd_sq)[cairo_idx])
paths_hex <- dijkstraFrom(gd_hex, getNodes(gd_hex)[cairo_idx])

dists_sq  <- gPath2dist(paths_sq)
dists_hex <- gPath2dist(paths_hex)

# ──────────────────────────────────────────────────────────────────
# 7. Tidy up into a comparison table
# ──────────────────────────────────────────────────────────────────

results <- data.frame(
  name      = targets$name,
  lon       = targets$lon,
  lat       = targets$lat,
  dist_sq   = dists_sq,
  dist_hex  = dists_hex,
  diff_km   = dists_hex - dists_sq,
  rel_diff  = (dists_hex - dists_sq) / dists_hex
)

# Sort by hexagonal distance for readability
results <- results[order(results$dist_hex), ]
print(results, row.names = FALSE)

# ──────────────────────────────────────────────────────────────────
# 8. Visualise: side-by-side maps of paths from Cairo
# ──────────────────────────────────────────────────────────────────

par(mfrow = c(2, 1), mar = c(2, 2, 3, 1))

plot(g_sq, reset = TRUE,
     main = "Square 1° grid — distances from Cairo")

plot(paths_sq, add = TRUE, lwd = 1.5, col = "darkblue")

plot(g_hex, reset = TRUE,
     main = "Hexagonal 110 km grid — distances from Cairo")

plot(paths_hex, add = TRUE, lwd = 1.5, col = "darkblue")

par(mfrow = c(1, 1))

# ──────────────────────────────────────────────────────────────────
# 9. Bias as a function of latitude
# ──────────────────────────────────────────────────────────────────
# The headline result for the paper: the discrepancy between square
# and hex distances grows with latitude

par(mar = c(4, 4, 3, 1))
plot(results$lat, results$rel_diff * 100,
     xlab = "Latitude (°N)",
     ylab = "(hex − square) / hex distance, %",
     main = "Square-grid distance bias as a function of latitude",
     pch  = 19, col = "darkred")
abline(h = 0, lty = 2, col = "grey50")
text(results$lat, results$rel_diff * 100,
     labels = results$name, pos = 4, cex = 0.6)

