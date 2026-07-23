##############################################################################
# Compare topoDistance to geoGraph on the Cairo-to-Eurasia route set
# across Eurasia, replicating the Di Santo et al. (2026) setup
##############################################################################

library(topoDistance)
library(elevatr)
library(sf)
library(terra)
library(geoGraph)

# ─── 1. Load the target points ───────────────────────────────────────

targets <- data.frame(
  name = c("Cairo (source)",
           "Istanbul", "Athens", "Rome", "Madrid", "Paris",
           "London", "Berlin", "Stockholm", "Helsinki",
           "Murmansk", "Tehran", "Tbilisi", "Baku",
           "Tashkent", "Almaty", "Novosibirsk", "Yakutsk",
           "Vladivostok", "Magadan", "Beijing", "Ulaanbaatar"),
  lon  = c(31.24, 28.98, 23.73, 12.50, -3.70, 2.35,
           -0.13, 13.40, 18.07, 24.94, 33.08, 51.39, 44.79, 49.87,
           69.27, 76.89, 82.93, 129.74, 131.89, 150.81, 116.41, 106.92),
  lat  = c(30.04, 41.01, 37.98, 41.90, 40.42, 48.86,
           51.51, 52.52, 59.33, 60.17, 68.97, 35.69, 41.72, 40.41,
           41.31, 43.24, 55.01, 62.04, 43.12, 59.56, 39.90, 47.89)
)

# ─── 2. Get a DEM covering the full Cairo-to-Magadan extent ──────────

# Compute a bounding box with margin
lon_range <- range(targets$lon) + c(-2, 2)
lat_range <- range(targets$lat) + c(-2, 2)

# Build sf polygon for elevatr
extent_sf <- sf::st_as_sf(
  data.frame(id = 1L),
  geometry = sf::st_as_sfc(
    sf::st_bbox(c(xmin = lon_range[1], xmax = lon_range[2],
                  ymin = lat_range[1], ymax = lat_range[2]),
                crs = 4326)
  )
)

# z = 4 gives ~6 km cells, which is a reasonable compromise between
# accuracy and computational cost for Eurasian-scale paths.
# For direct comparison to the Di Santo et al. analysis, z = 5 (~3 km)
# is more common — but at this extent, z = 5 gives ~100 million cells,
# which topoDist will handle very slowly. Start with z = 4.

dem <- elevatr::get_elev_raster(locations = extent_sf,
                                z         = 4,
                                clip      = "locations")

# ─── 3. Mask oceans ──────────────────────────────────────────────────

# topoDist doesn't distinguish land from sea. We need to set sea pixels
# to NA so the algorithm can't route through them.

dem_rast <- terra::rast(dem)
dem_land <- dem_rast
dem_land[dem_land <= 0] <- NA   # oceans (and below-sea-level land like
# the Dead Sea) become impassable
# This is the standard convention.

# Convert back to raster (topoDist expects a raster::RasterLayer, not
# terra::SpatRaster, in older versions)
dem_land_raster <- raster::raster(dem_land)

# ─── 4. Prepare points as a two-column matrix ────────────────────────

# topoDist wants a matrix of coordinates
points_matrix <- as.matrix(targets[, c("lon", "lat")])
rownames(points_matrix) <- targets$name

# ─── 5. Compute pairwise topographic distances ───────────────────────
# This is slow — expect 5-30 minutes for 22 points on a Eurasian DEM.
# The internal algorithm builds a lattice graph from every DEM pixel
# and runs Dijkstra from each source. Time scales roughly with
# (n_points × n_pixels).

topo_dists <- topoDistance::topoDist(
  DEM      = dem_land_raster,
  pts      = points_matrix,
  paths    = TRUE   # save memory; don't return path geometries
)

# Result is a matrix of distances in the DEM's linear units (metres for
# this projection). Rownames/colnames should match target names.

cat("topoDistance computation complete.\n")
print(round(topo_dists / 1000, 0))   # in km, rounded to whole numbers

# ─── 6. Extract distances from Cairo ─────────────────────────────────

# Cairo is the first row
cairo_topo <- topo_dists[1, ]

result <- data.frame(
  name          = targets$name,
  lat           = targets$lat,
  lon           = targets$lon,
  topo_dist_km  = cairo_topo / 1000
)

print(result)


##############################################################################
# 8. Compute geoGraph distances from Cairo for the same targets
##############################################################################

# ─── Build the same-name targets, with Cairo first ───────────────────

# Load the geoGraph graph you're using (worldgraph.40k, land-only)
g_geog <- dropCosts(worldgraph.40k)   # land-only, unweighted

# Build a gData object of the target points
target_gd <- new("gData", coords = targets[, c("lon", "lat")],
                 data = targets["name"], gGraph.name = "rawgraph.40k")

target_gd <- closestNode(target_gd, attr.name = "habitat", attr.value = "land")

# Verify all points connected (some coastal cities might need nudging)
isConnected(target_gd)

# ─── Compute least-cost paths from Cairo to all other targets ────────

cairo_node <- unlist(target_gd@nodes.id[1])   # first row of targets = Cairo

paths <- dijkstraFrom(target_gd, start = cairo_node)
dists_geograph <- gPath2dist(paths)

# Attach to the targets data frame
targets$geog_dist_raw <- dists_geograph

##############################################################################
# 9. Combine into a comparison table
##############################################################################

comparison <- data.frame(
  name         = targets$name,
  lat          = targets$lat,
  lon          = targets$lon,
  topo_m       = cairo_topo,                    # from your topoDistance run
  topo_km      = cairo_topo / 1000,
  geog_raw     = targets$geog_dist_raw          # geoGraph units, whatever they are
)

# Rank each distance from smallest to largest
comparison$rank_topo <- rank(comparison$topo_km)
comparison$rank_geog <- rank(comparison$geog_raw)

# Normalise both by their maximum, so each is on [0,1]
comparison$topo_norm <- comparison$topo_km / max(comparison$topo_km, na.rm = TRUE)
comparison$geog_norm <- comparison$geog_raw / max(comparison$geog_raw, na.rm = TRUE)

# Sort by topoDistance for readability
comparison <- comparison[order(comparison$topo_km), ]

print(comparison[, c("name", "lat", "topo_km", "geog_raw",
                     "topo_norm", "geog_norm", "rank_topo", "rank_geog")])

##############################################################################
# 10. Quantify agreement
##############################################################################

# Pearson correlation on raw values (scale-invariant)
r_pearson <- cor(comparison$topo_km, comparison$geog_raw, method = "pearson")

# Spearman correlation on ranks (checks ordering agreement)
r_spearman <- cor(comparison$topo_km, comparison$geog_raw, method = "spearman")

# Difference in normalised distances (target-by-target)
comparison$normalized_diff <- comparison$geog_norm - comparison$topo_norm

cat("Pearson r (topo vs geog):  ", round(r_pearson, 4), "\n")
cat("Spearman r (topo vs geog): ", round(r_spearman, 4), "\n")
cat("\nMean absolute normalised diff:",
    round(mean(abs(comparison$normalized_diff)), 4), "\n")
cat("Max absolute normalised diff:",
    round(max(abs(comparison$normalized_diff)), 4), "\n")

##############################################################################
# 11. Plot 1 — scatter of the two distances
##############################################################################

library(ggplot2)
library(ggrepel)

ggplot(comparison, aes(topo_km, geog_raw)) +
  geom_point(size = 2, colour = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, colour = "tomato", linewidth = 0.5) +
  ggrepel::geom_text_repel(aes(label = name), size = 3, max.overlaps = 20) +
  labs(x = "Topographic distance from Cairo (km, topoDistance)",
       y = "Least-cost distance from Cairo (raw units, geoGraph)",
       title = "geoGraph vs topoDistance: distances from Cairo across Eurasia",
       subtitle = sprintf("Pearson r = %.3f, Spearman r = %.3f",
                          r_pearson, r_spearman)) +
  theme_minimal()

##############################################################################
# 12. Plot 2 — normalised distances side by side, coloured by latitude
##############################################################################

# For each target, compare geog_norm vs topo_norm; residuals show
# the latitudinal bias

ggplot(comparison, aes(x = lat, y = normalized_diff)) +
  geom_point(size = 2, colour = "darkred") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_smooth(method = "lm", se = FALSE, colour = "steelblue", linewidth = 0.5) +
  ggrepel::geom_text_repel(aes(label = name), size = 3, max.overlaps = 20) +
  labs(x = "Target latitude (°N)",
       y = "Normalised distance difference (geoGraph − topoDist)",
       title = "Latitudinal bias in geoGraph vs topoDistance",
       subtitle = "Positive values: geoGraph gives longer relative distance") +
  theme_minimal()

##############################################################################
# 13. Statistical test: is there a latitude-dependent bias?
##############################################################################

# Regress the normalised difference on target latitude
lm_bias <- lm(normalized_diff ~ lat, data = comparison)
summary(lm_bias)

# If the slope is significantly non-zero, there's a systematic
# latitude-dependent bias.
cat("\nSlope of bias vs latitude:", 
    round(coef(lm_bias)["lat"], 5), "\n")
cat("p-value:",
    round(summary(lm_bias)$coefficients["lat", "Pr(>|t|)"], 4), "\n")