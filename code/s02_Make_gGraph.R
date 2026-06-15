# ──────────────────────────────────────────────────────────────────
# 0. Packages
# ──────────────────────────────────────────────────────────────────
library(geoGraph)
library(terra)
library(sf)

# ──────────────────────────────────────────────────────────────────
# 1. Inspect and (if necessary) reorient the DEM
# ──────────────────────────────────────────────────────────────────

dem_path <- "data/raw/original/scotese_paleodems/Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc_v2/Map48_PALEOMAP_1deg_Middle_Triassic_245Ma.nc"

dem <- terra::rast(dem_path)
print(dem)
plot(dem)

plot(dem, main = "Middle Triassic (~245 Ma) — Scotese & Wright (2018)")

terra::crs(dem) <- "EPSG:4326"

# Sanity check the elevation range
summary(values(dem))
hist(values(dem), breaks = 50,
     main = "Elevation distribution",
     xlab = "Elevation (m)")


# ──────────────────────────────────────────────────────────────────
# 2. Build a global hexagonal gGraph at ~110 km spacing
# ──────────────────────────────────────────────────────────────────

geo.box <- c(xmin = -180, xmax = 180, ymin = -90, ymax = 90)

# Spacing in km. 110 km ≈ 1° at the equator, so comparable
# in node density to the underlying 1° rectangular DEM.
g <- makeHexGrid(geo.box = geo.box, spacing = 173)

g
# Expect ~50,000–55,000 nodes for a 110 km global hex grid

# Quick plot to confirm the grid is sensible
plot(g, edge = TRUE, main = "Global hexagonal grid (110 km spacing)")

# ──────────────────────────────────────────────────────────────────
# 3. Assign elevation from the DEM to each node
# ──────────────────────────────────────────────────────────────────

# Mean elevation per hex cell — used to determine land vs sea
g <- assignByRaster(g, dem,
                    layer.name = "z",
                    fun = "max")

# Check that we got values for most nodes
elev <- getNodesAttr(g, attr.name = "z")
summary(elev)
sum(is.na(elev))   # how many nodes have no DEM coverage (should be ~0)

# ──────────────────────────────────────────────────────────────────
# 5. Set costs: land = 1, sea = impassable
# ──────────────────────────────────────────────────────────────────

# make a new nodes attribute "habitat" where all nodes with ele > 0 are 'land' and all others are 'sea'
habitat <- ifelse(is.na(elev) | elev <= 0, "sea",
                  ifelse(elev > 1000, "rugged", "land"))


# Add it to the gGraph
g <- setNodesAttr(g, attr.name = "habitat", values = habitat)

# Verify
table(habitat) 

colour_rules <- data.frame(
  habitat = c("land", "sea", "rugged"),
  color   = c("darkgreen",  "lightblue", "darkred"),
  stringsAsFactors = FALSE
)

g <- setColors(g, colour_rules)

cost_rules <- data.frame(
  habitat = c("land", "sea", "rugged"),
  cost    = c(1,      1e6, 7)     # sea cost so high it's effectively a barrier
)

g <- setCosts(g, attr.name = "habitat",
              cost.rules = cost_rules,
              method = "mean")

plot(g, shape = NULL)


