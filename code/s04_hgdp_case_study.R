#### hgdp case study ####
library(geoGraph)
library(readr)
library(vegan)

########################################
# 1. Load FST matrix and HGDP data
########################################

Fst_pairwise <- read_csv("data/raw/original/hgdp/Fst_pairwise.csv")
fst_mat <- as.matrix(as.data.frame(Fst_pairwise)[, -1])
rownames(fst_mat) <- Fst_pairwise[[1]]

hgdpData <- getData(hgdp)

########################################
# 2. Standardise population names and identify shared populations
########################################

# HGDP-side rename
hgdpData$pop_canonical <- hgdpData$Population
hgdpData$pop_canonical[hgdpData$Population == "Han (pooled)"] <- "Han"

# FST-side rename
rownames(fst_mat)[rownames(fst_mat) == "NAN-Melanesian"] <- "Melanesian"
colnames(fst_mat)[colnames(fst_mat) == "NAN-Melanesian"] <- "Melanesian"

# Intersection: 49 populations (drops Piapoco, Bergamo from FST;
# Italian, Colombian, TundraNentsi from HGDP)
shared_pops <- intersect(hgdpData$pop_canonical, rownames(fst_mat))

########################################
# 3. Subset both datasets to shared populations
########################################

fst_shared <- fst_mat[shared_pops, shared_pops]

keep <- hgdpData$pop_canonical %in% shared_pops
hgdp.sub <- hgdp[keep]

myGraph  <- dropCosts(worldgraph.40k)
hgdp.sub <- setGraph(hgdp.sub, "myGraph")

########################################
# 4. get unweighted geographic as well as geodetic distances for the shared populations
########################################

# get geodesic distances
coords  <- getCoords(hgdp.sub)
gc_dist <- sp::spDists(coords, longlat = TRUE)
rownames(gc_dist) <- colnames(gc_dist) <- hgdp.sub@data$Population

# get geographic distances
m <- dijkstraBetween(hgdp.sub)
geog_dist <- gPath2dist(m)

########################################
# 5. compare the geographic and geodetic distances with fst
########################################

mantel_r_gc <- vegan::mantel(as.dist(fst_shared),
                      as.dist(gc_dist),
                      method = "pearson")

mantel_r_geog <- vegan::mantel(as.dist(fst_shared),
                         as.dist(geog_dist),
                         method = "pearson")

mantel_r_gc #0.749
mantel_r_geog #0.8671


########################################
# 6. limit geographic paths to land
########################################

elevGraph <- readRDS("data/intermediate/elevation_graph/elevGraph.rds")

# Search over elevation cost coefficient
run_elev_coeff <- function(cost.coeff) {
  diff.cost <- function(x1, x2, cost.coeff) {
    1 + sqrt(abs(x1 - x2)) * cost.coeff
  }
  
  g <- setCosts(
    elevGraph,
    node.values = getNodesAttr(elevGraph)$elevation,
    method      = "function",
    FUN         = diff.cost,
    cost.coeff  = cost.coeff
  )
  
  assign(".gg_current_graph", g, envir = globalenv())
  on.exit(rm(".gg_current_graph", envir = globalenv()), add = TRUE)
  
  gd <- hgdp.sub
  gd@gGraph.name <- ".gg_current_graph"
  
  tryCatch({
    m  <- dijkstraBetween(gd)
    d  <- as.matrix(gPath2dist(m))
    
    # Attach population labels — d rows are in the order of populations
    # in hgdp.sub
    pops_order <- unique(getData(gd)$Population)
    pops_order[pops_order == "Han (pooled)"] <- "Han"
    
    stopifnot(length(pops_order) == nrow(d))
    rownames(d) <- colnames(d) <- pops_order
    
    # Now subset to fst_shared's order
    d <- d[rownames(fst_shared), colnames(fst_shared)]
    
    cor(as.numeric(as.dist(fst_shared)),
        as.numeric(as.dist(d)),
        method = "pearson")
  }, error = function(e) {
    message("Error at cost.coeff = ", cost.coeff, ": ", conditionMessage(e))
    NA_real_
  })
}

run_elev_coeff(0.01)

## Coarse log-spaced grid over ~5 orders of magnitude

coeffs <- 10^seq(-5, 5, length.out = 25)   # 0.00001 to 1

profile <- data.frame(
  cost.coeff = coeffs,
  mantel_r   = vapply(coeffs, run_elev_coeff, numeric(1))
)

print(profile)

## Quick plot

library(ggplot2)
ggplot(profile, aes(cost.coeff, mantel_r)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = mantel_r_geog$statistic, linetype = "dashed", colour = "grey40") +
  annotate("text", x = min(coeffs), y = mantel_r_geog$statistic,
           label = "unweighted (land-only)", hjust = 0, vjust = -0.3, size = 3) +
  scale_x_log10() +
  labs(x = "cost coefficient (log scale)",
       y = "Mantel r vs pairwise FST",
       title = "1 + sqrt(abs(x1 - x2)) * cost.coeff") +
  theme_minimal()

# get the best mantel r and the corresponding cost coefficient
best_row <- profile[which.max(profile$mantel_r), ]

########################################
# 7. compare the three distance matrices with fst
########################################

mantel_r_gc$statistic #0.749033
mantel_r_geog$statistic #0.8671144
best_row$mantel_r #0.8701551

# so on a worldwide scale, adding terrain heterogeneity to the geographic distance matrix 
# improves the correlation with FST only slightly. 
# The unweighted geographic distance matrix already captures most of the variation in FST.

