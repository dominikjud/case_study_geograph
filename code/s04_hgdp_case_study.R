#### hgdp case study ####

library(geoGraph)
library(readr)

########################################
# # 1. Load the HGDP data and Fst pairwise matrix
########################################

# fst
Fst_pairwise <- read_csv("data/raw/original/hgdp/Fst_pairwise.csv")
fst_wide <- as.data.frame(Fst_pairwise)
rownames(fst_wide) <- fst_wide[[1]]     
fst_wide <- fst_wide[, -1]
fst_mat  <- as.matrix(fst_wide)

# hgdp gData
hgdpData <- getData(hgdp)

########################################
# 2. check if the populations in the Fst_pairwise data are present in the hgdpData
########################################

populations_in_hgdp <- unique(hgdpData$Population)
populations_in_fst <- rownames(fst_mat)
setdiff <- setdiff(populations_in_fst, populations_in_hgdp)

# check how many overlapping populations there are
overlap_hgdp <- intersect(populations_in_fst, populations_in_hgdp)

########################################
# 3. Standardise population names and subset to shared populations
########################################

hgdpData$pop_canonical <- hgdpData$Population
hgdpData$pop_canonical[hgdpData$Population == "Han (pooled)"] <- "Han"

# For the FST matrix, we rename via the rownames/colnames
rownames(fst_mat)[rownames(fst_mat) == "NAN-Melanesian"] <- "Melanesian"
colnames(fst_mat)[colnames(fst_mat) == "NAN-Melanesian"] <- "Melanesian"

populations_in_fst_canon <- rownames(fst_mat)

shared_pops <- intersect(hgdpData$pop_canonical, populations_in_fst_canon)
length(shared_pops)   # should be 49

hgdp_shared <- hgdpData[hgdpData$pop_canonical %in% shared_pops, ]

fst_shared  <- fst_mat[shared_pops, shared_pops]  

pops <- getData(hgdp)$Population
pops_canonical <- pops
pops_canonical[pops == "Han (pooled)"] <- "Han"
keep <- pops_canonical %in% shared_pops
hgdp.sub <- hgdp[keep]

myGraph <- dropCosts(elevGraph)
hgdp.sub <- setGraph(hgdp.sub, "myGraph")

########################################
# 4. get geographic as well as geodetic distances for the shared populations
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

mantel_r_gc <- cor(as.numeric(as.dist(fst_shared)),
                   as.numeric(as.dist(gc_dist)),
                   method = "pearson")

mantel_r_geog <- cor(as.numeric(as.dist(fst_shared)),
                      as.numeric(as.dist(geog_dist)),
                      method = "pearson")

########################################
# 6. add the elevation data to the hgdp shared populations
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
  geom_hline(yintercept = mantel_r_geog, linetype = "dashed", colour = "grey40") +
  annotate("text", x = min(coeffs), y = mantel_r_geog,
           label = "unweighted (land-only)", hjust = 0, vjust = -0.3, size = 3) +
  scale_x_log10() +
  labs(x = "cost coefficient (log scale)",
       y = "Mantel r vs pairwise FST",
       title = "1 + abs(x1 - x2) * cost.coeff ") +
  theme_minimal()

## Refine around the peak with a finer grid or optim()

best_idx <- which.max(profile$mantel_r)
best_coarse <- profile$cost.coeff[best_idx]

cat("\nBest coarse coefficient:", best_coarse,
    "with Mantel r =", round(profile$mantel_r[best_idx], 4), "\n")

## Local refinement using optim() over a narrow log-window

nm <- optim(
  par     = log(best_coarse),
  fn      = function(log_c) -run_elev_coeff(exp(log_c)),
  method  = "Brent",              # 1D bounded search
  lower   = log(best_coarse) - 1, # one order of magnitude below
  upper   = log(best_coarse) + 1  # one order of magnitude above
)

best_coeff <- exp(nm$par)
best_r     <- -nm$value

cat("\nRefined optimum:\n",
    "  cost.coeff =", best_coeff, "\n",
    "  Mantel r   =", round(best_r, 4), "\n")



diff.cost <- function(x1, x2, cost.coeff) {
  1 + sqrt(x1 - x2) * cost.coeff        
}

elevGraph <- setCosts(
  elevGraph,
  node.values = getNodesAttr(elevGraph)$elevation,
  method      = "function",
  FUN         = diff.cost,
  cost.coeff  = 0.04504557
)

hgdpElevGraph <- setGraph(hgdp.sub, elevGraph)
isConnected(hgdpElevGraph)

m_elev <- dijkstraBetween(hgdpElevGraph)

elev_dist <- gPath2dist(m_elev)

mantel_r_elev <- cor(as.numeric(as.dist(fst_shared)),
                    as.numeric(as.dist(elev_dist)),
                    method = "pearson")

########################################
# 7. compare the three distance matrices with fst
########################################

mantel_r_gc #0.749033
mantel_r_geog #0.8671144
mantel_r_elev #0.8705
