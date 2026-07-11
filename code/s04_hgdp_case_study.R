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


