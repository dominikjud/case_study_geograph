# ──────────────────────────────────────────────────────────────────
# Download Scotese & Wright (2018) PaleoDEMs
# ──────────────────────────────────────────────────────────────────

library(utils)

# 1. Create a data directory
dir.create("data/scotese_paleodems", recursive = TRUE, showWarnings = FALSE)

# 2. Download the netCDF zip (~50–100 MB)
url <- paste0("https://www.earthbyte.org/webdav/ftp/Data_Collections/",
              "Scotese_Wright_2018_PaleoDEM/",
              "Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc.zip")

zip_path <- "data/scotese_paleodems/Scotese_Wright_2018_PaleoDEMs.zip"

# Increase timeout — the file is moderately large
options(timeout = 600)
download.file(url, destfile = zip_path, mode = "wb")

# 3. Unzip
unzip(zip_path, exdir = "data/scotese_paleodems/")

# 4. Inspect what's inside
list.files("data/scotese_paleodems/", pattern = "\\.nc$", recursive = TRUE)



library(terra)

# Example: load the Induan/Early Triassic DEM (~250 Ma)
dem_path <- "data/scotese_paleodems/Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc_v2/Map46.5_PALEOMAP_1deg_Late_Triassic_235Ma.nc"

dem <- terra::rast(dem_path)
print(dem)
plot(dem)

# Set CRS if not already set
if (is.na(terra::crs(dem)) || terra::crs(dem) == "") {
  terra::crs(dem) <- "EPSG:4326"
}
