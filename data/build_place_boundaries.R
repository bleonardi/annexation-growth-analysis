# Build data/city_place_boundaries.geojson
#
# For the subset of cities plotted on map.qmd (top ~100 by 2020 population,
# plus Columbus OH and Cincinnati OH forced in as the narrative anchor),
# pull the already-fetched tigris place polygons out of the reusable cache
# at /tmp/annex_scratch/raw_geo_pieces.rds (the same cache used to build
# data/city_centroid_distance.csv -- see build_centroid_distance.R) rather
# than hitting tigris::places() fresh.
#
# Output: data/city_place_boundaries.geojson (EPSG:4326), one polygon per
# selected place, with place_geoid/place_name/state_fips carried through for
# joining against city_pop_2000_2020.csv / city_centroid_distance.csv.

library(sf)
library(dplyr)
library(readr)

sf_use_s2(TRUE)

message("Loading city panel ...")
city <- read_csv("data/city_pop_2000_2020.csv", show_col_types = FALSE)

# Selection: top 100 by 2020 population, plus Columbus/Cincinnati OH forced in
top100 <- city %>%
  arrange(desc(population_2020)) %>%
  slice_head(n = 100) %>%
  pull(place_geoid)

anchors <- city %>%
  filter(place_name %in% c("Columbus", "Cincinnati"), state_name == "Ohio") %>%
  pull(place_geoid)

selected_geoids <- union(top100, anchors)
message("Selected ", length(selected_geoids), " place geoids (top 100 by pop2020 + Columbus/Cincinnati OH).")

message("Loading raw tract/place cache ...")
raw <- readRDS("/tmp/annex_scratch/raw_geo_pieces.rds")

state_fips_needed <- city %>%
  filter(place_geoid %in% selected_geoids) %>%
  pull(state_fips) %>%
  unique() %>%
  sort()

places_list <- lapply(state_fips_needed, function(sf_code) {
  piece <- raw[[sf_code]]
  if (is.null(piece)) return(NULL)
  piece$places %>%
    filter(GEOID %in% selected_geoids) %>%
    select(place_geoid = GEOID, place_name_geo = NAME, state_fips = STATEFP)
})

places_sel <- bind_rows(places_list) %>%
  st_as_sf() %>%
  st_transform(4326)

message("Matched ", nrow(places_sel), " of ", length(selected_geoids), " selected place polygons.")

st_write(places_sel, "data/city_place_boundaries.geojson", delete_dsn = TRUE, quiet = TRUE)
message("Wrote data/city_place_boundaries.geojson")
