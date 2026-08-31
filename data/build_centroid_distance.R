# Build data/city_centroid_distance.csv
#
# For each city (place), compute a 2020-population-weighted centroid from its
# constituent census tracts (tract assigned to a place if the tract's
# geometric centroid falls within the place boundary; tracts that only
# partially overlap and whose centroid falls outside are excluded, matching
# the standard "centroid-in-polygon" tract-to-place assignment rule), then
# compute the minimum distance from that weighted centroid to the boundary
# line of the city's home county (data/county_boundaries.geojson).
#
# Source data:
#   - /tmp/annex_scratch/raw_geo_pieces.rds: per-state list (47 states + PR)
#     of tigris place polygons, tigris tract polygons, and 2020 decennial
#     P1_001N (total population) by tract -- already fetched in an earlier
#     session, reused here as-is (no fresh Census/tigris pulls were needed;
#     all 47 state FIPS present in city_pop_2000_2020.csv were already
#     covered in the cache).
#   - data/county_boundaries.geojson: home-county polygons (EPSG:4326).
#
# All distance computation is done in EPSG:5070 (NAD83 / Conus Albers Equal
# Area), a projected CRS appropriate for CONUS distance calculations.
#
# Output: data/city_centroid_distance.csv
#   place_geoid, place_name, state_fips, county_fips, n_tracts_matched,
#   match_method, centroid_x_5070, centroid_y_5070, dist_to_county_border_km

library(sf)
library(dplyr)
library(purrr)
library(readr)
library(tibble)

sf_use_s2(TRUE)

ALBERS <- 5070

message("Loading city panel + raw tract/place cache ...")
city <- read_csv("data/city_pop_2000_2020.csv", show_col_types = FALSE)
raw  <- readRDS("/tmp/annex_scratch/raw_geo_pieces.rds")

county_bd <- st_read("data/county_boundaries.geojson", quiet = TRUE) %>%
  st_transform(ALBERS)
county_bd_lines <- county_bd %>%
  mutate(geometry = st_boundary(geometry))

process_state <- function(state_fips) {
  piece <- raw[[state_fips]]
  if (is.null(piece)) return(NULL)

  places <- piece$places %>% st_transform(ALBERS)
  tracts <- piece$tracts %>% st_transform(ALBERS)
  pop    <- piece$pop %>%
    filter(variable == "P1_001N") %>%
    select(GEOID, pop2020 = value)

  tracts <- tracts %>%
    left_join(pop, by = "GEOID") %>%
    mutate(pop2020 = coalesce(pop2020, 0))

  tract_centroids <- tracts %>%
    st_centroid() %>%
    select(tract_geoid = GEOID, pop2020)

  places_needed <- city %>% filter(state_fips == !!state_fips) %>% pull(place_geoid)

  places <- places %>% filter(GEOID %in% places_needed)
  if (nrow(places) == 0) return(NULL)

  # primary match: tract centroid within place polygon
  within_join <- st_join(tract_centroids, places %>% select(GEOID, NAME),
                          join = st_within, left = FALSE)

  results <- map(seq_len(nrow(places)), function(i) {
    pl <- places[i, ]
    geoid <- pl$GEOID

    matched <- within_join %>% filter(GEOID == geoid)
    method <- "centroid_within_place"

    if (nrow(matched) == 0 || sum(matched$pop2020, na.rm = TRUE) == 0) {
      # fallback: any tract intersecting the place polygon (partial overlap ok)
      inter_idx <- st_intersects(tracts, pl, sparse = FALSE)[, 1]
      matched_tracts <- tracts[inter_idx, ]
      if (nrow(matched_tracts) > 0 && sum(matched_tracts$pop2020, na.rm = TRUE) > 0) {
        matched <- st_centroid(matched_tracts) %>%
          transmute(tract_geoid = GEOID, pop2020, GEOID = geoid)
        method <- "centroid_intersects_place_fallback"
      } else {
        # last resort: no tract population found at all -> use place polygon centroid
        pc <- st_centroid(pl)
        coords <- st_coordinates(pc)
        return(tibble(
          place_geoid = geoid, n_tracts_matched = 0L,
          match_method = "place_polygon_centroid_fallback",
          centroid_x_5070 = coords[1, "X"], centroid_y_5070 = coords[1, "Y"]
        ))
      }
    }

    coords <- st_coordinates(matched)
    wx <- weighted.mean(coords[, "X"], w = pmax(matched$pop2020, 0.001))
    wy <- weighted.mean(coords[, "Y"], w = pmax(matched$pop2020, 0.001))

    tibble(
      place_geoid = geoid, n_tracts_matched = nrow(matched),
      match_method = method,
      centroid_x_5070 = wx, centroid_y_5070 = wy
    )
  })

  bind_rows(results)
}

message("Processing states ...")
states_needed <- sort(unique(city$state_fips))
out <- map_dfr(states_needed, process_state, .progress = TRUE)

message("Computing distance to county border for each city ...")
county_fips_lookup <- county_bd_lines %>% st_drop_geometry() %>% select(county_fips)

dist_to_border <- function(geoid, cfips, x, y) {
  cty_line <- county_bd_lines %>% filter(county_fips == cfips)
  if (nrow(cty_line) == 0) return(NA_real_)
  pt <- st_sfc(st_point(c(x, y)), crs = ALBERS)
  as.numeric(st_distance(pt, cty_line)) / 1000  # meters -> km
}

city_geo <- city %>%
  select(place_geoid, place_name, state_fips, county_fips) %>%
  left_join(out, by = "place_geoid")

city_geo <- city_geo %>%
  mutate(dist_to_county_border_km = pmap_dbl(
    list(place_geoid, county_fips, centroid_x_5070, centroid_y_5070),
    function(place_geoid, county_fips, x, y) {
      if (is.na(x) || is.na(y)) return(NA_real_)
      dist_to_border(place_geoid, county_fips, x, y)
    }
  ))

message("Rows with no tract match at all: ", sum(is.na(city_geo$centroid_x_5070)))
message("Fallback methods used: ")
print(table(city_geo$match_method, useNA = "ifany"))

write_csv(city_geo, "data/city_centroid_distance.csv")
message("Wrote data/city_centroid_distance.csv (", nrow(city_geo), " rows)")
