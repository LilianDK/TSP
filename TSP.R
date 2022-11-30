# Load packages ----------------------------------------------------------------
library(knitr)
library(dplyr)
library(ggplot2)
library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)
library(leaflet)
library(geosphere)
library(readr)

# Load data ----------------------------------------------------------------
locations <- read_delim("Santa_7_Distances.csv", 
                        delim = ";", escape_double = FALSE, col_types = cols(latitude = col_number(), 
                                                                             longitude = col_number()), trim_ws = TRUE)
name <- locations$name
popup <- locations$name

# Generate map 1 ---------------------------------------------------------------
leaflet(data = locations) %>% addTiles() %>%
  addMarkers(~longitude, ~latitude, popup = ~popup, label = ~popup)

# Data preparation for optimization --------------------------------------------
location_coords <- locations %>% select(longitude, latitude)

distance_matrix <- as.matrix(
  distm(location_coords, fun = distHaversine)
)/1000 #convert metres to kilometres

rownames(distance_matrix) <- locations$name
colnames(distance_matrix) <- locations$name

#specify the dimensions of the distance matrix
n <- nrow(locations)

#create a distance extraction function
dist_fun <- function(i, j) {
  vapply(seq_along(i), function(k) distance_matrix[i[k], j[k]], numeric(1L))
}

# Optimization model -----------------------------------------------------------
model <- MILPModel() %>%
  # binary variable is 1 if two nodes are connected and 0 otherwise
  add_variable(x[i, j], i = 1:n, j = 1:n, 
               type = "integer", lb = 0, ub = 1) %>%
  
  # a helper variable for the MTZ formulation of the tsp
  add_variable(u[i], i = 1:n, lb = 1, ub = n) %>% 
  
  # objective function set to minimize total travel distance
  set_objective(sum_expr(colwise(dist_fun(i, j)) * x[i, j], i = 1:n, j = 1:n), "min") %>%
  
  # elimination of circular motions
  set_bounds(x[i, i], ub = 0, i = 1:n) %>%
  
  # leave each city
  add_constraint(sum_expr(x[i, j], j = 1:n) == 1, i = 1:n) %>%
  
  # visit each city
  add_constraint(sum_expr(x[i, j], i = 1:n) == 1, j = 1:n) %>%
  
  # ensure no subtours (arc constraints)
  add_constraint(u[i] >= 2, i = 2:n) %>% 
  add_constraint(u[i] - u[j] + 1 <= (n - 1) * (1 - x[i, j]), i = 2:n, j = 2:n)
model

# Results output ---------------------------------------------------------------
result <- solve_model(model, with_ROI(solver = "glpk", verbose = TRUE))

result_val <- round(objective_value(result), 2)
result_val

solution <- get_solution(result, x[i, j]) %>% 
  filter(value > 0)

# Generate map 2 ---------------------------------------------------------------
paths <- select(solution, i, j) %>% 
  rename(from = i, to = j) %>% 
  mutate(trip_id = row_number()) %>% 
  inner_join(locations, by = c("from" = "id"))

paths_leaflet <- paths[1,]
paths_row <- paths[1,]

for (i in 1:n) {
  paths_row <- paths %>% filter(from == paths_row$to[1])
  
  paths_leaflet <- rbind(paths_leaflet, paths_row)
}

leaflet() %>% 
  addTiles() %>%
  addMarkers(data = locations, ~longitude, ~latitude) %>% 
  addPolylines(data = paths_leaflet, ~longitude, ~latitude, weight = 2)
