library(spatialRF)
library(kableExtra)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(randomForestExplainer)
library(pdp)


############ Load CSV file ##############
df_ground_photo <- read.csv("./CSV/ground_level_photo_extracted.csv") %>%
  mutate(
    geo_parsed = map(.geo, ~ fromJSON(.)$coordinates),
    x = map_dbl(geo_parsed, 1),
    y = map_dbl(geo_parsed, 2)
  ) %>%
  select(-geo_parsed) %>% 
  select(-c(system.index, .geo, vegType)) %>%
  na.omit() %>% 
  filter(Month >= 5) %>% 
  sample_frac(0.7)

head(df_ground_photo)

# Load a script
source("./Script/functions.R")

bands <- c("blue","green","red","nir","swir1","swir2")
vi <- c("bsi","evi","exg","msavi","ndmi","ndvi","osavi","satvi","savi","swirRatio")
minerals <- c("calcite","chlorite","dolomite","goethite","hematite","illite_muscovite","kaolinite") # excluding gypsum and vermiculite
soil <- c("soil_depth","clay","sand","silt","stone", "nitrogen","phosphorous","potassium","ph","soc")
topography <- c("elevation", "slope", "aspect")

predictors_1 <- c(vi)
predictors_all <- c(vi, bands, minerals, soil, topography)

distance_matrix <-data.frame(df_ground_photo$x, df_ground_photo$y) %>% 
  as.matrix() %>% 
  dist(method = 'euclidean') %>% 
  as.matrix()

#names of the response variable and the predictors
dependent.variable.name <- "BareGround"
predictor.variable.names <- predictors_all

#coordinates of the cases
xy <- df_ground_photo[, c("x", "y")]
colnames(xy) <- c("x","y")
#distance matrix
distance.matrix <- distance_matrix

#distance thresholds (same units as distance_matrix)
distance.thresholds <- c(0,0.01,0.05,0.1,0.15, 0.2, 0.5, 1)

#random seed for reproducibility
random.seed <- 1

preference.order <- c(
  "satvi",
  "swirRatio",
  "savi",
  "msavi"
)

predictor.variable.names <- spatialRF::auto_cor(
  x = df_ground_photo[, predictor.variable.names],
  cor.threshold = 0.8,
  preference.order = preference.order
) %>% 
  spatialRF::auto_vif(
    vif.threshold = 10,
    preference.order = preference.order
  )

predictor.variable.names$selected.variables

model.non.spatial <- spatialRF::rf(
  data = df_ground_photo,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,
  verbose = FALSE
)

spatialRF::plot_importance(
  model.non.spatial,
  verbose = FALSE
)

model.non.spatial <- spatialRF::rf_importance(
  model = model.non.spatial
)

model.non.spatial <- spatialRF::rf_evaluate(
  model = model.non.spatial,
  xy = xy,                  #data coordinates
  repetitions = 30,         #number of spatial folds
  training.fraction = 0.75, #training data fraction on each fold
  metrics = "r.squared",
  seed = random.seed,
  verbose = FALSE
)

names(model.non.spatial$evaluation)

pr <- df_ground_photo[, c("x", "y")]
pr$group.2 <- pr$group.1 <- "Training"
pr[model.non.spatial$evaluation$spatial.folds[[1]]$testing, "group.1"] <- "Testing"
pr[model.non.spatial$evaluation$spatial.folds[[25]]$testing, "group.2"] <- "Testing"
head(pr)

p1 <- ggplot2::ggplot() +
  #ggplot2::geom_sf(data = world, fill = "white") +
  ggplot2::geom_point(data = pr,
                      ggplot2::aes(
                        x = x,
                        y = y,
                        color = group.1
                      ),
                      size = 2
  ) +
  ggplot2::scale_color_viridis_d(
    direction = -1, 
    end = 0.5, 
    alpha = 0.8, 
    option = "F"
  ) +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "Group") +
  ggplot2::ggtitle("Spatial fold 1") + 
  ggplot2::theme(
    legend.position = "none", 
    plot.title = ggplot2::element_text(hjust = 0.5)
  ) + 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")

p2 <- ggplot2::ggplot() +
  #ggplot2::geom_sf(data = world, fill = "white") +
  ggplot2::geom_point(data = pr,
                      ggplot2::aes(
                        x = x,
                        y = y,
                        color = group.2
                      ),
                      size = 2
  ) +
  ggplot2::scale_color_viridis_d(
    direction = -1, 
    end = 0.5, 
    alpha = 0.8, 
    option = "F"
  ) +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "Group") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5)
  ) + 
  ggplot2::ggtitle("Spatial fold 25") + 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("")

p1 | p2

model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = random.seed
)

spatialRF::plot_moran(
  model.non.spatial, 
  verbose = FALSE
)

spatialRF::plot_moran(
  model.spatial, 
  verbose = FALSE
)

# indicating residuals are not correlated at any distance

p1 <- spatialRF::plot_importance(
  model.non.spatial, 
  verbose = FALSE) + 
  ggplot2::ggtitle("Non-spatial model") 

p2 <- spatialRF::plot_importance(
  model.spatial,
  verbose = FALSE) + 
  ggplot2::ggtitle("Spatial model")

p1 | p2 

comparison <- spatialRF::rf_compare(
  models = list(
    `Non-spatial` = model.non.spatial,
    `Spatial` = model.spatial
  ),
  xy = xy,
  repetitions = 30,
  training.fraction = 0.8,
  metrics = "r.squared",
  seed = random.seed
)
