library(sf)
library(tidyverse)
library(jsonlite)
library(ggplot2)
library(ggforce)
library(RColorBrewer)

#setwd("C:/Users/turat/Dropbox/01_UCT/04_DegradationMapping/07_Analysis/01_Estimate Cover from Ground Level Photo/30_PLSR_Model/BaregroundPrediction")
setwd("~/01_UCT/04_Degradation Mapping/07_Analysis/01_Estimate Cover from Ground Level Photo/30_PLSR_Model/BaregroundPrediction")


# Load necessary shapefile
study_area <- st_read("./SHP/StudyArea2024.shp") %>% 
  st_transform(crs = 32734)
south_africa <- st_read("./SHP/SouthernAfricanCountries.shp") %>% 
  filter(NAME=="South Africa")
cities <- st_read("./SHP/Cities.shp") %>% 
  filter(name %in% c("Upington","Kimberly","Beaufort West", "Graaf-Reinet","Springbok",
                     "Calvinia","Gqeberha"))

# ------------------------------------------------------------------------------
# Study Area Map

# Load ground-level photo metadata and convert to spatial data frame
sf_ground_photo <- read.csv("./CSV/ground_level_photo_extracted.csv") %>%
  mutate(
    # Parse JSON from .geo column to extract coordinates
    geo_parsed = map(.geo, ~ fromJSON(.)$coordinates),
    longitude = map_dbl(geo_parsed, 1),
    latitude = map_dbl(geo_parsed, 2)
  ) %>%
  dplyr::select(-geo_parsed, -system.index, -.geo, -vegType) %>%  # Drop unnecessary columns
  na.omit() %>%  # Remove rows with missing values (resulting in 1,017 rows)
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%  # Convert to sf object in WGS84
  st_transform(crs = 32734)  # Transform to UTM for spatial blocking

# add colour information to the dataframe
sf_base_map <- study_area %>% 
  mutate(Biome_BioRegion = case_when(
    T_BIOME == "Albany Thicket" ~ "Albany Thicket",
    T_BIOME == "Azonal Vegetation" ~ "Azonal Vegetation",
    T_BIOME == "Desert" ~ "Desert",
    T_BIOME == "Grassland" ~ "Grassland",
    TRUE ~ T_BIOREGIO
  )) %>%
  mutate(colour = case_when(
    Biome_BioRegion == "Albany Thicket" ~ "darkolivegreen",
    Biome_BioRegion == "Azonal Vegetation" ~ "blue1",
    Biome_BioRegion == "Desert" ~ "grey60",
    Biome_BioRegion == "Grassland" ~ "green2",
    Biome_BioRegion == "Bushmanland Bioregion" ~ "orange",
    Biome_BioRegion == "Lower Karoo Bioregion" ~ "brown",
    Biome_BioRegion == "Upper Karoo Bioregion" ~ "hotpink"
  )) %>%
  mutate(colour = factor(colour, levels = c("orange", "brown", "hotpink", "grey60", "green2", "darkolivegreen", "blue1")))


nama_karoo <- sf_base_map %>% 
  filter(T_BIOME=="Nama-Karoo")

gg_map_biome <- ggplot()+
  geom_sf(data=south_africa)+
  geom_sf(data=nama_karoo,fill="brown",color=NA, alpha=0.6) +
  geom_sf(data=cities)+
  coord_sf(xlim=c(17,28), ylim=c(-35,-26)) +
  theme_minimal()

gg_map_photosites <- ggplot() +
  geom_sf(data = sf_base_map, aes(fill=colour), color=NA, alpha=0.7)+
  scale_fill_identity()+
  geom_sf(data = sf_ground_photo,size=0.2) +
  coord_sf(xlim=c(17,28), ylim=c(-35,-26), crs=4326) +
  theme_minimal()

# export maps
ggsave(file="./Script/figures/nama-karoo.png",
       plot=gg_map_biome,
       device="png",
       width=1000, height=1000,
       units="px")

ggsave(file="./Script/figures/photo_sites.png",
       plot=gg_map_photosites,
       device="png",
       width=1000, height=1000,
       units="px")

# ------------------------------------------------------------------------------
# VIP of PLSR plot

# load pls models
best_non_growing <- readRDS("./Script/models/pls_non_growing.rds")
best_growing <- readRDS("./Script/models/pls_growing.rds")

# preparing a dataframe 
df_vip_non_growing <- VIP(best_non_growing$finalModel, best_non_growing$bestTune$ncomp) %>% 
  as.data.frame() %>% 
  rownames_to_column("vi") %>%
  rename(vip=".") %>% 
  mutate(colour = case_when(
    vi == "ndmi" ~ "blue3",
    vi == "bsi" ~ "red3",
    vi == "tc_bright" ~ "red3",
    vi == "swirRatio" ~ "orange",
    TRUE ~ "green4"
  )) %>% 
  mutate(vi=case_when(
    vi == "swirRatio" ~ "SWIR ratio",
    vi == "tc_bright" ~ "TC Brightness",
    vi == "tc_green" ~ "TC Greeness",
    TRUE ~ toupper(vi)
  )) %>% 
  mutate(vi = fct_reorder(vi, vip))

df_vip_growing <- VIP(best_growing$finalModel, best_growing$bestTune$ncomp) %>% 
  as.data.frame() %>% 
  rownames_to_column("vi") %>%
  rename(vip=".") %>% 
  mutate(colour = case_when(
    vi == "ndmi" ~ "blue3",
    vi == "bsi" ~ "red3",
    vi == "tc_bright" ~ "red3",
    vi == "swirRatio" ~ "orange",
    TRUE ~ "green4"
  )) %>% 
  mutate(vi=case_when(
    vi == "swirRatio" ~ "SWIR ratio",
    vi == "tc_bright" ~ "TC Brightness",
    TRUE ~ toupper(vi)
  )) %>% 
  mutate(vi = fct_reorder(vi, vip))

gg_vip_non_growing <- ggplot(data=df_vip_non_growing)+
  #  geom_hline(yintercept=c(1:11), linetype="dotted", alpha=0.5)+
  geom_segment(aes(x=vip, xend=0, y=vi, yend=vi, color=colour),
               linewidth=1.5, alpha=0.7)+
  geom_point(aes(x=vip,y=vi,color=colour),
             shape = 16,
             size=2.5,
             stroke=2)+
  scale_color_identity()+
  geom_vline(xintercept=1.0, linetype="dotted")+
  xlab("Variable Importance in Projection (VIP)")+
  ylab("Vegetation Indices")+
  xlim(c(0,1.5))+
  theme_bw()

gg_vip_growing <- ggplot(data=df_vip_growing)+
  #  geom_hline(yintercept=c(1:10), linetype="dotted", alpha=0.5)+
  geom_segment(aes(x=vip, xend=0, y=vi, yend=vi, color=colour),
               linewidth=1.5, alpha=0.7)+
  geom_point(aes(x=vip,y=vi,color=colour),
             shape = 16,
             size=2.5,
             stroke=2)+
  scale_color_identity()+
  geom_vline(xintercept=1.0, linetype="dotted")+
  xlab("Variable Importance in Projection (VIP)")+
  ylab("Vegetation Indices")+
  xlim(c(0,1.5))+
  theme_bw()

ggsave(file="./Script/figures/vip_non_growing.png",
       plot=gg_vip_non_growing,
       device="png",
       width=1500, height=800,
       units="px")

ggsave(file="./Script/figures/vip_growing.png",
       plot=gg_vip_growing,
       device="png",
       width=1500, height=800,
       units="px")

# ------------------------------------------------------------------------------
# PLSR Biplot

# Function to plot PLS biplot
plot_pls_biplot <- function(pls_model, comps = 1:2, label_scale = 1.2, score_scale=8) {
  
  # Extract scores
  S <- scores(pls_model)[, comps, drop = FALSE] %>% as.data.frame()
  colnames(S) <- c("Comp1", "Comp2")
  S$.outcome <- pls_model$model$.outcome
  
  # Correlation of predictors with scores
  df.cor <- model.matrix(pls_model) %>%
    as.data.frame() %>%
    cor(S[, c("Comp1", "Comp2")]) %>%
    as.data.frame() %>%
    mutate(variable = rownames(.)) %>%
    select(variable, axis_1 = "Comp1", axis_2 = "Comp2")
  
  # Correlation of Y
  Y.cor <- cor(pls_model$Yscores[, 1], S[, c("Comp1", "Comp2")])
  df.cor <- df.cor %>%
    bind_rows(data.frame(variable = "Y", axis_1 = Y.cor[1], axis_2 = Y.cor[2])) %>%
    mutate(axis_1 = axis_1*score_scale,
           axis_2 = axis_2*score_scale) %>% 
    mutate(label_x = axis_1 * label_scale,
           label_y = axis_2 * label_scale,
           color = ifelse(variable == "Y", "red", "blue"))
  
  # variance explained
  r2x <- explvar(pls_model)[comps]
  r2y_cum <- R2(pls_model)$val[comps+1]
  r2y_each <- c(r2y_cum[1], diff(r2y_cum))
  
  xlab_str <- sprintf("Score axis, Comp %d (R²X=%.1f%%, R²Y=%.1f%%)", 1, r2x[1], r2y_each[1]*100)
  ylab_str <- sprintf("Score axis, Comp %d (R²X=%.1f%%, R²Y=%.1f%%)", 2, r2x[2], r2y_each[2]*100)
  
  # Create plot
  p <- ggplot() +
    # Score points colored by outcome
    geom_point(data = S, aes(x = Comp1, y = Comp2, color = .outcome),
               size = 2, alpha = 0.6) +
    
    # Unit circle
    geom_circle(data = data.frame(x0 = 0, y0 = 0, r = score_scale),
                aes(x0 = x0, y0 = y0, r = r),
                inherit.aes = FALSE, color = "gray50") +
    
    # Variable arrows (fixed colors: blue or red)
    geom_segment(data = df.cor,
                 aes(x = 0, y = 0, xend = axis_1, yend = axis_2),
                 color = df.cor$color,
                 arrow = arrow(length = unit(0.2, "cm")),
                 linewidth = 1.5, alpha = 0.5) +
    
    # Labels for arrows
    geom_text(data = df.cor,
              aes(x = label_x, y = label_y, label = variable),
              color = df.cor$color, size = 4) +
    
    # Axes
    geom_hline(yintercept = 0, color = "gray80") +
    geom_vline(xintercept = 0, color = "gray80") +
    
    # Layout
    # 主軸・副軸
    scale_y_continuous(
      name = ylab_str,
      sec.axis = sec_axis(~ . /score_scale , name = "Correlation axis (scaled)"),
     limits = c(-(score_scale+1),score_scale+1)
    ) +
    scale_x_continuous(
      name = xlab_str,
     limits = c(-(score_scale+1),score_scale+1)
    ) +
    coord_fixed() +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)) +
    labs(x = paste0("Comp ", comps[1]),
         y = paste0("Comp ", comps[2]))
  
  # 連続変数ならグラデーションスケールを適用
  p <- p + scale_color_gradientn(
    colors = rev(brewer.pal(11, "BrBG")),
    limits = c(0,100),
    name = "Bare Ground (%)"
  )
  return(p)
}

gg_biplot_non_growing <- plot_pls_biplot(best_non_growing$finalModel)
gg_biplot_growing <- plot_pls_biplot(best_growing$finalModel)

ggsave(file="./Script/figures/biplot_non_growing.svg",
       plot=gg_biplot_non_growing,
       device="svg",
       width=2000, height=1500,
       units="px")

ggsave(file="./Script/figures/biplot_growing.svg",
       plot=gg_biplot_growing,
       device="svg",
       width=2000, height=1500,
       units="px")

# ------------------------------------------------------------------------------
# Validation plot

# ground level photo validation
valid_non_growing <- predict(best_non_growing, sf_non_growing_split$valid_df)
plot_obs_pred(sf_non_growing_split$valid_df$BareGround, valid_non_growing)

valid_growing <- predict(best_growing, sf_growing_log_split$valid_df)
plot_obs_pred(sf_growing_log_split$valid_df$BareGround, valid_growing)

png("./Script/figures/photo_validation_non_growing.png", width = 1500, height = 1600, res=300)  
plot_obs_pred(sf_non_growing_split$valid_df$BareGround, valid_non_growing)
dev.off() 

png("./Script/figures/photo_validation_growing.png", width = 1500, height = 1600, res=300)  
plot_obs_pred(sf_growing_split$valid_df$BareGround, valid_growing)
dev.off() 


# Transect Data Validation
df_transect <- read.csv("./CSV/transect_extracted.csv")

df_transect_non_growing <-   df_transect %>% 
  filter(Month >= 5) %>% 
  select(BareGround,all_of(vi),Notes) %>%
  na.omit()  # 52 rows

df_transect_growing <-  df_transect %>% 
  filter(Month < 5) %>% 
  select(BareGround,all_of(vi),Notes) %>%
  na.omit()  # 48 rows

# non-growing season
transect_non_growing <- predict(best_non_growing, df_transect_non_growing)
# plot_slope_flat(df_transect_non_growing, transect_non_growing)

# growing season
# log transformation
trans_vars <- func_log_transform(sf_growing, vi)[[2]] # "bsi","evi","msavi", "ndvi", "osavi", "satvi","savi","swirRatio","tc_bright","tc_green" 
# add_log <- func_log_transform(sf_growing, vi)[[4]] # 0

func_log_transform(sf_growing, vi)[[2]]
df_transect_growing_log <- df_transect_growing %>% 
  mutate(across(all_of(trans_vars), log))

transect_growing <- predict(best_growing, df_transect_growing_log)
plot_obs_pred(df_transect_growing$BareGround, transect_growing)

# export as png files
png("./Script/figures/transect_validation_non_growing.png", width = 1500, height = 1600, res=300)  
plot_slope_flat(df_transect_non_growing, transect_non_growing)
dev.off() 

# export as png files
png("./Script/figures/transect_validation_growing.png", width = 1500, height = 1600, res=300)  
plot_slope_flat(df_transect_growing, transect_growing)
dev.off() 
