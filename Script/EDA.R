# Load required libraries
library(tidyverse)
library(ggplot2)


############ Load and Clean CSV Data ############
# Read CSV and remove unnecessary columns and rows with NA values
df_ground_photo <- read.csv("../CSV/ground_level_photo_extracted.csv") %>%
  select(-c(system.index, .geo, vegType)) %>%
  na.omit()  # Resulting in 1,017 rows

############ Define Variable Groups #############
# Define variable groups for plotting
bands <- c("blue", "green", "red", "nir", "swir1", "swir2")
vi <- c("bsi", "evi", "exg", "msavi", "ndmi", "ndvi", "osavi", "satvi", "savi", "swirRatio","tc_bright","tc_green","tc_wet")
minerals <- c("calcite", "chlorite", "dolomite", "goethite", "hematite", "illite_muscovite", "kaolinite")
soil <- c("soil_depth","clay","sand","silt","stone", "nitrogen","phosphorous","potassium","ph","soc")
topography <- c("elevation", "slope", "aspect")

# Define color and order for bioregions
fill_order <- c("Bushmanland", "Upper Karoo", "Lower Karoo",
                "Desert", "Grassland", "Albany Thicket", "Azonal Vegetation")

############ Format Data for Plotting ############
df_plot <- df_ground_photo %>%
  select(all_of(vi), all_of(bands), all_of(minerals), all_of(soil), all_of(topography),
         BareGround, Month, BioRegion, Biome) %>%
  mutate(season = ifelse(Month >= 5, "Non-Growing Season", "Growing Season")) %>%
  mutate(Biome_BioRegion = case_when(
    Biome == "Albany Thicket" ~ "Albany Thicket",
    Biome == "Azonal Vegetation" ~ "Azonal Vegetation",
    Biome == "Desert" ~ "Desert",
    Biome == "Grassland" ~ "Grassland",
    TRUE ~ BioRegion
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

# Reshape data for each variable group
df_vi_plot <- df_plot %>%
  select(BareGround, all_of(vi), Biome_BioRegion, colour, season) %>%
  pivot_longer(cols = all_of(vi), names_to = "vegetation_indices", values_to = "value") %>%
  as.data.frame()

df_bands_plot <- df_plot %>%
  select(BareGround, all_of(bands), Biome_BioRegion, colour, season) %>%
  pivot_longer(cols = all_of(bands), names_to = "bands", values_to = "value") %>%
  as.data.frame()

df_minerals_plot <- df_plot %>%
  select(BareGround, all_of(minerals), Biome_BioRegion, colour, season) %>%
  pivot_longer(cols = all_of(minerals), names_to = "minerals", values_to = "value") %>%
  as.data.frame()

df_soil_plot <- df_plot %>%
  select(BareGround, all_of(soil), Biome_BioRegion, colour, season) %>%
  pivot_longer(cols = all_of(soil), names_to = "soil_properties", values_to = "value") %>%
  as.data.frame()

df_topography_plot <- df_plot %>%
  select(BareGround, all_of(topography), Biome_BioRegion, colour, season) %>%
  pivot_longer(cols = all_of(topography), names_to = "topography", values_to = "value") %>%
  as.data.frame()

############ Generate Plots ############
# Vegetation indices: Part 1 (first 5 indices)
p1 <- df_vi_plot %>%
  filter(vegetation_indices %in% vi[1:6]) %>%
  mutate(vegetation_indices = str_to_upper(vegetation_indices),
         vegetation_indices = case_when(
           vegetation_indices == "EXG" ~ "ExG",
           vegetation_indices == "SWIRRATIO" ~ "SWIR ratio",
           TRUE ~ vegetation_indices)) %>%
  ggplot(aes(x = BareGround, y = value, color = colour)) +
  geom_point(shape = 16, alpha = 0.4, size = 2) +
  scale_colour_identity(guide = "legend", labels = fill_order) +
  facet_grid(vegetation_indices ~ season, scales = "free") +
  labs(x = "Photo estimated fraction of bare ground (%)",
       y = "Variable measure",
       colour = "Biome/Bioregion") +
  theme_bw()

# Vegetation indices: Part 2 (next 5 indices)
p2 <- df_vi_plot %>%
  filter(vegetation_indices %in% vi[7:13]) %>%
  mutate(vegetation_indices = str_to_upper(vegetation_indices),
         vegetation_indices = case_when(
           vegetation_indices == "EXG" ~ "ExG",
           vegetation_indices == "SWIRRATIO" ~ "SWIR ratio",
           TRUE ~ vegetation_indices)) %>%
  ggplot(aes(x = BareGround, y = value, color = colour)) +
  geom_point(shape = 16, alpha = 0.4, size = 2) +
  scale_colour_identity(guide = "legend", labels = fill_order) +
  facet_grid(vegetation_indices ~ season, scales = "free") +
  labs(x = "Photo estimated fraction of bare ground (%)",
       y = "Variable measure",
       colour = "Biome/Bioregion") +
  theme_bw()

# Bands
p3 <- df_bands_plot %>%
  mutate(bands = str_to_upper(bands),
         bands = case_when(
           bands == "RED" ~ "Red",
           bands == "GREEN" ~ "Green",
           bands == "BLUE" ~"Blue",
           TRUE ~ bands)) %>%
  mutate(bands = factor(bands, levels = c("Blue","Green","Red","NIR","SWIR1","SWIR2"))) %>% 
  ggplot(aes(x = BareGround, y = value, color = colour)) +
  geom_point(shape = 16, alpha = 0.4, size = 2) +
  scale_colour_identity(guide = "legend", labels = fill_order) +
  facet_grid(bands ~ season, scales = "free") +
  labs(x = "Photo estimated fraction of bare ground (%)",
       y = "Variable measure",
       colour = "Biome/Bioregion") +
  theme_bw()

# Mineral indices
p4 <- df_minerals_plot %>%
  mutate(minerals = str_to_title(minerals)) %>%
  ggplot(aes(x = BareGround, y = value, color = colour)) +
  geom_point(shape = 16, alpha = 0.4, size = 2) +
  scale_colour_identity(guide = "legend", labels = fill_order) +
  facet_wrap(vars(minerals), scales = "free", ncol = 4) +
  labs(x = "Photo estimated fraction of bare ground (%)",
       y = "Variable measure",
       colour = "Biome/Bioregion") +
  theme_bw()

# Soil properties: Part 1 (first 5 indices)
p5 <- df_soil_plot %>%
  filter(soil_properties %in% soil[1:5]) %>%
  mutate(soil_properties = str_to_title(soil_properties),
         soil_properties = case_when(
           soil_properties == "Soil_depth" ~ "Soil Depth",
           TRUE ~ soil_properties)) %>%
  mutate(soil_properties = factor(soil_properties, levels = c("Soil Depth","Clay","Silt","Sand","Stone"))) %>% 
  ggplot(aes(x = BareGround, y = value, color = colour)) +
  geom_point(shape = 16, alpha = 0.4, size = 2) +
  scale_colour_identity(guide = "legend", labels = fill_order) +
  facet_grid(soil_properties ~ season, scales = "free") +
  labs(x = "Photo estimated fraction of bare ground (%)",
       y = "Variable measure",
       colour = "Biome/Bioregion") +
  theme_bw()

# Soil properties: Part 1 (next 5 indices)
p6 <- df_soil_plot %>%
  filter(soil_properties %in% soil[6:10]) %>%
  mutate(soil_properties = str_to_title(soil_properties),
         soil_properties = case_when(
           soil_properties == "Soc" ~ "SOC",
           soil_properties == "Ph" ~ "pH",
           TRUE ~ soil_properties)) %>%
  mutate(soil_properties = factor(soil_properties, levels = c("Nitrogen","Phosphorous","Potassium","pH","SOC"))) %>% 
  ggplot(aes(x = BareGround, y = value, color = colour)) +
  geom_point(shape = 16, alpha = 0.4, size = 2) +
  scale_colour_identity(guide = "legend", labels = fill_order) +
  facet_grid(soil_properties ~ season, scales = "free") +
  labs(x = "Photo estimated fraction of bare ground (%)",
       y = "Variable measure",
       colour = "Biome/Bioregion") +
  theme_bw()

# Topography variables
p7 <- df_topography_plot %>%
  mutate(topography = str_to_title(topography)) %>%
  ggplot(aes(x = BareGround, y = value, color = colour)) +
  geom_point(shape = 16, alpha = 0.4, size = 2) +
  scale_colour_identity(guide = "legend", labels = fill_order) +
  facet_wrap(vars(topography), scales = "free") +
  labs(x = "Photo estimated fraction of bare ground (%)",
       y = "Variable measure",
       colour = "Biome/Bioregion") +
  theme_bw()

############ Save Plots as PNG Files ############
ggsave("./script/figures/vi_cor_1.png", plot = p1, width = 500, height = 700, units = "px", dpi = 300, scale = 4)
ggsave("./script/figures/vi_cor_2.png", plot = p2, width = 500, height = 700, units = "px", dpi = 300, scale = 4)
ggsave("./script/figures/bands_cor.png", plot = p3, width = 500, height = 700, units = "px", dpi = 300, scale = 4)
ggsave("./script/figures/minerals_cor.png", plot = p4, width = 500, height = 200, units = "px", dpi = 300, scale = 6)
ggsave("./script/figures/soil_cor_1.png", plot = p5, width = 500, height = 700, units = "px", dpi = 300, scale = 4)
ggsave("./script/figures/soil_cor_2.png", plot = p6, width = 500, height = 700, units = "px", dpi = 300, scale = 4)
ggsave("./script/figures/topo_cor.png", plot = p7, width = 500, height = 150, units = "px", dpi = 300, scale = 5)


############ Biplot for minerals and Vegetation indices ############
df_biplot <- df_plot %>%
  select(BareGround, all_of(vi), all_of(minerals), Biome_BioRegion, colour, season) %>%
  pivot_longer(cols = all_of(vi), names_to = "vegetation_indices", values_to = "value") %>%
  as.data.frame()

df_biplot %>% 
  filter(vegetation_indices %in% vi[1:5]) %>%
  # filter(vegetation_indices %in% vi[6:10]) %>%
  ggplot(aes(x = illite_muscovite, y = value, color = BareGround)) +
  geom_point(shape = 16, alpha = 0.4, size = 2) +
  scale_colour_gradient2(low="red",mid="yellow", high="darkgreen", midpoint =50) +
  facet_grid(vegetation_indices~season, scales="free") +
  labs(x = "Variable measure",
       y = "NDVI",
       colour = "Fraction of Bare Ground (%)") +
  theme_bw()

############ Correlation matrix for VIs ############
library(corrplot)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
df_plot %>% 
  select(all_of(vi)) %>% 
  cor() %>% 
  corrplot(type = "upper",
           method = "color",
           col=col(200),
           addCoef.col = "black", # Add coefficient of correlation
           number.cex = 0.8, # font size
           number.font = 1,
           tl.col = "black", 
           tl.srt = 45,
           diag=FALSE)
