###############################################
### 1. Load libraries
###############################################
library(tidyverse)
library(igraph)
library(ggraph)
library(readr)
library(vegan)
library(cluster)
library(countrycode)
library(geosphere)
library(patchwork)
library(maps)

###############################################
### 2. Load metadata
###############################################
genomes <- read_csv("/home/usuario/virgibacillus_analysis/onlyvirgi.csv") %>%
  rename(
    gcf_raw = gcf,
    species = species,
    origin = geographic_origin,
    environment = environment_type
  ) %>%
  mutate(gcf = str_extract(gcf_raw, "GCF_[0-9]+\\.[0-9]+"))

###############################################
### 3. Load BiG-SCAPE network
###############################################
edges <- read_tsv("mix_onlyvirgi.network") %>%
  rename(
    from = `Clustername 1`,
    to   = `Clustername 2`,
    distance = `Raw distance`,
    jaccard = `Jaccard index`,
    bgc_class = `Combined group`
  )

nodes <- read_tsv("mix_onlyvirgi_Network_Annotations_Full.tsv") %>%
  rename(
    bgc = BGC,
    product = `Product Prediction`,
    class = `BiG-SCAPE class`,
    organism = Organism
  ) %>%
  mutate(gcf = str_extract(bgc, "GCF_[0-9]+\\.[0-9]+")) %>%
  left_join(genomes, by = "gcf")

###############################################
### 4. Plot BGC network
###############################################
g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

ggraph(g, layout = "fr") +
  geom_edge_link(aes(alpha = 1 - distance), color = "grey70") +
  geom_node_point(aes(color = species, shape = class), size = 4) +
  theme_void() +
  theme(legend.position = "right")

###############################################
### 5. Build GCF matrix
###############################################
clusters <- read_tsv("mix_onlyvirgi_clustering.tsv") %>%
  rename(bgc = `#BGC Name`, family = `Family Number`) %>%
  mutate(gcf = str_extract(bgc, "GCF_[0-9]+\\.[0-9]+"))

mat <- clusters %>%
  select(gcf, family) %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = family,
    values_from = value,
    values_fill = 0
  )

gcf_matrix <- mat %>%
  column_to_rownames("gcf") %>%
  as.matrix()

###############################################
### 6. NMDS
###############################################
bray_dist <- vegdist(gcf_matrix, method = "bray")

nmds <- metaMDS(gcf_matrix, distance = "bray", k = 2, trymax = 100)

nmds_points <- as.data.frame(nmds$points) %>%
  mutate(gcf = rownames(.)) %>%
  left_join(genomes, by = "gcf")

# NMDS plots
p_species <- ggplot(nmds_points, aes(MDS1, MDS2, color = species)) +
  geom_point(size = 4) + theme_minimal()

p_environment <- ggplot(nmds_points, aes(MDS1, MDS2, color = environment)) +
  geom_point(size = 4) + theme_minimal()

p_origin <- ggplot(nmds_points, aes(MDS1, MDS2, color = origin)) +
  geom_point(size = 4) + theme_minimal()

p_species / p_environment / p_origin

###############################################
### 7. PERMANOVA + dispersion
###############################################
adonis2(bray_dist ~ environment, data = nmds_points)
anova(betadisper(bray_dist, nmds_points$environment))

###############################################
### 8. Geographic distances
###############################################
data(world.cities)

country_coords <- world.cities %>%
  group_by(country.etc) %>%
  slice_max(pop) %>%
  ungroup() %>%
  mutate(country_std = case_when(
    country.etc == "Korea, South" ~ "South Korea",
    TRUE ~ country.etc
  )) %>%
  select(country_std, lat, lon = long)

# Clean country names
nmds_geo <- nmds_points %>%
  mutate(origin = str_trim(origin),
         origin = case_when(
           origin == "Antartic" ~ "Antarctica",
           origin == "Micronesia" ~ "Micronesia, Federated States of",
           origin == "NP" ~ NA_character_,
           TRUE ~ origin
         ),
         origin_std = countrycode(origin, "country.name", "country.name"),
         origin_std2 = case_when(
           origin_std == "South Korea" ~ "Korea South",
           origin_std == "United States" ~ "USA",
           TRUE ~ origin_std
         )) %>%
  left_join(country_coords, by = c("origin_std2" = "country_std"))

geo_valid <- nmds_geo %>% filter(!is.na(lat), !is.na(lon))

geo_dist <- distm(geo_valid[, c("lon", "lat")], fun = distHaversine) / 1000

###############################################
### 9. Align matrices
###############################################
bray_mat <- as.matrix(bray_dist)
common_ids <- intersect(rownames(bray_mat), geo_valid$gcf)

bray_sub <- bray_mat[common_ids, common_ids]
geo_valid_sub <- geo_valid[match(common_ids, geo_valid$gcf), ]

geo_dist <- distm(geo_valid_sub[, c("lon", "lat")], fun = distHaversine) / 1000

###############################################
### 10. Environment distance (Gower)
###############################################
env_dist <- daisy(geo_valid_sub["environment"], metric = "gower")
env_dist <- as.matrix(env_dist)

###############################################
### 11. Mantel + partial Mantel
###############################################
mantel(bray_sub, geo_dist, permutations = 9999)
mantel.partial(bray_sub, geo_dist, env_dist, permutations = 9999)
mantel.partial(bray_sub, env_dist, geo_dist, permutations = 9999)

###############################################
### 12. Variance partitioning
###############################################
geo_df <- as.data.frame(cmdscale(as.dist(geo_dist), k = 2))
env_df <- as.data.frame(cmdscale(as.dist(env_dist), k = 2))

vp <- varpart(bray_sub, geo_df, env_df)
plot(vp, bg = c("skyblue", "lightgreen"), Xnames = c("Geography", "Environment"))
vp$part$indfract

###############################################
### 13. Distance–decay
###############################################
bray_vec <- bray_sub[lower.tri(bray_sub)]
geo_vec  <- geo_dist[lower.tri(geo_dist)]

df_decay <- data.frame(
  geographic_distance_km = geo_vec,
  gcf_dissimilarity = bray_vec
)

ggplot(df_decay, aes(geographic_distance_km, gcf_dissimilarity)) +
  geom_point(alpha = 0.2, color = "steelblue") +
  geom_smooth(method = "loess", color = "darkred") +
  theme_minimal(base_size = 14)

###############################################
### 14. Alpha diversity (environment)
###############################################
alpha_richness <- rowSums(gcf_matrix)

alpha_env <- data.frame(
  gcf = rownames(gcf_matrix),
  richness = alpha_richness
) %>%
  left_join(genomes %>% select(gcf, environment), by = "gcf")

ggplot(alpha_env, aes(environment, richness, fill = environment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

kruskal.test(richness ~ environment, data = alpha_env)

###############################################
### 15. Alpha diversity (country)
###############################################
alpha_country <- data.frame(
  gcf = rownames(gcf_matrix),
  richness = alpha_richness
) %>%
  left_join(genomes %>% select(gcf, origin), by = "gcf")

ggplot(alpha_country, aes(origin, richness, fill = origin)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

kruskal.test(richness ~ origin, data = alpha_country)
