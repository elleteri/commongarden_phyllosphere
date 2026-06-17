# Orchard Common Garden - Microbial & Chemistry analysis, and visualization 2026
# Install and load necessary packages ####
# library(ANCOMBC)
library(effects)
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(vegan)
library(pairwiseAdonis)
library(lme4)
library(MASS)
library(qiime2R)
library(ggthemes)
library(sf)
library(terra)
library(geosphere)
library(metacoder)

# Read in cleaned data ####
## METADATA
mdITS_OCG <- read.csv("data_csv/metadata_OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T) #238 of 24 vars

mdITS_OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant", "Ecoregion")] <- lapply(mdITS_OCG[, c("Ploidy", "Subspecies", "Subsp_ploidy", "Year", "Plant", "Ecoregion")], as.factor)

mdITS_OCG <- droplevels(mdITS_OCG) 

## ASV DATA
asvITS_OCG <- read.csv("data_csv/asvITS.OCG.csv",head=T, row.names = 1, check.names = F,stringsAsFactors = T)
summary(rowSums(asvITS_OCG)) #0
summary(colSums(asvITS_OCG)) #1

### Remove duplicates from ASV
asvITS_OCG <- asvITS_OCG[!rownames(asvITS_OCG) %in% c('CAT.2.9_2012v1', 'CAV.2.7_2012v2',
                                                      'NVT.2.9_2012v2', 'ORT.2.10_2012v2', 
                                                      'WAT.1.4_2012v2', 'WAT.1.9_2012v2',
                                                      'WAT.2.8_2012v1'), ]

## Remove negative control
asvITS_OCG <- asvITS_OCG[!(row.names(asvITS_OCG) == "NEG_8-28-21"),]

## Remove MTW.3.7.R_2012
asvITS_OCG <- asvITS_OCG[!(row.names(asvITS_OCG) == "MTW.3.7.R_2012"),] 

## Remove consistent outlier
asvITS_OCG <- asvITS_OCG[!(row.names(asvITS_OCG) == "UTW.1.4_2021"),] 

rownames(asvITS_OCG) <- sub("v[12]$", "", rownames(asvITS_OCG)) # remove v1 and v2 from sample names to match metadata, which has already been cleaned to remove v1 and v2

#histogram of asv values
# hist(rowSums(asvITS_OCG), main="Histogram of ASV counts per sample", xlab="ASV counts", breaks=20)

asvITS_OCG <- subset(asvITS_OCG, row.names(asvITS_OCG) %in% row.names(mdITS_OCG)) ## 189
mdITS_OCG <- subset(mdITS_OCG, row.names(mdITS_OCG) %in% row.names(asvITS_OCG)) ## 189

row.names(asvITS_OCG) == row.names(mdITS_OCG) # sanity check:TRUE

mdITS_OCG <- droplevels(mdITS_OCG)

# Filtering
asvITS_OCG <- asvITS_OCG[, colSums(asvITS_OCG) >= 10] # step one: remove asvs that have < 10 seqs across all samples. 
asvITS_OCG <- asvITS_OCG[rowSums(asvITS_OCG) > 0, ] # step two: remove samples that have zero counds after above filtering
asvITS_OCG <- asvITS_OCG[rowSums(asvITS_OCG) >= 500, ] # step three: remove samples that have < 500 seqs after above filtering. 146 samples remain.
asvITS_OCG <- asvITS_OCG[, colSums(asvITS_OCG) > 0] # step four: remove any asvs that are now zero after filtering 146 of 1627 asvs.

mdITS_OCG <- subset(mdITS_OCG, row.names(mdITS_OCG) %in% row.names(asvITS_OCG)) ## 146

row.names(asvITS_OCG) == row.names(mdITS_OCG) # sanity check:TRUE

## TAXONOMY
tax.ITS <- read.csv("data_csv/taxonomy/taxonomy_updated.csv", head=T, row.names = 1, check.names = F) 

## Read in height data for 2012 and 2019 
# read in Bryce's height data 
b_height_data_2012 <- read.csv("data_csv/Orch_vol_cal9-14.csv", row.names = 1) # 468 of 3 vars

# convert row names to a column
b_height_data_2012 <- b_height_data_2012 %>%
  rownames_to_column(var = "Garden Plant ID")

# subset md to 2012 
mdITS_2012_h <- subset(mdITS_OCG, Year == "2012") #97 of 24 vars

mdITS_2012_h$`Garden Plant ID` <- as.character(mdITS_2012_h$`Garden Plant ID`)

# add height data to md for 2012 bu garden plant ID
mdITS_2012_h <- mdITS_2012_h %>%
  left_join(
    b_height_data_2012,
    by = "Garden Plant ID"
  )

mdITS_2012_h$Height4.12 <- as.numeric(as.character(mdITS_2012_h$Height4.12))
rownames(mdITS_2012_h) <- mdITS_2012_h$Description

# # 2015
# canopy_2015 <-st_read("data_csv/shp_files/orchard_digitized_2015_v1.shp")
# # plot(canopy_2015$geometry)
# 
# chm_2015 <- rast("data_csv/shp_files/orchard2015_chm_IDTM_clip2_3cm.tif")
# # plot(chm_2015)
# # plot(st_geometry(canopy_2015), add = TRUE)
# 
# canopy_2015 <- st_set_crs(canopy_2015, st_crs(chm_2015))
# 
# st_crs(canopy_2015) # check the coordinate reference system of the shapefile
# st_crs(chm_2015) # check the coordinate reference system of the raster
# 
# canopy_vect_2015 <- vect(canopy_2015)
# 
# extracted_2015 <- extract(chm_2015, canopy_vect_2015, fun = NULL, ID = TRUE)
# 
# # Map the real plant IDs onto the row index
# extracted_2015$plant_id <- canopy_2015$Plant_ID[extracted_2015$ID]
# str(extracted_2015)
# 
# # remove row with plant id 1 
# extracted_2015 <- extracted_2015[extracted_2015$plant_id != 1, ]
# 
# height_stats_2015 <- extracted_2015 %>%
#   dplyr::group_by(plant_id) %>%
#   dplyr::summarise(
#     mean_height_2015 = mean(orchard2015_chm_IDTM_clip2_3cm, na.rm = TRUE),
#     max_height_2015  = max(orchard2015_chm_IDTM_clip2_3cm, na.rm = TRUE)
#   )

# canopy height data read in from 2019
canopy_2019 <-st_read("data_csv/shp_files/orchard_digitized_2019_v1_wgs84utm.shp")
# Drop Z coordinates from the geometry
canopy_2019 <- st_zm(canopy_2019, drop = TRUE, what = "ZM")
# plot(canopy_2019$geometry)
st_crs(canopy_2019)

canopy_2019 <- st_set_crs(canopy_2019, 32611)  # EPSG for WGS84 UTM Zone 11N

chm_2019 <- rast("data_csv/shp_files/orchard2019_chm_IDTM_clip2_3cm.tif")
# plot(chm_2019)
# plot(st_geometry(canopy_2019), add = TRUE)

canopy_2019 <- st_transform(canopy_2019, st_crs(chm_2019))

st_crs(canopy_2019) == st_crs(chm_2019) # # check the coordinate reference system of the shapefile  against the coordinate reference system of the raster

canopy_vect_2019 <- vect(canopy_2019)

extracted_2019 <- extract(chm_2019, canopy_vect_2019, fun = NULL, ID = TRUE)

# Map the real plant IDs onto the row index
extracted_2019$plant_id <- canopy_2019$Tag[extracted_2019$ID]
str(extracted_2019)

height_stats_2019 <- extracted_2019 %>%
  dplyr::group_by(plant_id) %>%
  dplyr::summarise(
    mean_height_2019 = mean(orchard2019_chm_IDTM_clip2_3cm, na.rm = TRUE),
    max_height_2019  = max(orchard2019_chm_IDTM_clip2_3cm, na.rm = TRUE)
  )

# subset md to 2021
mdITS_2021_h <- subset(mdITS_OCG, Year == "2021") # 

# add height data to mdITS_2021_h
mdITS_2021_h <- mdITS_2021_h %>%
  left_join(height_stats_2019 %>% 
              dplyr::rename(mean_height_2019 = mean_height_2019,
                            max_height_2019 = max_height_2019),
            by = c("Garden Plant ID" = "plant_id")) # 49 of 27 

mdITS_2021_h <- mdITS_2021_h[!is.na(mdITS_2021_h$max_height_2019), ] # 47 
rownames(mdITS_2021_h) <- mdITS_2021_h$Description

# combine md: now it includes the height data 
# mdITS_OCG_h <- mdITS_OCG %>%
#   left_join(height_stats_2015 %>% 
#               dplyr::rename(mean_height_2015 = mean_height_2015,
#                             max_height_2015 = max_height_2015),
#             by = c("Garden Plant ID" = "plant_id")) %>%
#   left_join(height_stats_2019 %>%
#               dplyr::rename(mean_height_2019 = mean_height_2019,
#                             max_height_2019 = max_height_2019),
#             by = c("Garden Plant ID" = "plant_id")) # 146 

# Read in GC chemistry 
OCG_GC_2012 <- read.csv("data_csv/OCG_GC_2012_thresholded.csv", head=T, row.names = 1, check.names = F, stringsAsFactors = F) # 154 of 51 vars
OCG_GC_2021 <- read.csv("data_csv/OCG_GC_2021_thresholded.csv", head=T, row.names = 1, check.names = F, stringsAsFactors = F) # 70 of 45 vars

# Alpha diversity asv level####
## Rarefying
set.seed(12)
asvITS_OCG <- as.matrix(asvITS_OCG)  # convert if needed
asvITS_OCG.r <- suppressWarnings(rrarefy(asvITS_OCG, 500)) ## rarefy: Warning message
asvITS_OCG.r <- asvITS_OCG.r[,colSums(asvITS_OCG.r) > 0] # each sample needs at least x seqs.
summary(colSums(asvITS_OCG.r)) #1
summary(rowSums(asvITS_OCG.r)) #500

asvITS_OCG.richness <- specnumber(asvITS_OCG.r)
mdITS_OCG <- cbind(mdITS_OCG, richness = asvITS_OCG.richness)

hist(asvITS_OCG.richness, main="Histogram of ASV richness", xlab="ASV richness", breaks=20) # check distribution for model family choice
glm.OCG <- glm.nb(richness ~ Year, data=mdITS_OCG)
summary(glm.OCG) 
sum(residuals(glm.OCG, type = "pearson")^2) / df.residual(glm.OCG) # check for overdispersion 

# seperate md for 2012 and 2021 to run separate models for each year
# 2012 
asvITS.2012_h.r <- subset(asvITS_OCG.r, row.names(asvITS_OCG.r) %in% row.names(mdITS_2012_h)) #97 
mdITS_2012_h <- subset(mdITS_2012_h, row.names(mdITS_2012_h) %in% row.names(asvITS_OCG.r)) #97 
asvITS_OCG_2012.richness <- specnumber(asvITS.2012_h.r)
mdITS_2012_h <- cbind(mdITS_2012_h, richness_2012 = asvITS_OCG_2012.richness)
mdITS_2012_h <- droplevels(mdITS_2012_h)

glm.OCG_2012 <- glm.nb(
  richness_2012 ~ Subsp_ploidy/Site + Height4.12,
  data = mdITS_2012_h
)
summary(glm.OCG_2012)

# figure 1A - ASV richness by subspecies ploidy for 2012
ggplot(mdITS_2012_h, aes(x=Subsp_ploidy, y=richness_2012)) +
  geom_violin(aes(fill = Subsp_ploidy)) + 
  geom_jitter() +
  scale_fill_manual(values = c("#784116", "#9CB6A9", "#B0E2FF", "#A3CC51", "#BF7417")) +
  labs(x = "Subspecies:cytotype", y = "ASV richness") +
  theme_classic() +
  theme(legend.position = "none")

# 2021
asvITS.2021_h.r <- subset(asvITS_OCG.r, row.names(asvITS_OCG.r) %in% row.names(mdITS_2021_h)) #47 
asvITS_OCG_2021.richness <- specnumber(asvITS.2021_h.r)
mdITS_2021_h <- cbind(mdITS_2021_h, richness_2021 = asvITS_OCG_2021.richness)

glm.OCG_2021 <- glm.nb(
  richness_2021 ~ Subsp_ploidy/Site + max_height_2019,
  data = mdITS_2021_h
)
summary(glm.OCG_2021)

# supplementary figure - ASV richness by subspecies ploidy for 2021
ggplot(mdITS_2021_h, aes(x=Subsp_ploidy, y=richness_2021)) +
  geom_violin(aes(fill = Subsp_ploidy)) +
  geom_jitter() +
  scale_fill_manual(values = c("#784116", "#9CB6A9", "#A3CC51", "#BF7417")) +
  labs(x = "Subspecies:cytotype", y = "ASV richness") +
  theme_classic() +
  theme(legend.position = "none")

# Betadispersion ####
asvITS_OCG_2012.dispersion <- betadisper(
  vegdist(asvITS.2012_h.r, method = "bray"),
  factor(mdITS_2012_h$Subsp_ploidy, labels = c("T_2n", "T_4n", "V_2n", "V_4n", "W_4n"))
)
anova(asvITS_OCG_2012.dispersion) # sig
TukeyHSD(asvITS_OCG_2012.dispersion)
plot(asvITS_OCG_2012.dispersion)

asvITS_OCG_2021.dispersion <- betadisper(
  vegdist(asvITS.2021_h.r, method = "bray"),
  factor(mdITS_2021_h$Subsp_ploidy, labels = c("T_2n", "T_4n", "V_4n", "W_4n"))
)
anova(asvITS_OCG_2021.dispersion) # sig 
TukeyHSD(asvITS_OCG_2021.dispersion)
plot(asvITS_OCG_2021.dispersion)

#Beta diversity - NMDS plots#### 
set.seed(23)
asvITS_OCG.nmds <- metaMDS(asvITS_OCG.r, trymax = 500) ### Solution reached
# save(asvITS_OCG.nmds, file = "nmds/asvITS_OCG_nmds.rda") #save the nmds so you won't need to run it again
# load("nmds/asvITS_OCG_nmds.rda") #load it to use in code anytime after the initial run

ordiplot(asvITS_OCG.nmds, type = "t",display = "sites",cex = .6)
rownames(asvITS_OCG.nmds$points) == rownames(mdITS_OCG)

### PERMANOVA for year ####
asvITS_OCG.subsp_yr <- adonis2(asvITS_OCG.r ~ mdITS_OCG$Year)
asvITS_OCG.subsp_yr

# ## 2012 asv NMDS ####
# NMDS for height 
set.seed(84)
asvITS.2012_h.nmds <- suppressWarnings(metaMDS(asvITS.2012_h.r, trymax = 500)) ###solution reached! warning message
# save(asvITS.2012.nmds, file = "nmds/asvITS.2012.nmds.rda")
# load("nmds/asvITS.2012.nmds.rda")

ordiplot(asvITS.2012_h.nmds, type = "t",display = "sites",cex = .6)
rownames(asvITS.2012_h.nmds$points) == rownames(mdITS_2012_h)

#SUBSP PLOIDY - Figure 1B
plot(asvITS.2012_h.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",
     col= c("#784116", "#9CB6A9", "#B0E2FF", "#A3CC51", "#BF7417")[mdITS_2012_h$Subsp_ploidy],
     pch=c(17))
legend("topleft",
       legend=c("T_2n","T_4n","V_2n","V_4n","W_4n"),
       col= c("#784116", "#9CB6A9", "#B0E2FF", "#A3CC51", "#BF7417"),
       pch=17,
       cex=0.8,
       bty = "n")
ordispider(asvITS.2012_h.nmds,groups = mdITS_2012_h$Subsp_ploidy, show.groups = "T_2n", col = "#784101")
ordispider(asvITS.2012_h.nmds,groups = mdITS_2012_h$Subsp_ploidy, show.groups = "V_4n", col = "#A3CC51")
ordispider(asvITS.2012_h.nmds,groups = mdITS_2012_h$Subsp_ploidy, show.groups = "W_4n", col = "#BF7417")
ordispider(asvITS.2012_h.nmds,groups = mdITS_2012_h$Subsp_ploidy, show.groups = "T_4n", col = "#9CB6A9")
ordispider(asvITS.2012_h.nmds,groups = mdITS_2012_h$Subsp_ploidy, show.groups = "V_2n", col = "#B0E2FF")

### 2012 PERMANOVA and adonis for subspecies ploidy, and height ####
asvITS.2012.height <- adonis2(
  asvITS.2012_h.r ~ Subsp_ploidy + Height4.12,
  data = mdITS_2012_h,
  permutations = 999,
  strata = mdITS_2012_h$Site,
  by = "margin"
)

#pairwiseadonis
asvITS.2012.subsp.pw <- pairwise.adonis(asvITS.2012_h.r, mdITS_2012_h$Subsp_ploidy, p.adjust.m = 'holm')

# ## 2021 asv NMDS ####
# NMDS for height 
set.seed(93)
asvITS.2021_h.nmds <- suppressWarnings(metaMDS(asvITS.2021_h.r, trymax = 500)) ###solution reached! 
# save(asvITS.2021.nmds, file = "nmds/asvITS.2021.nmds.rda")
# load("nmds/asvITS.2021.nmds.rda")

ordiplot(asvITS.2021_h.nmds, type = "t",display = "sites",cex = .6)
rownames(asvITS.2021_h.nmds$points) == rownames(mdITS_2021_h)

# drop subspecies ploidy levels from 2021 metadata since we lose a subsp_pl level
mdITS_2021_h <- droplevels(mdITS_2021_h)

#SUBSPECIES PLOIDY
plot(asvITS.2021_h.nmds$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2",,
     col= c("#784116","#9CB6A9",'#A3CC51',"#BF7417")[mdITS_2021_h$Subsp_ploidy],
     pch=19)
legend("topright",
       legend=c("T_2n","T_4n","V_4n","W_4n"),
       col= c("#784116","#9CB6A9",'#A3CC51',"#BF7417"),
       pch=19,
       cex=0.8,
       bty = "n")
ordispider(asvITS.2021_h.nmds,groups = mdITS_2021_h$Subsp_ploidy, show.groups = "T_2n", col = "#784116")
ordispider(asvITS.2021_h.nmds,groups = mdITS_2021_h$Subsp_ploidy, show.groups = "T_4n", col = "#9CB6A9")
ordispider(asvITS.2021_h.nmds,groups = mdITS_2021_h$Subsp_ploidy, show.groups = "V_4n", col = "#A3CC51")
ordispider(asvITS.2021_h.nmds,groups = mdITS_2021_h$Subsp_ploidy, show.groups = "W_4n", col = "#BF7417")
#there is no V_2n in 2021

### 2021 PERMANOVA and adonis for subspecies ploidy, and height ####
asvITS.2021.height <- adonis2(asvITS.2021_h.r ~ mdITS_2021_h$Subsp_ploidy  + mdITS_2021_h$max_height_2019, strata = mdITS_2021_h$Site, by = "margin")

#Bar chart of ASV level ####
asvITSt <- t(asvITS_OCG.r)
asvITSto <- asvITSt[order(rowSums(asvITSt),decreasing = T),]
asvITStop <- asvITSto[1:59,]
unknown <- colSums(asvITSto[60:nrow(asvITSto),])
asvITStopo <- rbind(asvITStop,unknown)

barplot(asvITStopo,legend.text=F,axes=F,cex.names = .5,las=2, args.legend = list(x = "topleft", bty = "n", inset=c(-0.15, 0)))
barplot(asvITStopo,legend.text=T,axes=F,cex.names = .3,las=2, args.legend = list(x = "topleft", cex = 0.3, bty = "n", inset=c(-0.11, 0)))

#Bar chart of 2012 and 2021 Genus level ####
tax.ITS.p <- cbind(Feature.ID=rownames(tax.ITS),tax.ITS)
tax.ITS.p <- parse_taxonomy(tax.ITS.p)
tax.ITS.OCG.p <- subset(tax.ITS.p, rownames(tax.ITS.p) %in% rownames(t(asvITS_OCG.r)))

#Function 
prep_asv_barplot <- function(asv_mat, top_frac = 0.5) {
  
  asv_t <- t(asv_mat)
  asv_sums <- rowSums(asv_t)
  
  ord <- order(asv_sums, decreasing = TRUE)
  asv_t <- asv_t[ord, ]
  asv_sums <- asv_sums[ord]
  
  cum_prop <- cumsum(asv_sums) / sum(asv_sums)
  N <- which(cum_prop >= top_frac)[1]
  
  top_asvs <- asv_t[1:N, ]
  other <- colSums(asv_t[(N + 1):nrow(asv_t), ])
  
  rbind(top_asvs, Other = other)
}

asv_bar_2012 <- prep_asv_barplot(asvITS.2012_h.r, top_frac = 0.5) # 49 and 103
asv_bar_2021 <- prep_asv_barplot(asvITS.2021_h.r, top_frac = 0.5) # 24 and 49

reshape_asv_bar <- function(asv_bar, metadata, year_label,
                            id_col, group_col, plant_col) {
  
  as.data.frame(t(asv_bar)) %>%
    rownames_to_column("SampleID") %>%
    left_join(
      metadata %>%
        dplyr::select(all_of(c(id_col, group_col, plant_col))) %>%
        dplyr::rename(
          SampleID = all_of(id_col),
          Subsp_ploidy = all_of(group_col),
          Plant = all_of(plant_col)
        ),
      by = "SampleID"
    ) %>%
    pivot_longer(
      -c(SampleID, Subsp_ploidy, Plant),
      names_to = "ASV",
      values_to = "Abundance"
    ) %>%
    mutate(Year = year_label)
}

# make Description column row names in my metadata sets 
rownames(mdITS_2012_h) <- mdITS_2012_h$Description
rownames(mdITS_2021_h) <- mdITS_2021_h$Description

# Adjust id_col and group_col to match your actual metadata column names
df_2012 <- reshape_asv_bar(asv_bar_2012, mdITS_2012_h, "2012", 
                           id_col = "Description",      # your actual ID column name
                           group_col = "Subsp_ploidy",
                           plant_col = "Plant") # your actual group column name

df_2021 <- reshape_asv_bar(asv_bar_2021, mdITS_2021_h, "2021",
                           id_col = "Description",
                           group_col = "Subsp_ploidy",
                           plant_col = "Plant")

df_all <- bind_rows(df_2012, df_2021)

# See which SampleIDs are coming up as NA
df_all %>% 
  filter(is.na(Subsp_ploidy)) %>% 
  pull(SampleID) %>% 
  unique()

df_all <- df_all %>%
  mutate(Subsp_ploidy = case_when(
    SampleID == "ORTV.2.4_2012" ~ "T_2n",
    TRUE ~ Subsp_ploidy                      # keep everything else as is
  ))

tax_labels <- tax.ITS.OCG.p %>%
  rownames_to_column("ASV") %>%
  mutate(tax_label = case_when(
    !is.na(Genus)  ~ Genus,
    !is.na(Family) ~ paste0(Family, " (family)"),
    !is.na(Order)  ~ paste0(Order, " (order)"),
    !is.na(Class)  ~ paste0(Class, " (class)"),
    TRUE           ~ "Unknown"
  )) %>%
  dplyr::select(ASV, tax_label)

# Rejoin to df_all
df_all <- df_all %>%
  left_join(tax_labels, by = "ASV") %>%
  mutate(tax_label = ifelse(ASV == "Other", "Other", tax_label))

df_all <- df_all %>%
  mutate(tax_label = case_when(
    tax_label == "Unknown" ~ "Other",      # merge Unknown into Other
    ASV == "Other" ~ "Other",
    TRUE ~ tax_label
  ))

# all_asvs <- unique(df_all$ASV)
# n_total <- length(all_asvs)
# asv_colors <- setNames(colorRampPalette(customcol)(n_total), all_asvs)

customcol <- c("gray", "#c969a1", "#859b6c", "#ce4441","#62929a","#ee8577", "#004163", "#eb7926", "#b695bc", "#ffbb44")

all_labels <- unique(df_all$tax_label)
asv_colors <- setNames(colorRampPalette(customcol)(length(all_labels)), all_labels)
asv_colors[grep("Other", names(asv_colors))] <- "grey70"

ggplot(df_all, aes(x = SampleID, y = Abundance, fill = tax_label)) +
  geom_bar(stat = "identity", position = "fill", width = 0.9) +
  scale_fill_manual(values = asv_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_grid(Year ~ Subsp_ploidy,
             scales = "free_x",
             space = "free_x") + 
  labs(y = "Relative Abundance",
       x = NULL,
       fill = "Genus") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
    axis.text.y = element_text(size = 7),
    strip.text.x = element_text(angle = 0, size = 7, face = "bold"),
    strip.text.y = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = NA),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.3, "cm"),
    panel.spacing = unit(0.2, "lines")
  ) +
  guides(fill = guide_legend(ncol = 2))

# Calculate total abundance per tax_label for ordering
df_all <- df_all %>%
  mutate(tax_label = as.character(tax_label))

label_summary <- df_all %>%
  filter(tax_label != "Other") %>%
  group_by(tax_label) %>%
  summarise(total = sum(Abundance)) %>%
  arrange(desc(total))

label_order <- label_summary$tax_label

# Put Other at the top (it will appear at bottom of stacked bar)
label_order <- c("Other", label_order)

# Apply factor ordering
df_all <- df_all %>%
  mutate(tax_label = factor(tax_label, levels = label_order))

# Update colors to match new order
asv_colors <- setNames(colorRampPalette(customcol)(length(label_order)), label_order)
asv_colors["Other"] <- "grey70"

# Replot
barchart_figure <- ggplot(df_all, aes(x = Plant, y = Abundance, fill = tax_label)) +
  geom_bar(stat = "identity", position = "fill", width = 0.9) +
  scale_fill_manual(values = customcol) +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_grid(Year ~ Subsp_ploidy,
             scales = "free_x",
             space = "free_x") +
  labs(y = "Relative Abundance",
       x = NULL,
       fill = "Genus") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
    axis.text.y = element_text(size = 7),
    strip.text.x = element_text(angle = 0, size = 7, face = "bold"),
    strip.text.y = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = NA),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.3, "cm"),
    panel.spacing = unit(0.2, "lines")
  ) +
  guides(fill = guide_legend(ncol = 2))

#METACODER ####
# 2012 SUBSPECIES PLOIDY HEATMAP
# transpose 
# rownames(mdITS.2012)[rownames(mdITS.2012) == "V413_2012"] <- "ORTV.2.4_2012"
asvITS.2012 <- subset(asvITS_OCG, row.names(asvITS_OCG) %in% row.names(mdITS.2012)) #97 
asvITS.2012_t <- t(asvITS.2012)
tax.ITS_2012 <- subset(tax.ITS, rownames(tax.ITS) %in% rownames(asvITS.2012_t))
#order row names
tax.ITS_2012 <- tax.ITS_2012[order(row.names(tax.ITS_2012)), ]
asvITS.2012_t <- asvITS.2012_t[order(row.names(asvITS.2012_t)), ]
row.names(asvITS.2012_t) == row.names(tax.ITS_2012) # sanity check:TRUE

asvITS.2012_t <- as.data.frame(asvITS.2012_t)# return asv table to data frame format after transposing

#Create taxmap object for use with Metacoder functions #
##Create data frame matching taxonomy information with ASV sequences ###
totMC_OCG_2012 <- asvITS.2012_t # new data frame to hold combined asv + taxonomy information
totMC_OCG_2012$Taxon <- tax.ITS_2012$Taxon # append taxonomy information to ASV data frame
totMC_OCG_2012$Taxon <- gsub("k__Eukaryota", "k__Fungi", totMC_OCG_2012$Taxon)
tax.ITS_2012$Taxon <- gsub("k__Eukaryota", "k__Fungi", tax.ITS_2012$Taxon)
row.names(tax.ITS_2012) == row.names(totMC_OCG_2012)# sanity check for ASVs: TRUE
tax.ITS_2012$Taxon == totMC_OCG_2012$Taxon # sanity check for taxonomy: TRUE

totMC_OCG_2012 <- totMC_OCG_2012 %>%
  mutate(Taxon = sapply(Taxon, function(x) {
    parts <- strsplit(x, ";")[[1]]
    # Keep only parts where after __ there is no k__ p__ c__ etc
    # i.e. remove parts where the value after __ starts with another rank prefix
    parts <- parts[grepl("^[kpcofgs]__", parts)]  # must start with a rank prefix
    parts <- parts[!grepl("^[kpcofgs]__[kpcofgs]__", parts)]  # remove ones where value is another rank
    parts <- parts[!grepl("^[kpcofgs]__$", parts)]  # remove empty ranks
    # Also remove parts where value after __ contains k__Fungi (whole string repeated)
    parts <- parts[!grepl("__.*k__", parts)]
    paste(parts, collapse = ";")
  }))

## Create taxmap object #
obj <- parse_tax_data(totMC_OCG_2012,
                      class_cols = "Taxon",                   # name of column that contains input taxon data
                      class_sep = ";",                        # character that separates taxon data
                      class_regex = "^(.+)__(.+)$",           # regex to identify taxon entries
                      class_key = c(tax_rank = "info",        # this is the key that labels each column pulled from the regex, since we had two sections for each identifier we need two columns
                                    tax_name = "taxon_name"))

print(obj) 
print(obj$data$tax_data)
print(obj$data$class_data)
obj$data$class_data <- NULL     # class_data is repetitive/unnecessary
names(obj$data) <- "asv_counts" # rename "data" to something more relevant
print(obj)

#DIFFERENTIAL HEAT TREES FOR OCG
##Adjust taxmap so we can investigate taxa instead of ASVs #
#Now we need to calculate abundances based on taxon not ASV
obj$data$tax_abund <- metacoder::calc_taxon_abund(obj, "asv_counts",
                                                  cols = row.names(mdITS.2012))
print(obj$data$tax_abund)
##**Notice** ###
# Metacoder did not play well with separate naming conventions. As such, all heatmaps are each written into the same places: "diff_table" and "diff_heattree_color"
# To look at different maps, run code starting from the beginning of each section: "Heatmap matrix by _____"

#HEATMAP MATRIX BY 2012 SUBSPECIES PLOIDY 
obj$data$diff_table <- suppressWarnings(
  metacoder::compare_groups(obj, data = "tax_abund",
                            cols = row.names(mdITS.2012),
                            groups = mdITS.2012$Subsp_ploidy))
print(obj$data$diff_table)
obj <- mutate_obs(obj, "diff_table",
                  wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))

#lets look at the p-values and see if there is any significance
range(obj$data$diff_table$wilcox_p_value, finite = TRUE)
# [1] 0.003266539 1.00000000
# the lower range is significant

## Focus only on significant taxa
obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0

print(obj$data$diff_table)

sig_taxa <- obj$data$diff_table %>%
  filter(wilcox_p_value <= 0.05) %>%
  arrange(wilcox_p_value) %>%
  mutate(taxon_name = taxon_names(obj)[taxon_id])  # add readable taxon names

sig_taxa_clean <- sig_taxa %>%
  dplyr::select(taxon_name, treatment_1, treatment_2, 
                log2_median_ratio, wilcox_p_value) %>%
  arrange(wilcox_p_value)
print(sig_taxa_clean)


set.seed(1)
diff_heattree_color <- metacoder::heat_tree_matrix(obj, data = "diff_table",
                                                   node_size = n_obs,
                                                   node_label = taxon_names,
                                                   node_color = log2_median_ratio,
                                                   node_color_range = diverging_palette(),
                                                   node_color_trans = "linear",
                                                   node_color_interval = c(-3, 3),
                                                   edge_color_interval = c(-3, 3),
                                                   node_size_axis_label = "Number of ASVs",
                                                   node_color_axis_label = "Log2 ratio median proportions",
                                                   layout = "davidson-harel",
                                                   initial_layout = "reingold-tilford")

## This plot takes a while to load
print(diff_heattree_color) ## Show taxonomic heat tree
jpeg(filename="metacoder_heatmap_2012_subsp_ploidy.jpeg", width= 823, height= 546, units = "px", quality = 100, res = 500)

# 2021 SUBSPECIES PLOIDY HEATMAP
# transpose
asvITS.2021 <- subset(asvITS_OCG, row.names(asvITS_OCG) %in% row.names(mdITS.2021)) #97
asvITS.2021_t <- t(asvITS.2021)
tax.ITS_2021 <- subset(tax.ITS, rownames(tax.ITS) %in% rownames(asvITS.2021_t))
#order row names
tax.ITS_2021 <- tax.ITS_2021[order(row.names(tax.ITS_2021)), ]
asvITS.2021_t <- asvITS.2021_t[order(row.names(asvITS.2021_t)), ]
row.names(asvITS.2021_t) == row.names(tax.ITS_2021) # sanity check:TRUE

asvITS.2021_t <- as.data.frame(asvITS.2021_t)# return asv table to data frame format after transposing

#Create taxmap object for use with Metacoder functions #
##Create data frame matching taxonomy information with ASV sequences ###
totMC_OCG_2021 <- asvITS.2021_t # new data frame to hold combined asv + taxonomy information
totMC_OCG_2021$Taxon <- tax.ITS_2021$Taxon # append taxonomy information to ASV data frame
totMC_OCG_2021$Taxon <- gsub("k__Eukaryota", "k__Fungi", totMC_OCG_2021$Taxon)
totMC_OCG_2021$Taxon <- gsub("k__Eukaryota", "k__Fungi", totMC_OCG_2021$Taxon)
row.names(totMC_OCG_2021) == row.names(totMC_OCG_2021)# sanity check for ASVs: TRUE
totMC_OCG_2021$Taxon == totMC_OCG_2021$Taxon # sanity check for taxonomy: TRUE

totMC_OCG_2021 <- totMC_OCG_2021 %>%
  mutate(Taxon = sapply(Taxon, function(x) {
    parts <- strsplit(x, ";")[[1]]
    # Keep only parts where after __ there is no k__ p__ c__ etc
    # i.e. remove parts where the value after __ starts with another rank prefix
    parts <- parts[grepl("^[kpcofgs]__", parts)]  # must start with a rank prefix
    parts <- parts[!grepl("^[kpcofgs]__[kpcofgs]__", parts)]  # remove ones where value is another rank
    parts <- parts[!grepl("^[kpcofgs]__$", parts)]  # remove empty ranks
    # Also remove parts where value after __ contains k__Fungi (whole string repeated)
    parts <- parts[!grepl("__.*k__", parts)]
    paste(parts, collapse = ";")
  }))

row.names(tax.ITS_2021) == row.names(totMC_OCG_2021)# sanity check for ASVs: TRUE
tax.ITS_2021$Taxon == totMC_OCG_2021$Taxon # sanity check for taxonomy: TRUE

## Create taxmap object #
obj <- parse_tax_data(totMC_OCG_2021,
                      class_cols = "Taxon",                   # name of column that contains input taxon data
                      class_sep = ";",                        # character that separates taxon data
                      class_regex = "^(.+)__(.+)$",           # regex to identify taxon entries
                      class_key = c(tax_rank = "info",        # this is the key that labels each column pulled from the regex, since we had two sections for each identifier we need two columns
                                    tax_name = "taxon_name"))

print(obj)
print(obj$data$tax_data)
print(obj$data$class_data)
obj$data$class_data <- NULL     # class_data is repetitive/unnecessary
names(obj$data) <- "asv_counts" # rename "data" to something more relevant
print(obj)

#DIFFERENTIAL HEAT TREES FOR OCG
##Adjust taxmap so we can investigate taxa instead of ASVs #
#Now we need to calculate abundances based on taxon not ASV
obj$data$tax_abund <- metacoder::calc_taxon_abund(obj, "asv_counts",
                                                  cols = row.names(mdITS.2021))
print(obj$data$tax_abund)
##**Notice** ###
# Metacoder did not play well with separate naming conventions. As such, all heatmaps are each written into the same places: "diff_table" and "diff_heattree_color"
# To look at different maps, run code starting from the beginning of each section: "Heatmap matrix by _____"

#HEATMAP MATRIX BY SUBSPECIES PLOIDY
obj$data$diff_table <- suppressWarnings(
  metacoder::compare_groups(obj, data = "tax_abund",
                            cols = row.names(mdITS.2021),
                            groups = mdITS.2021$Subsp_ploidy))
print(obj$data$diff_table)
obj <- mutate_obs(obj, "diff_table",
                  wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))

range(obj$data$diff_table$wilcox_p_value, finite = TRUE)

# Procrustes for GC and fungal data####
# 2012 procrustes for GC and fungi
asvITS_OCG_GC_2012.r <- subset(asvITS.2012_h.r, row.names(asvITS.2012_h.r) %in% row.names(OCG_GC_2012)) 
OCG_GC_2012_asv <- subset(OCG_GC_2012, row.names(OCG_GC_2012) %in% row.names(asvITS_OCG_GC_2012.r)) # 96

# order the row names
OCG_GC_2012_asv <- OCG_GC_2012_asv[order(row.names(OCG_GC_2012_asv)), ]
asvITS_OCG_GC_2012.r <- asvITS_OCG_GC_2012.r[order(row.names(asvITS_OCG_GC_2012.r)), ]

row.names(OCG_GC_2012_asv) == row.names(asvITS_OCG_GC_2012.r) # sanity check: TRUE

# make nas zeros
OCG_GC_2012_asv <- scale(OCG_GC_2012_asv)
OCG_GC_2012_asv[is.na(OCG_GC_2012_asv)] <- 0

procrustes_nmds_fungi_2012 <- vegdist(asvITS_OCG_GC_2012.r, method = "bray")
procrustes_gc_2012 <- vegdist(OCG_GC_2012_asv, method = "euclidean")
pro_gc_fun_2012 <- protest(procrustes_nmds_fungi_2012, procrustes_gc_2012, permutations = 999)

# 2021 procrustes for GC 
asvITS_OCG_GC_2021.r <- subset(asvITS.2021_h.r, row.names(asvITS.2021_h.r) %in% row.names(OCG_GC_2021)) 
OCG_GC_2021_asv <- subset(OCG_GC_2021, row.names(OCG_GC_2021) %in% row.names(asvITS_OCG_GC_2021.r)) 

# order the row names
OCG_GC_2021_asv <- OCG_GC_2021_asv[order(row.names(OCG_GC_2021_asv)), ]
asvITS_OCG_GC_2021.r <- asvITS_OCG_GC_2021.r[order(row.names(asvITS_OCG_GC_2021.r)), ]

row.names(OCG_GC_2021_asv) == row.names(asvITS_OCG_GC_2021.r) # sanity check: TRUE

# make nas zeros
OCG_GC_2021_asv <- scale(OCG_GC_2021_asv)
OCG_GC_2021_asv[is.na(OCG_GC_2021_asv)] <- 0

set.seed(95)
procrustes_nmds_fungi_2021 <- vegdist(asvITS_OCG_GC_2021.r, method = "bray")
procrustes_gc_2021 <- vegdist(OCG_GC_2021_asv, method = "euclidean")
pro_gc_fun_2021 <- protest(procrustes_gc_2021, procrustes_nmds_fungi_2021, permutations = 999)

# Procrustes between spatial data in garden and fungal data ####
# 2012
row.names(mdITS_2012_h) == row.names(asvITS.2012_h.r) # TRUE
set.seed(852)
spatial_dist_2012 <- vegdist(mdITS_2012_h[, c("x", "y")], method = "euclidean")
fungi_spatial_dist_2012 <- vegdist(asvITS.2012_h.r, method = "bray")
pro_spatial_fungi_2012 <- protest(spatial_dist_2012, fungi_spatial_dist_2012, permutations = 999) 

# 2021
row.names(mdITS_2021_h) == row.names(asvITS.2021_h.r) # TRUE
set.seed(473)
spatial_dist_2021 <- vegdist(mdITS_2021_h[, c("x", "y")], method = "euclidean")
fungi_spatial_dist_2021 <- vegdist(asvITS.2021_h.r, method = "bray")
pro_spatial_fungi_2021 <- protest(spatial_dist_2021, fungi_spatial_dist_2021, permutations = 999) 

# Procrustes between location and fungal data ####
# 2012
geo_dist_2012 <- distm(mdITS_2012_h[, c("Longitude", "Latitude")], fun = distHaversine)
geo_dist_2012 <- as.dist(geo_dist_2012)

set.seed(47)
fungi_location_dist_2012 <- vegdist(asvITS.2012_h.r, method = "bray")
pro_location_fungi_2012 <- protest(geo_dist_2012, fungi_location_dist_2012, permutations = 999) 

# 2021
geo_dist_2021 <- distm(mdITS_2021_h[, c("Longitude", "Latitude")], fun = distHaversine)
geo_dist_2021 <- as.dist(geo_dist_2021)

set.seed(20)
fungi_location_dist_2021 <- vegdist(asvITS.2021_h.r, method = "bray")
pro_location_fungi_2021 <- protest(geo_dist_2021, fungi_location_dist_2021, permutations = 999) 




