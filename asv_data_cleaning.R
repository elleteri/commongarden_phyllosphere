#Cleaning to just Orchard Common Gardens samples in the metadata and rhe asv data
setwd("/Users/ellehorwath/Documents/Orchard_Common_Garden/commongarden")

#Packages
cran_packages <- c("dplyr", "readr", "tidyr", "tidyverse")

# Function to install and load CRAN packages
install_and_load_cran <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Install and load all CRAN packages
sapply(cran_packages, install_and_load_cran)

#ASV data read in 
asvITS<- read.csv("data_csv/asv-table-dada2-ITS-sagebrush.csv",head=T,row.names=1, check.names = F) #5983 obs of 463 variable
asvITS<- asvITS[,order(colnames(asvITS))] # order samples alphabetically
summary(rowSums(asvITS)) #1.0
summary(colSums(asvITS)) #0

#METADATA read in 
mdITS <- read.csv("data_csv/Sagebrush2021_Mapping_both_4-12-22.csv", head=T, row.names = 1, check.names = F,stringsAsFactors = T) #505 obs of 21 variables.
mdITS <- mdITS[order(row.names(mdITS)),]
mdITS <- subset(mdITS, row.names(mdITS) %in% colnames(asvITS)) #463 of 21 variables
colnames(asvITS) == row.names(mdITS) # sanity check true.

mdITS$Description <- sub("_2020$", "_2021", mdITS$Description) # Replace "2020" with "2021" at the end of the strings

##TAXONOMY read in
tax.ITS <- read.csv("~/Documents/Orchard_Common_Garden/Shared_OCG_Code/data_csv/taxonomy.csv", head=T, row.names = 1, check.names = F) #5983 obs of 2 variables
row.names(asvITS) == row.names(tax.ITS) #TRUE

#Cleaning data#
asvITS <- asvITS[rowSums(asvITS) > 0,]
summary(rowSums(asvITS))
summary(colSums(asvITS)) 
asvITS.t <- t(asvITS) # transpose rows and columns
asvITS.t <- asvITS.t[order(row.names(asvITS.t)),] # order samples alphabetically
asvITS.t <- asvITS.t[,order(colnames(asvITS.t))] # order asvs alphabetically

summary(rowSums(asvITS.t))
summary(colSums(asvITS.t)) 

mdITS <- subset(mdITS, row.names(mdITS) %in% row.names(asvITS.t)) 

asvITS.t2 <- asvITS.t[!(row.names(asvITS.t) %in% c("NEG_9-30-21","AH1919","AHM20207","AHM20125")),] #outliers that need to be removed

# "UTWV.2.7_2012","COVW.2.2_2012","IDV.5.1_2012","CAV.1.7_2012", "UTWV.2.8_2012","WAW.1.4_2012", 'COVW.2.4_2012', 'IDV.2.8_2012'

asvITS.t2 <- asvITS.t2[,colSums(asvITS.t2) > 0]

summary(rowSums(asvITS.t2)) #0
summary(colSums(asvITS.t2)) #1.0

mdITS2 <- subset(mdITS, row.names(mdITS) %in% row.names(asvITS.t2)) 

####OCG subsetting, beta diversity across states, subspecies, ploidy and year ####
asvITS.OCG <- subset(asvITS.t2, mdITS2$Project=="OCG")
asvITS.OCG <- asvITS.OCG[,colSums(asvITS.OCG) > 0]

summary(rowSums(asvITS.OCG))
summary(colSums(asvITS.OCG)) 

asvITS.OCG <- asvITS.OCG[,colSums(asvITS.OCG) > 0]
summary(rowSums(asvITS.OCG))
summary(colSums(asvITS.OCG))

mdITS.OCG <- subset(mdITS, row.names(mdITS) %in% row.names(asvITS.OCG)) ##198 of 21 var

#Write csv for cleaned metadata and asv table
write.csv(asvITS.OCG, file = "data_csv/asvITS.OCG.csv") #206
write.csv(mdITS.OCG, file = "data_csv/metadataITS_OCG.csv") #206 of 21

#Clear Global Environment 
rm(list = ls())

