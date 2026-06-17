Common Garden Phyllosphere Project

This repository contains the script for analyzing, and visualizing the fungal phyllosphere data from big sagebrush leaves from a commmon garden in Orchard, ID.

The goal of this study is to investigate what host genetic factors affect the fungal phyllosphere in a common garden setting with co-occurring subspecies. We sampled sagebrush leaves from over 70 plants in both 2012 and 2021. For leaf microbial analyses we use ITS sequencing.

Script:

1. fungal_asv_ocg_analysis.R: We start with assessing differences between our two time points data was collected. We evauluate ASV richness across subspecies-cytotype groups and maximum height for each year. Following this is beta diversity assessment across subspecies-cytotype groups and height using PERMANOVAS and NMDS plots. We use taxa barplots at the genus level for taxa differences between sampling years and subspecies-cytotype groups. Then we used metacoder to evauluate whether there werre any differential taxa between subspecies-cytoype groups. We use procrustes to test for correlation between our fungal community composition and seed-source population location, spatial location in the common garden, and phytochemistry using gas chromatography data from leaves.

Folder:
1. data_csv

   - contains all of the csv files needed for the above script (ASV table, metadata, GC data, height data for 2012 data)
        * height data for 2021 comes from a canopy height model from Olsoy et al work (see manuscript for details)
   
     a. taxonomy
     
       - contains the fasta files pulling from remote BLAST and known taxa from previous artemisia work
       - taxonomy.R: script showing the taxonomy assignment that was done using NCBI BLAST remotely in R
