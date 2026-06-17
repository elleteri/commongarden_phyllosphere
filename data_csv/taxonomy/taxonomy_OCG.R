# new taxonomy for OCG - UNITE database has been done in qiime, but we will need to check our built fasta files for our known taxa 

Sys.setenv(PATH = paste(
  Sys.getenv("PATH"),
  "/Users/ellehorwath/Documents/Artificial_Leaves/artificial_leaf/blast/bin",
  sep = ":"
))

system("blastn -version")

# ============================================================
# STEP 1: LOAD BLAST DATABASES
# system("makeblastdb -in ./taxonomy/Dry_Creek_Phyllosphere_ITS_cultures.fasta -dbtype nucl -out ./taxonomy/Dry_Creek_Phyllosphere_ITS_cultures")
bl_ITS <- blast(db = "./taxonomy/Dry_Creek_Phyllosphere_ITS_cultures", type = "blastn")

# ============================================================
# STEP 2: LOAD REPRESENTATIVE SEQUENCES
seq_ITS <- readDNAStringSet("taxonomy/its-dna-sequences.fasta")
taxMC <- read.csv("taxonomy/taxonomy.csv", head=T, row.names = 1, check.names = F)

# ============================================================
# STEP 3: BLAST TAXONOMY FUNCTION (Chadwick's approach) 
# Count resolved levels
taxMC$n_levels_resolved <- sapply(taxMC$Taxon, function(x) {
  parts <- strsplit(x, ";")[[1]]
  sum(nchar(gsub("^[kpcofgs]__", "", parts)) > 0)
})

# Pull out kingdom-only ASVs
kingdom_only <- taxMC %>%
  filter(n_levels_resolved <= 1)

seq_reblast <- seq_ITS[names(seq_ITS) %in% rownames(kingdom_only)]

# ============================================================
# STEP 6: RUN BLAST FOR BOTH DATASETS
# (these will take a while - run overnight if needed)
local_hits <- c()

for (a in 1:length(seq_reblast@ranges@NAMES)) {
  print(paste("BLASTing ASV", a, "of", length(seq_reblast@ranges@NAMES)))
  
  tax1 <- predict(bl_ITS, seq_reblast[a], BLAST_args = "-perc_identity 90")
  
  # Handle completely empty dataframe AND NA bitscore
  if (nrow(tax1) == 0 || is.na(tax1$bitscore[1])) {
    local_hits[[a]] <- data.frame(
      asv      = names(seq_reblast[a]),
      has_hit  = FALSE,
      sseqid   = NA,
      pident   = NA,
      bitscore = NA
    )
  } else {
    local_hits[[a]] <- data.frame(
      asv      = names(seq_reblast[a]),
      has_hit  = TRUE,
      sseqid   = as.character(tax1$sseqid[1]),
      pident   = tax1$pident[1],
      bitscore = tax1$bitscore[1]
    )
  }
}

local_hits_df <- do.call(rbind, local_hits)
hits_only <- local_hits_df %>% filter(has_hit == TRUE)

# Extract taxonomy string after the | character
hits_only <- hits_only %>%
  mutate(
    tax_string   = gsub(".*\\|", "", sseqid),
    superkingdom = ifelse(grepl("k__", tax_string),
                          gsub(".*k__([^;]+).*", "\\1", tax_string), NA),
    phylum       = ifelse(grepl("p__", tax_string),
                          gsub(".*p__([^;]+).*", "\\1", tax_string), NA),
    class        = ifelse(grepl("c__", tax_string),
                          gsub(".*c__([^;]+).*", "\\1", tax_string), NA),
    order        = ifelse(grepl("o__", tax_string),
                          gsub(".*o__([^;]+).*", "\\1", tax_string), NA),
    family       = ifelse(grepl("f__", tax_string),
                          gsub(".*f__([^;]+).*", "\\1", tax_string), NA),
    genus        = ifelse(grepl("g__", tax_string),
                          gsub(".*g__([^;]+).*", "\\1", tax_string), NA),
    species      = ifelse(grepl("s__", tax_string),
                          gsub(".*s__([^;]+).*", "\\1", tax_string), NA)
  )

# Reconstruct Taxon string in UNITE/QIIME2 format
hits_only <- hits_only %>%
  mutate(Taxon = paste0(
    "k__", ifelse(is.na(superkingdom), "", superkingdom), ";",
    "p__", ifelse(is.na(phylum), "", phylum), ";",
    "c__", ifelse(is.na(class), "", class), ";",
    "o__", ifelse(is.na(order), "", order), ";",
    "f__", ifelse(is.na(family), "", family), ";",
    "g__", ifelse(is.na(genus), "", genus), ";",
    "s__", ifelse(is.na(species), "", species)
  ),
  Confidence = pident  # no confidence score for BLAST hits
  ) %>%
  select(asv, Taxon, Confidence)

tax_updated <- taxMC
tax_updated[hits_only$asv, "Taxon"] <- hits_only$Taxon
tax_updated[hits_only$asv, "Confidence"] <- hits_only$Confidence

taxMC[hits_only$asv[1:5], ]   # before
tax_updated[hits_only$asv[1:5], ]    # after

write.csv(tax_updated, "taxonomy_updated.csv")
