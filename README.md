# alpha-diversity
### 1. Install.packages ####
library(phyloseq)
library(dplyr)
library(ggplot2)
library(vegan)
library(ape)
library(Biostrings)

setwd("D:/Marine_Iguanas_Project/MARINE_IGUANAS/SECOND PAPER 2025/rbcL")
list.files()

### convert csv with sequences into fasta if not done already
library(tidyverse)
library(readr)
list.files()
csv = read_csv("rbcL_296zOTUs_sequences&taxonomy.csv")
writeLines(paste0(">", csv$zOTU, "\n", csv$sequence), "D:/Marine_Iguanas_Project/MARINE_IGUANAS/SECOND PAPER 2025/rbcL/rbcL_zOTU_sequences&taxonomy.fasta")

### 2. Import data ####
otu_mat    <- read.table("OTUtable694x296_rbcL_15082025.txt", 
                         header = TRUE, row.names = 1, sep = "\t", check.names=FALSE)
dim(otu_mat)
meta_df   <- read.table("metadata_rbcL_15082025.txt", 
                        header = TRUE, row.names = 1, sep = "\t", check.names=FALSE)
dim(meta_df)
tax_df      <- read.table("taxonomy_rbcL_296zOTUs_15082025.txt",   # your taxonomy file
                          header = TRUE, row.names = 1, sep = "\t", check.names=FALSE)
dim(tax_df)

seqs <- readDNAStringSet("rbcL_zOTU_sequences.fasta", format="fasta")

### Convert to phyloseq components ####
OTU <- otu_table(as.matrix(otu_mat), taxa_are_rows = TRUE)
SAM <- sample_data(meta_df)
TAX <- tax_table(as.matrix(tax_df))
REFSEQ <- refseq(seqs)

# Create physeq, first without sequences to merge samples
physeq <- phyloseq(OTU, TAX, SAM)
physeq

### 3. Identify and remove low-replicate Locations ####
# Build a small df of sample → Location
sample_info <- data.frame(
  SampleID = sample_names(physeq),
  Location = sample_data(physeq)$Location,
  stringsAsFactors = FALSE
)
# Count samples per Location
loc_counts <- sample_info %>%
  count(Location, name = "n_samples")

# Which Locations have ≥5 samples?
keep_locs <- loc_counts %>%
  filter(n_samples >= 5) %>%
  pull(Location)

# Prune phyloseq to keep only those samples
physeq_filt <- subset_samples(physeq, Location %in% keep_locs)

# (Optional) drop any taxa that now have zero counts
physeq <- prune_taxa(taxa_sums(physeq_filt) > 0, physeq_filt)
physeq

### 4. PLOT TREEEEEEEEEEEEEEE####
###Merge all samples from each island together

phy_island <- merge_samples(physeq, "Island", fun = sum)
sample_data(phy_island) <- data.frame(
  Island = sample_names(phy_island),
  row.names = sample_names(phy_island)
)
phy_island

### Align reference sequences #####

# Incorporate sequences in physeq element
physeq_refseq = merge_phyloseq(phy_island, REFSEQ)
physeq_refseq

# Load DECIPHER
if (!requireNamespace("DECIPHER", quietly=TRUE)) {
  BiocManager::install("DECIPHER")
}
library(DECIPHER)

# Extract sequences and align
seqs <- refseq(physeq_refseq) #changed from: seqs <- refseq(physeq)
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

# Quick check
alignment[1:5]  # print first 5 aligned sequences

### Infer a phylogenetic tree with phangorn ####
# 4a. Load phangorn
if (!requireNamespace("phangorn", quietly=TRUE)) {
  install.packages("phangorn")
}
library(phangorn)

# 4b. Convert the DECIPHER alignment into a phyDat object
phang.align <- phyDat(as.matrix(alignment), type = "DNA")

# 4c. Estimate a substitution model distance matrix (GTR by default)
dm <- dist.ml(phang.align)

# 4d. Build a neighbor-joining (NJ) tree
treeNJ <- NJ(dm)

# 4e. Optionally midpoint‐root the tree
treeNJ <- midpoint(treeNJ)

# 4f. (Optional) perform a quick ML optimization
fit  <- pml(treeNJ, data = phang.align)
fitGTR <- update(fit, k = 4, inv = 0.2)
fitGTR <- optim.pml(fitGTR,
                    model = "GTR",
                    optInv = TRUE,
                    optGamma = TRUE,
                    rearrangement = "stochastic")

# 4g. Extract the optimized tree
final_tree <- fitGTR$tree

# 4h. Merge into your phyloseq object
physeq_tree <- merge_phyloseq(phy_island, phy_tree(final_tree))

# 4i. Quick plot
plot_tree(physeq_tree,
          ladderize  = "left",
          label.tips = "Family",    # or "taxa_names" for OTU IDs
          color = "Island",
          size = "abundance",
          title = "rbcL OTU Phylogeny (NJ + ML optimized)")

# 4j. collapse OTUs at family/genus/sp level
phy_Fam <- tax_glom(phy_island, taxrank = "Family")
phy_Fam <- merge_phyloseq(phy_Fam, phy_tree(phy_tree(physeq_tree)))
plot_tree(phy_Fam,
          nodelabf = nodeplotblank,
          ladderize   = "left",
          label.tips  = "Family",
          color       = "Island",
          size = "abundance",
          title = "rbcL OTU Family tree with samples by island (NJ)")

#or Genus
phy_genus <- tax_glom(phy_island, taxrank = "Genus")
phy_genus <- merge_phyloseq(phy_genus, phy_tree(phy_tree(physeq_tree)))
plot_tree(phy_genus,
          nodelabf = nodeplotblank,
          ladderize   = TRUE,
          label.tips  = "Genus",
          color       = "Island",
          size = "abundance",
          title = "rbcL OTU Genus tree with samples by island (NJ)")

plot_tree(phy_genus,
          nodelabf = nodeplotblank,
          ladderize   = TRUE,
          label.tips  = "Genus",
          size = "abundance",
          title = "rbcL OTU Genus tree with samples by island (NJ)",
          color = "Island") +
  scale_color_manual(name="Island", values = island_colors)

### 5. Visualize Descriptive stats per group####
# Compute total reads & number of samples per sample
df_stats <- data.frame(
  SampleID     = sample_names(physeq),
  TotalReads   = sample_sums(physeq),
  Island       = sample_data(physeq)$Island,
  Location     = sample_data(physeq)$Location,
  Subspecies   = sample_data(physeq)$Marine_iguana_subspecies
)
# Summaries by Island, Location, Subspecies
desc_island <- df_stats %>%
  group_by(Island) %>%
  summarise(
    n_samples   = n(),
    mean_reads  = mean(TotalReads),
    median_reads= median(TotalReads),
    sd_reads    = sd(TotalReads),
    min_reads   = min(TotalReads),
    max_reads   = max(TotalReads)
  )

desc_loc <- df_stats %>%
  group_by(Location) %>%
  summarise(
    n_samples   = n(),
    mean_reads  = mean(TotalReads),
    median_reads= median(TotalReads),
    sd_reads    = sd(TotalReads),
    min_reads   = min(TotalReads),
    max_reads   = max(TotalReads)
  )

desc_subsp <- df_stats %>%
  group_by(Subspecies) %>%
  summarise(
    n_samples   = n(),
    mean_reads  = mean(TotalReads),
    median_reads= median(TotalReads),
    sd_reads    = sd(TotalReads),
    min_reads   = min(TotalReads),
    max_reads   = max(TotalReads)
  )
# View tables
print(desc_island)
print(desc_loc)
print(desc_subsp)

### 6. Alpha-diversity patterns by group ####
# Estimate richness/diversity
diversity_df <- estimate_richness(physeq, 
                                  measures = c("Observed", "Shannon", "Chao1"))
# confirm the three columns are there
head(diversity_df)
diversity_df$SampleID <- rownames(diversity_df)

# Merge with sample metadata
div_stats <- diversity_df %>%
  left_join(df_stats, by = "SampleID")

# your fixed palette for islands
island_colors <- c(
  "PINTA"         = "#CC0033",
  "MARCHENA"      = "#E7298A",
  "GENOVESA"      = "#FFC0CB",
  "ISABELA"       = "#4DAF4A",
  "FERNANDINA"    = "#1B9E77",
  "SANTIAGO"      = "#33FFCC",
  "SANTA_CRUZ"    = "#6B8E23",
  "SANTA_FE"      = "#FFCC33",
  "SAN_CRISTOBAL" = "#0066cc",
  "ESPANOLA"      = "#FF7F00",
  "FLOREANA"      = "#A65628",
  "DARWIN"        = "#999999",
  "WOLF"          = "#984EA3"
)
# your fixed palette for subspecies
subspecies_colors <- c(
  "A_c_sielmanni"     = "#E41A1C",
  "A_c_hayampi"       = "#E7298A",
  "A_c_nanus"         = "#FFC0CB",
  "A_c_cristatus"     = "#1B9E77",
  "A_c_wikelskii"     = "#00FFFF",
  "A_c_hassi"         = "#6B8E23",
  "A_c_trillmichi"    = "#FFFF33",
  "A_c_godzilla"     = "#0066CC",
  "A_c_mertensi"     = "#66CCFF", 
  "A_c_venustissimus" = "#FF7F00",
  "A_c_jeffreysi"    = "#999999"
)

# example plot
p <- plot_richness(physeq, 
                   x ="Island",
                   measures = c("Chao1", "Shannon"),
                   color = "Island") +
  scale_color_manual(name="Island", values = island_colors) +
  theme_light()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
print(p)

#islands+subspecies
p <- plot_richness(physeq, 
                   x ="Island",
                   measures = c("Observed", "Chao1", "Shannon"),
                   color = "Marine_iguana_subspecies") +
  scale_color_manual(name="Marine_iguana_subspecies", values = subspecies_colors) +
  theme_light()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)

#basico
plot_richness(physeq, 
              x ="Island",
              measures = c("Chao1", "Shannon"), 
              color= "Island"+
                theme_minimal())

# visualize these results with ggplot
ggplot(div_stats, aes(x = Island, y = Shannon, fill = Island)) +
  geom_boxplot() +
  theme_bw() +
  labs(
    title = "Shannon Diversity across Islands",
    x = "Island",
    y = "Shannon index"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Boxplots of Shannon by Island, Location, Subspecies
p1 <- ggplot(div_stats, aes(x = Island, y = Shannon)) +
  geom_boxplot() + theme_bw() + 
  labs(title = "Shannon Diversity by Island") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p1)

p2 <- ggplot(div_stats, aes(x = Location, y = Shannon)) +
  geom_boxplot() + theme_bw() + 
  labs(title = "Shannon Diversity by Location") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p2)

# 4. b. Make Island a factor with levels matching your palette names
div_stats$Island <- factor(div_stats$Island, levels = names(island_colors))

# (optional) quick mismatch check:
#setdiff(unique(div_stats$Island), names(island_colors))

p2 <- ggplot(div_stats, aes(x = Location, y = Shannon, fill = Island)) +
geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~ Island, scales = "free_x", nrow = 3) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    legend.position = "none"
  ) +
  labs(title = "Shannon Diversity by Location")

# Con los colores correctos de islas
p2 <- ggplot(div_stats, aes(x = Location, y = Shannon)) +
  geom_boxplot(aes(fill = Island), outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +   # points stay neutral
  facet_wrap(~ Island, scales = "free_x", nrow = 3) +
  scale_fill_manual(values = island_colors, drop = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    legend.position = "none"
  ) +
  labs(title = "Shannon Diversity by Location")
print(p2)

# all boxes the same size
  div_stats$Island <- factor(div_stats$Island, levels = names(island_colors))

p2 <- ggplot(div_stats, aes(x = Location, y = Shannon)) +
  geom_boxplot(aes(fill = Island), outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  facet_grid(~ Island, scales = "free_x", space = "free_x") +  # <- key change
  scale_fill_manual(values = island_colors, guide = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  labs(title = "Shannon Diversity by Location")

p3 <- ggplot(div_stats, aes(x = Subspecies, y = Shannon)) +
  geom_boxplot() + theme_bw() + 
  labs(title = "Shannon Diversity by Subspecies") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print plots
print(p1)
print(p2)
print(p3)

#test if its the same with vegan
H_est <- estimate_richness(physeq, measures = "Shannon")[,1]
mat <- as(otu_table(physeq), "matrix")
if (taxa_are_rows(physeq)) mat <- t(mat)
H_veg <- diversity(mat, index = "shannon")
# should be TRUE
all.equal(unname(H_est), unname(H_veg), tolerance = 1e-12)

library(phyloseq); library(vegan)

er <- estimate_richness(physeq, measures = "Shannon")
H  <- setNames(er$Shannon, rownames(er))          # Shannon per sample
D  <- sample_sums(physeq)[names(H)]               # depth aligned to H

cor.test(H, D, method = "spearman", exact = FALSE)




###aaa
#check sequencing depth is not biasing shannon
depth <- sample_sums(physeq)
cor.test(H_est, depth, method = "spearman")  # ideally weak / non-sig
> library(dplyr); library(ggplot2)
er <- estimate_richness(physeq, "Shannon"); H <- setNames(er$Shannon, rownames(er))
H  <- setNames(er$Shannon, rownames(er))          # Shannon per sample
D  <- sample_sums(physeq)[names(H)]               # depth aligned to H

cor.test(H, D, method = "spearman", exact = FALSE)
###bbb


### 7. Plot ABUNDANCES per taxa per island ####
p <- plot_bar(physeq, x= "Island", 
              fill = "Genus")+
  theme_minimal()

# Remove the black border
p$layers[[1]]$aes_params$colour <- NA
# Print the plot
p

# 7.a Plot RELATIVE ABUNDANCES per taxa per island
# merge your raw counts by Island
physeq_isl <- merge_samples(physeq, "Island")
physeq_isl

physeq_isl_rel <- transform_sample_counts(physeq_isl, function(x) x / sum(x))
round(sample_sums(physeq_isl_rel), 3)

# build the barplot of island‐level relative abundances
p_isl_rel <- plot_bar(physeq_isl_rel,
                      x      = "Sample",    # note: after merge, the sample names ARE the islands
                      fill   = "Genus")+
  theme_minimal()

# remove the black borders again
p_isl_rel$layers[[1]]$aes_params$colour <- NA

# relabel the axes
p_isl_rel +
  xlab("Island") +
  ylab("Relative Abundance")

#6.B.a Plot RELATIVE ABUNDANCES per taxa per location
# merge your raw counts by Island
physeq_loc <- merge_samples(physeq, "Location")
physeq_loc
physeq_loc_rel <- transform_sample_counts(physeq_loc, function(x) x / sum(x))
round(sample_sums(physeq_loc_rel), 3)

# build the barplot of island‐level relative abundances
p_loc_rel <- plot_bar(physeq_loc_rel,
                      x      = "Sample",    # note: after merge, the sample names ARE the islands
                      fill   = "Genus")+
  theme_minimal()

# remove the black borders again
p_loc_rel$layers[[1]]$aes_params$colour <- NA
p_loc_rel

#####**** per loc, per island, same width
## *Collapse zOTUs to Genus
ps_genus <- tax_glom(physeq, taxrank = "Genus", NArm = FALSE)

## ** Merge samples by Location (one bar per Location)
ps_loc <- merge_samples(ps_genus, "Location")

## ***Re-attach Island info to the merged object
loc_to_island <- unique(meta_df[, c("Location","Island")])
rownames(loc_to_island) <- loc_to_island$Location
sample_data(ps_loc) <- sample_data(loc_to_island[sample_names(ps_loc), , drop = FALSE])

## ****Convert to relative abundance per Location
ps_loc_rel <- transform_sample_counts(ps_loc, function(x) if (sum(x) > 0) x/sum(x) else x)

## /Long format for plotting
df_plot <- psmelt(ps_loc_rel)   # has columns: Sample(=Location), Island, Genus, Abundance

## /* Order locations so islands stay together (like your boxplot)
df_plot <- df_plot[order(df_plot$Island, df_plot$Sample), ]
df_plot$Sample <- factor(df_plot$Sample, levels = unique(df_plot$Sample))
# (optional) keep your custom island order if you defined island_colors earlier
if (exists("island_colors")) {
  df_plot$Island <- factor(df_plot$Island, levels = names(island_colors))
}

## /** Plot: stacked bars by Genus, facets per Island
library(ggplot2)
p_genus <- ggplot(df_plot, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", width = 0.6) +
  facet_grid(~ Island, scales = "free_x", space = "free_x") +  # <<< same-width bars, island blocks sized by #locations
  labs(title = "Genus relative abundance by Location",
       x = "Location", y = "Relative abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        panel.spacing = unit(0.8, "lines"))
print(p_genus)

#*****

##### 10. Helper to run a Kruskal–Wallis test for a metric ~ group
kw_test <- function(df, metric, group){
  f <- reformulate(group, response = metric)
  out <- kruskal.test(f, data = df)
  data.frame(
    metric   = metric,
    group_by = group,
    statistic= unname(out$statistic),
    df       = unname(out$parameter),
    p_value  = out$p.value,
    stringsAsFactors = FALSE
  )
}
# Run for both metrics across the three grouping factors
kw_results <- dplyr::bind_rows(
  kw_test(div_stats, "Shannon",  "Island"),
  kw_test(div_stats, "Observed", "Island"),
  kw_test(div_stats, "Shannon",  "Location"),
  kw_test(div_stats, "Observed", "Location"),
  kw_test(div_stats, "Shannon",  "Subspecies"),
  kw_test(div_stats, "Observed", "Subspecies")
)
print(kw_results)

# ── Pairwise Wilcoxon (BH) helper ─────────────────────────────────────────────
pairwise_wilcox_df <- function(df, metric, group){
  d <- df %>%
    dplyr::filter(!is.na(.data[[metric]]), !is.na(.data[[group]]))
  
  pw <- pairwise.wilcox.test(
    x      = d[[metric]],
    g      = d[[group]],
    p.adjust.method = "BH",
    exact  = FALSE
  )
  # Convert upper-tri p-value matrix to tidy data frame
  m <- pw$p.value
  out <- as.data.frame(as.table(m), stringsAsFactors = FALSE)
  colnames(out) <- c("group1","group2","p_adj")
  out <- out %>%
    dplyr::filter(!is.na(p_adj)) %>%
    dplyr::mutate(
      metric   = metric,
      group_by = group,
      sig = dplyr::case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ ""
      )
    ) %>%
    dplyr::arrange(p_adj)
  out
}

# Run for both metrics across Island, Location, Subspecies
pw_results <- dplyr::bind_rows(
  pairwise_wilcox_df(div_stats, "Shannon",  "Island"),
  pairwise_wilcox_df(div_stats, "Observed", "Island"),
  pairwise_wilcox_df(div_stats, "Shannon",  "Location"),
  pairwise_wilcox_df(div_stats, "Observed", "Location"),
  pairwise_wilcox_df(div_stats, "Shannon",  "Subspecies"),
  pairwise_wilcox_df(div_stats, "Observed", "Subspecies")
)

# View the top contrasts (smallest adjusted p-values first)
head(pw_results, 20)
print(pw_results, 40)
