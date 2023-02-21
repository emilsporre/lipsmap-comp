library(ggbiplot)
library(tidyverse)
library(rlist)
library(MSstats)
library(RColorBrewer)
library(shades)
library(gridExtra)
library(stringr)
library(NormalyzerDE)
library(ggrepel)

#############################################################################################
### define input
#############################################################################################

pep_path <- list.files("/home/user/lipsmap-comp/data/pcc6803/raw_data/20210113/", pattern = ".txt")
pep_path <- lapply(pep_path, function(x){
  x <- paste("/home/user/lipsmap-comp/data/pcc6803/raw_data/20210113/", x, sep = "")
})
uniprot_path <- "/home/user/lipsmap-comp/data/pcc6803/annotations/uniprot_pcc6803_20220401.tab"
layout_path <- c("/home/user/lipsmap-comp/data/pcc6803/layout/layout_pcc6803_20210113.csv")
min_pep <- 3
mode <- "group" ### either "group" or "treatm", decides if the samples are grouped by group or treatment
save_dir <- "/home/user/lipsmap-comp/results/"
organism <- "pcc6803"
#############################################################################################
### Load data
#############################################################################################
layout <- read_csv(layout_path)
colnames(layout) <- c("well", "name")
metabolites <- unique(substr(layout$name, 1, nchar(layout$name)-4))
pep <- lapply(pep_path, read_tsv)

#### count the number of detected peptides in each sample ####
pep_count <- lapply(pep, function(x){
  well <- colnames(x)[4:ncol(x)]
  count <- colSums(x[,4:ncol(x)] != 0)
  df <- data.frame(well, count)
  df <- mutate(df, well = str_extract(well, "[A-Z][0-9]+\\.mzML$") %>% str_remove(".mzML"))
  df <- left_join(df, layout, by = "well")
  return(df)
})

names(pep_count) <- lapply(pep_count, function(x){
  substr(x$name, 1, nchar(x$name)-4)[1]
})

#### perform vsn normalization and transform to long format ####
pep <- lapply(pep, function(x){  
  mtx <- as.matrix(x[,4:ncol(x)])
  mtx[mtx == 0] <- NA
  mtx <- performVSNNormalization(mtx)
  mtx <- 2^mtx
  df <- x
  df[,4:ncol(df)] <- mtx
  df <- pivot_longer(df, cols = ends_with("mzML"), names_to = "Run", values_to = "Intensity")
  return(df)
})
pep <- bind_rows(pep)
uniprot <- read_tsv(uniprot_path)
date <- substring(pep$Run[1], 4,11)

#############################################################################################
### Work up data
#############################################################################################

# Remove decoy peptides in the result and report
nbr_decoy = length(filter(pep, grepl("DECOY", Protein))$Protein)

pep <- pep %>% filter(!(grepl("DECOY", Protein)))

# Rename samples
pep <- mutate(pep, well = str_extract(Run, "[A-Z][0-9]+\\.mzML$") %>% str_remove(".mzML"))
pep <- left_join(pep, layout, by = "well")
pep <- pep %>% select(-Run) %>% rename(Run = name)

# Add variables describing the condition and replicate of each sample/run.
pep <- pep %>% separate(col = Run, into = c("group", "treatm", "BioReplicate"), sep = "_", remove = F) %>% 
  mutate(BioReplicate = as.numeric(BioReplicate))

# Change all bioreplicates to 1 to facilitate group comparison
pep <- pep %>% mutate(BioReplicate = 1)

# Remove numFragments and Protein variable
pep <- pep %>% select(-numFragments) %>% rename(PeptideSequence = Peptide)

# Add condition variable
pep <- pep %>% unite("Condition", group:treatm, sep = "_", remove = F)

# Remove peptides not detected in min_pep replicates
pep_list <- pep %>% group_by(PeptideSequence, Condition) %>% group_split
pep_list <- list.filter(pep_list, length(which(Intensity != 0)) >= (min_pep)) %>% bind_rows()
pep_list <- pep_list %>% group_by(PeptideSequence, group) %>% group_split
pep_list <- list.filter(pep_list, length(Intensity) == 12)
pep_list <- bind_rows(pep_list)

if (mode == "group") {
  pep_list <- pep_list %>% group_by(group) %>% group_split()
} else {
  pep_list <- pep_list %>% group_by(treatm) %>% group_split()
}

# Add peptide IDs
pep_list <- lapply(pep_list, function(x){
  df <-x
  df$ProteinName <- paste("p_", as.numeric(as.factor(df$Protein)), sep = "")
  return(df)
}) 

# Create conversion file to be able to later combine results with uniprot annotation. Also handles ambiguous peptides
conv_pep <- lapply(pep_list, function(x){
  df <- x
  df <- df[,c(2,10)] 
  df <- df %>% rename(Peptide_ID = ProteinName, Entry = Protein)
  df <- unique(df)
  amb <- str_split(df$Entry, ";")
  df$Entry <- lapply(amb, function(x){
    df <- x
    df <- gsub("^.*\\|\\s*(.+?)\\s*\\|.*$", "\\1", df)
    df <- paste(df, collapse = ".;.")
    return(df)
  })
  return(df)
})

# Add variables PrecursorCharge, FragmentIon, ProductCharge and assign value 0 to all entries 
# in the column. Add variable IsotopeLabelType and assign “L” (light) to all entries. Order 
# columns correctly.
pep_list <- lapply(pep_list, function(tbl) {
  tbl %>%
    mutate(PrecursorCharge = 0, FragmentIon = 0, ProductCharge = 0, IsotopeLabelType = "L") %>%
    select(-PeptideSequence:-IsotopeLabelType, -Protein, ProteinName, PeptideSequence, PrecursorCharge,
           FragmentIon, ProductCharge, IsotopeLabelType, Condition, BioReplicate, Run,
           Intensity)
}) 

# Randomize observations
pep_list <- lapply(pep_list, sample_frac)

#############################################################################################
### Process data with MSstats
#############################################################################################
quant_pep <- lapply(pep_list, function(tbl) {
  quant <- dataProcess(
    tbl,
    logTrans = 2,
    normalization = FALSE,
    featureSubset = "all",
    censoredInt = "NA")
})

names(quant_pep) <- lapply(quant_pep, function(x){
  df <-x$FeatureLevelData
  name <- substr(df$GROUP, 1, nchar(as.vector(df$GROUP))-2)[1]
  return(name)
})

#####################################
# Determine the number of comparisons to be done for each group

nbr_exp <- lapply(quant_pep, function(x){
  length(levels(x$ProteinLevelData$GROUP)) -1
})

# Create a comparison matrix for MSstats groupComparison()

comp_mtx <- lapply(nbr_exp, function(x){
  cbind(matrix(rep(-1, x)), diag(x))
})

comp_pep <- vector(mode = "list", length = length(quant_pep))

# Run MSstats statistical comparison, merge output with uniprot annotation and annotate significant peptides
# Save 
for (q in 1:length(quant_pep)) {
  colnames(comp_mtx[[q]]) <- levels(quant_pep[[q]]$FeatureLevelData$GROUP)
  rownames(comp_mtx[[q]]) <- levels(quant_pep[[q]]$FeatureLevelData$GROUP)[-1]
  comp_pep[[q]] <- groupComparison(contrast.matrix = comp_mtx[[q]], data = quant_pep[[q]])
  
  comp_pep[[q]]$ComparisonResult <- comp_pep[[q]]$ComparisonResult %>% rename(Peptide_ID = Protein)
  comp_pep[[q]]$ComparisonResult <- left_join(comp_pep[[q]]$ComparisonResult, conv_pep[[q]])
  comp_pep[[q]]$ComparisonResult <- separate_rows(comp_pep[[q]]$ComparisonResult, Entry, sep = ";")
  comp_pep[[q]]$ComparisonResult$ambig <- ifelse(grepl("\\.", comp_pep[[q]]$ComparisonResult$Entry), "yes", "no")
  comp_pep[[q]]$ComparisonResult$Entry <- gsub("\\.", "", comp_pep[[q]]$ComparisonResult$Entry)
  comp_pep[[q]] <- left_join(comp_pep[[q]]$ComparisonResult, uniprot, by = "Entry")
  
  comp_pep[[q]] <- comp_pep[[q]] %>% mutate(Sign = case_when(adj.pvalue < 0.01 ~ "Sign", TRUE ~ "Unsign")) 
  comp_pep[[q]] <- comp_pep[[q]][!is.infinite(comp_pep[[q]]$log2FC),]
  write_tsv(comp_pep[[q]], paste(save_dir, "comparison_results/pepcomp_", organism, "_", date, "_", substr(comp_pep[[q]]$Label[1],1,nchar(comp_pep[[q]]$Label[1])-2), sep = ""))
  write_tsv(quant_pep[[q]]$FeatureLevelData, paste(save_dir, "processed_data/pepproc_", organism, "_", date, "_", substr(comp_pep[[q]]$Label[1],1,nchar(comp_pep[[q]]$Label[1])-2), sep = ""))
}

names(comp_pep) <- lapply(comp_pep, function(x){
  df <-x
  name <- substr(df$Label, 1, nchar(as.vector(df$Label))-2)[1]
  return(name)
})

######### P-value distribution plot #########
pval_dist <- lapply(metabolites, function(x){
  g <- ggplot(comp_pep[[x]], aes(x = pvalue)) + geom_histogram(bins = 100) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2)) + theme_bw()
  return(g)
})

######### QQ plot #########
qq <- lapply(metabolites, function(x){
  df <- quant_pep[[x]]$FeatureLevelData
  g <- ggplot(df, aes(sample = ABUNDANCE, colour = GROUP)) + stat_qq(size = 0.5) + stat_qq_line() + theme_bw()
  return(g)
})

######### Peptide Count plot #########
peptide_count <- vector(mode = "list", length = length(metabolites))

peptide_count <- lapply(metabolites, function(x){
  df <- pep_count[[x]]
  g <- ggplot(data = df, aes(x = name, y = count)) + geom_bar(stat = "identity", width = 0.8, color = "black") + theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + theme(aspect.ratio = 9/16)
})

############### PCA plot ###########################
pca_plot <- lapply(metabolites, function(x){
  proc <- quant_pep[[x]]$FeatureLevelData
  proc <- proc %>% mutate(INTENSITY = replace(INTENSITY, is.na(INTENSITY), 0))
  proc <- rename(proc, Run = originalRUN, Protein_ID = PROTEIN) %>%  select(-PEPTIDE) %>%
    separate(col = Run,
             into = c("Exp", "Conc", "Repl"),
             sep = "_", remove = F)
  proc$Peptide_ID <- paste("p_", as.numeric(as.factor(proc$FEATURE)), sep = "") 
  proc_wide <- proc %>% select(Peptide_ID, Run, INTENSITY) %>%
    arrange(Run) %>% pivot_wider(names_from = Run, values_from = INTENSITY)
  pep_before <- nrow(proc_wide)
  proc_after <- proc_wide[apply(proc_wide!=0, 1, all),]
  pep_after <- nrow(proc_after)
  cat("", sep = "\n\n")
  print(paste(pep_before-pep_after, " out of ",  pep_before, " peptides (",
              round(100*(pep_before-pep_after)/pep_before, 1),
              "%) were not detected across all MS runs.",
              sep = ""))
  rm(pep_after, pep_before)
  
  pca_matrix <- proc_wide %>% column_to_rownames("Peptide_ID") %>% t
  pca <- prcomp(pca_matrix, center = T, scale = F)
  
  # PCA biplot
  groups <- proc %>% select(Run, Exp, Conc, Repl) %>% distinct %>%
    unite("Condition", Exp:Conc, sep = "_", remove = F) %>% arrange(Condition)
  nbr_colors <- groups$Exp %>% unique %>% length
  nbr_shades <- groups$Conc %>% unique %>% length
  palette <- brewer.pal(nbr_colors,"Set1")[-6]
  shade_high <- saturation(palette, delta(+10))
  
  n = 1
  shades <- list(shade_high)
  while(n < (nbr_shades)) {
    shade_n <- saturation(shades[[length(shades)]], delta(-0.4))
    shades[[n+1]] <- shade_n
    n=n+1
  }
  palette <- shades[[length(shades)]]
  for(i in (length(shades)-1):1) {
    palette <- append(palette, shades[[i]])
  }
  palette <- matrix(palette, nrow = nbr_shades, byrow = T)
  palette <- c(palette)

  g <- ggbiplot(pca, choices = c(1, 2), var.axes = F, labels = rownames(pca_matrix),
                groups = groups$Condition, ellipse = TRUE) + 
    scale_color_manual(values = palette) +
    theme_bw()
  return(g)
})

# Combine and save as a grid of quality control plots
for (q in 1:length(metabolites)) {
  qc_plot <- grid.arrange(pca_plot[[q]], qq[[q]], pval_dist[[q]], peptide_count[[q]], nrow = 2)
  ggsave(filename = paste(save_dir, "/plots/qc_plots/qc_", organism, "_", date, "_", substr(comp_pep[[q]]$Label[1],1,nchar(comp_pep[[q]]$Label[1])-2), ".png", sep = ""),
         plot = qc_plot, units = "mm", width = 400, height = 230, dpi = 600)
}

# Volcano plots with labels
comp_pep <- lapply(comp_pep, function(x){
  df <- x %>% separate(col = Label, into = c("Molecule", "Conc"), sep = "_", remove = F)
  return(df)
  })

g <- lapply(metabolites, function(x){
  conc_labs <- c("1" = "Low Concentration", "2" = "High concentration")
  g <- ggplot(comp_pep[[x]], aes(x = log2FC, y = -log10(adj.pvalue), label = `Gene names  (primary )`)) +
    geom_point(mapping = aes(color = Sign), size = 1, alpha = 0.5) +
    geom_text_repel(data = subset(comp_pep[[x]], adj.pvalue < 0.01),
                    size = 2,
                    color =  "darkcyan",
                    box.padding = 0.5,
                    point.padding = 0.2,
                    force = 10,
                    segment.size = 0.1,
                    segment.color = "black",
                    segment.alpha = 0.3,
                    max.iter = 1000,
                    max.overlaps = 50) +
    scale_color_manual(values = c("red", "black"), name = element_blank()) +
    geom_vline(xintercept = c(-1, 1),
               linetype="dashed", size = 0.2 , alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05),
               linetype="dashed", size = 0.2 , alpha = 0.5) +
    facet_grid(Conc~Molecule, labeller = labeller(Conc = conc_labs)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          strip.background = element_blank()) +
    xlab("Log2(FC)") +
    ylab("-Log10(q-value)")
  
  ggsave(filename = paste(save_dir, "plots/volcano/volcano_", x, "_", organism, "_", date, ".png", sep = ""),
         plot = g, units = "mm", width = 120, height = 200, dpi = 300)
})

