# Run before app.R

# Clear workspace
rm(list = ls())

library(shiny)
library(shinyBS)
library(shinyWidgets)
library(tidyverse)
library(stringr)
library(Matrix)
library(Seurat)
library(metacell)
library(scater)
library(Biostrings) 
library(igraph)
library(plotly)
library(data.tree)
library(collapsibleTree)
library(phylocanvas)
library(msaR)
library(ape)
library(shinycssloaders)
library(shinyalert)
library(shinyhelper)

#setwd("/home/flamanna/LampreyBrainApps/")
source(file = "../../Utils/functions.R")


# Set dir to Metacell
base_dir <- "Data/Cell_trees/"
sample <- "am_brain_liger_k75"
out_path <- str_c(base_dir, sample, "/")
assay = "RNA"

# Load mc_hc and mc_sup
load(file = str_c(out_path, "mc_hc.Rdata"))
load(file =  str_c(out_path, "mc_sup.Rdata"))

# Load Seurat object
in_path <-  "Data/Seurat_objects/"
chromium_seurat <- readRDS(file = str_c(in_path, sample, ".rds"))
Idents(chromium_seurat) <- "cell_types_def"

# Get meta data data frame from Seurat object
meta_data <- chromium_seurat[[]] %>% 
  as_tibble()

# Load and tidy up the log fold change table
#----------------------------------------------------------------------------------------------------------
# Load df table
markers <- read_tsv(str_c("Data/Seurat_markers/", sample, "_cell_types_def.tsv"),
                    col_types = cols(
                      p_val = col_double(),
                      avg_logFC = col_double(),
                      pct.1 = col_double(),
                      pct.2 = col_double(),
                      p_val_adj = col_double(),
                      cluster = col_character(),
                      gene = col_character()
                    )) %>% 
  dplyr::mutate(p_val = round(p_val, digits = 2),
                avg_logFC = round(avg_logFC, digits = 2),
                pct.1 = round(pct.1, digits = 2),
                pct.2 = round(pct.2, digits = 2),
                p_val_adj = round(p_val_adj, digits = 2))

markers_top <- markers %>%
  group_by(cluster) %>% 
  top_n(n = 100,
        wt = avg_logFC)

# Load gene categories
tf <- scan(file = "Data/Lookup_tables/tf_lamprey_complete_list.txt", 
           what = "character")
cf <- scan(file = "Data/Lookup_tables/cf_lamprey_complete_list.txt", 
           what = "character")
nm <- scan(file = "Data/Lookup_tables/nm_lamprey_complete_list.txt",
           what = "character")
nr <- scan(file = "Data/Lookup_tables/nr_lamprey_complete_list.txt",
           what = "character")
nt <- scan(file = "Data/Lookup_tables/nt_lamprey_complete_list.txt",
           what = "character")
np <- scan(file = "Data/Lookup_tables/np_lamprey_complete_list.txt",
           what = "character")
npr <- scan(file = "Data/Lookup_tables/npr_lamprey_complete_list.txt",
            what = "character")

markers_top <- markers_top %>% 
  dplyr::mutate(gene_category = ifelse(gene %in% tf, "transcription_factor", 
                                       ifelse(gene %in% cf, "transcription_cofactor", 
                                              ifelse(gene %in% nm, "neurotransmitter_metabolism", 
                                                     ifelse(gene %in% nr, "neurotransmitter_receptor", 
                                                            ifelse(gene %in% nt, "neurotransmitter_transport", 
                                                                   ifelse(gene %in% np, "neuropeptide", 
                                                                          ifelse(gene %in% npr, "neuropeptide_receptor", NA))))))))

markers_top$gene_category <- factor(markers_top$gene_category)

meta_data_summ <- meta_data %>% 
  group_by(cell_types_def) %>% 
  dplyr::summarise(mean_umis = mean(nCount_RNA),
                   n_cells = n(),
                   family = sample(Level2, 1),
                   description = sample(Description, 1),
                   putative_location = sample(Putative_location, 1)) %>% 
  dplyr::mutate(mean_umis = round(mean_umis, digits = 0))

names(meta_data_summ) <- c("cluster", "mean_umis", "n_cells", "family", "description", "putative_location")

markers_top <- markers_top %>% 
  dplyr::inner_join(meta_data_summ, by = "cluster") %>% 
  dplyr::select(cluster, description, family, putative_location, n_cells, pct.1, pct.2, mean_umis, avg_logFC, p_val_adj, gene, gene_category)

markers_all <- read_tsv(str_c("Data/Seurat_markers/", sample, "_all", "_cell_types_def.tsv"),
                        col_types = cols(
                          p_val = col_double(),
                          avg_logFC = col_double(),
                          pct.1 = col_double(),
                          pct.2 = col_double(),
                          p_val_adj = col_double(),
                          cluster = col_character(),
                          gene = col_character()
                        )) %>% 
  dplyr::mutate(p_val = round(p_val, digits = 2),
                avg_logFC = round(avg_logFC, digits = 2),
                pct.1 = round(pct.1, digits = 2),
                pct.2 = round(pct.2, digits = 2),
                p_val_adj = round(p_val_adj, digits = 2))
#----------------------------------------------------------------------------------------------------------

# Convert hclust to dendrogram 
dend <- as.dendrogram(tree_hc)

# Convert dendrogram to data.tree
dtree <- as.Node(dend)

# Get number of cells per leaf 
ord <- names(dtree$Get(attribute = "names", filterFun = isLeaf))
mc_summ <- meta_data_summ$n_cells
names(mc_summ) <- meta_data_summ$cluster
mc_summ <- mc_summ[base::order(factor(names(mc_summ), levels = ord))]
mc_summ <- as.numeric(mc_summ)

# Get cluster location 
mc_loc <- meta_data_summ$putative_location
names(mc_loc) <- meta_data_summ$cluster
mc_loc <- mc_loc[base::order(factor(names(mc_loc), levels = ord))]
mc_loc <- as.character(mc_loc)
mc_loc[is.na(mc_loc)] <- "undef"

# Get cluster description 
mc_desc <- meta_data_summ$description
names(mc_desc) <- meta_data_summ$cluster
mc_desc <- mc_desc[base::order(factor(names(mc_desc), levels = ord))]
mc_desc <- as.character(mc_desc)
mc_desc[is.na(mc_desc)] <- "undef"

# Set number of cells at leaves
dtree$Set(n_cells = mc_summ, filterFun = isLeaf)

# Add number of cells at each node
dtree$Do(function(node) node$n_cells <- data.tree::Aggregate(node, attribute = "n_cells", aggFun = sum), traversal = "post-order")

# Set cluster location at leaves
dtree$Set(location = mc_loc, filterFun = isLeaf)
dtree$Do(function(node) node$location <- data.tree::Aggregate(node, 
                                                              attribute = "location", 
                                                              aggFun = assign), 
         traversal = "post-order",
         filterFun = isLeaf)

# Set cluster location at nodes
dtree$Set(location = "undef", filterFun = isNotLeaf)

# Set cluster description at leaves
dtree$Set(description = mc_desc, filterFun = isLeaf)
dtree$Do(function(node) node$description <- data.tree::Aggregate(node, 
                                                                 attribute = "description", 
                                                                 aggFun = assign), 
         traversal = "post-order",
         filterFun = isLeaf)

# Set cluster description at nodes
dtree$Set(description = "undef", filterFun = isNotLeaf)

# Setup tooltip
dtree_names <- dtree$Get("name", traversal = "post-order", filterFun = isLeaf)
dtree_descriptions <- dtree$Get("description", traversal = "post-order", filterFun = isLeaf)
dtree_locations <- dtree$Get("location", traversal = "post-order", filterFun = isLeaf)
dtree_n_cells <- dtree$Get("n_cells", traversal = "post-order", filterFun = isLeaf)

tooltip <- str_c("<font size='+1'>",
                 "<b>Cluster:</b> ", dtree_names, 
                 "<br>",
                 "<b>Description:</b> ", dtree_descriptions,
                 "<br>",
                 "<b>Location:</b> ", dtree_locations,
                 "<br>",
                 "<b>N. of cells:</b> ", dtree_n_cells,
                 "<br>",
                 "</font>")

# Set tooltip at leaves
dtree$Set(tooltip = tooltip, filterFun = isLeaf)
dtree$Do(function(node) node$tooltip <- data.tree::Aggregate(node, 
                                                             attribute = "tooltip", 
                                                             aggFun = assign), 
         traversal = "post-order",
         filterFun = isLeaf)

# Set tooltip at nodes
dtree$Do(function(node) node$tooltip <- data.tree::Aggregate(node, 
                                                             attribute = "n_cells", 
                                                             aggFun = sum), 
         traversal = "post-order",
         filterFun = isNotLeaf)

# Set leaves and nodes colors
dtree_cols <- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(length(mc_summ))
names(dtree_cols) <- ord
dtree_cols <- dtree_cols[base::order(factor(names(dtree_cols), levels = ord))]
dtree$Set(colors = dtree_cols, filterFun = isLeaf)
dtree$Set(colors = "white", filterFun = isNotLeaf)

# Prepare data frame for plotting clustering res values
res <- c(str_subset(string = names(meta_data), pattern = "cell_types_def"))
meta_data_res <- chromium_seurat[[]] %>% 
  rownames_to_column(var = "barcode") %>% 
  as_tibble() %>% 
  tidyr::gather(key =  "res", 
                value = "res_value", 
                res)

# Load alignments and (resolved) gene trees
aln_list <- readRDS(file = "Data/Gene_info/aln_list.rds")
aln_list_replace()

tree_list <- readRDS(file = "Data/Gene_info/tree_list.rds")
tree_list_replace()

r_tree_list <- readRDS(file = "Data/Gene_info/r_tree_list.rds")
r_tree_list_replace()

# msa_list <- list(aln_list, alg_list)
phylo_list <- list(r_tree_list, tree_list)

# Load gene to transcript table
gene_trans <- read_tsv(file = "Data/Gene_info/gene_trans.tsv") %>% 
  map_dfr(str_replace, "PMZ_", "PMZ-")

# Load SwissProt results
sprot <- read_tsv("Data/Lookup_tables/sprot.tsv") %>%
  dplyr::mutate(gene_id = str_replace(gene_id, "PMZ_", "PMZ-"),
                transcript_id = str_replace(transcript_id, "PMZ_", "PMZ-")) %>% 
  right_join(dplyr::select(gene_trans, gene_id))

# Load orthologs

# Sea squirt
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ciona <- read_tsv("Data/Orthologues/Orthologues_C_intestinalis/C_intestinalis__v__P_marinus.tsv")
c_ens_gene <- read_tsv("Data/Genomes/sea_squirt/ensembl_id_gene_name.tsv")
ciona_ens <- strsplit(ciona$C_intestinalis, split = ", ")
ciona_ens <- data.frame(Orthogroup = rep(ciona$Orthogroup, sapply(ciona_ens, length)), C_intestinalis = unlist(ciona_ens))
ciona_ens <- inner_join(ciona_ens, c_ens_gene, by = c("C_intestinalis" = "ensembl_transcript_id")) %>%
  dplyr::rename("ensembl_transcript_id" = C_intestinalis) %>% 
  as_tibble()

ciona_ort <- strsplit(ciona$P_marinus, split = ", ")
ciona_ort <- data.frame(Orthogroup = rep(ciona$Orthogroup, sapply(ciona_ort, length)), P_marinus = unlist(ciona_ort)) %>%
  as_tibble() %>% 
  dplyr::mutate(P_marinus = str_replace(P_marinus, "PMZ_", "PMZ-")) %>% 
  dplyr::rename("transcript_id" = P_marinus) %>% 
  inner_join(gene_trans, by = "transcript_id")

ciona_ortho <- inner_join(ciona_ens, ciona_ort)
ciona_ens_url <- str_replace(ciona_ortho$ensembl_gene_id, "^ENS", "http://www.ensembl.org/id/ENS")

ciona_ortho <- ciona_ortho %>% 
  dplyr::mutate(ensembl_gene_id = str_c("<a href='", ciona_ens_url, "' target='_blank'>", ensembl_gene_id, "</a>"))
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Hagfish
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
hagfish <- read_tsv("Data/Orthologues/Orthologues_E_burgeri/E_burgeri__v__P_marinus.tsv")
h_ens_gene <- read_tsv("Data/Genomes/hagfish/ensembl_id_gene_name.tsv")
hagfish_ens <- strsplit(hagfish$E_burgeri, split = ", ")
hagfish_ens <- data.frame(Orthogroup = rep(hagfish$Orthogroup, sapply(hagfish_ens, length)), E_burgeri = unlist(hagfish_ens))
hagfish_ens <- inner_join(hagfish_ens, h_ens_gene, by = c("E_burgeri" = "ensembl_transcript_id")) %>% 
  dplyr::rename("ensembl_transcript_id" = E_burgeri) %>% 
  as_tibble()

hagfish_ort <- strsplit(hagfish$P_marinus, split = ", ")
hagfish_ort <- data.frame(Orthogroup = rep(hagfish$Orthogroup, sapply(hagfish_ort, length)), P_marinus = unlist(hagfish_ort)) %>%
  as_tibble() %>% 
  dplyr::mutate(P_marinus = str_replace(P_marinus, "PMZ_", "PMZ-")) %>% 
  dplyr::rename("transcript_id" = P_marinus) %>% 
  inner_join(gene_trans, by = "transcript_id")

hagfish_ortho <- inner_join(hagfish_ens, hagfish_ort)
hagfish_ens_url <- str_replace(hagfish_ortho$ensembl_gene_id, "^ENS", "http://www.ensembl.org/id/ENS")

hagfish_ortho <- hagfish_ortho %>% 
  dplyr::mutate(ensembl_gene_id = str_c("<a href='", hagfish_ens_url, "' target='_blank'>", ensembl_gene_id, "</a>"))
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Elephant shark
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
shark <- read_tsv("Data/Orthologues/Orthologues_C_milii/C_milii__v__P_marinus.tsv")
h_ens_gene <- read_tsv("Data/Genomes/elephant_shark/ensembl_id_gene_name.tsv")
shark_ens <- strsplit(shark$C_milii, split = ", ")
shark_ens <- data.frame(Orthogroup = rep(shark$Orthogroup, sapply(shark_ens, length)), C_milii = unlist(shark_ens))
shark_ens <- inner_join(shark_ens, h_ens_gene, by = c("C_milii" = "ensembl_transcript_id")) %>% 
  dplyr::rename("ensembl_transcript_id" = C_milii) %>% 
  as_tibble()

shark_ort <- strsplit(shark$P_marinus, split = ", ")
shark_ort <- data.frame(Orthogroup = rep(shark$Orthogroup, sapply(shark_ort, length)), P_marinus = unlist(shark_ort)) %>%
  as_tibble() %>% 
  dplyr::mutate(P_marinus = str_replace(P_marinus, "PMZ_", "PMZ-")) %>% 
  dplyr::rename("transcript_id" = P_marinus) %>% 
  inner_join(gene_trans, by = "transcript_id")

shark_ortho <- inner_join(shark_ens, shark_ort)
shark_ens_url <- str_replace(shark_ortho$ensembl_gene_id, "^ENS", "http://www.ensembl.org/id/ENS")

shark_ortho <- shark_ortho %>% 
  dplyr::mutate(ensembl_gene_id = str_c("<a href='", shark_ens_url, "' target='_blank'>", ensembl_gene_id, "</a>"))
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Spotted gar
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
gar <- read_tsv("Data/Orthologues/Orthologues_L_oculatus/L_oculatus__v__P_marinus.tsv")
h_ens_gene <- read_tsv("Data/Genomes/spotted_gar/ensembl_id_gene_name.tsv")
gar_ens <- strsplit(gar$L_oculatus, split = ", ")
gar_ens <- data.frame(Orthogroup = rep(gar$Orthogroup, sapply(gar_ens, length)), L_oculatus = unlist(gar_ens))
gar_ens <- inner_join(gar_ens, h_ens_gene, by = c("L_oculatus" = "ensembl_transcript_id")) %>% 
  dplyr::rename("ensembl_transcript_id" = L_oculatus) %>% 
  as_tibble()

gar_ort <- strsplit(gar$P_marinus, split = ", ")
gar_ort <- data.frame(Orthogroup = rep(gar$Orthogroup, sapply(gar_ort, length)), P_marinus = unlist(gar_ort)) %>%
  as_tibble() %>% 
  dplyr::mutate(P_marinus = str_replace(P_marinus, "PMZ_", "PMZ-")) %>% 
  dplyr::rename("transcript_id" = P_marinus) %>% 
  inner_join(gene_trans, by = "transcript_id")

gar_ortho <- inner_join(gar_ens, gar_ort)
gar_ens_url <- str_replace(gar_ortho$ensembl_gene_id, "^ENS", "http://www.ensembl.org/id/ENS")

gar_ortho <- gar_ortho %>% 
  dplyr::mutate(ensembl_gene_id = str_c("<a href='", gar_ens_url, "' target='_blank'>", ensembl_gene_id, "</a>"))
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Zebrafish
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
zebrafish <- read_tsv("Data/Orthologues/Orthologues_D_rerio/D_rerio__v__P_marinus.tsv")
h_ens_gene <- read_tsv("Data/Genomes/zebrafish/ensembl_id_gene_name.tsv")
zebrafish_ens <- strsplit(zebrafish$D_rerio, split = ", ")
zebrafish_ens <- data.frame(Orthogroup = rep(zebrafish$Orthogroup, sapply(zebrafish_ens, length)), D_rerio = unlist(zebrafish_ens))
zebrafish_ens <- inner_join(zebrafish_ens, h_ens_gene, by = c("D_rerio" = "ensembl_transcript_id")) %>% 
  dplyr::rename("ensembl_transcript_id" = D_rerio) %>% 
  as_tibble()

zebrafish_ort <- strsplit(zebrafish$P_marinus, split = ", ")
zebrafish_ort <- data.frame(Orthogroup = rep(zebrafish$Orthogroup, sapply(zebrafish_ort, length)), P_marinus = unlist(zebrafish_ort)) %>%
  as_tibble() %>% 
  dplyr::mutate(P_marinus = str_replace(P_marinus, "PMZ_", "PMZ-")) %>% 
  dplyr::rename("transcript_id" = P_marinus) %>% 
  inner_join(gene_trans, by = "transcript_id")

zebrafish_ortho <- inner_join(zebrafish_ens, zebrafish_ort)
zebrafish_ens_url <- str_replace(zebrafish_ortho$ensembl_gene_id, "^ENS", "http://www.ensembl.org/id/ENS")

zebrafish_ortho <- zebrafish_ortho %>% 
  dplyr::mutate(ensembl_gene_id = str_c("<a href='", zebrafish_ens_url, "' target='_blank'>", ensembl_gene_id, "</a>"))
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Coelacanth
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
coelacanth <- read_tsv("Data/Orthologues/Orthologues_L_chalumnae/L_chalumnae__v__P_marinus.tsv")
h_ens_gene <- read_tsv("Data/Genomes/coelacanth/ensembl_id_gene_name.tsv")
coelacanth_ens <- strsplit(coelacanth$L_chalumnae, split = ", ")
coelacanth_ens <- data.frame(Orthogroup = rep(coelacanth$Orthogroup, sapply(coelacanth_ens, length)), L_chalumnae = unlist(coelacanth_ens))
coelacanth_ens <- inner_join(coelacanth_ens, h_ens_gene, by = c("L_chalumnae" = "ensembl_transcript_id")) %>% 
  dplyr::rename("ensembl_transcript_id" = L_chalumnae) %>% 
  as_tibble()

coelacanth_ort <- strsplit(coelacanth$P_marinus, split = ", ")
coelacanth_ort <- data.frame(Orthogroup = rep(coelacanth$Orthogroup, sapply(coelacanth_ort, length)), P_marinus = unlist(coelacanth_ort)) %>%
  as_tibble() %>% 
  dplyr::mutate(P_marinus = str_replace(P_marinus, "PMZ_", "PMZ-")) %>% 
  dplyr::rename("transcript_id" = P_marinus) %>% 
  inner_join(gene_trans, by = "transcript_id")

coelacanth_ortho <- inner_join(coelacanth_ens, coelacanth_ort)
coelacanth_ens_url <- str_replace(coelacanth_ortho$ensembl_gene_id, "^ENS", "http://www.ensembl.org/id/ENS")

coelacanth_ortho <- coelacanth_ortho %>% 
  dplyr::mutate(ensembl_gene_id = str_c("<a href='", coelacanth_ens_url, "' target='_blank'>", ensembl_gene_id, "</a>"))
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Xenopus
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
xenopus <- read_tsv("Data/Orthologues/Orthologues_X_tropicalis/X_tropicalis__v__P_marinus.tsv")
h_ens_gene <- read_tsv("Data/Genomes/xenopus/ensembl_id_gene_name.tsv")
xenopus_ens <- strsplit(xenopus$X_tropicalis, split = ", ")
xenopus_ens <- data.frame(Orthogroup = rep(xenopus$Orthogroup, sapply(xenopus_ens, length)), X_tropicalis = unlist(xenopus_ens))
xenopus_ens <- inner_join(xenopus_ens, h_ens_gene, by = c("X_tropicalis" = "ensembl_transcript_id")) %>% 
  dplyr::rename("ensembl_transcript_id" = X_tropicalis) %>% 
  as_tibble()

xenopus_ort <- strsplit(xenopus$P_marinus, split = ", ")
xenopus_ort <- data.frame(Orthogroup = rep(xenopus$Orthogroup, sapply(xenopus_ort, length)), P_marinus = unlist(xenopus_ort)) %>%
  as_tibble() %>% 
  dplyr::mutate(P_marinus = str_replace(P_marinus, "PMZ_", "PMZ-")) %>% 
  dplyr::rename("transcript_id" = P_marinus) %>% 
  inner_join(gene_trans, by = "transcript_id")

xenopus_ortho <- inner_join(xenopus_ens, xenopus_ort)
xenopus_ens_url <- str_replace(xenopus_ortho$ensembl_gene_id, "^ENS", "http://www.ensembl.org/id/ENS")

xenopus_ortho <- xenopus_ortho %>% 
  dplyr::mutate(ensembl_gene_id = str_c("<a href='", xenopus_ens_url, "' target='_blank'>", ensembl_gene_id, "</a>"))
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Chicken
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
chicken <- read_tsv("Data/Orthologues/Orthologues_G_gallus/G_gallus__v__P_marinus.tsv")
h_ens_gene <- read_tsv("Data/Genomes/chicken/ensembl_id_gene_name.tsv")
chicken_ens <- strsplit(chicken$G_gallus, split = ", ")
chicken_ens <- data.frame(Orthogroup = rep(chicken$Orthogroup, sapply(chicken_ens, length)), G_gallus = unlist(chicken_ens))
chicken_ens <- inner_join(chicken_ens, h_ens_gene, by = c("G_gallus" = "ensembl_transcript_id")) %>% 
  dplyr::rename("ensembl_transcript_id" = G_gallus) %>% 
  as_tibble()

chicken_ort <- strsplit(chicken$P_marinus, split = ", ")
chicken_ort <- data.frame(Orthogroup = rep(chicken$Orthogroup, sapply(chicken_ort, length)), P_marinus = unlist(chicken_ort)) %>%
  as_tibble() %>% 
  dplyr::mutate(P_marinus = str_replace(P_marinus, "PMZ_", "PMZ-")) %>% 
  dplyr::rename("transcript_id" = P_marinus) %>% 
  inner_join(gene_trans, by = "transcript_id")

chicken_ortho <- inner_join(chicken_ens, chicken_ort)
chicken_ens_url <- str_replace(chicken_ortho$ensembl_gene_id, "^ENS", "http://www.ensembl.org/id/ENS")

chicken_ortho <- chicken_ortho %>% 
  dplyr::mutate(ensembl_gene_id = str_c("<a href='", chicken_ens_url, "' target='_blank'>", ensembl_gene_id, "</a>"))
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Mouse
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mouse <- read_tsv("Data/Orthologues/Orthologues_M_musculus/M_musculus__v__P_marinus.tsv")
h_ens_gene <- read_tsv("Data/Genomes/mouse/ensembl_id_gene_name.tsv")
mouse_ens <- strsplit(mouse$M_musculus, split = ", ")
mouse_ens <- data.frame(Orthogroup = rep(mouse$Orthogroup, sapply(mouse_ens, length)), M_musculus = unlist(mouse_ens))
mouse_ens <- inner_join(mouse_ens, h_ens_gene, by = c("M_musculus" = "ensembl_transcript_id")) %>% 
  dplyr::rename("ensembl_transcript_id" = M_musculus) %>% 
  as_tibble()

mouse_ort <- strsplit(mouse$P_marinus, split = ", ")
mouse_ort <- data.frame(Orthogroup = rep(mouse$Orthogroup, sapply(mouse_ort, length)), P_marinus = unlist(mouse_ort)) %>%
  as_tibble() %>% 
  dplyr::mutate(P_marinus = str_replace(P_marinus, "PMZ_", "PMZ-")) %>% 
  dplyr::rename("transcript_id" = P_marinus) %>% 
  inner_join(gene_trans, by = "transcript_id")

mouse_ortho <- inner_join(mouse_ens, mouse_ort)
mouse_ens_url <- str_replace(mouse_ortho$ensembl_gene_id, "^ENS", "http://www.ensembl.org/id/ENS")
mouse_linnarson <- str_c("http://mousebrain.org/genes/", mouse_ortho$external_gene_name, ".html")

mouse_ortho <- mouse_ortho %>% 
  dplyr::mutate(ensembl_gene_id = str_c("<a href='", mouse_ens_url, "' target='_blank'>", ensembl_gene_id, "</a>")) %>% 
  dplyr::mutate(external_gene_name = str_c("<a href='", mouse_linnarson, "' target='_blank'>", external_gene_name, "</a>"))
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Human
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
human <- read_tsv("Data/Orthologues/Orthologues_H_sapiens/H_sapiens__v__P_marinus.tsv")
h_ens_gene <- read_tsv("Data/Genomes/human/ensembl_id_gene_name.tsv")
human_ens <- strsplit(human$H_sapiens, split = ", ")
human_ens <- data.frame(Orthogroup = rep(human$Orthogroup, sapply(human_ens, length)), H_sapiens = unlist(human_ens))
human_ens <- inner_join(human_ens, h_ens_gene, by = c("H_sapiens" = "ensembl_transcript_id")) %>% 
  dplyr::rename("ensembl_transcript_id" = H_sapiens) %>% 
  as_tibble()

human_ort <- strsplit(human$P_marinus, split = ", ")
human_ort <- data.frame(Orthogroup = rep(human$Orthogroup, sapply(human_ort, length)), P_marinus = unlist(human_ort)) %>%
  as_tibble() %>% 
  dplyr::mutate(P_marinus = str_replace(P_marinus, "PMZ_", "PMZ-")) %>% 
  dplyr::rename("transcript_id" = P_marinus) %>% 
  inner_join(gene_trans, by = "transcript_id")

human_ortho <- inner_join(human_ens, human_ort)
human_ens_url <- str_replace(human_ortho$ensembl_gene_id, "^ENS", "http://www.ensembl.org/id/ENS")

human_ortho <- human_ortho %>% 
  dplyr::mutate(ensembl_gene_id = str_c("<a href='", human_ens_url, "' target='_blank'>", ensembl_gene_id, "</a>"))
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
meta_data <- chromium_seurat[[]] %>% 
  dplyr::select(tSNE2d_1,
                tSNE2d_2,
                UMAP2d_cos_1,
                UMAP2d_cos_2)

Idents(chromium_seurat) <- factor(Idents(chromium_seurat), 
                                  levels = ord, 
                                  ordered = TRUE)

chromium_seurat_sce <- as.SingleCellExperiment(chromium_seurat)

# Remove unused objects
rm(ciona,
   c_ens_gene,
   ciona_ens,
   ciona_ort,
   ciona_ens_url,
   hagfish,
   h_ens_gene,
   hagfish_ens,
   hagfish_ort,
   hagfish_ens_url,
   shark,
   shark_ens,
   shark_ort,
   shark_ens_url,
   gar,
   gar_ens,
   gar_ort,
   gar_ens_url,
   zebrafish,
   zebrafish_ens,
   zebrafish_ort,
   zebrafish_ens_url,
   coelacanth,
   coelacanth_ens,
   coelacanth_ort,
   coelacanth_ens_url,
   xenopus,
   xenopus_ens,
   xenopus_ort,
   xenopus_ens_url,
   chicken,
   chicken_ens,
   chicken_ort,
   chicken_ens_url,
   mouse,
   mouse_ens,
   mouse_ort,
   mouse_ens_url,
   mouse_linnarson,
   human,
   human_ens,
   human_ort,
   human_ens_url,
   tree_list,
   r_tree_list,
   tree_hc,
   markers,
   dend,
   dtree_descriptions,
   dtree_locations,
   dtree_n_cells,
   dtree_names,
   mc_desc,
   mc_loc,
   mc_summ,
   nm,
   np,
   npr,
   nr,
   nt,
   ord,
   tf,
   tooltip)

save.image(file = "Apps/am_brain/prepared_data.Rdata")