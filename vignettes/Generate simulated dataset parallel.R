library(umap)
library(ggplot2)
library(phangorn)
library(DCLEAR)
library(TedSim)
library(TreeDist)
# library(slingshot)
library(readxl)
library(Seurat)
library(plyr)
library(DrImpute)
library(Rtsne)
library(SummarizedExperiment)
library(ape)
library(phytools)
library(plyr)



source("/home/krushna/Documents/lineage-method/TedSim/TedSim/R/GE_LB_sim.R")
source("/home/krushna/Documents/lineage-method/TedSim/TedSim/R/GeneExpression_sim.R")
source("/home/krushna/Documents/lineage-method/TedSim/TedSim/R/LineageBarcode_sim.R")
source("/home/krushna/Documents/lineage-method/TedSim/TedSim/R/SimulateCIFs.R")
source("/home/krushna/Documents/lineage-method/TedSim/TedSim/R/Utilities.R")
load("/home/krushna/Documents/lineage-method/TedSim/TedSim/R/sysdata.rda")


# Mutation rate: µ = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4].
# Dropout: pd = 0 or 1.
# We run 10 instances for each combination of the variables.

# Number of cells (ncells)	1024
# Number of genes (ngenes)	500
# Number of Identity Vectors (NIV )	30
# State Identity Vector stepsize (step)	0.5
# Maximum number of state shifts for one division (max walk)	5
# Identity Vector center (starting value for diff-IF)	1
# Number of diff-Identity Vectors (Ndif f )	20
# nondiff-SIV standard deviation (σ)	0.5
# Probability of nonzero gene effect (ge prob)	0.3
# Probability of outlier gene (prob hge)	0.03
# Mean of capture efficiency α (alpha mean)	0.1
# Standard deviation of capture efficiency α (alpha sd)	0.1
# Number of target sites (Nchar)	16

generate_date <- function(mu, p_d, unif_on, seed){
  phyla <- read.tree(text='((t1:4, t2:4, t3:4, t4:4, t5:4, t6:4):2);')
  ncells <- 1024
  ngenes <- 500
  n_cif <- 30
  cif_step <- 0.5
  max_walk <- 5
  n_diff <- 20
  N_char <- 16

  p_a <- 0.8 # defaults


  # mu <- 0.05 #  µ = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4].
  # p_d <- 1 # # Dropout: pd = 0 or 1.
  # unif_on = FALSE

  N_nodes <- 2*ncells-2

  # seed = 1


  set.seed(seed)

  returnlist <- SIFGenerate(phyla,n_diff,step = cif_step)
  cifs <- SimulateCIFs(ncells,phyla,p_a = p_a,n_CIF = n_cif,n_diff = n_diff,step = cif_step,p_d = p_d, Sigma = 0.5, N_char = N_char, max_walk = max_walk, SIF_res = returnlist, unif_on = unif_on)
  #
  #
  #We only need the leaf cells for experiments
  cif_leaves <- lapply(c(1:3),function(parami){
    cif_leaves_all <- cifs[[1]][[parami]][c(1:ncells),]
    return(cif_leaves_all)
  })
  cif_res <- list(cif_leaves,cifs[[2]])
  states <- cifs[[2]]
  states <- states[1:N_nodes,]
  states_leaves <- states[1:ncells,]
  muts <- cifs[[7]]
  rownames(muts) <- paste("cell",states[,4],sep = "_")
  muts_leaves <- muts[1:ncells,]
  muts_nd <- cifs[[8]]
  rownames(muts_nd) <- paste("cell",states[,4],sep = "_")
  muts_leaves_nd <- muts_nd[1:ncells,]

  tree_ct <- cifs[[4]]
  tip_label <- c()
  for (node in tree_ct$tip.label){
    tip <- paste('cell_',substr(node,2,nchar(node)),sep = '')
    tip_label <- c(tip_label,tip)
  }
  tree_ct$tip.label <- tip_label


  #simulate true counts
  true_counts_res <- CIF2Truecounts(ngenes = 500,ncif = n_cif,ge_prob = 0.3,ncells = ncells, cif_res = cif_res)

  data(gene_len_pool)
  gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)

  
  observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]], meta_cell=true_counts_res[[3]], protocol="UMI", alpha_mean=0.2, alpha_sd=0.05, gene_len=gene_len, depth_mean=1e5, depth_sd=3e3)

  # removing compleate mutation if a selected gene is not captured
  # muts <- CaptureDrop(observed_counts[[1]],muts)
  #
  profile_res <- Generate_profile_multi(observed_counts,muts,states_leaves)
  profile_out <- profile_res[[1]]
  top_genes <- profile_res[[2]]

  folder <- paste0("./Datasets/",ncells,"cells-",mu,"mu-", p_d, "pd-", as.integer(unif_on), "unif-", seed,"seed/")
  if (dir.exists(folder)){
    print("Folder already exists")
  } else {
    dir.create(folder)
    gene_expression_dir <- paste0(folder,"counts.csv")
    cell_meta_dir <-  paste0(folder,"cell_meta.csv")
    tree_gt_dir <-  paste0(folder,"tree_gt_bin.nw")
    character_matrix_dir <- paste0(folder,"character_matrix.txt")
    character_matrix_dir_nd <- paste0(folder,"character_matrix_non_drop.txt")
    combined_profile_dir <- paste0("profile_c.txt")
    top_genes_dir <- paste0(folder,"top_genes.txt")

    write.csv(observed_counts[[1]],gene_expression_dir,row.names = FALSE)
    write.csv(states_leaves,cell_meta_dir)
    write.tree(tree_ct,tree_gt_dir)
    write.table(muts_leaves,character_matrix_dir)
    write.table(muts_leaves_nd,character_matrix_dir_nd)
    write.table(profile_out,file = combined_profile_dir,row.names = FALSE,sep = "\t", quote = FALSE)
    write.table(subset(top_genes,select = 1),top_genes_dir,row.names = FALSE,sep = "\t", quote = FALSE)
  }
}


library(parallel)

mu <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
p_d <- c(1,0)
unif_on <- c(TRUE, FALSE)
seed <- c(0,1,12,125,1234,10,20,30,40,50)

# Generate all possible combinations
combinations <- expand.grid(mu, p_d, unif_on, seed)
colnames(combinations) <- c("mu", "p_d", "unif_on", "seed")

mcret <- mclapply(1:nrow(combinations), function(i) {
  generate_date(mu = combinations[i,1], p_d = combinations[i,2],
                unif_on = combinations[i,3], seed = combinations[i,4])
}, mc.cores = 8)

# mcret <- mclapply(1:2, function(i) {
#   generate_date(mu = combinations[i,1], p_d = combinations[i,2],
#                 unif_on = combinations[i,3], seed = combinations[i,4])
# }, mc.cores = 8)
