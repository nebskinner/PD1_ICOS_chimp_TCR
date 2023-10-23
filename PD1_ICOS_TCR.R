##### Analysis of TCR sequences extracted from bulk RNA seq data

#### Load libraries
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(stringi)
library(ggpubr)
library(ggtext)
library(rstatix)
library(gridExtra)
library(vegan)

library(ff)
library(RecordLinkage)
library(RcppAlgos)
library(ggseqlogo)

#### Functions for use later in the analysis
### functions for formatting
make_italics <- function(x) {
  as.expression(lapply(x, function(y) bquote(italic(.(y)))))
}

make_bold <- function(x) {
  as.expression(lapply(x, function(y) bquote(bold(.(y)))))
}

### function to find amino acid usage at each position of a sequence:
AA_usage <- function (AA_vec, AA_length) {
  AA_use_df <- data.frame()
  for(i in 1:AA_length){
    AAs <- paste0(substring(AA_vec,i,i), collapse="")
    A_count <- str_count(AAs, "A")
    R_count <- str_count(AAs, "R")
    N_count <- str_count(AAs, "N")
    D_count <- str_count(AAs, "D")
    C_count <- str_count(AAs, "C")
    E_count <- str_count(AAs, "E")
    Q_count <- str_count(AAs, "Q")
    G_count <- str_count(AAs, "G")
    H_count <- str_count(AAs, "H")
    I_count <- str_count(AAs, "I")
    L_count <- str_count(AAs, "L")
    K_count <- str_count(AAs, "K")
    M_count <- str_count(AAs, "M")
    F_count <- str_count(AAs, "F")
    P_count <- str_count(AAs, "P")
    S_count <- str_count(AAs, "S")
    T_count <- str_count(AAs, "T")
    W_count <- str_count(AAs, "W")
    Y_count <- str_count(AAs, "Y")
    V_count <- str_count(AAs, "V")
    use_freq<- c(A_count, R_count, N_count, D_count, C_count, E_count, Q_count, G_count, H_count, I_count, L_count, K_count, M_count, F_count, P_count, S_count, T_count, W_count, Y_count, V_count)
    AA_use_df <- rbind.data.frame(AA_use_df, use_freq)
  }
  
  AA_freq <- t(AA_use_df)
  rownames(AA_freq) <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  colnames(AA_freq) <- paste0("Pos", 1:AA_length)
  
  return(AA_freq)
}


#### Load data
### subject descriptives
chimp_attr <- data.frame(samp = c("s001", "s002", "s003", "s004", "s006", "s007", "s008", "s009", "s010", "s011", "s012", "s013", "s014", "s015"),
                   subj = c(rep("4x0395", 3), rep("4x0339", 2), rep("4x0526", 3), rep("4x0405", 3), rep("4x0312", 3)),
                   week = c(rep("wk11", 3), rep("wk8", 11)),
                   PD1_ICOS = c("high", "intermediate", "negative", "high", "negative", "high", "intermediate", "negative", "high", "intermediate", "negative", "high", "intermediate", "negative"),
                   samp2 = c("S47", "S48", "S49", "S50", "S52", "S53", "S54", "S55", "S56", "S57", "S58", "S59", "S60", "S61"))
### TRB data 
chimp_TRB_pre <- data.frame()
for(i in 1:nrow(chimp_attr)){
  clonotype_file <- paste("/Users/nes002/Desktop/chimp_TCR/Bulk_TCR_seq/MIXCR/TRB/p21002-", chimp_attr$samp[i], chimp_attr$subj[i], chimp_attr$week[i], "PD1-ICOS", chimp_attr$PD1_ICOS[i], chimp_attr$samp2[i], "clones_TRB.tsv", sep = "_")
  clonotype_file <- sub("-_", "-", clonotype_file)
  clonotype_file <- sub("_clones", ".clones", clonotype_file) 
  
  df <- read.delim(clonotype_file,  header=TRUE, sep="\t") %>%
    mutate(subj = chimp_attr$subj[i]) %>%
    mutate(week = chimp_attr$week[i]) %>%
    mutate(PD1_ICOS = chimp_attr$PD1_ICOS[i])
  
  chimp_TRB_pre <- rbind.data.frame(chimp_TRB_pre, df, stringsAsFactors = FALSE)
}

chimp_TRB <- chimp_TRB_pre %>%
  select(clone_ID = cloneId, read_count = readCount, read_fraction = readFraction, sequence = targetSequences, v_gene = allVHitsWithScore, j_gene = allJHitsWithScore, c_gene = allCHitsWithScore, cdr3_aa = aaSeqCDR3, cdr3_nt = nSeqCDR3, subj, week, PD1_ICOS) %>%
  mutate(cdr3_len = nchar(cdr3_aa)) 

# formatting V, J, and C genes
chimp_TRB$v_gene <- gsub("\\,.*", "", chimp_TRB$v_gene)
chimp_TRB$v_gene <- gsub("\\*.*", "", chimp_TRB$v_gene)

chimp_TRB$j_gene <- gsub("\\,.*", "", chimp_TRB$j_gene)
chimp_TRB$j_gene <- gsub("\\*.*", "", chimp_TRB$j_gene)

chimp_TRB$c_gene <- gsub("\\,.*", "", chimp_TRB$c_gene)
chimp_TRB$c_gene <- gsub("\\*.*", "", chimp_TRB$c_gene)


### TRA data 
chimp_TRA_pre <- data.frame()
for(i in 1:nrow(chimp_attr)){
  clonotype_file <- paste("/Users/nes002/Desktop/chimp_TCR/Bulk_TCR_seq/MIXCR/TRA/p21002-", chimp_attr$samp[i], chimp_attr$subj[i], chimp_attr$week[i], "PD1-ICOS", chimp_attr$PD1_ICOS[i], chimp_attr$samp2[i], "clones_TRAD.tsv", sep = "_")
  clonotype_file <- sub("-_", "-", clonotype_file)
  clonotype_file <- sub("_clones", ".clones", clonotype_file) 
  
  df <- read.delim(clonotype_file,  header=TRUE, sep="\t") %>%
    mutate(subj = chimp_attr$subj[i]) %>%
    mutate(week = chimp_attr$week[i]) %>%
    mutate(PD1_ICOS = chimp_attr$PD1_ICOS[i])
  
  chimp_TRA_pre <- rbind.data.frame(chimp_TRA_pre, df, stringsAsFactors = FALSE)
}

chimp_TRA <- chimp_TRA_pre %>%
  select(clone_ID = cloneId, read_count = readCount, read_fraction = readFraction, sequence = targetSequences, v_gene = allVHitsWithScore, j_gene = allJHitsWithScore, c_gene = allCHitsWithScore, cdr3_aa = aaSeqCDR3, cdr3_nt = nSeqCDR3, subj, week, PD1_ICOS) %>%
  filter(grepl("TRA", v_gene)) %>%
  filter(grepl("TRA", j_gene)) %>%
  mutate(cdr3_len = nchar(cdr3_aa))

# formatting V, J, and C genes
chimp_TRA$v_gene <- gsub("\\,.*", "", chimp_TRA$v_gene)
chimp_TRA$v_gene <- gsub("\\*.*", "", chimp_TRA$v_gene)

chimp_TRA$j_gene <- gsub("\\,.*", "", chimp_TRA$j_gene)
chimp_TRA$j_gene <- gsub("\\*.*", "", chimp_TRA$j_gene)

chimp_TRA$c_gene <- gsub("\\,.*", "", chimp_TRA$c_gene)
chimp_TRA$c_gene <- gsub("\\*.*", "", chimp_TRA$c_gene)


### quantify number of reads per subject/timepoint
num_reads_TRB <- chimp_TRB %>%
  group_by(subj, week, PD1_ICOS) %>%
  mutate(num_reads = (sum(read_count))) %>%
  distinct(subj,subj, week, PD1_ICOS, num_reads)

num_reads_TRA <- chimp_TRA %>%
  group_by(subj, week, PD1_ICOS) %>%
  mutate(num_reads = (sum(read_count))) %>%
  distinct(subj,subj, week, PD1_ICOS, num_reads)


#### IGHV gene usage
### TRB
chimp_TRBVJ_pre <- chimp_TRB %>%
  group_by(v_gene, PD1_ICOS) %>%
  mutate(v_count = n()) %>%
  group_by(PD1_ICOS) %>%
  mutate(v_total = n()) %>%
  mutate(v_prop = v_count/v_total) %>%
  mutate(vj_gene = paste(v_gene, j_gene)) %>%
  group_by(vj_gene, PD1_ICOS) %>%
  mutate(vj_count = n()) %>%
  mutate(vj_prop = vj_count/v_total) %>%
  mutate(v_gene = gsub("TRBV", "", v_gene)) %>%
  ungroup()

chimp_TRBV <- chimp_TRBVJ_pre %>%
  distinct(v_gene, PD1_ICOS, v_count, v_total, v_prop) %>%
  complete(v_gene, PD1_ICOS, fill = list(v_count = 0, v_prop = 0)) %>%
  group_by(PD1_ICOS) %>%
  mutate(v_total = sum(v_count))

### Statistics
## list of all TRBV genes
listTRBV <- unique(chimp_TRBV$v_gene)

## fisher testing
chimp_TRBV_stat <- data.frame()
for(i in 1:length(listTRBV)){
  VG <- filter(chimp_TRBV, v_gene == listTRBV[i])
  VG_neg <- VG %>%
    filter(PD1_ICOS == "negative") 
  VG_int <- VG %>%
    filter(PD1_ICOS == "intermediate")
  VG_high <- VG %>%
    filter(PD1_ICOS == "high")
  
  high_neg <- fisher.test(matrix(c(VG_high$v_count, VG_neg$v_count, VG_high$v_total-VG_high$v_count, VG_neg$v_total-VG_neg$v_count), ncol=2))$p.value
  
  high_int <- fisher.test(matrix(c(VG_high$v_count, VG_int$v_count, VG_high$v_total-VG_high$v_count, VG_int$v_total-VG_int$v_count), ncol=2))$p.value
  
  int_neg <- fisher.test(matrix(c(VG_int$v_count, VG_neg$v_count, VG_int$v_total-VG_int$v_count, VG_neg$v_total-VG_neg$v_count), ncol=2))$p.value
  
  VG_df <- cbind.data.frame(v_gene = listTRBV[i], high_neg, high_int, int_neg)
  chimp_TRBV_stat <- rbind.data.frame(chimp_TRBV_stat, VG_df, stringsAsFactors = FALSE)
}

## correcting P-values for multiple comparisons
chimp_TRBV_stat_corr <-data.frame(chimp_TRBV_stat$v_gene, matrix(p.adjust(as.vector(as.matrix(chimp_TRBV_stat[,-1])), method="BH"),ncol=3))
colnames(chimp_TRBV_stat_corr) <- colnames(chimp_TRBV_stat)

## format statistical comparisons and add p signficance asterisks
chimp_TRBV_sig <- chimp_TRBV_stat_corr %>%
  filter(high_neg < 0.05 | high_int < 0.05 | int_neg < 0.05) %>%
  gather(key = group, value = p.adj, -v_gene) %>%
  mutate(group1 = gsub("_.*", "", group)) %>%
  mutate(group2 = gsub(".*_", "", group)) %>%
  mutate(group1 = gsub("neg", "negative", group1)) %>%
  mutate(group1 = gsub("int", "intermediate", group1)) %>%
  mutate(group2 = gsub("neg", "negative", group2)) %>%
  mutate(group2 = gsub("int", "intermediate", group2)) %>%
  mutate(.y. = "v_prop") %>%
  mutate(p.adj.signif = ifelse(p.adj < 0.0001, "****",
                           ifelse(p.adj < 0.001, "***",
                                  ifelse(p.adj < 0.01, "**",
                                         ifelse(p.adj < 0.05, "*", "ns"))))) %>%
  mutate(numbering = sub("[0-9]*-", "", v_gene)) %>%
  mutate(fam = sub("-.*", "", v_gene)) %>%
  group_by(fam) %>%
  arrange(as.numeric(fam), as.numeric(numbering)) %>%
  ungroup() %>%
  select(v_gene, group1, group2, p.adj, p.adj.signif)
  
chimp_TRBV_sig$v_gene <- factor(chimp_TRBV_sig$v_gene, levels = unique(chimp_TRBV_sig$v_gene))


### Plots
## heatmaps
# format data
chimp_TRBV_hm_pre <-chimp_TRBV[,c(1:2,5)] %>%
  spread(key = v_gene, value = v_prop)

chimp_TRBV_hm_pre$PD1_ICOS <- factor(chimp_TRBV_hm_pre$PD1_ICOS, levels=c("high", "intermediate", "negative"))

chimp_TRBV_hm <- chimp_TRBV_hm_pre[order(chimp_TRBV_hm_pre$PD1_ICOS),]

chimp_TRBV_hm <- chimp_TRBV_hm[,-1]

# reorder TRBV genes to be in numerical order rather than alphanumeric (i.e., 1, 2, 3... rather than 1, 10, 11)
reorder_chimp_TRBV <- data.frame(vname = colnames(chimp_TRBV_hm)) %>%
  mutate(numbering = sub("[0-9]*-", "", vname)) %>%
  mutate(fam = sub("-.*", "", vname)) %>%
  group_by(fam) %>%
  arrange(as.numeric(fam), as.numeric(numbering))

chimp_TRBV_hm_arr <- t(chimp_TRBV_hm)[reorder_chimp_TRBV$vname,] 

# define heatmap color palette
my_palette <-colorRampPalette(c("#f8fafc", "#ebf1f5", "#c4d4e2", "#b1c6d9", "#9db8cf", "#89aac5", "#769cbc", "#628db2", "#4f7fa9", "#2E5984"))(100)

# plot
pheatmap(chimp_TRBV_hm_arr, scale="none", cluster_rows = FALSE, cluster_cols=FALSE, border_color = NA, legend = TRUE, show_rownames = TRUE, show_colnames = TRUE, cellheight = 13, cellwidth = 25, fontsize = 13, fontsize_row = 14, fontsize_col = 12, angle_col = "45", labels_row = make_italics(rownames(chimp_TRBV_hm_arr)), labels_col = make_bold(c("High", "Int", "Low")), main = "", color= my_palette)

### barplot (only showing genes with statistically significant differences)
## format data for plotting
chimp_TRBV_subset <- chimp_TRBVJ_pre %>%
  group_by(subj, v_gene, PD1_ICOS) %>%
  mutate(v_count_subj = n()) %>%
  group_by(subj, PD1_ICOS) %>%
  mutate(v_total_subj = n()) %>%
  mutate(v_prop_subj = v_count_subj/v_total_subj) %>%
  filter(v_gene %in% chimp_TRBV_sig$v_gene) %>%
  distinct(v_gene, subj, PD1_ICOS, v_prop, v_prop_subj) %>%
  mutate(v_prop_subj_round = round(v_prop_subj, digits =3)) %>%
  group_by(v_gene, PD1_ICOS, v_prop_subj_round) %>%
  mutate(overlap = n()) %>%
  mutate(overlap = ifelse(overlap > 1, "yes", "no")) %>%
  mutate(numbering = sub("[0-9]*-", "", v_gene)) %>%
  mutate(fam = sub("-.*", "", v_gene)) %>%
  group_by(fam) %>%
  arrange(as.numeric(fam), as.numeric(numbering)) %>%
  ungroup() %>%
  select(-fam, -numbering)

chimp_TRBV_subset$v_gene <- factor(chimp_TRBV_subset$v_gene, levels = unique(chimp_TRBV_subset$v_gene))

## plot
chimp_TRBV_sig_p <- ggplot() +
  geom_bar(data = chimp_TRBV_subset, aes(x = PD1_ICOS, y = v_prop, fill = PD1_ICOS), position="dodge", stat="identity", size=1, alpha = 0.8) +
  geom_point(data = chimp_TRBV_subset[chimp_TRBV_subset$overlap == "no",], aes(x = PD1_ICOS, y = v_prop_subj), color = "black", size=1.5, shape = 19) +
  geom_point(data = chimp_TRBV_subset[chimp_TRBV_subset$overlap == "yes",], aes(x = PD1_ICOS, y = v_prop_subj), color = "black", size=1.5, shape = 19, position = position_dodge(width = 0.9)) +
  facet_wrap(. ~ v_gene, nrow = 2, strip.position = "top") +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#999999"), labels = c("High", "Intermediate", "Negative")) +
  #scale_color_manual(values = c("#225EA8", "#1D91C0", "#41B6C4", "#7FCDBB", "#C7E9B4"))+
  scale_y_continuous(breaks = seq(from = 0, to = 0.15, by = .05), limits = c(0,0.153), expand = c(0,0))+
  guides(fill = guide_legend(order =1),
         color = guide_legend(override.aes = list(size=2.45, order =2)))+
  stat_pvalue_manual(chimp_TRBV_sig, inherit.aes = FALSE, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(0.13, 0.12, 0.12, 0.12, 0.12, 0.13, 0.12, 0.13, 0.12, 0.13, 0.12, 0.13, 0.12, 0.12, 0.12, 0.12, 0.13, 0.12, 0.14, 0.13, 0.12, 0.12, 0.12, 0.13, 0.12), tip.length = 0, vjust = 0.6, size=6, bracket.size = 0.4) +
  theme(axis.text.y=element_text(color="black", size=14),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title=element_text(size=16),
        legend.position = "none",
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 14), 
        strip.text = element_text(face="bold.italic", size=13, color = "black"),
        strip.background = element_rect(fill="#F5F5F5", color="black", linewidth = 0.75),
        strip.placement = "inside")+
  labs(x= "",
       y= expression("Proportion"~italic(TRBV)~ " Usage"),
       color = "Subject",
       fill = "PD-1/ICOS Expression")
chimp_TRBV_sig_p 


#### TRA
chimp_TRAVJ_pre <- chimp_TRA %>%
  group_by(v_gene, PD1_ICOS) %>%
  mutate(v_count = n()) %>%
  group_by(PD1_ICOS) %>%
  mutate(v_total = n()) %>%
  mutate(v_prop = v_count/v_total) %>%
  mutate(vj_gene = paste(v_gene, j_gene)) %>%
  group_by(vj_gene, PD1_ICOS) %>%
  mutate(vj_count = n()) %>%
  mutate(vj_prop = vj_count/v_total) %>%
  mutate(v_gene = gsub("TRAV", "", v_gene)) %>%
  mutate(v_gene = gsub("DV[0-9]", "", v_gene)) %>%
  ungroup()

chimp_TRAV <- chimp_TRAVJ_pre %>%
  distinct(v_gene, PD1_ICOS, v_count, v_total, v_prop) %>%
  complete(v_gene, PD1_ICOS, fill = list(v_count = 0, v_prop = 0)) %>%
  group_by(PD1_ICOS) %>%
  mutate(v_total = sum(v_count)) 

### Statistics
## list of all TRAV genes
listTRAV <- unique(chimp_TRAV$v_gene)

## fisher testing
chimp_TRAV_stat <- data.frame()
for(i in 1:length(listTRAV)){
  VG <- filter(chimp_TRAV, v_gene == listTRAV[i])
  VG_neg <- VG %>%
    filter(PD1_ICOS == "negative") 
  VG_int <- VG %>%
    filter(PD1_ICOS == "intermediate")
  VG_high <- VG %>%
    filter(PD1_ICOS == "high")
  
  high_neg <- fisher.test(matrix(c(VG_high$v_count, VG_neg$v_count, VG_high$v_total-VG_high$v_count, VG_neg$v_total-VG_neg$v_count), ncol=2))$p.value
  
  high_int <- fisher.test(matrix(c(VG_high$v_count, VG_int$v_count, VG_high$v_total-VG_high$v_count, VG_int$v_total-VG_int$v_count), ncol=2))$p.value
  
  int_neg <- fisher.test(matrix(c(VG_int$v_count, VG_neg$v_count, VG_int$v_total-VG_int$v_count, VG_neg$v_total-VG_neg$v_count), ncol=2))$p.value
  
  VG_df <- cbind.data.frame(v_gene = listTRAV[i], high_neg, high_int, int_neg)
  chimp_TRAV_stat <- rbind.data.frame(chimp_TRAV_stat, VG_df, stringsAsFactors = FALSE)
}

## correcting P-values for multiple comparisons
chimp_TRAV_stat_corr <-data.frame(chimp_TRAV_stat$v_gene, matrix(p.adjust(as.vector(as.matrix(chimp_TRAV_stat[,-1])), method="BH"),ncol=3))
colnames(chimp_TRAV_stat_corr) <- colnames(chimp_TRAV_stat)

## format statistical comparisons and add p signficance asterisks
chimp_TRAV_sig <- chimp_TRAV_stat_corr %>%
  filter(high_neg < 0.05 | high_int < 0.05 | int_neg < 0.05) %>%
  gather(key = group, value = p.adj, -v_gene) %>%
  mutate(group1 = gsub("_.*", "", group)) %>%
  mutate(group2 = gsub(".*_", "", group)) %>%
  mutate(group1 = gsub("neg", "negative", group1)) %>%
  mutate(group1 = gsub("int", "intermediate", group1)) %>%
  mutate(group2 = gsub("neg", "negative", group2)) %>%
  mutate(group2 = gsub("int", "intermediate", group2)) %>%
  mutate(.y. = "v_prop") %>%
  mutate(p.adj.signif = ifelse(p.adj < 0.0001, "****",
                           ifelse(p.adj < 0.001, "***",
                                  ifelse(p.adj < 0.01, "**",
                                         ifelse(p.adj < 0.05, "*", "ns"))))) %>%
  mutate(numbering = sub("[0-9]*-", "", v_gene)) %>%
  mutate(fam = sub("-.*", "", v_gene)) %>%
  group_by(fam) %>%
  arrange(as.numeric(fam), as.numeric(numbering)) %>%
  ungroup() %>%
  select(v_gene, group1, group2, p.adj, p.adj.signif)
  
chimp_TRAV_sig$v_gene <- factor(chimp_TRAV_sig$v_gene, levels = unique(chimp_TRAV_sig$v_gene))


### Plots
## heatmaps
# format data
chimp_TRAV_hm_pre <-chimp_TRAV[,c(1:2,5)] %>%
  spread(key = v_gene, value = v_prop)

chimp_TRAV_hm_pre$PD1_ICOS <- factor(chimp_TRAV_hm_pre$PD1_ICOS, levels=c("high", "intermediate", "negative"))

chimp_TRAV_hm <- chimp_TRAV_hm_pre[order(chimp_TRAV_hm_pre$PD1_ICOS),]

chimp_TRAV_hm <- chimp_TRAV_hm[,-1]

# reorder TRBV genes to be in numerical order rather than alphanumeric (i.e., 1, 2, 3... rather than 1, 10, 11)
reorder_chimp_TRAV <- data.frame(vname = colnames(chimp_TRAV_hm)) %>%
  mutate(numbering = sub("[0-9]*-", "", vname)) %>%
  mutate(fam = sub("-.*", "", vname)) %>%
  group_by(fam) %>%
  arrange(as.numeric(fam), as.numeric(numbering))

chimp_TRAV_hm_arr <- t(chimp_TRAV_hm)[reorder_chimp_TRAV$vname,] 

# plot
pheatmap(chimp_TRAV_hm_arr, scale="none", cluster_rows = FALSE, cluster_cols=FALSE, border_color = NA, legend = TRUE, show_rownames = TRUE, show_colnames = TRUE, cellheight = 13, cellwidth = 25, fontsize = 13, fontsize_row = 14, fontsize_col = 12, angle_col = "45", labels_row = make_italics(rownames(chimp_TRAV_hm_arr)), labels_col = make_bold(c("High", "Int", "Low")), main = "", color= my_palette)


### barplot (only showing genes with statistically significant differences)
## format data for plotting
chimp_TRAV_subset <- chimp_TRAVJ_pre %>%
  group_by(subj, v_gene, PD1_ICOS) %>%
  mutate(v_count_subj = n()) %>%
  group_by(subj, PD1_ICOS) %>%
  mutate(v_total_subj = n()) %>%
  mutate(v_prop_subj = v_count_subj/v_total_subj) %>%
  filter(v_gene %in% chimp_TRAV_sig$v_gene) %>%
  distinct(v_gene, subj, PD1_ICOS, v_prop, v_prop_subj) %>%
  mutate(v_prop_subj_round = round(v_prop_subj, digits =3)) %>%
  group_by(v_gene, PD1_ICOS, v_prop_subj_round) %>%
  mutate(overlap = n()) %>%
  mutate(overlap = ifelse(overlap > 1, "yes", "no")) %>%
  mutate(numbering = sub("[0-9]*-", "", v_gene)) %>%
  mutate(fam = sub("-.*", "", v_gene)) %>%
  group_by(fam) %>%
  arrange(as.numeric(fam), as.numeric(numbering)) %>%
  ungroup() %>%
  select(-fam, -numbering)

chimp_TRAV_subset$v_gene <- factor(chimp_TRAV_subset$v_gene, levels = unique(chimp_TRAV_subset$v_gene))

# plot
chimp_TRAV_sig_p <- ggplot() +
  geom_bar(data = chimp_TRAV_subset, aes(x = PD1_ICOS, y = v_prop, fill = PD1_ICOS), position="dodge", stat="identity", size=1, alpha = 0.8) +
  geom_point(data = chimp_TRAV_subset[chimp_TRAV_subset$overlap == "no",], aes(x = PD1_ICOS, y = v_prop_subj), color = "black", size=1.5, shape = 19) +
  geom_point(data = chimp_TRAV_subset[chimp_TRAV_subset$overlap == "yes",], aes(x = PD1_ICOS, y = v_prop_subj), color = "black", size=1.5, shape = 19, position = position_dodge(width = 0.9)) +
  facet_wrap(. ~ v_gene, nrow = 1, strip.position = "top") +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#999999"), labels = c("High", "Intermediate", "Negative")) +
 scale_y_continuous(breaks = seq(from = 0, to = 0.10, by = .02), limits = c(0,0.10), expand = c(0,0))+
  guides(fill = guide_legend(order =1),
         color = guide_legend(override.aes = list(size=2.45, order =2)))+
  stat_pvalue_manual(chimp_TRAV_sig, inherit.aes = FALSE, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(0.09, 0.085, 0.09), tip.length = 0, vjust = 0.6, size=6, bracket.size = 0.4) +
  theme(axis.text.y=element_text(color="black", size=14),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title=element_text(size=16),
        legend.position = "none",
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 14), 
        strip.text = element_text(face="bold.italic", size=13, color = "black"),
        strip.background = element_rect(fill="#F5F5F5", color="black", linewidth = 0.75),
        strip.placement = "inside")+
  labs(x= "",
       y= expression("Proportion"~italic(TRAV)~ " Usage"),
       color = "Subject",
       fill = "PD-1/ICOS Expression")
chimp_TRAV_sig_p


#### Identification of public clonotypes (same V gene, J gene, CDRH3 sequence present in 2 members of the same group)
### TRB
pub_clon_TRB <- chimp_TRB %>%
  distinct(subj, PD1_ICOS, cdr3_aa, v_gene, j_gene) %>%
  group_by(cdr3_aa, v_gene, j_gene) %>%
  mutate(num_clon = n()) %>%
  filter(num_clon > 1) %>%
  mutate(num_subj = n_distinct(subj)) %>%
  filter(num_subj > 1) %>%
  summarise(group = paste(PD1_ICOS, collapse = ", "),
            subject = paste(subj, collapse = ", ")) %>%
  filter(group == "high, high" | group == "high, high, high" | group == "intermediate, intermediate" | group == "negative, negative")

## find proportions
# total num clonotypes per group
TRB_total_clono <- chimp_TRB %>%
  distinct(subj, PD1_ICOS, cdr3_aa, v_gene, j_gene) %>%
  group_by(PD1_ICOS) %>%
  mutate(num_clon = n()) %>%
  ungroup() %>%
  distinct(PD1_ICOS, num_clon)

# proprotions
pub_clon_prop_TRB <- pub_clon_TRB %>%
  mutate(group = gsub("\\,.*", "", group)) %>%
  group_by(group) %>%
  mutate(num_pub = n()) %>%
  distinct(group, num_pub) %>%
  ungroup %>%
  mutate(total = TRB_total_clono$num_clon) %>%
  mutate(num_nonpub = total - num_pub) %>%
  mutate(prop = num_pub/total) 

## Statistics
# make matrices for fisher test
prop_pub_TRB_ni <- matrix(c(pub_clon_prop_TRB[pub_clon_prop_TRB$group == "negative",]$num_pub,
                        pub_clon_prop_TRB[pub_clon_prop_TRB$group == "intermediate",]$num_pub,
                        pub_clon_prop_TRB[pub_clon_prop_TRB$group == "negative",]$num_nonpub,
                        pub_clon_prop_TRB[pub_clon_prop_TRB$group == "intermediate",]$num_nonpub), ncol=2)
prop_pub_TRB_nh <- matrix(c(pub_clon_prop_TRB[pub_clon_prop_TRB$group == "negative",]$num_pub,
                        pub_clon_prop_TRB[pub_clon_prop_TRB$group == "high",]$num_pub,
                        pub_clon_prop_TRB[pub_clon_prop_TRB$group == "negative",]$num_nonpub,
                        pub_clon_prop_TRB[pub_clon_prop_TRB$group == "high",]$num_nonpub), ncol=2)
prop_pub_TRB_ih <- matrix(c(pub_clon_prop_TRB[pub_clon_prop_TRB$group == "intermediate",]$num_pub,
                        pub_clon_prop_TRB[pub_clon_prop_TRB$group == "high",]$num_pub,
                        pub_clon_prop_TRB[pub_clon_prop_TRB$group == "intermediate",]$num_nonpub,
                        pub_clon_prop_TRB[pub_clon_prop_TRB$group == "high",]$num_nonpub), ncol=2)

# data frame with adjusted p values and significance indicated 
pub_clon_stat_TRB <- data.frame(.y. = "prop", 
                                group1 = c("negative", "negative", "intermediate"), 
                                group2 = c("intermediate", "high", "high"), 
                                p.adj = p.adjust(c(fisher.test(prop_pub_TRB_ni)$p.value, fisher.test(prop_pub_TRB_nh)$p.value, fisher.test(prop_pub_TRB_ih)$p.value), "BH"),
                                p.adj.signif = c("ns", "****", "****"))

# format TRB chain
pub_clon_prop_TRB$chain <- c(paste0("TCR ", "*\U03B2*"))

## Plot
pub_TRB_p <- ggplot() +
  geom_bar(data = pub_clon_prop_TRB, aes(x = group, y = prop*100, fill = group), position="dodge", stat="identity", size=2) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#999999")) +
  scale_x_discrete(breaks = c("high", "intermediate", "negative"), labels = c("High", "Int.", "Neg.")) +
  facet_grid(. ~ chain)+
  stat_pvalue_manual(pub_clon_stat_TRB, inherit.aes = FALSE, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(0.75, 0.79), tip.length = 0, vjust = 0.75, size=5, bracket.size = 0.4) +
  theme(axis.text=element_text(color="black", size=12),
        panel.background = element_rect(fill = "white"),
        #panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title=element_text(size=14),
        legend.position = "none",
        strip.text = element_markdown(face="bold", size=16, color = "black"),
        strip.background = element_rect(fill="#F5F5F5", color="lightgray", linewidth = 0.25))+
  labs(x= "PD-1/ICOS Expression",
       y= "% Public Clonotypes")
pub_TRB_p

### TRA
pub_clon_TRA <- chimp_TRA %>%
  distinct(subj, PD1_ICOS, cdr3_aa, v_gene, j_gene) %>%
  group_by(cdr3_aa, v_gene, j_gene) %>%
  mutate(num_clon = n()) %>%
  filter(num_clon > 1) %>%
  mutate(num_subj = n_distinct(subj)) %>%
  filter(num_subj > 1) %>%
  summarise(group = paste(PD1_ICOS, collapse = ", "),
            subject = paste(subj, collapse = ", ")) %>%
  filter(group == "high, high" | group == "intermediate, intermediate" | group == "negative, negative")

# find proportions
# total num clonotypes per group
TRA_total_clono <- chimp_TRA %>%
  distinct(subj, PD1_ICOS, cdr3_aa, v_gene, j_gene) %>%
  group_by(PD1_ICOS) %>%
  mutate(num_clon = n()) %>%
  ungroup() %>%
  distinct(PD1_ICOS, num_clon)

# proportions
pub_clon_prop_TRA <- pub_clon_TRA %>%
  mutate(group = gsub("\\,.*", "", group)) %>%
  group_by(group) %>%
  mutate(num_pub = n()) %>%
  distinct(group, num_pub) %>%
  ungroup %>%
  mutate(total = TRA_total_clono$num_clon) %>%
  mutate(num_nonpub = total - num_pub) %>%
  mutate(prop = num_pub/total) 

## Statistics
# make matrices for fisher testing
prop_pub_TRA_ni <- matrix(c(pub_clon_prop_TRA[pub_clon_prop_TRA$group == "negative",]$num_pub,
                        pub_clon_prop_TRA[pub_clon_prop_TRA$group == "intermediate",]$num_pub,
                        pub_clon_prop_TRA[pub_clon_prop_TRA$group == "negative",]$num_nonpub,
                        pub_clon_prop_TRA[pub_clon_prop_TRA$group == "intermediate",]$num_nonpub), ncol=2)
prop_pub_TRA_nh <- matrix(c(pub_clon_prop_TRA[pub_clon_prop_TRA$group == "negative",]$num_pub,
                        pub_clon_prop_TRA[pub_clon_prop_TRA$group == "high",]$num_pub,
                        pub_clon_prop_TRA[pub_clon_prop_TRA$group == "negative",]$num_nonpub,
                        pub_clon_prop_TRA[pub_clon_prop_TRA$group == "high",]$num_nonpub), ncol=2)
prop_pub_TRA_ih <- matrix(c(pub_clon_prop_TRA[pub_clon_prop_TRA$group == "intermediate",]$num_pub,
                        pub_clon_prop_TRA[pub_clon_prop_TRA$group == "high",]$num_pub,
                        pub_clon_prop_TRA[pub_clon_prop_TRA$group == "intermediate",]$num_nonpub,
                        pub_clon_prop_TRA[pub_clon_prop_TRA$group == "high",]$num_nonpub), ncol=2)

# data frame with adjusted p values and significance indicated 
pub_clon_stat_TRA <- data.frame(.y. = "prop", 
                                group1 = c("negative", "negative", "intermediate"), 
                                group2 = c("intermediate", "high", "high"), 
                                p.adj = p.adjust(c(fisher.test(prop_pub_TRA_ni)$p.value, fisher.test(prop_pub_TRA_nh)$p.value, fisher.test(prop_pub_TRA_ih)$p.value), "BH"),
                                p.adj.signif = c("ns", "ns", "ns"))

# format TRA chain
pub_clon_prop_TRA$chain <- c(paste0("TCR ", "*\U03B1*"))

## Plot
pub_TRA_p <- ggplot() +
  geom_bar(data = pub_clon_prop_TRA, aes(x = group, y = prop*100, fill = group), position="dodge", stat="identity", size=2) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#999999")) +
  scale_x_discrete(breaks = c("high", "intermediate", "negative"), labels = c("High", "Int.", "Neg.")) +
  facet_grid(. ~ chain)+
  theme(axis.text=element_text(color="black", size=12),
        panel.background = element_rect(fill = "white"),
        #panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title=element_text(size=14),
        legend.position = "none",
        strip.text = element_markdown(face="bold", size=16, color = "black"),
        strip.background = element_rect(fill="#F5F5F5", color="lightgray", linewidth = 0.25))+
  labs(x= "PD-1/ICOS Expression",
       y= "% Public Clonotypes")
pub_TRA_p

## Plot TRB and TRA together

# make combined data frame
pub_clon_prop_TRBA <- rbind.data.frame(pub_clon_prop_TRB, pub_clon_prop_TRA, stringsAsFactors = FALSE)
pub_clon_stat_TRBA <- rbind.data.frame(mutate(pub_clon_stat_TRB, chain = "TCR *β*"), mutate(pub_clon_stat_TRA, chain = "TCR *α*"), stringsAsFactors = FALSE)

pub_clon_prop_TRBA$chain <- factor(pub_clon_prop_TRBA$chain, levels = c("TCR *β*", "TCR *α*"))
pub_clon_stat_TRBA$chain <- factor(pub_clon_stat_TRBA$chain, levels = c("TCR *β*", "TCR *α*"))

# plot
pub_TRBA_p <- ggplot() +
  geom_bar(data = pub_clon_prop_TRBA, aes(x = group, y = prop*100, fill = group), position="dodge", stat="identity", size=2) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#999999")) +
  scale_x_discrete(breaks = c("high", "intermediate", "negative"), labels = c("High", "Int.", "Neg.")) +
  stat_pvalue_manual(pub_clon_stat_TRBA, inherit.aes = FALSE, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(0.75, 0.79), tip.length = 0, vjust = 0.75, size=5, bracket.size = 0.4) +
  facet_wrap(. ~ chain, scales = "free_y")+
  theme(axis.text=element_text(color="black", size=14),
        panel.background = element_rect(fill = "white"),
        #panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title=element_text(size=16),
        legend.position = "none",
        strip.text = element_markdown(face="bold", size=16, color = "black"),
        strip.background = element_rect(fill="#F5F5F5", color="lightgray", linewidth = 0.25))+
  labs(x= "PD-1/ICOS Expression",
       y= "% Public Clonotypes")
pub_TRBA_p

##### Alignment of CDRH3s
#### PD1/ICOS high public clonotype
pub_high_TRB <- pub_clon_TRB %>%
  filter(grepl("high", group)) %>%
  mutate(cdr3_len = nchar(cdr3_aa)) 

# write.table(pub_high_TRB, "~/Desktop/chimp_TCR/Bulk_TCR_seq/pub_high_TRB.txt", quote=FALSE, sep="\t", row.names = FALSE)
# perform alignment using algorithm (written in R script run on the computing cluster) and download

### read in alignment
pub_high_TRB_align <- read.csv("~/Desktop/chimp_TCR/Bulk_TCR_seq/pub_high_TRB_alignment.csv") %>%
  select(-X)
  

#### PD1/ICOS negative clonotypes 
neg_TRB <- chimp_TRB %>%
  filter(PD1_ICOS == "negative") %>%
  filter(cdr3_len <= max(pub_high_TRB$cdr3_len) & cdr3_len >= min(pub_high_TRB$cdr3_len))

# write.table(neg_TRB, "~/Desktop/neg_TRB", quote=FALSE, sep="\t", row.names = FALSE)
# perform alignment using algorithm (written in R script run on the computing cluster) and download

### read in alignment
neg_TRB_align <- read.csv("~/Desktop/chimp_TCR/Bulk_TCR_seq/neg_TRB_alignment.csv") %>%
  select(-X)
  

#### Determine amino acid usage for alignments
# The goal is to identify amino acids enriched in given positions in public PD1/ICOS clonotypes. Will quantify amino acids used by public PD1/ICOS high at each position and determine if usage is statistically significantly higher than for PD1/ICOS negative usage at that position

### Determine usage
pub_high_cdr3_usage <- AA_usage(pub_high_TRB_align$Seq, nchar(pub_high_TRB_align$Seq[1]))
neg_cdr3_usage <- AA_usage(neg_TRB_align$Seq, nchar(neg_TRB_align$Seq[1]))

### Statistics 
pub_high_usage_stat_pre <- data.frame()
for(i in 1:nchar(pub_high_TRB_align$Seq[1])){
  # filter to only include analysis of amino acids present in public PD1/ICOS high clonotypes
  pub_filt <- pub_high_cdr3_usage[,i]
  pub_filt2 <- pub_filt[pub_filt!=0]
  neg_filt <- neg_cdr3_usage[,i] 
  neg_filt2 <- neg_filt[names(neg_filt) %in% names(pub_filt2)]

  amino_acids <- names(pub_filt2)
  
  # fisher testing of proportion of AA usage at each position 
  pub_neg_fish <- data.frame()
  for(j in 1:length(amino_acids)){
    AA <- amino_acids[j]
    AA_count_pub <- pub_filt[names(pub_filt) == AA]
    AA_count_neg <- neg_filt[names(neg_filt) == AA]
    AA_total_pub <- sum(pub_filt)
    AA_total_neg <- sum(neg_filt)
    fish_mat <-matrix(c(AA_count_pub, AA_count_neg, AA_total_pub-AA_count_pub, AA_total_neg-AA_count_neg), ncol = 2)
    fish_test <- cbind.data.frame(AA = AA, 
                                  comp_prop = AA_count_pub/AA_total_pub, 
                                  neg_prop = AA_count_neg/AA_total_neg,
                                  P_value = fisher.test(fish_mat)$p.value,
                                  comp = "pub")  
    pub_neg_fish <- rbind.data.frame(pub_neg_fish, fish_test, stringsAsFactors = FALSE)
  }
  
  pub_neg_fish2 <- pub_neg_fish %>%
    mutate(Position = i)
  
  pub_high_usage_stat_pre <- rbind.data.frame(pub_high_usage_stat_pre , pub_neg_fish2, stringsAsFactors = FALSE)
}

## Correct p-values for multiple comparisons and find fold-change
# convert IMGT numbering to Kabat numbering
kabat_cdr3_TRB <- data.frame(Position = c(1:17), Kabat = c("92", "93", "94", "95", "96", "97", "98", "99", "100", "100a", "100b", "100c", "100d", "100e", "101", "102", "103"))

# for fold-change calculation, have to replace 0s to avoid attempting to divide by 0. We will replace 0 with 1/2 of the smallest non-zero number 
cdr3_align_zero_corr_IGH <- min(
  min(pub_high_usage_stat_pre$comp_prop[pub_high_usage_stat_pre$comp_prop != 0])/2,
  min(pub_high_usage_stat_pre$neg_prop[pub_high_usage_stat_pre$neg_prop != 0])/2)

# replace 0s, calculate fold change, adjust P values
TRB_AA_usage_stat <- pub_high_usage_stat_pre %>%
  left_join(kabat_cdr3_TRB, by= "Position") %>%
  mutate(comp_prop_adj = ifelse(comp_prop == 0, cdr3_align_zero_corr_IGH, comp_prop)) %>%
  mutate(neg_prop_adj = ifelse(neg_prop == 0, cdr3_align_zero_corr_IGH, neg_prop)) %>%
  mutate(FC = comp_prop_adj/neg_prop_adj) %>%
  mutate(P_value_adj = p.adjust(P_value, "fdr"))

# determine significance (fold change > 1.6, P value < 0.05)
TRB_AA_usage_sig <- TRB_AA_usage_stat  %>%
  mutate(log2FC = log2(FC)) %>%
  filter(log2FC >= 0.6) %>%
  filter(P_value_adj < 0.05) %>%
  group_by(AA, Position) %>%
  mutate(num_use = n()) %>%
  arrange(desc(num_use), Position, AA) %>%
  ungroup()

TRB_AA_usage_sig$Kabat <- factor(TRB_AA_usage_sig$Kabat, levels = c("92", "93", "94", "95", "96", "97", "98", "99", "100", "100a", "100b", "100c", "100d", "100e", "101", "102", "103"))

TRB_AA_usage_sig_arr <- TRB_AA_usage_sig %>%
  group_by(AA, Kabat) %>%
  summarise(comp = str_c(comp, collapse="_")) %>%
  arrange(Kabat, AA)


#### Motif analysis: identified 2 possible motifs based on the enriched AA analysis above (YRGxAT and Q.GQ). Will search for these motifs in PD1/ICOS high clonotypes to see if they are enriched

### YRG.AT
motif_TRB <- chimp_TRB %>%
  mutate(motif = ifelse(grepl("YRG.AT", cdr3_aa), "yes", "no")) %>%
  group_by(PD1_ICOS) %>%
  mutate(total = n()) %>%
  group_by(PD1_ICOS, motif) %>%
  mutate(count = n()) %>%
  distinct(PD1_ICOS, motif, count, total) %>%
  ungroup() %>%
  complete(PD1_ICOS, motif, fill = list(count = 0)) %>%
  group_by(PD1_ICOS) %>%
  fill(total, .direction = "updown") %>%
  ungroup() %>%
  mutate(prop = count/total)


### Statistics
## fisher testing
motif_TRB_high_neg <- fisher.test(matrix(c(motif_TRB[motif_TRB$PD1_ICOS == "high" & motif_TRB$motif == "yes",]$count,
                                 motif_TRB[motif_TRB$PD1_ICOS == "negative" & motif_TRB$motif == "yes",]$count,
                                 motif_TRB[motif_TRB$PD1_ICOS == "high" & motif_TRB$motif == "no",]$count,
                                 motif_TRB[motif_TRB$PD1_ICOS == "negative" & motif_TRB$motif == "no",]$count),
                                 ncol=2))$p.value
motif_TRB_high_int <- fisher.test(matrix(c(motif_TRB[motif_TRB$PD1_ICOS == "high" & motif_TRB$motif == "yes",]$count,
                                 motif_TRB[motif_TRB$PD1_ICOS == "intermediate" & motif_TRB$motif == "yes",]$count,
                                 motif_TRB[motif_TRB$PD1_ICOS == "high" & motif_TRB$motif == "no",]$count,
                                 motif_TRB[motif_TRB$PD1_ICOS == "intermediate" & motif_TRB$motif == "no",]$count),
                                 ncol=2))$p.value
motif_TRB_int_neg <- fisher.test(matrix(c(motif_TRB[motif_TRB$PD1_ICOS == "intermediate" & motif_TRB$motif == "yes",]$count,
                                 motif_TRB[motif_TRB$PD1_ICOS == "negative" & motif_TRB$motif == "yes",]$count,
                                 motif_TRB[motif_TRB$PD1_ICOS == "intermediate" & motif_TRB$motif == "no",]$count,
                                 motif_TRB[motif_TRB$PD1_ICOS == "negative" & motif_TRB$motif == "no",]$count),
                                 ncol=2))$p.value

## correct p values for multiple comparisons and add significance 
motif_TRB_stats <- data.frame(.y. = rep("motif", 3),
                              group1 = c("high", "high", "intermediate"),
                              group2 = c("negative", "intermediate", "negative"),
                              p.adj = p.adjust(c(motif_TRB_high_neg, motif_TRB_high_int, motif_TRB_int_neg), "BH"),
                              p.adj.signif = c("****", "****", "ns"))                       
                        
### Plot
## proportion of each group using the motif
motif_TRB_yes <- motif_TRB[motif_TRB$motif == "yes",]
motif_TRB_yes$PD1_ICOS <- factor(motif_TRB_yes$PD1_ICOS, levels = c("high", "intermediate", "negative"))

# plot
motif_TRB_p <- ggplot() +
  geom_bar(data = motif_TRB_yes, aes(x = PD1_ICOS, y = prop*100, fill = PD1_ICOS), position="dodge", stat="identity", size=2) +
  geom_text(data = motif_TRB_yes, aes(x = PD1_ICOS, y = prop*100, label = paste0("n = ",count)), position=position_dodge(width=0.9), vjust=-0.35, size=4, fontface = "bold.italic")+
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#999999")) +
  scale_x_discrete(breaks = c("high", "intermediate", "negative"), labels = c("High", "Int.", "Neg.")) +
  stat_pvalue_manual(motif_TRB_stats, inherit.aes = FALSE, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(0.81, 0.85), tip.length = 0, vjust = 0.75, size=5, bracket.size = 0.4) +
  theme(axis.text=element_text(color="black", size=12),
        panel.background = element_rect(fill = "white"),
        #panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title=element_text(size=14),
        legend.position = "none",
        strip.text = element_markdown(face="bold", size=16, color = "black"),
        strip.background = element_rect(fill="#F5F5F5", color="lightgray", linewidth = 0.25))+
  labs(x= "PD-1/ICOS Expression",
       y= "% YRGxAT Usage")
motif_TRB_p


### Q.GQ
motif_TRB2 <- chimp_TRB %>%
  mutate(motif = ifelse(grepl("Q.GQ", cdr3_aa), "yes", "no")) %>%
  group_by(PD1_ICOS) %>%
  mutate(total = n()) %>%
  group_by(PD1_ICOS, motif) %>%
  mutate(count = n()) %>%
  distinct(PD1_ICOS, motif, count, total) %>%
  ungroup() %>%
  complete(PD1_ICOS, motif, fill = list(count = 0)) %>%
  group_by(PD1_ICOS) %>%
  fill(total, .direction = "updown") %>%
  ungroup() %>%
  mutate(prop = count/total)


### Statistics
## fisher testing
motif_TRB2_high_neg <- fisher.test(matrix(c(motif_TRB2[motif_TRB2$PD1_ICOS == "high" & motif_TRB2$motif == "yes",]$count,
                                 motif_TRB2[motif_TRB2$PD1_ICOS == "negative" & motif_TRB2$motif == "yes",]$count,
                                 motif_TRB2[motif_TRB2$PD1_ICOS == "high" & motif_TRB2$motif == "no",]$count,
                                 motif_TRB2[motif_TRB2$PD1_ICOS == "negative" & motif_TRB2$motif == "no",]$count),
                                 ncol=2))$p.value
motif_TRB2_high_int <- fisher.test(matrix(c(motif_TRB2[motif_TRB2$PD1_ICOS == "high" & motif_TRB2$motif == "yes",]$count,
                                 motif_TRB2[motif_TRB2$PD1_ICOS == "intermediate" & motif_TRB2$motif == "yes",]$count,
                                 motif_TRB2[motif_TRB2$PD1_ICOS == "high" & motif_TRB2$motif == "no",]$count,
                                 motif_TRB2[motif_TRB2$PD1_ICOS == "intermediate" & motif_TRB2$motif == "no",]$count),
                                 ncol=2))$p.value
motif_TRB2_int_neg <- fisher.test(matrix(c(motif_TRB2[motif_TRB2$PD1_ICOS == "intermediate" & motif_TRB2$motif == "yes",]$count,
                                 motif_TRB2[motif_TRB2$PD1_ICOS == "negative" & motif_TRB2$motif == "yes",]$count,
                                 motif_TRB2[motif_TRB2$PD1_ICOS == "intermediate" & motif_TRB2$motif == "no",]$count,
                                 motif_TRB2[motif_TRB2$PD1_ICOS == "negative" & motif_TRB2$motif == "no",]$count),
                                 ncol=2))$p.value

## correct p values for multiple comparisons and add significance 
motif_TRB2_stats <- data.frame(.y. = rep("motif", 3),
                              group1 = c("high", "high", "intermediate"),
                              group2 = c("negative", "intermediate", "negative"),
                              p.adj = p.adjust(c(motif_TRB2_high_neg, motif_TRB2_high_int, motif_TRB2_int_neg), "BH"),
                              p.adj.signif = c("**", "**", "ns"))                       
                        
### Plot
## proportion of each group using the motif
motif_TRB2_yes <- motif_TRB2[motif_TRB2$motif == "yes",]
motif_TRB2_yes$PD1_ICOS <- factor(motif_TRB2_yes$PD1_ICOS, levels = c("high", "intermediate", "negative"))

# plot
motif_TRB2_p <- ggplot() +
  geom_bar(data = motif_TRB2_yes, aes(x = PD1_ICOS, y = prop*100, fill = PD1_ICOS), position="dodge", stat="identity", size=2) +
  geom_text(data = motif_TRB2_yes, aes(x = PD1_ICOS, y = prop*100, label = paste0("n = ",count)), position=position_dodge(width=0.9), vjust=-0.35, size=4, fontface = "bold.italic")+
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#999999")) +
  scale_x_discrete(breaks = c("high", "intermediate", "negative"), labels = c("High", "Int.", "Neg.")) +
  stat_pvalue_manual(motif_TRB2_stats, inherit.aes = FALSE, hide.ns = TRUE, label = "{p.adj.signif}", y.position = c(1.12, 1.15), tip.length = 0, vjust = 0.75, size=5, bracket.size = 0.4) +
  theme(axis.text=element_text(color="black", size=12),
        panel.background = element_rect(fill = "white"),
        #panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title=element_text(size=14),
        legend.position = "none",
        strip.text = element_markdown(face="bold", size=16, color = "black"),
        strip.background = element_rect(fill="#F5F5F5", color="lightgray", linewidth = 0.25))+
  labs(x= "PD-1/ICOS Expression",
       y= "% QxGQ Usage")
motif_TRB2_p
