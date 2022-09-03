library("file2meco")
library("microeco")
library("tidyverse")

# importing the dataset
dataset=qiime2meco(ASV_data="table.qza", phylo_tree="rep_seqs_aligned_masked_tree_rooted.qza", 
            sample_data="metadata.tsv", taxonomy_data="classification.qza")

# make the taxonomic information unified, very important
dataset$tax_table=tidy_taxonomy(dataset$tax_table)

# renaming groups
levels(dataset$sample_table$treatment.group)=c("FH","FW","NP","Control","Inositol","Folic","NC")

# selecting only the relevant groups
dataset$sample_table = subset(dataset$sample_table, treatment.group == "Control"  |
                                treatment.group =="Inositol" | treatment.group =="Folic")

# removing ASVs unassigned at kingdom level and archaea
dataset$tax_table =subset(dataset$tax_table,Kingdom == "k__Bacteria")

#make data consistent
dataset$tidy_dataset()

########################
# ABUNDANCE ANALYSIS
########################

# Some sequences unassigned at phylum level with a very low abundance are creating
# an "Others" category. It can be removed:
dataset$tax_table=dataset$tax_table = subset(dataset$tax_table,Phylum != "p__")
#make data consistent again
dataset$tidy_dataset()

# calculating abundances at each taxonomic rank (stored in $taxa_abund)
dataset$cal_abund()

# PHYLA

# create object with transformed abundance data
# select the 10 most abundant phyla (in this case it doesn't make a difference)
t0phyla_c=trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10,show = 0)


# create a plot with all samples
g0phyla_c=t0phyla_c$plot_bar(color_values = RColorBrewer::brewer.pal(8, "Set1"),
                             others_color = "grey70", facet = "treatment.group", 
                             xtext_keep = FALSE, legend_text_italic = FALSE)

# same as before, but with the samples grouped
t0phyla_g = trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10, groupmean = "treatment.group")
g0phyla_g = t0phyla_g$plot_bar(color_values = RColorBrewer::brewer.pal(8, "Set1"),
                               others_color = "grey70", legend_text_italic = FALSE) +  
  theme_classic() + theme(axis.title.y = element_text(size = 18))


# GENERA

#selecting top 40
t0genera_40 = trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 40)
g0genera_heatmap=t0genera_40$plot_heatmap(facet = "treatment.group", xtext_keep = FALSE, 
                                          withmargin = FALSE,ytext_size = 7)

#selecting top 15
t0genera_15_g = trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 15,
                                 groupmean = "treatment.group")
g0genera_bar=t0genera_15_g$plot_bar(legend_text_italic = FALSE, 
                                    color_values = c(RColorBrewer::brewer.pal(12, "Paired"),
                                                     RColorBrewer::brewer.pal(8, "Set2")))+  
  theme_classic() + theme(axis.title.y = element_text(size = 18))+guides(fill=guide_legend(ncol=2,title="Genus"))

# VENN PLOT

# creating a new dataset with merged samples into groups
t0merged = dataset$merge_samples(use_group = "treatment.group")

# making a venn plot. absolute numbers are ASVs
t0venn = trans_venn$new(t0merged, ratio = "seqratio")
g0venn=t0venn$plot_venn()

####################
# DIVERSITY ANALYSIS
####################

# ALPHA DIVERSITY

# creating transformed alpha diversity object
t0alpha = trans_alpha$new(dataset = dataset, group = "treatment.group")

# test the difference between groups (several methods: KW, wilcox, anova)

# wilcoxon pairwise comparison
t0alpha$cal_diff(method = "wilcox")

# select only significant comparisons
t0alpha$res_diff = base::subset(t0alpha$res_diff, Significance != "ns")

# Observed richness (ASV count)
g0alpha_observed=t0alpha$plot_alpha(measure="Observed",pair_compare = TRUE)

# for some reason the significant layers are displaced. package bug? fix:
g0alpha_observed[["layers"]][[3]][["stat_params"]][["xmin"]]=g0alpha_observed[["layers"]][[3]][["stat_params"]][["xmin"]]-3
g0alpha_observed[["layers"]][[3]][["stat_params"]][["xmax"]]=g0alpha_observed[["layers"]][[3]][["stat_params"]][["xmax"]]-3


# Chao1 index
g0alpha_chao1=t0alpha$plot_alpha(measure="Chao1",pair_compare = TRUE)

# sig label fix:
g0alpha_chao1[["layers"]][[3]][["stat_params"]][["xmin"]]=g0alpha_chao1[["layers"]][[3]][["stat_params"]][["xmin"]]-3
g0alpha_chao1[["layers"]][[3]][["stat_params"]][["xmax"]]=g0alpha_chao1[["layers"]][[3]][["stat_params"]][["xmax"]]-3


# Shannon index
g0alpha_shannon=t0alpha$plot_alpha(measure="Shannon",pair_compare = TRUE)

# sig label fix (No significant comparisons):
g0alpha_shannon[["layers"]][[3]]=NULL


# BETA DIVERSITY

# creating transformed beta diversity object
library(GUniFrac)
dataset$cal_betadiv(unifrac = TRUE)

t0beta = trans_beta$new(dataset = dataset, group = "treatment.group",measure = "wei_unifrac")

# Pouse PCoA as an example, PCA or NMDS is also available
t0beta$cal_ordination(ordination = "PCoA")

# plotting
g0beta_PCoA=t0beta$plot_ordination(plot_color = "treatment.group", plot_shape = "treatment.group", 
                       plot_type = c("point", "ellipse"),ellipse_chull_fill = FALSE)+
  theme_bw()

# calculate and plot sample distances WITHIN groups
t0beta$cal_group_distance()
g0beta_distance_w=t0beta$plot_group_distance(distance_pair_stat = TRUE,hide_ns=T)

# calculate and plot sample distances BETWEEN groups
t0beta$cal_group_distance(within_group = FALSE)
g0beta_distance_b=t0beta$plot_group_distance(distance_pair_stat = TRUE,hide_ns=T)


# clustering (not done here)

# library(ggdendro)
# g0beta_dendro=t0beta$plot_clustering(group = "treatment.group")

# PERTMANOVA for all groups when manova_all = TRUE
t0beta$cal_manova(manova_all = TRUE)
g0beta_MANOVA_all=t0beta$res_manova
knitr::kable(g0beta_MANOVA_all)
# PERMANOVA for each paired groups
t0beta$cal_manova(manova_all = FALSE)
g0beta_MANOVA_pairwise=t0beta$res_manova
knitr::kable(g0beta_MANOVA_pairwise)

########################
# DIFFERENTIAL ANALYSIS
########################

# METASTAT

# calculate differential test results (stored in $res_diff)
t0diff_meta = trans_diff$new(dataset = dataset, method = "metastat", group = "treatment.group", taxa_level = "Genus")

# PC vs Inositol
g2diff_meta=t0diff_meta$plot_diff_abund(use_number = 1:20, select_group = "Control - Inositol")

# PC vs FA
g3diff_meta=t0diff_meta$plot_diff_abund(use_number = 1:20, select_group = "Control - Folic")

# WILCOXON

# Phylum
# calculate differential test results (stored in $res_diff)
t0diff_wilcoxp = trans_diff$new(dataset = dataset, method = "wilcox", 
                                group = "treatment.group", taxa_level = "Phylum")


# Selecting the significant comparisons
t0diff_wilcoxp$res_diff = subset(t0diff_wilcoxp$res_diff, Significance == "***" | 
                                  Significance == "**" | Significance == "*")

g0diff_wilcoxp=t0diff_wilcoxp$plot_diff_abund(add_sig = T, keep_prefix = F,
                                              add_sig_label = "Significance",                                              
                                              group_order = c("Control","Inositol","Folic"))

# Genera (same as before)
# calculate differential test results (stored in $res_diff)
t0diff_wilcoxg = trans_diff$new(dataset = dataset, method = "wilcox", 
                                group = "treatment.group", taxa_level = "Genus")


# Selecting the significant comparisons
t0diff_wilcoxg$res_diff = subset(t0diff_wilcoxg$res_diff, Significance == "***" | 
                                  Significance == "**" | Significance == "*")

g0diff_wilcoxg=t0diff_wilcoxg$plot_diff_abund(add_sig = T, keep_prefix = F,
                                              add_sig_label = "Significance",
                                              group_order = c("Control","Inositol","Folic"))

# Selecting only the three comparisons of interest
#t0diff_wilcoxg$res_diff = subset(t0diff_wilcoxg$res_diff, Comparison== "NP - CP" |
#                                   Comparison== "Control - Folic" |  Comparison== "Control - Inositol")

#g0diff_wilcoxg2=t0diff_wilcoxg$plot_diff_abund(add_sig = T, keep_prefix = F,
#                                               add_sig_label = "Significance",
#                                               group_order = c("Control","Inositol","Folic"))


# LEFSE

# calculate differential test results (stored in $res_diff)
t0diff_lefse = trans_diff$new(dataset = dataset, method = "lefse", group = "treatment.group", alpha = 0.05, lefse_subgroup = NULL)

# Threshold is used for the LDA score selection
t0diff_lefse$plot_diff_bar(threshold = 4)

# select the 20 taxa with the highest LDA (log10)
t0diff_lefse$plot_diff_bar(use_number = 1:30, width = 0.8, group_order = c("Control", "Inositol", "Folic"))

# plot the relative abundances
g0diff_lefse_bar=t0diff_lefse$plot_diff_abund(use_number = 1:30, group_order = c("Control", "Inositol", "Folic"))

# clade_label_level 5 represent phylum level in this analysis
g0diff_lefse_cladogram=t0diff_lefse$plot_diff_cladogram(use_taxa_num = 100, 
                                                        use_feature_num = 50, 
                                                        clade_label_level = 5, 
                                                        group_order = c("Control", "Inositol", "Folic"))
g0diff_lefse_cladogram

# LDA score is a measure of the effect size: how different the groups are based on
# that taxonomic grouping


#############
# NETWORK CORRELATION ANALYSIS
############
# igraph package is needed fo the creation of networks
# WGCNA package is needed for the calculation of weighted correlation
# library("igraph")
# library("WGCNA")
# t0network = trans_network$new(dataset = dataset, cor_method = "spearman", 
#                              use_WGCNA_pearson_spearman = TRUE, filter_thres = 0.0001)

# t0network$cal_network(COR_p_thres = 0.01, COR_optimization = TRUE)

# t0network$cal_module(method = "cluster_fast_greedy")

# saving network as 
# install.packages("rgexf")
# library(rgexf)
# t0network$save_network()

# FUNCTIONAL ANALYSIS

t0func = trans_func$new(dataset)

# anotate function from database
t0func$cal_spe_func(prok_database = "FolicPROTAX")

#calculate percentages (doesnt take into account abundances by default)
t0func$cal_spe_func_perc(abundance_weighted = TRUE)

# construct a network for the example
fnetwork = trans_network$new(dataset = dataset, cal_cor = "base", taxa_level = "OTU", 
                             filter_thres = 0.0001, cor_method = "spearman")

fnetwork$cal_network(p_thres = 0.01, COR_cut = 0.7)
fnetwork$cal_module()

# convert module info to microtable object
fobj = fnetwork$trans_comm(use_col = "module")
fobj_func = trans_func$new(fobj)
fobj_func$cal_spe_func(prok_database = "FolicPROTAX")
fobj_func$cal_spe_func_perc(abundance_weighted = TRUE)
fobj_func$plot_spe_func_perc()
fobj_func$func_group_list


####################
#FUNCTIONAL ANALYSIS
#####################

# picrust2 was used to obtain the pathways from kegg and metacyc

# MetaCyc pathway
pathway_table = read.delim("path_abun_unstrat_1.tsv", row.names = 1)

sampletable=dataset$sample_table
rownames(sampletable)=gsub("-",".",rownames(dataset$sample_table))

data("MetaCyc_pathway_map")
picrust = microtable$new(otu_table = pathway_table, tax_table = MetaCyc_pathway_map,
                         sample_table = sampletable)
picrust$tidy_dataset()
picrust$cal_abund()

# transformando a un objeto trans_abund
t0pic = trans_abund$new(picrust, taxrank = "pathway", groupmean = "treatment.group")
t0pic$plot_bar(legend_text_italic = FALSE)

# Making statistical test
t0picd <- trans_diff$new(dataset = picrust, method = "wilcox", group = "treatment.group", alpha = 0.05, lefse_subgroup = NULL)
g0picd1=t0picd$plot_diff_bar()


# Making statistycal tests: Pairwise

# Comparison 2: pregnant vs pregnant+ chiro
t0metacycd3=trans_diff$new(dataset = picrust, method = "wilcox", group = "treatment.group", 
                        alpha = 0.05, lefse_subgroup = NULL)
t0metacycd3$res_diff = subset(t0metacycd3$res_diff, Comparison== "Control - Inositol")
t0metacycd3$res_abund = subset(t0metacycd3$res_abund, Group== "Control" | Group=="Inositol")
t0metacycd3$res_diff = subset(t0metacycd3$res_diff, Significance == "***" | 
                             Significance == "**" | Significance == "*")
g0metad3=t0metacycd3$plot_diff_abund(add_sig = T, keep_prefix = F,add_sig_label = "Significance")+
  theme(axis.title.x = element_text(size=10))+scale_fill_manual(values=c("#D95F02","#7570B3"))+
  scale_color_manual(values=c("#D95F02","#7570B3"))

# Comparison 3: pregnant vs pregnant folic
t0metacycd4=trans_diff$new(dataset = picrust, method = "wilcox", group = "treatment.group", 
                        alpha = 0.05, lefse_subgroup = NULL)
t0metacycd4$res_diff = subset(t0metacycd4$res_diff, Comparison== "Control - Folic")
t0metacycd4$res_abund = subset(t0metacycd4$res_abund, Group== "Control" | Group=="Folic")
t0metacycd4$res_diff = subset(t0metacycd4$res_diff, Significance == "***" | 
                             Significance == "**" | Significance == "*")
g0metad4=t0metacycd4$plot_diff_abund(add_sig = T, keep_prefix = F,add_sig_label = "Significance")+
  theme(axis.title.x = element_text(size=10))+scale_fill_manual(values=c("#D95F02","#E7298A"))+
  scale_color_manual(values=c("#D95F02","#E7298A"))



# USING KEGG
kegg_table = read.delim("pred_metagenome_unstrat_descrip.tsv", row.names = 1)

# changing row name format to match the tax table
rownames(kegg_table)=gsub("K", "ko",rownames(kegg_table))

#importing tax table
data("Tax4Fun2_KEGG")
kegg = microtable$new(otu_table = kegg_table[-c(1:2)], tax_table = Tax4Fun2_KEGG$ptw_desc,
                      sample_table = sampletable)
kegg$tidy_dataset()

# calculating abundance
kegg$cal_abund()

# transforming into an abundance object. level set to 3 (lowest)
t0kegg = trans_abund$new(kegg, taxrank = "Level.3", groupmean = "treatment.group")
t0kegg$plot_bar(legend_text_italic = FALSE)

# Making statistycal tests: Pairwise

# Comparison 2: pregnant vs pregnant+ chiro
t0keggd3=trans_diff$new(dataset = kegg, method = "wilcox", group = "treatment.group", 
                        alpha = 0.05, lefse_subgroup = NULL)
t0keggd3$res_diff = subset(t0keggd3$res_diff, Comparison== "Control - Inositol")
t0keggd3$res_abund = subset(t0keggd3$res_abund, Group== "Control" | Group=="Inositol")
t0keggd3$res_diff = subset(t0keggd3$res_diff, Significance == "***" | 
                             Significance == "**" | Significance == "*")
g0keggd3=t0keggd3$plot_diff_abund(add_sig = T, keep_prefix = F,add_sig_label = "Significance")+
  theme(axis.title.x = element_text(size=10))+scale_fill_manual(values=c("#D95F02","#7570B3"))+
  scale_color_manual(values=c("#D95F02","#7570B3"))

# Comparison 3: pregnant vs pregnant folic
t0keggd4=trans_diff$new(dataset = kegg, method = "wilcox", group = "treatment.group", 
                        alpha = 0.05, lefse_subgroup = NULL)
t0keggd4$res_diff = subset(t0keggd4$res_diff, Comparison== "Control - Folic")
t0keggd4$res_abund = subset(t0keggd4$res_abund, Group== "Control" | Group=="Folic")
t0keggd4$res_diff = subset(t0keggd4$res_diff, Significance == "***" | 
                             Significance == "**" | Significance == "*")
g0keggd4=t0keggd4$plot_diff_abund(add_sig = T, keep_prefix = F,add_sig_label = "Significance")+
  theme(axis.title.x = element_text(size=10))+scale_fill_manual(values=c("#D95F02","#E7298A"))+
  scale_color_manual(values=c("#D95F02","#E7298A"))
###############################################################################

###################
# FECES AND OTHER COMPARISONS
###################

# repeating the first part

# importing the dataset

feces=qiime2meco(ASV_data="table.qza", phylo_tree="rep_seqs_aligned_masked_tree_rooted.qza", 
                 sample_data="metadata.tsv", taxonomy_data="classification.qza")

# make the taxonomic information unified, very important
feces$tax_table=tidy_taxonomy(feces$tax_table)

# renaming groups
levels(feces$sample_table$treatment.group)=c("FH","FW","NP","CP","Inositol","Folic","NC")
feces$sample_table$treatment.group=relevel(feces$sample_table$treatment.group, "FW")

# selecting only the relevant groups
feces$sample_table = subset(feces$sample_table, treatment.group == "FW" | treatment.group =="FH" | 
                              treatment.group =="NP" | treatment.group =="CP")

# removing ASVs unassigned at kingdom level and archaea
feces$tax_table =subset(feces$tax_table,Kingdom == "k__Bacteria")

#make data consistent
feces$tidy_dataset()


########################
# ABUNDANCE ANALYSIS
########################

# Some sequences unassigned at phylum level with a very low abundance are creating
# an "Others" category. It can be removed:
feces$tax_table=feces$tax_table = subset(feces$tax_table,Phylum != "p__")
#make data consistent again
feces$tidy_dataset()

# calculating abundances at each taxonomic rank (stored in $taxa_abund)
feces$cal_abund()

# PHYLA

# create object with transformed abundance data
# select the 10 most abundant phyla (in this case it doesn't make a difference)
f0phyla_c=trans_abund$new(dataset = feces, taxrank = "Phylum", ntaxa = 10,show = 0)


# create a plot with all samples
gf0phyla_c=f0phyla_c$plot_bar(color_values = RColorBrewer::brewer.pal(8, "Set1"),
                             others_color = "grey70", facet = "treatment.group", 
                             xtext_keep = FALSE, legend_text_italic = FALSE)

# same as before, but with the samples grouped
f0phyla_g = trans_abund$new(dataset = feces, taxrank = "Phylum", ntaxa = 10, groupmean = "treatment.group")
gf0phyla_g = f0phyla_g$plot_bar(color_values = RColorBrewer::brewer.pal(8, "Set1"),
                               others_color = "grey70", legend_text_italic = FALSE) +  
  theme_classic() + theme(axis.title.y = element_text(size = 18))

gf0phyla_g$facet
# GENERA

#selecting top 40
f0genera_40 = trans_abund$new(dataset = feces, taxrank = "Genus", ntaxa = 40)
gf0genera_heatmap=f0genera_40$plot_heatmap(facet = "treatment.group", xtext_keep = FALSE, withmargin = FALSE)

#selecting top 10
f0genera_15_g = trans_abund$new(dataset = feces, taxrank = "Genus", ntaxa = 15,
                                groupmean = "treatment.group")
gf0genera_bar=f0genera_15_g$plot_bar(legend_text_italic = FALSE, 
                                    color_values = c(RColorBrewer::brewer.pal(12, "Paired"),
                                                     RColorBrewer::brewer.pal(8, "Set2")))+  
  theme_classic() + theme(axis.title.y = element_text(size = 18))+guides(fill=guide_legend(ncol=2,title="Genus"))

# VENN PLOT

# creating a new dataset with merged samples into groups
f0merged = feces$merge_samples(use_group = "treatment.group")

# making a venn plot. absolute numbers are ASVs
f0venn = trans_venn$new(f0merged, ratio = "seqratio")
gf0venn=f0venn$plot_venn()

####################
# DIVERSITY ANALYSIS
####################

# ALPHA DIVERSITY

# creating transformed alpha diversity object
f0alpha = trans_alpha$new(dataset = feces, group = "treatment.group")

# test the difference between groups (several methods: KW, wilcox, anova)

# wilcoxon pairwise comparison
f0alpha$cal_diff(method = "wilcox")

# select only significant comparisons
f0alpha$res_diff = base::subset(f0alpha$res_diff, Significance != "ns")

# Observed richness (ASV count)
gf0alpha_observed=f0alpha$plot_alpha(measure="Observed",pair_compare = TRUE)

# No significant, removing layer
gf0alpha_observed[["layers"]][[3]]=NULL

# Chao1 index
gf0alpha_chao1=f0alpha$plot_alpha(measure="Chao1",pair_compare = TRUE)

# sig label fix:
gf0alpha_chao1[["layers"]][[3]]=NULL

# Shannon index
gf0alpha_shannon=f0alpha$plot_alpha(measure="Shannon",pair_compare = TRUE)

# sig label fix:
gf0alpha_shannon[["layers"]][[3]]=NULL


# BETA DIVERSITY

# creating transformed beta diversity object
library(GUniFrac)
feces$cal_betadiv(unifrac = TRUE)

f0beta = trans_beta$new(dataset = feces, group = "treatment.group",measure = "wei_unifrac")

# Pouse PCoA as an example, PCA or NMDS is also available
f0beta$cal_ordination(ordination = "PCoA")

# plotting
gf0beta_PCoA=f0beta$plot_ordination(plot_color = "treatment.group", plot_shape = "treatment.group", 
                                   plot_type = c("point", "ellipse"),ellipse_chull_fill = FALSE)+
  theme_bw()+theme(legend.position = "right")+labs(col="Group",shape="Group")

# calculate and plot sample distances WITHIN groups
f0beta$cal_group_distance()
gf0beta_distance_w=f0beta$plot_group_distance(distance_pair_stat = TRUE,hide_ns=T)

# calculate and plot sample distances BETWEEN groups
f0beta$cal_group_distance(within_group = FALSE)
gf0beta_distance_b=f0beta$plot_group_distance(distance_pair_stat = TRUE,hide_ns=T)

# PERTMANOVA for all groups when manova_all = TRUE
f0beta$cal_manova(manova_all = TRUE)
gf0beta_MANOVA_all=f0beta$res_manova
knitr::kable(gf0beta_MANOVA_all)

# PERMANOVA for each paired groups
f0beta$cal_manova(manova_all = FALSE)
gf0beta_MANOVA_pairwise=f0beta$res_manova
knitr::kable(gf0beta_MANOVA_pairwise)

###########
# DIFF
###########
# WILCOXON

# Phylum
# calculate differential test results (stored in $res_diff)
f0diff_wilcoxp = trans_diff$new(dataset = feces, method = "wilcox", 
                                group = "treatment.group", taxa_level = "Phylum")


# Selecting the significant comparisons
f0diff_wilcoxp$res_diff = subset(f0diff_wilcoxp$res_diff, Significance == "***" | 
                                   Significance == "**" | Significance == "*")

gf0diff_wilcoxp=f0diff_wilcoxp$plot_diff_abund(add_sig = T, keep_prefix = F,
                                              add_sig_label = "Significance",                                              
                                              group_order = c("FW","FH","NP","CP"))

# Selecting only the two comparisons of interest
f0diff_wilcoxp$res_diff = subset(f0diff_wilcoxp$res_diff, Comparison== "NP - FH" |
                                   Comparison== "NP - PF" |  Comparison== "FH - FW" |  
                                   Comparison== "NP - CP")

gf0diff_wilcoxp2=f0diff_wilcoxp$plot_diff_abund(add_sig = T, keep_prefix = F,
                                                add_sig_label = "Significance",
                                                group_order = c("FW","FH","NP","CP"))

# Genera (same as before)
# calculate differential test results (stored in $res_diff)
f0diff_wilcoxg = trans_diff$new(dataset = feces, method = "wilcox", 
                                group = "treatment.group", taxa_level = "Genus")


# Selecting the significant comparisons
f0diff_wilcoxg$res_diff = subset(f0diff_wilcoxg$res_diff, Significance == "***" | 
                                   Significance == "**" | Significance == "*")

gf0diff_wilcoxg=f0diff_wilcoxg$plot_diff_abund(add_sig = T, keep_prefix = F,
                                              add_sig_label = "Significance",
                                              group_order = c("FW","FH","NP","CP"))

# Selecting only the two comparisons of interest
f0diff_wilcoxg$res_diff = subset(f0diff_wilcoxg$res_diff, Comparison== "NP - FH" |
                                   Comparison== "NP - PF" |  Comparison== "FH - FW" |  
                                   Comparison== "NP - CP")

gf0diff_wilcoxg2=f0diff_wilcoxg$plot_diff_abund(add_sig = T, keep_prefix = F,
                                               add_sig_label = "Significance",
                                               group_order = c("FW","FH","NP","CP"))

############
#IMAGE SAVING
############
library(ggpubr)

figure2=ggarrange(ggarrange(ggarrange(gf0phyla_g,gf0genera_bar,labels=c("A","B"),widths=c(1,1.6)),
                            ggarrange(gf0diff_wilcoxp2,gf0diff_wilcoxg2,labels=c("C","D")),nrow=2))

ggsave("Images/figure2.png",figure2,width=3200, height = 2700, units="px")


figure3=ggarrange(ggarrange(gf0alpha_observed,gf0alpha_chao1,
                            gf0alpha_shannon, labels=c("A","B","C"),nrow=1),
                  gf0beta_PCoA,ncol=1, widths=c(1,1),labels=c("","D"))
ggsave("Images/figure3.png",figure3,width=3000, height = 3500, units="px")  

figure4=ggarrange(g0phyla_g,g0genera_bar,g0diff_wilcoxp,g0diff_wilcoxg,labels = "AUTO",
          heights=c(1,1),widths=c(1,1.5))

ggsave("Images/figure4.png",figure4,width=3700, height = 2700, units="px")


figure5=ggarrange(ggarrange(g0alpha_observed,g0alpha_chao1,g0alpha_shannon,labels = "AUTO",ncol=3),
                  g0beta_PCoA,nrow = 2,labels = c("","D"))

ggsave("Images/figure5.png",figure5,width=3000, height = 3500, units="px")


figure6=ggarrange(g0metad3,g0metad4,g0keggd3,g0keggd4,labels = "AUTO",nrow=4,
                  align = "hv",heights=c(1.2,2,2,3))

ggsave("Images/figure6.png",figure6,width=2500, height = 2500, units="px")

              

################################
#Individual plots
# Abundances
ggsave("Images/Phyla.png",g0phyla_g, width=2000, height=2000, units="px")
ggsave("Images/Genera.png",g0genera_bar, width=2000, height=2000, units="px")
ggsave("Images/Heatmap.png",g0genera_heatmap, width=2000, height=2000, units="px")
ggsave("Images/Venn.png",g0venn, width=2000, height=2000, units="px")

# Diversity
ggsave("Images/ADobserved.png",g0alpha_observed, width=2000, height=2000, units="px")
ggsave("Images/ADChao1.png",g0alpha_chao1, width=2000, height=2000, units="px")
ggsave("Images/ADShannon.png",g0alpha_shannon, width=2000, height=2000, units="px")

ggsave("Images/BDPCoA.png",g0beta_PCoA, width=2000, height=2000, units="px")
ggsave("Images/BDWithin.png",g0beta_distance_w, width=2000, height=2000, units="px")
ggsave("Images/BDBetween.png",g0beta_distance_b, width=2000, height=2000, units="px")

# Differential
ggsave("Images/DiffP.png",g0diff_wilcoxp, width=2000, height=2000, units="px")
ggsave("Images/Diffg.png",g0diff_wilcoxg, width=2000, height=2000, units="px")
ggsave("Images/Diffg2.png",g0diff_wilcoxg2, width=2000, height=2000, units="px")

# Functional
ggsave("Images/Pathways.png",g0picd, width=2000, height=2000, units="px")

