library("file2meco")
library("microeco")
library("tidyverse")

###############################################################################

###########################################
# Priliminary comparisons: feces and others
###########################################

# importing the dataset

feces=qiime2meco(ASV_data="table.qza", phylo_tree="rep_seqs_aligned_masked_tree_rooted.qza", 
                 sample_data="metadata.tsv", taxonomy_data="classification.qza")

# make the taxonomic information unified, very important
feces$tax_table=tidy_taxonomy(feces$tax_table)

# renaming groups
levels(feces$sample_table$treatment.group)=c("LP_NP_NS_F","WT_NP_NS_F","LP_NP_NS_C","LP_Pr_NS_C","LP_Pr_DQ_C","LP_Pr_FA_C","NC")
feces$sample_table$treatment.group=relevel(feces$sample_table$treatment.group, "WT_NP_NS_F")

# selecting only the relevant groups
feces$sample_table = subset(feces$sample_table, treatment.group == "WT_NP_NS_F" | treatment.group =="LP_NP_NS_F" | 
                              treatment.group =="LP_NP_NS_C" | treatment.group =="LP_Pr_NS_C")

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
  theme_classic() + theme(axis.title.y = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust=1))


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
  theme_classic() + theme(axis.title.y = element_text(size = 18),
                          axis.text.x = element_text(angle = 45, hjust=1))+
  guides(fill=guide_legend(ncol=2,title="Genus"))

# VENN PLOT

# creating a new dataset with merged samples into groups
f0merged = feces$merge_samples(use_group = "treatment.group")

# making a venn plot. absolute numbers are ASVs
f0venn = trans_venn$new(f0merged, ratio = "seqratio")
gf0venn=f0venn$plot_venn()



#######################
# DIFFERENTIAL ANALYSIS
#######################

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
                                               group_order = c("WT_NP_NS_F","LP_NP_NS_F","LP_NP_NS_C","LP_Pr_NS_C"))+
  theme(legend.position = "bottom")+guides(fill=guide_legend(nrow=2,byrow=T,reverse = T))

# Selecting only the two comparisons of interest
f0diff_wilcoxp$res_diff = subset(f0diff_wilcoxp$res_diff, Comparison== "LP_NP_NS_C - LP_NP_NS_F" |
                                   Comparison== "LP_NP_NS_C - PF" |  Comparison== "LP_NP_NS_F - WT_NP_NS_F" |  
                                   Comparison== "LP_NP_NS_C - LP_Pr_NS_C")

gf0diff_wilcoxp2=f0diff_wilcoxp$plot_diff_abund(add_sig = T, keep_prefix = F,
                                                add_sig_label = "Significance",
                                                group_order = c("WT_NP_NS_F","LP_NP_NS_F","LP_NP_NS_C","LP_Pr_NS_C"))+
  theme(legend.position = "bottom")+guides(fill=guide_legend(nrow=2,byrow=T,reverse = T))

# Genera (same as before)
# calculate differential test results (stored in $res_diff)
f0diff_wilcoxg = trans_diff$new(dataset = feces, method = "wilcox", 
                                group = "treatment.group", taxa_level = "Genus")


# Selecting the significant comparisons
f0diff_wilcoxg$res_diff = subset(f0diff_wilcoxg$res_diff, Significance == "***" | 
                                   Significance == "**" | Significance == "*")

gf0diff_wilcoxg=f0diff_wilcoxg$plot_diff_abund(add_sig = T, keep_prefix = F,
                                               add_sig_label = "Significance",
                                               group_order = c("WT_NP_NS_F","LP_NP_NS_F",
                                                               "LP_NP_NS_C","LP_Pr_NS_C"))+
  theme(legend.position = "bottom")+guides(fill=guide_legend(nrow=2,byrow=T,reverse = T))

# Selecting only the two comparisons of interest
f0diff_wilcoxg$res_diff = subset(f0diff_wilcoxg$res_diff, Comparison== "LP_NP_NS_C - LP_NP_NS_F" |
                                   Comparison== "LP_NP_NS_C - PF" |  Comparison== "LP_NP_NS_F - WT_NP_NS_F" |  
                                   Comparison== "LP_NP_NS_C - LP_Pr_NS_C")

gf0diff_wilcoxg2=f0diff_wilcoxg$plot_diff_abund(add_sig = T, keep_prefix = F,
                                                add_sig_label = "Significance",
                                                group_order = c("WT_NP_NS_F","LP_NP_NS_F","LP_NP_NS_C","LP_Pr_NS_C"))+
  theme(legend.position = "bottom")+guides(fill=guide_legend(nrow=2,byrow=T,reverse = T))


####################
# DIVERSITY ANALYSIS
####################

# ALPHA DIVERSITY

# Rarefaction is performed to have the same sample count
# Check the lowest ammount of sequences and make that the sample size
feces$sample_sums() %>% range
feces$rarefy_samples(sample.size =12662)

# creating transformed alpha diversity object
f0alpha = trans_alpha$new(dataset = feces, group = "treatment.group")

# test the difference between groups (several methods: KW, wilcox, anova)

# wilcoxon pairwise comparison
f0alpha$cal_diff(method = "wilcox")

# select only significant comparisons
f0alpha$res_diff = base::subset(f0alpha$res_diff, Significance != "ns")

# Observed richness (ASV count)
gf0alpha_observed=f0alpha$plot_alpha(measure="Observed",pair_compare = TRUE, xtext_angle = 45)

# No significant, removing layer
gf0alpha_observed[["layers"]][[3]]=NULL

# Chao1 index
gf0alpha_chao1=f0alpha$plot_alpha(measure="Chao1",pair_compare = TRUE, xtext_angle = 45)

# sig label fix:
gf0alpha_chao1[["layers"]][[3]]=NULL

# Shannon index
gf0alpha_shannon=f0alpha$plot_alpha(measure="Shannon",pair_compare = TRUE, xtext_angle = 45)

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

#######################################################################
# MAIN ANALYSIS (SUPPLEMENTATION)
#######################################################################

# importing the dataset
dataset=qiime2meco(ASV_data="table.qza", phylo_tree="rep_seqs_aligned_masked_tree_rooted.qza", 
            sample_data="metadata.tsv", taxonomy_data="classification.qza")

# make the taxonomic information unified, very important
dataset$tax_table=tidy_taxonomy(dataset$tax_table)

# renaming groups
levels(dataset$sample_table$treatment.group)=c("LP_NP_NS_F","WT_NP_NS_F","LP_NP_NS_C",
                                               "LP_Pr_NS_C","LP_Pr_DQ_C","LP_Pr_FA_C","NC")

# reordering factor levels
dataset$sample_table$treatment.group=factor(dataset$sample_table$treatment.group, 
                                            levels=c("LP_NP_NS_F","WT_NP_NS_F","LP_NP_NS_C",
                                                     "LP_Pr_NS_C","LP_Pr_FA_C","LP_Pr_DQ_C","NC"))

# selecting only the relevant groups
dataset$sample_table = subset(dataset$sample_table, treatment.group == "LP_Pr_NS_C"  |
                                treatment.group =="LP_Pr_FA_C" | treatment.group =="LP_Pr_DQ_C")

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
                             xtext_keep = FALSE, legend_text_italic = FALSE,
                             order_facet = c("LP_Pr_NS_C","LP_Pr_FA_C","LP_Pr_DQ_C"))

# same as before, but with the samples grouped
t0phyla_g = trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10, groupmean = "treatment.group")
g0phyla_g = t0phyla_g$plot_bar(color_values = RColorBrewer::brewer.pal(8, "Set1"),
                               others_color = "grey70", legend_text_italic = FALSE,
                               order_x = c("LP_Pr_NS_C","LP_Pr_FA_C","LP_Pr_DQ_C")) +  
  theme_classic() + theme(axis.title.y = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1))


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
                                                     RColorBrewer::brewer.pal(8, "Set2")),
                                    order_x = c("LP_Pr_NS_C","LP_Pr_FA_C","LP_Pr_DQ_C"))+  
  theme_classic() + theme(axis.title.y = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1))+guides(fill=guide_legend(ncol=2,title="Genus"))

# VENN PLOT

# creating a new dataset with merged samples into groups
t0merged = dataset$merge_samples(use_group = "treatment.group")

# making a venn plot. absolute numbers are ASVs
t0venn = trans_venn$new(t0merged, ratio = "seqratio")
g0venn=t0venn$plot_venn()

########################
# DIFFERENTIAL ANALYSIS
########################

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
                                              group_order = c("LP_Pr_NS_C","LP_Pr_FA_C","LP_Pr_DQ_C"))+
  theme(legend.position = "bottom")+guides(fill=guide_legend(nrow=2,byrow=T,reverse = T))

#Selecting only the two comparisons of interest
t0diff_wilcoxp$res_diff = subset(t0diff_wilcoxp$res_diff, 
                                 Comparison== "LP_Pr_NS_C - LP_Pr_FA_C" |  
                                   Comparison== "LP_Pr_NS_C - LP_Pr_DQ_C")

g0diff_wilcoxp2=t0diff_wilcoxp$plot_diff_abund(add_sig = T, keep_prefix = F,
                                               add_sig_label = "Significance",
                                               group_order = c("LP_Pr_NS_C","LP_Pr_FA_C","LP_Pr_DQ_C"))+
  theme(legend.position = "bottom")+guides(fill=guide_legend(nrow=2,byrow=T,reverse = T))

# Genera (same as before)
# calculate differential test results (stored in $res_diff)
t0diff_wilcoxg = trans_diff$new(dataset = dataset, method = "wilcox", 
                                group = "treatment.group", taxa_level = "Genus")


# Selecting the significant comparisons
t0diff_wilcoxg$res_diff = subset(t0diff_wilcoxg$res_diff, Significance == "***" | 
                                   Significance == "**" | Significance == "*")

g0diff_wilcoxg=t0diff_wilcoxg$plot_diff_abund(add_sig = T, keep_prefix = F,
                                              add_sig_label = "Significance",
                                              group_order = c("LP_Pr_NS_C","LP_Pr_FA_C","LP_Pr_DQ_C"))+
  theme(legend.position = "bottom")+guides(fill=guide_legend(nrow=2,byrow=T,reverse = T))

#Selecting only the two comparisons of interest
t0diff_wilcoxg$res_diff = subset(t0diff_wilcoxg$res_diff, 
                                 Comparison== "LP_Pr_NS_C - LP_Pr_FA_C" |  
                                   Comparison== "LP_Pr_NS_C - LP_Pr_DQ_C")

g0diff_wilcoxg2=t0diff_wilcoxg$plot_diff_abund(add_sig = T, keep_prefix = F,
                                               add_sig_label = "Significance",
                                               group_order = c("LP_Pr_NS_C","LP_Pr_FA_C","LP_Pr_DQ_C"))+
  theme(legend.position = "bottom")+guides(fill=guide_legend(nrow=2,byrow=T,reverse = T))



####################
#FUNCTIONAL ANALYSIS
#####################

# picrust2 was used to obtain the pathways from kegg and metacyc

# Loading the MetaCyc pathways obtained in picrust
pathway_table = read.delim("path_abun_unstrat_1.tsv", row.names = 1)

# Rewriting names to match
sampletable=dataset$sample_table
rownames(sampletable)=gsub("-",".",rownames(dataset$sample_table))

# Loading MetaCyc data
data("MetaCyc_pathway_map")
picrust = microtable$new(otu_table = pathway_table, tax_table = MetaCyc_pathway_map,
                         sample_table = sampletable)
picrust$tidy_dataset()
picrust$cal_abund()

# making a trans_abund object
t0pic = trans_abund$new(picrust, taxrank = "pathway", groupmean = "treatment.group")
t0pic$plot_bar(legend_text_italic = FALSE)

# Making statistical test
t0picd <- trans_diff$new(dataset = picrust, method = "wilcox", group = "treatment.group", alpha = 0.05, lefse_subgroup = NULL)
g0picd1=t0picd$plot_diff_bar()


# Making statistycal tests: Pairwise


# Comparison 1: pregnant vs pregnant folic
t0metacycd4=trans_diff$new(dataset = picrust, method = "wilcox", group = "treatment.group", 
                           alpha = 0.05, lefse_subgroup = NULL)
t0metacycd4$res_diff = subset(t0metacycd4$res_diff, Comparison== "LP_Pr_NS_C - LP_Pr_FA_C")
t0metacycd4$res_abund = subset(t0metacycd4$res_abund, Group== "LP_Pr_NS_C" | Group=="LP_Pr_FA_C")
t0metacycd4$res_diff = subset(t0metacycd4$res_diff, Significance == "***" | 
                                Significance == "**" | Significance == "*")
g0metad4=t0metacycd4$plot_diff_abund(add_sig = T, keep_prefix = F,add_sig_label = "Significance")+
  theme(axis.title.x = element_text(size=10))+scale_fill_manual(values=c("#D95F02","#1B9E77"))+
  scale_color_manual(values=c("#D95F02","#1B9E77"))

# Comparison 2: pregnant vs pregnant+ chiro
t0metacycd3=trans_diff$new(dataset = picrust, method = "wilcox", group = "treatment.group", 
                           alpha = 0.05, lefse_subgroup = NULL)
t0metacycd3$res_diff = subset(t0metacycd3$res_diff, Comparison== "LP_Pr_NS_C - LP_Pr_DQ_C")
t0metacycd3$res_abund = subset(t0metacycd3$res_abund, Group== "LP_Pr_NS_C" | Group=="LP_Pr_DQ_C")
t0metacycd3$res_diff = subset(t0metacycd3$res_diff, Significance == "***" | 
                                Significance == "**" | Significance == "*")
g0metad3=t0metacycd3$plot_diff_abund(add_sig = T, keep_prefix = F,add_sig_label = "Significance")+
  theme(axis.title.x = element_text(size=10))+scale_fill_manual(values=c("#7570B3","#1B9E77"))+
  scale_color_manual(values=c("#7570B3","#1B9E77"))



# USING KEGG

# Loading KEGG pathways obtained in picrust2
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

# Comparison 1: pregnant vs pregnant folic
t0keggd4=trans_diff$new(dataset = kegg, method = "wilcox", group = "treatment.group", 
                        alpha = 0.05, lefse_subgroup = NULL)
t0keggd4$res_diff = subset(t0keggd4$res_diff, Comparison== "LP_Pr_NS_C - LP_Pr_FA_C")
t0keggd4$res_abund = subset(t0keggd4$res_abund, Group== "LP_Pr_NS_C" | Group=="LP_Pr_FA_C")
t0keggd4$res_diff = subset(t0keggd4$res_diff, Significance == "***" | 
                             Significance == "**" | Significance == "*")
g0keggd4=t0keggd4$plot_diff_abund(add_sig = T, keep_prefix = F,add_sig_label = "Significance")+
  theme(axis.title.x = element_text(size=10))+scale_fill_manual(values=c("#D95F02","#1B9E77"))+
  scale_color_manual(values=c("#D95F02","#1B9E77"))

# Comparison 2: pregnant vs pregnant+ chiro
t0keggd3=trans_diff$new(dataset = kegg, method = "wilcox", group = "treatment.group", 
                        alpha = 0.05, lefse_subgroup = NULL)
t0keggd3$res_diff = subset(t0keggd3$res_diff, Comparison== "LP_Pr_NS_C - LP_Pr_DQ_C")
t0keggd3$res_abund = subset(t0keggd3$res_abund, Group== "LP_Pr_NS_C" | Group=="LP_Pr_DQ_C")
t0keggd3$res_diff = subset(t0keggd3$res_diff, Significance == "***" | 
                             Significance == "**" | Significance == "*")
g0keggd3=t0keggd3$plot_diff_abund(add_sig = T, keep_prefix = F,add_sig_label = "Significance")+
  theme(axis.title.x = element_text(size=10))+scale_fill_manual(values=c("#7570B3","#1B9E77"))+
  scale_color_manual(values=c("#7570B3","#1B9E77"))


####################
# DIVERSITY ANALYSIS
####################

# ALPHA DIVERSITY

# Rarefaction is performed to have the same sample count
# Check the lowest ammount of sequences and make that the sample size
dataset$sample_sums() %>% range
dataset$rarefy_samples(sample.size =19287)

# creating transformed alpha diversity object
t0alpha = trans_alpha$new(dataset = dataset, group = "treatment.group")

# test the difference between groups (several methods: KW, wilcox, anova)

# wilcoxon pairwise comparison
t0alpha$cal_diff(method = "wilcox")

# select only significant comparisons
t0alpha$res_diff = base::subset(t0alpha$res_diff, Significance != "ns")

# Observed richness (ASV count)
g0alpha_observed=t0alpha$plot_alpha(measure="Observed",pair_compare = TRUE, xtext_angle = 45)

# for some reason the significant layers are displaced. package bug? fix:
g0alpha_observed[["layers"]][[3]][["stat_params"]][["xmin"]]=g0alpha_observed[["layers"]][[3]][["stat_params"]][["xmin"]]-3
g0alpha_observed[["layers"]][[3]][["stat_params"]][["xmax"]]=g0alpha_observed[["layers"]][[3]][["stat_params"]][["xmax"]]-3


# Chao1 index
g0alpha_chao1=t0alpha$plot_alpha(measure="Chao1",pair_compare = TRUE, xtext_angle = 45)

# sig label fix:
g0alpha_chao1[["layers"]][[3]][["stat_params"]][["xmin"]]=g0alpha_chao1[["layers"]][[3]][["stat_params"]][["xmin"]]-3
g0alpha_chao1[["layers"]][[3]][["stat_params"]][["xmax"]]=g0alpha_chao1[["layers"]][[3]][["stat_params"]][["xmax"]]-3


# Shannon index
g0alpha_shannon=t0alpha$plot_alpha(measure="Shannon",pair_compare = TRUE, xtext_angle = 45)

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
  theme_bw()+guides(fill=guide_legend())+labs(col="Group",shape="Group")

# calculate and plot sample distances WITHIN groups
t0beta$cal_group_distance()
g0beta_distance_w=t0beta$plot_group_distance(distance_pair_stat = TRUE,hide_ns=T)

# calculate and plot sample distances BETWEEN groups
t0beta$cal_group_distance(within_group = FALSE)
g0beta_distance_b=t0beta$plot_group_distance(distance_pair_stat = TRUE,hide_ns=T)

# Statistical tests for beta diversity
# PERMANOVA for all groups when manova_all = TRUE
t0beta$cal_manova(manova_all = TRUE)
g0beta_MANOVA_all=t0beta$res_manova
knitr::kable(g0beta_MANOVA_all)

# PERMANOVA for each paired groups
t0beta$cal_manova(manova_all = FALSE)
g0beta_MANOVA_pairwise=t0beta$res_manova
knitr::kable(g0beta_MANOVA_pairwise)


###############################################################################

############
#IMAGE SAVING
############

# Using ggpubr package to plot the images together with labels

library(ggpubr)

figure2=ggarrange(ggarrange(gf0phyla_g,gf0genera_bar,labels=c("A","B"),widths=c(1,1.8)),
                  ggarrange(gf0diff_wilcoxp2,gf0diff_wilcoxg2,labels=c("C","D"),widths = c(1.1,1.4)),nrow=2)

ggsave("Images/figure2.png",figure2,width=3000, height = 2500, units="px")


figure3=ggarrange(ggarrange(gf0alpha_observed,gf0alpha_chao1,
                            gf0alpha_shannon, labels=c("A","B","C"),nrow=1),
                  gf0beta_PCoA,ncol=1, widths=c(1,1),labels=c("","D"))
ggsave("Images/figure3.png",figure3,width=3000, height = 3500, units="px")  

figure4=ggarrange(ggarrange(g0phyla_g,g0genera_bar,labels=c("A","B"),widths=c(1,1.6)),
                  ggarrange(g0diff_wilcoxp2,g0diff_wilcoxg2,labels = c("C","D"),widths=c(1.1,1.6)),nrow=2)

ggsave("Images/figure4.png",figure4,width=3000, height = 2700, units="px")


figure5=ggarrange(ggarrange(g0alpha_observed,g0alpha_chao1,g0alpha_shannon,labels = "AUTO",ncol=3),
                  g0beta_PCoA,nrow = 2,labels = c("","D"))

ggsave("Images/figure5.png",figure5,width=3000, height = 3500, units="px")


figure6=ggarrange(g0metad4,g0keggd4,g0metad3,g0keggd3,labels = "AUTO",nrow=4,
                  align = "hv",heights=c(2,3,1.2,2))

ggsave("Images/figure6.png",figure6,width=2700, height = 2500, units="px")

              
