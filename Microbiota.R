
#Import from QIIME2
library(qiime2R)
qiime=qza_to_phyloseq("../data/table.qza", tree="..data/rep_seqs_aligned_masked_tree_rooted.qza", 
                      metadata="../data/metadata0.tsv", taxonomy="../data/classification.qza")

# Removing the negative control and the unassigned ASVs at phylum level
library(phyloseq)
qiimed=prune_samples(sample_names(qiime)!="CN-16S",qiime)
qiimed=subset_taxa(qiimed,Phylum!="NA")

# Rarefaction curve
library(vegan)
rarecurve=rarecurve(t(otu_table(qiimed)), step=50, cex=0.5)

library(tidyverse)
rare=map_dfr(rarecurve,bind_rows) %>% 
  bind_cols(sample=rownames(qiimed@sam_data),.) %>%
  bind_cols(Group=qiimed@sam_data$treatment.group,.) %>%
  pivot_longer(-c(sample,Group))%>%
  drop_na() %>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name) 

p0rare=ggplot(data=rare,aes(x=n_seqs,y=value,group=sample,colour=Group))+
  geom_line(size=1.2)+facet_wrap(~Group,nrow=2)+
  labs(x= "Sequencing depth",y="Number of ASVs")+theme_bw()+
  theme(strip.background = element_blank()) 

# RELATIVE ABUNDANCE PLOTS

# Calculating relative abundance
qiimerel=transform_sample_counts(qiimed, function(x) x/sum(x))

#######################
# PHYLA COMPARISON
#######################

# COMPLETE-ALL SAMPLES
p0phylac=plot_bar(qiimerel, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(qiimerel@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",nrow=2) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_brewer(palette="Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank(),panel.grid.major = element_blank(),
        strip.background = element_blank())+ylab("Relative abundance") #Width 856

# GROUPED
p0phylag=plot_bar(transform_sample_counts(qiimerel, function(x) x/9), x="treatment.group",fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")+
  guides(color=guide_legend(ncol=1))+
  scale_color_brewer(palette="Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank(),panel.grid.major = element_blank(),
        strip.background = element_blank())+ylab("Relative abundance") 

#Plotting phyla side by side
qiime_phyla_abundance=psmelt(tax_glom(qiimerel,taxrank="Phylum"))[,c(1,2,3,5,14)]

p0phylas=ggplot(qiime_phyla_abundance,aes(x=reorder(Phylum, Abundance, FUN = median,decreasing=T), 
                                          y=Abundance,fill=treatment.group)) +
  geom_boxplot()+ theme_bw()+theme(axis.title.x = element_blank())+ 
  labs(col='Group') +ylab("Relative abundance")+
  geom_point(position = position_dodge(width=0.75),color="black")


#######################
# TOP 10 GENUS COMPARISON
#######################

qiimetop=tax_glom(qiimerel, taxrank = "Genus")

leasttaxa=intersect(intersect(intersect(intersect(intersect(names(sort(taxa_sums(subset_samples(qiimetop,treatment.group=="Feces WT")), 
                                                                       decreasing = TRUE))[-c(1:10)],
                                                            names(sort(taxa_sums(subset_samples(qiimetop,treatment.group=="Feces HET")), 
                                                                       decreasing = TRUE))[-c(1:10)]),
                                                  names(sort(taxa_sums(subset_samples(qiimetop,treatment.group=="HET-non pregnant")), 
                                                             decreasing = TRUE))[-c(1:10)]),
                                        names(sort(taxa_sums(subset_samples(qiimetop,treatment.group=="HET-pregnant")), 
                                                   decreasing = TRUE))[-c(1:10)]),
                              names(sort(taxa_sums(subset_samples(qiimetop,treatment.group=="HET-pregnant chiro")), 
                                         decreasing = TRUE))[-c(1:10)]),
                    names(sort(taxa_sums(subset_samples(qiimetop,treatment.group=="HET-pregnant folic")), 
                               decreasing = TRUE))[-c(1:10)])

# Merging all other taxa

qiimetop@tax_table[leasttaxa,"Genus"]=c("0ther")

# COMPLETE
p0topc=plot_bar(transform_sample_counts(qiimetop, function(x) x/sum(x)), fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(qiimetop@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=3) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(16))+
  scale_fill_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(16))+
  theme_bw()+theme(axis.title.x = element_blank(), panel.background = element_blank(),
                   axis.text.x = element_text(angle = 90), panel.grid.major = element_blank(),
                   strip.background = element_blank())+ylab("Relative abundance")

# GROUPED
p0topg=plot_bar(transform_sample_counts(qiimetop, function(x) x/(sum(x)*9)), x="treatment.group",fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(qiimetop@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=6) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(16))+
  scale_fill_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(16))+
  theme(axis.title.x = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(),strip.background = element_blank(),
        axis.text.x = element_text(angle = 0))+ylab("Relative abundance")

# plotting all 
qiime_top_abundance=psmelt(tax_glom(qiimetop,taxrank="Genus"))[,c(1,2,3,5,18)]

p0tops=ggplot(qiime_top_abundance,aes(x=reorder(Genus, Abundance, FUN = median,decreasing=T), 
                                      y=Abundance,fill=treatment.group)) +
  geom_boxplot()+ theme_bw()+theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 90))+ 
  labs(col='Group') +ylab("Relative abundance")+
  geom_point(position = position_dodge(width=0.75),color="black")


#################
#Alpha diversity
#################
library(ggpubr)
p0alpha=plot_richness(qiimed, x="treatment.group", measures=c("Observed","Chao1","Shannon"), 
                      color="treatment.group")+ theme(legend.position="none")+
  geom_boxplot()+ geom_point(color="black") + theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
        strip.background = element_blank(),axis.text.x = element_text(angle = 90,hjust=1))+
  stat_compare_means(aes(label = paste0(..p.format.., "\n", ..p.signif..)),
                     vjust=1,hide.ns = T,show.legend=FALSE,
                     comparisons = list(c("Feces HET", "Feces WT"),
                                        c("Feces HET", "HET-non pregnant"),
                                        c("HET-non pregnant", "HET-pregnant"),
                                        c("HET-pregnant", "HET-pregnant folic"),
                                        c("HET-pregnant", "HET-pregnant chiro")))



#################
# Beta diversity
#################
qiimedist = phyloseq::distance(qiimed, method="unifrac", weighted=T)

qiimeperm=adonis2(qiimedist ~ sample_data(qiimed)$treatment.group,data = data.frame(sample_data(qiimed)))

p0beta=plot_ordination(qiimed, ordinate(qiimed, method="PCoA", distance=qiimedist), 
                       color="treatment.group") +theme_bw()+ theme(aspect.ratio=1)+ 
  stat_ellipse(aes(group=treatment.group))+  stat_ellipse(aes(group=treatment.group))+
  annotate("text", x = -Inf, y = Inf, hjust = -3.2, vjust = 1.1, #Adjust to display in plot
           label= paste0("PERMANOVA", "\n", "p-value: ",qiimeperm$`Pr(>F)`[1]))

################################################################################


###################
# COMPARISON 1 FECES HET VS FECES WT
###################

# Subsetting
sub1=subset_samples(qiimed,treatment.group=="Feces WT" | treatment.group=="Feces HET")
sub1rel=subset_samples(qiimerel,treatment.group=="Feces WT" | treatment.group=="Feces HET")

# Comparing phyla abundance & statistical tests
sub1pphyla=ggplot(qiime_phyla_abundance[which(qiime_phyla_abundance$treatment.group==c("Feces WT", "Feces HET")),],
                  aes(x=reorder(Phylum, Abundance, FUN = median,decreasing=T), y=Abundance, color=treatment.group)) +
  geom_boxplot()+ theme_bw()+theme(axis.title.x = element_blank(),
                                   axis.text.x = element_text(angle = 90,hjust=1,vjust=0.4))+ 
  labs(col='Group') +ylab("Relative abundance (%)")+
  stat_compare_means(aes(label = paste0(..p.format.., "\n", ..p.signif..)),
                     vjust=1,hide.ns = T,show.legend=FALSE)+ylim(0,1)
# Some phyla were undetected in some samples


# Beta diversity
sub1dist = phyloseq::distance(sub1, method="unifrac", weighted=T)

sub1perm=adonis2(sub1dist ~ sample_data(sub1)$treatment.group,data = data.frame(sample_data(sub1)))

sub1pbeta=plot_ordination(sub1, ordinate(sub1, method="PCoA", distance=sub1dist), 
                          color="treatment.group") +theme_bw()+ theme(aspect.ratio=1)+ 
  stat_ellipse(aes(group=treatment.group))+  stat_ellipse(aes(group=treatment.group))+
  annotate("text", x = -Inf, y = Inf, hjust = -3.2, vjust = 1.1, #Adjust to display in plot
           label= paste0("PERMANOVA", "\n", "p-value: ",sub1perm$`Pr(>F)`[1]))

# Calculating top 10 ASVs

sub1top=tax_glom(sub1rel, taxrank = "Genus")

sub1leasttaxa=intersect(names(sort(taxa_sums(subset_samples(sub1top,treatment.group=="Feces WT")), 
                                   decreasing = TRUE))[-c(1:10)],
                        names(sort(taxa_sums(subset_samples(sub1top,treatment.group=="Feces HET")), 
                                   decreasing = TRUE))[-c(1:10)])
# Merging all other taxa

sub1top@tax_table[sub1leasttaxa,"Genus"]=c("0ther")

length(unique(sub1rel@tax_table[,"Genus"]))

# COMPLETE
sub1ptopc=plot_bar(transform_sample_counts(sub1top, function(x) x/sum(x)), fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(sub1top@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(12))+
  scale_fill_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(12))+
  theme_bw()+theme(axis.title.x = element_blank(), panel.background = element_blank(),
                   axis.text.x = element_text(angle = 90), panel.grid.major = element_blank(),
                   strip.background = element_blank())+ylab("Relative abundance")

# GROUPED
sub1ptopg=plot_bar(transform_sample_counts(sub1top, function(x) x/(sum(x)*9)), x="treatment.group",fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(sub1top@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(12))+
  scale_fill_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(12))+
  theme(axis.title.x = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(),strip.background = element_blank(),
        axis.text.x = element_text(angle = 0))+ylab("Relative abundance")

# Comparing genera abundance & statistical tests
sub1ptops=ggplot(qiime_top_abundance[which(qiime_top_abundance$treatment.group==c("Feces WT", "Feces HET")),],
                 aes(x=reorder(Genus, Abundance, FUN = median,decreasing=T), y=Abundance, color=treatment.group)) +
  geom_boxplot()+ theme_bw()+theme(axis.title.x = element_blank(),
                                   axis.text.x = element_text(angle = 90,hjust=1,vjust=0.4))+ 
  labs(col='Group') +ylab("Relative abundance")+
  stat_compare_means(aes(label = paste0(..p.format.., "\n", ..p.signif..)),
                     vjust=1,hide.ns = T,show.legend=FALSE)+ylim(0,1)





# Building a Heat tree
library(metacoder)
sub1metac=parse_phyloseq(sub1, class_regex = "(.*)", class_key = "taxon_name")

# Calculate relative abundance and per-taxon
sub1metac$data$rel_abd <- calc_obs_props(sub1metac, "otu_table", other_cols = T)
sub1metac$data$tax_rel_abd <- calc_taxon_abund(sub1metac, "rel_abd")

sub1metac$data$diff_table <- compare_groups(sub1metac,
                                            data = "tax_rel_abd",
                                            cols = sub1metac$data$sample_data$sample_id, # What columns of sample data to use
                                            groups = sub1metac$data$sample_data$treatment.group) # What category each sample is assigned to

#Correction for multiple tests
sub1metac$data$diff_table$wilcox_p_value <- p.adjust(sub1metac$data$diff_table$wilcox_p_value,
                                                     method = "fdr")

range(sub1metac$data$diff_table$log2_median_ratio, finite = TRUE) 

#Set non-significant values to zero
sub1metac$data$diff_table$log2_median_ratio[sub1metac$data$diff_table$wilcox_p_value > 0.05] <- 0

# plot heat tree matrix
set.seed(1)
sub1pht=heat_tree(sub1metac,
                  node_label = taxon_names,
                  node_size = n_obs,
                  node_color = log2_median_ratio, 
                  node_color_interval = c(-3, 3),
                  node_color_range = diverging_palette(),
                  node_size_axis_label = "ASV count",
                  node_color_axis_label = "Log2 ratio of median counts",
                  layout = "davidson-harel", 
                  initial_layout = "reingold-tilford") 

ggsave("../images/sub1pphyla.png",sub1pphyla)
ggsave("../images/sub1ptopc.png",sub1ptopc)
ggsave("../images/sub1ptopg.png",sub1ptopg)
ggsave("../images/sub1ptops.png",sub1ptops)
ggsave("../images/sub1pbeta.png",sub1pbeta)
ggsave("../images/sub1pht.png",sub1pht)

####################################################

###################
# COMPARISON 2 FECES HET VS HET-non pregnant
###################

# Subsetting
sub2=subset_samples(qiimed,treatment.group=="Feces HET" | treatment.group=="HET-non pregnant")
sub2rel=subset_samples(qiimerel,treatment.group=="Feces HET" | treatment.group=="HET-non pregnant")

# Comparing phyla abundance & statistical tests
sub2pphyla=ggplot(qiime_phyla_abundance[which(qiime_phyla_abundance$treatment.group==c("Feces HET", "HET-non pregnant")),],
                  aes(x=reorder(Phylum, Abundance, FUN = median,decreasing=T), y=Abundance, color=treatment.group)) +
  geom_boxplot()+ theme_bw()+theme(axis.title.x = element_blank(),
                                   axis.text.x = element_text(angle = 90,hjust=1,vjust=0.4))+ 
  labs(col='Group') +ylab("Relative abundance (%)")+
  stat_compare_means(aes(label = paste0(..p.format.., "\n", ..p.signif..)),
                     vjust=1,hide.ns = T,show.legend=FALSE)+ylim(0,1)
# Some phyla were undetected in some samples


# Beta diversity
sub2dist = phyloseq::distance(sub2, method="unifrac", weighted=T)

sub2perm=adonis2(sub2dist ~ sample_data(sub2)$treatment.group,data = data.frame(sample_data(sub2)))

sub2pbeta=plot_ordination(sub2, ordinate(sub2, method="PCoA", distance=sub2dist), 
                          color="treatment.group") +theme_bw()+ theme(aspect.ratio=1)+ 
  stat_ellipse(aes(group=treatment.group))+  stat_ellipse(aes(group=treatment.group))+
  annotate("text", x = -Inf, y = Inf, hjust = -3.2, vjust = 1.1, #Adjust to display in plot
           label= paste0("PERMANOVA", "\n", "p-value: ",sub2perm$`Pr(>F)`[1]))

# Calculating top 10 ASVs

sub2top=tax_glom(sub2rel, taxrank = "Genus")

sub2leasttaxa=intersect(names(sort(taxa_sums(subset_samples(sub2top,treatment.group=="Feces HET")), 
                                   decreasing = TRUE))[-c(1:10)],
                        names(sort(taxa_sums(subset_samples(sub2top,treatment.group=="HET-non pregnant")), 
                                   decreasing = TRUE))[-c(1:10)])
# Merging all other taxa

sub2top@tax_table[sub2leasttaxa,"Genus"]=c("0ther")

# COMPLETE
sub2ptopc=plot_bar(transform_sample_counts(sub2top, function(x) x/sum(x)), fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(sub2top@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(13))+
  scale_fill_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(13))+
  theme_bw()+theme(axis.title.x = element_blank(), panel.background = element_blank(),
                   axis.text.x = element_text(angle = 90), panel.grid.major = element_blank(),
                   strip.background = element_blank())+ylab("Relative abundance")

# GROUPED
sub2ptopg=plot_bar(transform_sample_counts(sub2top, function(x) x/(sum(x)*9)), x="treatment.group",fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(sub2top@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(13))+
  scale_fill_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(13))+
  theme(axis.title.x = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(),strip.background = element_blank(),
        axis.text.x = element_text(angle = 0))+ylab("Relative abundance")

# Comparing genera abundance & statistical tests
sub2ptops=ggplot(qiime_top_abundance[which(qiime_top_abundance$treatment.group==c("Feces HET", "HET-non pregnant")),],
                 aes(x=reorder(Genus, Abundance, FUN = median,decreasing=T), y=Abundance, color=treatment.group)) +
  geom_boxplot()+ theme_bw()+theme(axis.title.x = element_blank(),
                                   axis.text.x = element_text(angle = 90,hjust=1,vjust=0.4))+ 
  labs(col='Group') +ylab("Relative abundance")+
  stat_compare_means(aes(label = paste0(..p.format.., "\n", ..p.signif..)),
                     vjust=1,hide.ns = T,show.legend=FALSE)+ylim(0,1)





# Building a Heat tree
library(metacoder)
sub2metac=parse_phyloseq(sub2, class_regex = "(.*)", class_key = "taxon_name")

# Calculate relative abundance and per-taxon
sub2metac$data$rel_abd <- calc_obs_props(sub2metac, "otu_table", other_cols = T)
sub2metac$data$tax_rel_abd <- calc_taxon_abund(sub2metac, "rel_abd")

sub2metac$data$diff_table <- compare_groups(sub2metac,
                                            data = "tax_rel_abd",
                                            cols = sub2metac$data$sample_data$sample_id, # What columns of sample data to use
                                            groups = sub2metac$data$sample_data$treatment.group) # What category each sample is assigned to

#Correction for multiple tests
sub2metac$data$diff_table$wilcox_p_value <- p.adjust(sub2metac$data$diff_table$wilcox_p_value,
                                                     method = "fdr")

range(sub2metac$data$diff_table$log2_median_ratio, finite = TRUE) 

#Set non-significant values to zero
sub2metac$data$diff_table$log2_median_ratio[sub2metac$data$diff_table$wilcox_p_value > 0.05] <- 0

# plot heat tree matrix
set.seed(1)
sub2pht=heat_tree(sub2metac,
                  node_label = taxon_names,
                  node_size = n_obs,
                  node_color = log2_median_ratio, 
                  node_color_interval = c(-3, 3),
                  node_color_range = diverging_palette(),
                  node_size_axis_label = "ASV count",
                  node_color_axis_label = "Log2 ratio of median counts",
                  layout = "davidson-harel", 
                  initial_layout = "reingold-tilford") 

ggsave("../images/sub2pphyla.png",sub2pphyla)
ggsave("../images/sub2ptopc.png",sub2ptopc)
ggsave("../images/sub2ptopg.png",sub2ptopg)
ggsave("../images/sub2ptops.png",sub2ptops)
ggsave("../images/sub2pbeta.png",sub2pbeta)
ggsave("../images/sub2pht.png",sub2pht)


####################################################

###################
# COMPARISON 3 HET-non pregnant VS HET-pregnant
###################

# Subsetting
sub3=subset_samples(qiimed,treatment.group=="HET-non pregnant" | treatment.group=="HET-pregnant")
sub3rel=subset_samples(qiimerel,treatment.group=="HET-non pregnant" | treatment.group=="HET-pregnant")

# Comparing phyla abundance & statistical tests
sub3pphyla=ggplot(qiime_phyla_abundance[which(qiime_phyla_abundance$treatment.group==c("HET-non pregnant", "HET-pregnant")),],
                  aes(x=reorder(Phylum, Abundance, FUN = median,decreasing=T), y=Abundance, color=treatment.group)) +
  geom_boxplot()+ theme_bw()+theme(axis.title.x = element_blank(),
                                   axis.text.x = element_text(angle = 90,hjust=1,vjust=0.4))+ 
  labs(col='Group') +ylab("Relative abundance (%)")+
  stat_compare_means(aes(label = paste0(..p.format.., "\n", ..p.signif..)),
                     vjust=1,hide.ns = T,show.legend=FALSE)+ylim(0,1)
# Some phyla were undetected in some samples


# Beta diversity
sub3dist = phyloseq::distance(sub3, method="unifrac", weighted=T)

sub3perm=adonis2(sub3dist ~ sample_data(sub3)$treatment.group,data = data.frame(sample_data(sub3)))

sub3pbeta=plot_ordination(sub3, ordinate(sub3, method="PCoA", distance=sub3dist), 
                          color="treatment.group") +theme_bw()+ theme(aspect.ratio=1)+ 
  stat_ellipse(aes(group=treatment.group))+  stat_ellipse(aes(group=treatment.group))+
  annotate("text", x = -Inf, y = Inf, hjust = -3.2, vjust = 1.1, #Adjust to display in plot
           label= paste0("PERMANOVA", "\n", "p-value: ",sub3perm$`Pr(>F)`[1]))

# Calculating top 10 ASVs

sub3top=tax_glom(sub3rel, taxrank = "Genus")

sub3leasttaxa=intersect(names(sort(taxa_sums(subset_samples(sub3top,treatment.group=="HET-non pregnant")), 
                                   decreasing = TRUE))[-c(1:10)],
                        names(sort(taxa_sums(subset_samples(sub3top,treatment.group=="HET-pregnant")), 
                                   decreasing = TRUE))[-c(1:10)])
# Merging all other taxa

sub3top@tax_table[sub3leasttaxa,"Genus"]=c("0ther")

# COMPLETE
sub3ptopc=plot_bar(transform_sample_counts(sub3top, function(x) x/sum(x)), fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(sub3top@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(13))+
  scale_fill_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(13))+
  theme_bw()+theme(axis.title.x = element_blank(), panel.background = element_blank(),
                   axis.text.x = element_text(angle = 90), panel.grid.major = element_blank(),
                   strip.background = element_blank())+ylab("Relative abundance")

# GROUPED
sub3ptopg=plot_bar(transform_sample_counts(sub3top, function(x) x/(sum(x)*9)), x="treatment.group",fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(sub3top@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(13))+
  scale_fill_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(13))+
  theme(axis.title.x = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(),strip.background = element_blank(),
        axis.text.x = element_text(angle = 0))+ylab("Relative abundance")

# Comparing genera abundance & statistical tests
sub3ptops=ggplot(qiime_top_abundance[which(qiime_top_abundance$treatment.group==c("HET-non pregnant", "HET-pregnant")),],
                 aes(x=reorder(Genus, Abundance, FUN = median,decreasing=T), y=Abundance, color=treatment.group)) +
  geom_boxplot()+ theme_bw()+theme(axis.title.x = element_blank(),
                                   axis.text.x = element_text(angle = 90,hjust=1,vjust=0.4))+ 
  labs(col='Group') +ylab("Relative abundance")+
  stat_compare_means(aes(label = paste0(..p.format.., "\n", ..p.signif..)),
                     vjust=1,hide.ns = T,show.legend=FALSE)+ylim(0,1)





# Building a Heat tree
library(metacoder)
sub3metac=parse_phyloseq(sub3, class_regex = "(.*)", class_key = "taxon_name")

# Calculate relative abundance and per-taxon
sub3metac$data$rel_abd <- calc_obs_props(sub3metac, "otu_table", other_cols = T)
sub3metac$data$tax_rel_abd <- calc_taxon_abund(sub3metac, "rel_abd")

sub3metac$data$diff_table <- compare_groups(sub3metac,
                                            data = "tax_rel_abd",
                                            cols = sub3metac$data$sample_data$sample_id, # What columns of sample data to use
                                            groups = sub3metac$data$sample_data$treatment.group) # What category each sample is assigned to

#Correction for multiple tests
sub3metac$data$diff_table$wilcox_p_value <- p.adjust(sub3metac$data$diff_table$wilcox_p_value,
                                                     method = "fdr")

range(sub3metac$data$diff_table$log2_median_ratio, finite = TRUE) 

#Set non-significant values to zero
sub3metac$data$diff_table$log2_median_ratio[sub3metac$data$diff_table$wilcox_p_value > 0.05] <- 0

# plot heat tree matrix
set.seed(1)
sub3pht=heat_tree(sub3metac,
                  node_label = taxon_names,
                  node_size = n_obs,
                  node_color = log2_median_ratio, 
                  node_color_interval = c(-3, 3),
                  node_color_range = diverging_palette(),
                  node_size_axis_label = "ASV count",
                  node_color_axis_label = "Log2 ratio of median counts",
                  layout = "davidson-harel", 
                  initial_layout = "reingold-tilford") 

ggsave("../images/sub3pphyla.png",sub3pphyla)
ggsave("../images/sub3ptopc.png",sub3ptopc)
ggsave("../images/sub3ptopg.png",sub3ptopg)
ggsave("../images/sub3ptops.png",sub3ptops)
ggsave("../images/sub3pbeta.png",sub3pbeta)
ggsave("../images/sub3pht.png",sub3pht)


####################################################

###################
# COMPARISON 4 HET-pregnant VS HET-pregnant chiro
###################

# Subsetting
sub4=subset_samples(qiimed,treatment.group=="HET-pregnant" | treatment.group=="HET-pregnant chiro")
sub4rel=subset_samples(qiimerel,treatment.group=="HET-pregnant" | treatment.group=="HET-pregnant chiro")

# Comparing phyla abundance & statistical tests
sub4pphyla=ggplot(qiime_phyla_abundance[which(qiime_phyla_abundance$treatment.group==c("HET-pregnant", "HET-pregnant chiro")),],
                  aes(x=reorder(Phylum, Abundance, FUN = median,decreasing=T), y=Abundance, color=treatment.group)) +
  geom_boxplot()+ theme_bw()+theme(axis.title.x = element_blank(),
                                   axis.text.x = element_text(angle = 90,hjust=1,vjust=0.4))+ 
  labs(col='Group') +ylab("Relative abundance (%)")+
  stat_compare_means(aes(label = paste0(..p.format.., "\n", ..p.signif..)),
                     vjust=1,hide.ns = T,show.legend=FALSE)+ylim(0,1)
# Some phyla were undetected in some samples


# Beta diversity
sub4dist = phyloseq::distance(sub4, method="unifrac", weighted=T)

sub4perm=adonis2(sub4dist ~ sample_data(sub4)$treatment.group,data = data.frame(sample_data(sub4)))

sub4pbeta=plot_ordination(sub4, ordinate(sub4, method="PCoA", distance=sub4dist), 
                          color="treatment.group") +theme_bw()+ theme(aspect.ratio=1)+ 
  stat_ellipse(aes(group=treatment.group))+  stat_ellipse(aes(group=treatment.group))+
  annotate("text", x = -Inf, y = Inf, hjust = -3.2, vjust = 1.1, #Adjust to display in plot
           label= paste0("PERMANOVA", "\n", "p-value: ",sub4perm$`Pr(>F)`[1]))

# Calculating top 10 ASVs

sub4top=tax_glom(sub4rel, taxrank = "Genus")

sub4leasttaxa=intersect(names(sort(taxa_sums(subset_samples(sub4top,treatment.group=="HET-pregnant")), 
                                   decreasing = TRUE))[-c(1:10)],
                        names(sort(taxa_sums(subset_samples(sub4top,treatment.group=="HET-pregnant chiro")), 
                                   decreasing = TRUE))[-c(1:10)])
# Merging all other taxa

sub4top@tax_table[sub4leasttaxa,"Genus"]=c("0ther")

# COMPLETE
sub4ptopc=plot_bar(transform_sample_counts(sub4top, function(x) x/sum(x)), fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(sub4top@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(13))+
  scale_fill_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(13))+
  theme_bw()+theme(axis.title.x = element_blank(), panel.background = element_blank(),
                   axis.text.x = element_text(angle = 90), panel.grid.major = element_blank(),
                   strip.background = element_blank())+ylab("Relative abundance")

# GROUPED
sub4ptopg=plot_bar(transform_sample_counts(sub4top, function(x) x/(sum(x)*9)), x="treatment.group",fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(sub4top@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(13))+
  scale_fill_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(13))+
  theme(axis.title.x = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(),strip.background = element_blank(),
        axis.text.x = element_text(angle = 0))+ylab("Relative abundance")

# Comparing genera abundance & statistical tests
sub4ptops=ggplot(qiime_top_abundance[which(qiime_top_abundance$treatment.group==c("HET-pregnant", "HET-pregnant chiro")),],
                 aes(x=reorder(Genus, Abundance, FUN = median,decreasing=T), y=Abundance, color=treatment.group)) +
  geom_boxplot()+ theme_bw()+theme(axis.title.x = element_blank(),
                                   axis.text.x = element_text(angle = 90,hjust=1,vjust=0.4))+ 
  labs(col='Group') +ylab("Relative abundance")+
  stat_compare_means(aes(label = paste0(..p.format.., "\n", ..p.signif..)),
                     vjust=1,hide.ns = T,show.legend=FALSE)+ylim(0,1)





# Building a Heat tree
library(metacoder)
sub4metac=parse_phyloseq(sub4, class_regex = "(.*)", class_key = "taxon_name")

# Calculate relative abundance and per-taxon
sub4metac$data$rel_abd <- calc_obs_props(sub4metac, "otu_table", other_cols = T)
sub4metac$data$tax_rel_abd <- calc_taxon_abund(sub4metac, "rel_abd")

sub4metac$data$diff_table <- compare_groups(sub4metac,
                                            data = "tax_rel_abd",
                                            cols = sub4metac$data$sample_data$sample_id, # What columns of sample data to use
                                            groups = sub4metac$data$sample_data$treatment.group) # What category each sample is assigned to

#Correction for multiple tests
sub4metac$data$diff_table$wilcox_p_value <- p.adjust(sub4metac$data$diff_table$wilcox_p_value,
                                                     method = "fdr")

range(sub4metac$data$diff_table$log2_median_ratio, finite = TRUE) 

#Set non-significant values to zero
sub4metac$data$diff_table$log2_median_ratio[sub4metac$data$diff_table$wilcox_p_value > 0.05] <- 0

# plot heat tree matrix
set.seed(1)
sub4pht=heat_tree(sub4metac,
                  node_label = taxon_names,
                  node_size = n_obs,
                  node_color = log2_median_ratio, 
                  node_color_interval = c(-3, 3),
                  node_color_range = diverging_palette(),
                  node_size_axis_label = "ASV count",
                  node_color_axis_label = "Log2 ratio of median counts",
                  layout = "davidson-harel", 
                  initial_layout = "reingold-tilford") 

ggsave("../images/sub4pphyla.png",sub4pphyla)
ggsave("../images/sub4ptopc.png",sub4ptopc)
ggsave("../images/sub4ptopg.png",sub4ptopg)
ggsave("../images/sub4ptops.png",sub4ptops)
ggsave("../images/sub4pbeta.png",sub4pbeta)
ggsave("../images/sub4pht.png",sub4pht)

####################################################

###################
# COMPARISON 5 HET-pregnant VS HET-pregnant folic
###################

# Subsetting
sub5=subset_samples(qiimed,treatment.group=="HET-pregnant" | treatment.group=="HET-pregnant folic")
sub5rel=subset_samples(qiimerel,treatment.group=="HET-pregnant" | treatment.group=="HET-pregnant folic")

# Comparing phyla abundance & statistical tests
sub5pphyla=ggplot(qiime_phyla_abundance[which(qiime_phyla_abundance$treatment.group==c("HET-pregnant", "HET-pregnant folic")),],
                  aes(x=reorder(Phylum, Abundance, FUN = median,decreasing=T), y=Abundance, color=treatment.group)) +
  geom_boxplot()+ theme_bw()+theme(axis.title.x = element_blank(),
                                   axis.text.x = element_text(angle = 90,hjust=1,vjust=0.4))+ 
  labs(col='Group') +ylab("Relative abundance (%)")+
  stat_compare_means(aes(label = paste0(..p.format.., "\n", ..p.signif..)),
                     vjust=1,hide.ns = T,show.legend=FALSE)+ylim(0,1)
# Some phyla were undetected in some samples


# Beta diversity
sub5dist = phyloseq::distance(sub5, method="unifrac", weighted=T)

sub5perm=adonis2(sub5dist ~ sample_data(sub5)$treatment.group,data = data.frame(sample_data(sub5)))

sub5pbeta=plot_ordination(sub5, ordinate(sub5, method="PCoA", distance=sub5dist), 
                          color="treatment.group") +theme_bw()+ theme(aspect.ratio=1)+ 
  stat_ellipse(aes(group=treatment.group))+  stat_ellipse(aes(group=treatment.group))+
  annotate("text", x = -Inf, y = Inf, hjust = -3.2, vjust = 1.1, #Adjust to display in plot
           label= paste0("PERMANOVA", "\n", "p-value: ",sub5perm$`Pr(>F)`[1]))

# Calculating top 10 ASVs

sub5top=tax_glom(sub5rel, taxrank = "Genus")

sub5leasttaxa=intersect(names(sort(taxa_sums(subset_samples(sub5top,treatment.group=="HET-pregnant")), 
                                   decreasing = TRUE))[-c(1:10)],
                        names(sort(taxa_sums(subset_samples(sub5top,treatment.group=="HET-pregnant folic")), 
                                   decreasing = TRUE))[-c(1:10)])
# Merging all other taxa

sub5top@tax_table[sub5leasttaxa,"Genus"]=c("0ther")

# COMPLETE
sub5ptopc=plot_bar(transform_sample_counts(sub5top, function(x) x/sum(x)), fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(sub5top@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(13))+
  scale_fill_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(13))+
  theme_bw()+theme(axis.title.x = element_blank(), panel.background = element_blank(),
                   axis.text.x = element_text(angle = 90), panel.grid.major = element_blank(),
                   strip.background = element_blank())+ylab("Relative abundance")

# GROUPED
sub5ptopg=plot_bar(transform_sample_counts(sub5top, function(x) x/(sum(x)*9)), x="treatment.group",fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(sub5top@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(13))+
  scale_fill_manual(values = colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")))(13))+
  theme(axis.title.x = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(),strip.background = element_blank(),
        axis.text.x = element_text(angle = 0))+ylab("Relative abundance")

# Comparing genera abundance & statistical tests
sub5ptops=ggplot(qiime_top_abundance[which(qiime_top_abundance$treatment.group==c("HET-pregnant", "HET-pregnant folic")),],
                 aes(x=reorder(Genus, Abundance, FUN = median,decreasing=T), y=Abundance, color=treatment.group)) +
  geom_boxplot()+ theme_bw()+theme(axis.title.x = element_blank(),
                                   axis.text.x = element_text(angle = 90,hjust=1,vjust=0.4))+ 
  labs(col='Group') +ylab("Relative abundance")+
  stat_compare_means(aes(label = paste0(..p.format.., "\n", ..p.signif..)),
                     vjust=1,hide.ns = T,show.legend=FALSE)+ylim(0,1)





# Building a Heat tree
library(metacoder)
sub5metac=parse_phyloseq(sub5, class_regex = "(.*)", class_key = "taxon_name")

# Calculate relative abundance and per-taxon
sub5metac$data$rel_abd <- calc_obs_props(sub5metac, "otu_table", other_cols = T)
sub5metac$data$tax_rel_abd <- calc_taxon_abund(sub5metac, "rel_abd")

sub5metac$data$diff_table <- compare_groups(sub5metac,
                                            data = "tax_rel_abd",
                                            cols = sub5metac$data$sample_data$sample_id, # What columns of sample data to use
                                            groups = sub5metac$data$sample_data$treatment.group) # What category each sample is assigned to

#Correction for multiple tests
sub5metac$data$diff_table$wilcox_p_value <- p.adjust(sub5metac$data$diff_table$wilcox_p_value,
                                                     method = "fdr")

range(sub5metac$data$diff_table$log2_median_ratio, finite = TRUE) 

#Set non-significant values to zero
sub5metac$data$diff_table$log2_median_ratio[sub5metac$data$diff_table$wilcox_p_value > 0.05] <- 0

# plot heat tree matrix
set.seed(1)
sub5pht=heat_tree(sub5metac,
                  node_label = taxon_names,
                  node_size = n_obs,
                  node_color = log2_median_ratio, 
                  node_color_interval = c(-3, 3),
                  node_color_range = diverging_palette(),
                  node_size_axis_label = "ASV count",
                  node_color_axis_label = "Log2 ratio of median counts",
                  layout = "davidson-harel", 
                  initial_layout = "reingold-tilford") 

ggsave("../images/sub5pphyla.png",sub5pphyla)
ggsave("../images/sub5ptopc.png",sub5ptopc)
ggsave("../images/sub5ptopg.png",sub5ptopg)
ggsave("../images/sub5ptops.png",sub5ptops)
ggsave("../images/sub5pbeta.png",sub5pbeta)
ggsave("../images/sub5pht.png",sub5pht)