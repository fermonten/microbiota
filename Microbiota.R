library(phyloseq)
library(ggplot2)
library(qiime2R)
library(vegan)
library(RColorBrewer)
library(tidyverse)
library(MicrobiotaProcess)
library(VennDiagram)
library(metacoder)
library(phyloseqCompanion)
library(ggpubr)

# MOTHUR
# Import mothur data
#mothur_data <- import_mothur(mothur_shared_file = sharedfile,
#mothur_constaxonomy_file = taxfile,
#mothur_tree_file = treefile)

#QIIME
qiime=qza_to_phyloseq("table.qza", tree="rep_seqs_aligned_masked_tree_rooted.qza", 
                      metadata="metadata0.tsv", taxonomy="classification.qza")

# Removing the negative control and the unassigned ASVs at phylum level
qiimed=prune_samples(sample_names(qiime)!="CN-16S",qiime)
qiimed=subset_taxa(qiimed,Phylum!="NA")

ciego=subset_samples(qiimed, sample.type=="ciego")
heces=subset_samples(qiimed, sample.type=="heces")


# Common ASVs Venn Diagram
venn.diagram(get_vennlist(ciego,factorNames="treatment.group"),filename="vennmicro.png",
             fill = c("yellow", "green", "blue", "red"),width=4000,height=4000)
venn.diagram(get_vennlist(heces,factorNames="treatment.group"),filename="vennmicroh.png",
             fill = c("blue", "red"),width=4000,height=4000,ext.pos=c(340))
venn.diagram(get_vennlist(comp,factorNames="treatment.group"),filename="vennmicroc.png",
             fill = c("blue", "yellow"),width=4000,height=4000,ext.pos=c(340))

# Extra plots not used

#plot_heatmap(qiimed)
#plot_tree(qiimed)

# Rarefaction curve
rarecurve=rarecurve(t(otu_table(qiimed)), step=50, cex=0.5)
rare=map_dfr(rarecurve,bind_rows) %>% 
  bind_cols(sample=rownames(qiimed@sam_data),.) %>%
  bind_cols(Group=qiimed@sam_data$treatment.group,.) %>%
  pivot_longer(-c(sample,Group))%>%
  drop_na() %>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name) 


ggplot(data=rare,aes(x=n_seqs,y=value,group=sample,colour=Group))+
  geom_line(size=1.2)+facet_wrap(~Group,nrow=2)+
  labs(x= "Sequencing depth",y="Number of ASVs")+theme_bw()+
  theme(strip.background = element_blank()) 

# Gráficos de abundancia de secuencias

# Calculando la abundancia relativa
qiimerel=transform_sample_counts(qiimed, function(x) x/sum(x))
ciegorel=transform_sample_counts(ciego, function(x) x/sum(x))
hecesrel=transform_sample_counts(heces, function(x) x/sum(x))

# Generating relative abundance table
qiimegroups=merge_samples(qiimerel,group = "treatment.group")
abundancetable=psmelt(tax_glom(qiimegroups,taxrank="Phylum"))
abutable=abundancetable[,c(2,3,14)]
abutable=cbind(abutable[which(abutable$Sample=="HET-pregnant chiro"),][,-c(1,3)],
               abutable[which(abutable$Sample=="HET-pregnant folic"),][,-c(1,3)],
               abutable[which(abutable$Sample=="HET-pregnant"),][,-c(1,3)],
               abutable[which(abutable$Sample=="HET-non pregnant"),][,-c(1,3)],
               abutable[which(abutable$Sample=="Feces HET"),][,-c(1,3)],
               abutable[which(abutable$Sample=="Feces WT"),][,-c(1,3)]) 
rownames(abutable)=unique(abundancetable[,14])
colnames(abutable)=unique(abundancetable[,2])      
abutable=abutable/9

write.table(abutable,file="relabundance.csv")

# Representando en gráficos de barras, por filo

# COMPLETE-ALL SAMPLES
#A=
plot_bar(qiimerel, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(qiimerel@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_brewer(palette="Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme(legend.position = c(0.7,0.1), legend.justification = c(1, 0),
        axis.title.x = element_blank(),axis.title.y = element_blank(),
        panel.background = element_blank(),panel.grid.major = element_blank(),
        strip.background = element_blank()) 

# GROUPED
plot_bar(transform_sample_counts(qiimerel, function(x) x/9), x="treatment.group",fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")+
  guides(color=guide_legend(ncol=1))+
  scale_color_brewer(palette="Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank(),panel.grid.major = element_blank(),
        strip.background = element_blank())+ylab("Relative abundance") 

#Plotting phyla side by side
qiime_phyla_abundance=psmelt(tax_glom(qiimerel,taxrank="Phylum"))[,c(1,2,3,5,14)]

ggplot(qiime_phyla_abundance,aes(x=reorder(Phylum, Abundance, FUN = median,decreasing=T), 
                                 y=Abundance,fill=treatment.group)) +
  geom_boxplot()+ theme_bw()+theme(axis.title.x = element_blank())+ 
  labs(col='Group') +ylab("Relative abundance")+
  geom_point(position = position_dodge(width=0.75),color="black")


######################
# COMPARISON 1 FECES HET VS FECES WT
###################
phyla1=ggplot(qiime_phyla_abundance[which(qiime_phyla_abundance$treatment.group==c("Feces WT", "Feces HET")),],aes(x=reorder(Phylum, Abundance, FUN = median,decreasing=T), y=Abundance, color=treatment.group)) +
  geom_boxplot()+ theme_bw()+theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 90,hjust=1))+ 
  labs(col='Group') +ylab("Relative abundance (%)")+stat_compare_means(aes(label = paste0(..p.format.., "\n", ..p.signif..)),vjust=1,hide.ns = T,show.legend=FALSE)+ylim(0,1)
# Some phyla were undetected in some samples

subset1=subset_samples(qiimed,treatment.group=="Feces WT" | treatment.group=="Feces HET")
subset1rel=transform_sample_counts(subset1, function(x) x/sum(x))

# Beta diversity
dist1 = phyloseq::distance(subset1, method="unifrac", weighted=T)

perm1=adonis2(dist1 ~ sample_data(subset1)$treatment.group,data = data.frame(sample_data(subset1)))

PCOA1=plot_ordination(subset1, ordinate(subset1, method="PCoA", distance=dist1), 
                      color="treatment.group") +theme_bw()+ theme(aspect.ratio=1)+ 
  stat_ellipse(aes(group=treatment.group))+  stat_ellipse(aes(group=treatment.group))+
  annotate("text", x = -Inf, y = Inf, hjust = -3.2, vjust = 1.1, 
           label= paste0("PERMANOVA", "\n", "p-value: ",perm1$`Pr(>F)`[1]))


###################
# COMPARISON 2 FECES HET VS CECUM HET
###################
phyla2=ggplot(qiime_phyla_abundance[which(qiime_phyla_abundance$treatment.group==c("Feces HET", "HET-non pregnant")),],aes(x=reorder(Phylum, Abundance, FUN = median,decreasing=T), y=Abundance, color=treatment.group)) +
  geom_boxplot()+ theme_bw()+theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 90,hjust=1))+ 
  labs(col='Group') +ylab("Relative abundance (%)")+stat_compare_means(aes(label = paste0(..p.format.., "\n", ..p.signif..)),vjust=1,hide.ns = T,show.legend=FALSE)
# Some phyla were undetected in some samples

subset2=subset_samples(qiimed,treatment.group=="Feces HET" | treatment.group=="HET-non pregnant")
subset2rel=transform_sample_counts(subset2, function(x) x/sum(x))

# Beta diversity
dist2 = phyloseq::distance(subset2, method="unifrac", weighted=T)

perm2=adonis2(dist2 ~ sample_data(subset2)$treatment.group,data = data.frame(sample_data(subset2)))
#SIGNIFICANT

PCOA2=plot_ordination(subset2, ordinate(subset2, method="PCoA", distance=dist2), 
                      color="treatment.group") +theme_bw()+ theme(aspect.ratio=1)+ 
  stat_ellipse(aes(group=treatment.group))+  stat_ellipse(aes(group=treatment.group))+
  annotate("text", x = -Inf, y = Inf, hjust = -3.2, vjust = 1.1, 
           label= paste0("PERMANOVA", "\n", "p-value: ",perm2$`Pr(>F)`[1]))



###################
# COMPARISON 3 CECUM HET (NON PREGNANT) vs CECUM HET (PREGNANT)
###################

phyla3=ggplot(qiime_phyla_abundance[which(qiime_phyla_abundance$treatment.group==c("HET-non pregnant","HET-pregnant")),],aes(x=reorder(Phylum, Abundance, FUN = median,decreasing=T), y=Abundance, color=treatment.group)) +
  geom_boxplot()+ theme_bw()+theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 90,hjust=1))+ 
  labs(col='Group') +ylab("Relative abundance (%)")+stat_compare_means(aes(label = paste0(..p.format.., "\n", ..p.signif..)),vjust=1,hide.ns = T,show.legend=FALSE)+ylim(0,1)
# Some phyla were undetected in some samples

subset3=subset_samples(qiimed,treatment.group=="HET-non pregnant" | treatment.group=="HET-pregnant")
subset3rel=transform_sample_counts(subset3, function(x) x/sum(x))

# Beta diversity
dist3 = phyloseq::distance(subset3, method="unifrac", weighted=T)

perm3=adonis2(dist3 ~ sample_data(subset3)$treatment.group,data = data.frame(sample_data(subset3)))

PCOA3=plot_ordination(subset3, ordinate(subset3, method="PCoA", distance=dist3), 
                      color="treatment.group") +theme_bw()+ theme(aspect.ratio=1)+ 
  stat_ellipse(aes(group=treatment.group))+  stat_ellipse(aes(group=treatment.group))+
  annotate("text", x = -Inf, y = Inf, hjust = -3.2, vjust = 1.1, 
           label= paste0("PERMANOVA", "\n", "p-value: ",perm3$`Pr(>F)`[1]))

###################
# COMPARISON 4 CECUM HET (PREGNANT) vs CECUM HET (PREGNANT+D-CHIRO INOSITOL)
###################

phyla4=ggplot(qiime_phyla_abundance[which(qiime_phyla_abundance$treatment.group==c("HET-pregnant","HET-pregnant chiro")),],aes(x=reorder(Phylum, Abundance, FUN = median,decreasing=T), y=Abundance, color=treatment.group)) +
  geom_boxplot()+ theme_bw()+theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 90,hjust=1))+ 
  labs(col='Group') +ylab("Relative abundance (%)")+stat_compare_means(aes(label = paste0(..p.format.., "\n", ..p.signif..)),vjust=1,hide.ns = T,show.legend=FALSE)+ylim(0,1)
# Some phyla were undetected in some samples

subset4=subset_samples(qiimed,treatment.group=="HET-pregnant" | treatment.group=="HET-pregnant chiro")
subset4rel=transform_sample_counts(subset4, function(x) x/sum(x))

# Beta diversity
dist4 = phyloseq::distance(subset4, method="unifrac", weighted=T)

perm4=adonis2(dist4 ~ sample_data(subset4)$treatment.group,data = data.frame(sample_data(subset4)))
#SIGNIFICANT

PCOA4=plot_ordination(subset4, ordinate(subset4, method="PCoA", distance=dist4), 
                      color="treatment.group") +theme_bw()+ theme(aspect.ratio=1)+ 
  stat_ellipse(aes(group=treatment.group))+  stat_ellipse(aes(group=treatment.group))+
  annotate("text", x = -Inf, y = Inf, hjust = -3.2, vjust = 1.1, 
           label= paste0("PERMANOVA", "\n", "p-value: ",perm4$`Pr(>F)`[1]))


###################
# COMPARISON 5 CECUM HET (PREGNANT) vs CECUM HET (PREGNANT+FOLIC)
###################

phyla5=ggplot(qiime_phyla_abundance[which(qiime_phyla_abundance$treatment.group==c("HET-pregnant","HET-pregnant folic")),],aes(x=reorder(Phylum, Abundance, FUN = median,decreasing=T), y=Abundance, color=treatment.group)) +
  geom_boxplot()+ theme_bw()+theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 90,hjust=1))+ 
  labs(col='Group') +ylab("Relative abundance (%)")+stat_compare_means(aes(label = paste0(..p.format.., "\n", ..p.signif..)),vjust=1,hide.ns = T,show.legend=FALSE)+ylim(0,1)
# Some phyla were undetected in some samples

subset5=subset_samples(qiimed,treatment.group=="HET-pregnant" | treatment.group=="HET-pregnant chiro")
subset5rel=transform_sample_counts(subset5, function(x) x/sum(x))

# Beta diversity
dist5 = phyloseq::distance(subset5, method="unifrac", weighted=T)

perm5=adonis2(dist5 ~ sample_data(subset5)$treatment.group,data = data.frame(sample_data(subset5)))
#SIGNIFICANT

PCOA5=plot_ordination(subset5, ordinate(subset5, method="PCoA", distance=dist5), 
                      color="treatment.group") +theme_bw()+ theme(aspect.ratio=1)+ 
  stat_ellipse(aes(group=treatment.group))+ 
  stat_ellipse(aes(group=treatment.group))+
  annotate("text", x = -Inf, y = Inf, hjust = -3.2, vjust = 1.1, 
           label= paste0("PERMANOVA", "\n", "p-value: ",perm5$`Pr(>F)`[1]))




##########################################################

names(sort(taxa_sums(enterotype), TRUE)[1:10]) 

# Calculating top 10 ASVs

#1

qiime10_1 = prune_taxa(names(sort(taxa_sums(tax_glom(subset1rel, taxrank = "Genus")),
                                  decreasing = TRUE))[1:10], tax_glom(subset1rel, taxrank = "Genus"))

plot_bar(qiime10_1, fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(qiime10_1@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(brewer.pal(9,"Set1"))(10))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(9,"Set1"))(10))+
  theme_bw()+ theme(panel.background = element_blank(),panel.grid.major = element_blank(),
                    strip.background = element_blank(),axis.text.x = element_text(angle = 90)) 

#2
qiime10_2 = prune_taxa(names(sort(taxa_sums(tax_glom(subset2rel, taxrank = "Genus")),
                                  decreasing = TRUE))[1:10], tax_glom(subset2rel, taxrank = "Genus"))
# COMPLETE
plot_bar(qiime10_2, fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(qiime10_2@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(brewer.pal(9,"Set1"))(10))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(9,"Set1"))(10))+
  theme_bw()+ theme(panel.background = element_blank(),panel.grid.major = element_blank(),
                    strip.background = element_blank(),axis.text.x = element_text(angle = 90)) 

# GROUPED
plot_bar(qiime10_2, x="treatment.group", fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(qiime10_2@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(brewer.pal(9,"Set1"))(10))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(9,"Set1"))(10))+
  theme_bw()+ theme(panel.background = element_blank(),panel.grid.major = element_blank(),
                    strip.background = element_blank(),axis.text.x = element_text(angle = 90)) 

#3

qiime10_3 = prune_taxa(names(sort(taxa_sums(tax_glom(subset3rel, taxrank = "Genus")),
                                  decreasing = TRUE))[1:10], tax_glom(subset3rel, taxrank = "Genus"))
# COMPLETE
plot_bar(qiime10_3, fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(qiime10_3@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(brewer.pal(9,"Set1"))(10))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(9,"Set1"))(10))+
  theme_bw()+ theme(panel.background = element_blank(),panel.grid.major = element_blank(),
                    strip.background = element_blank(),axis.text.x = element_text(angle = 90)) 

# GROUPED
plot_bar(qiime10_3, x="treatment.group",fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(qiime10_3@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(brewer.pal(9,"Set1"))(10))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(9,"Set1"))(10))+
  theme_bw()+ theme(panel.background = element_blank(),panel.grid.major = element_blank(),
                    strip.background = element_blank(),axis.text.x = element_text(angle = 90)) 


#4
qiime10_4 = prune_taxa(names(sort(taxa_sums(tax_glom(subset4rel, taxrank = "Genus")),
                                  decreasing = TRUE))[1:10], tax_glom(subset4rel, taxrank = "Genus"))
# COMPLETE
plot_bar(qiime10_4, fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(qiime10_4@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(brewer.pal(9,"Set1"))(10))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(9,"Set1"))(10))+
  theme_bw()+ theme(panel.background = element_blank(),panel.grid.major = element_blank(),
                    strip.background = element_blank(),axis.text.x = element_text(angle = 90)) 
# GROUPED
plot_bar(qiime10_4, x="treatment.group",fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(qiime10_4@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(brewer.pal(9,"Set1"))(10))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(9,"Set1"))(10))+
  theme_bw()+ theme(panel.background = element_blank(),panel.grid.major = element_blank(),
                    strip.background = element_blank(),axis.text.x = element_text(angle = 90)) 

#5
qiime10_5 = prune_taxa(names(sort(taxa_sums(tax_glom(subset5rel, taxrank = "Genus")),
                                  decreasing = TRUE))[1:10], tax_glom(subset5rel, taxrank = "Genus"))
# COMPLETE
plot_bar(qiime10_5, fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(qiime10_5@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(brewer.pal(9,"Set1"))(10))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(9,"Set1"))(10))+
  theme_bw()+ theme(panel.background = element_blank(),panel.grid.major = element_blank(),
                    strip.background = element_blank(),axis.text.x = element_text(angle = 90)) 
# GROUPED
plot_bar(qiime10_5, x="treatment.group",fill = "Genus")+ 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_wrap(~factor(treatment.group,levels=c(sort(unique(qiime10_5@sam_data$treatment.group),
                                                   decreasing = T))),scales = "free",ncol=4) + 
  guides(color=guide_legend(ncol=1))+
  scale_color_manual(values = colorRampPalette(brewer.pal(9,"Set1"))(10))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(9,"Set1"))(10))+
  theme_bw()+ theme(panel.background = element_blank(),panel.grid.major = element_blank(),
                    strip.background = element_blank(),axis.text.x = element_text(angle = 90)) 

# Heat tree
metac=parse_phyloseq(qiimed, class_regex = "(.*)", class_key = "taxon_name")

# Calculate relative abundance and per-taxon
metac$data$rel_abd <- calc_obs_props(metac, "otu_table", other_cols = T)
metac$data$tax_rel_abd <- calc_taxon_abund(metac, "rel_abd")

metac$data$diff_table <- compare_groups(metac,
                                        data = "tax_rel_abd",
                                        cols = metac$data$sample_data$sample_id, # What columns of sample data to use
                                        groups = metac$data$sample_data$treatment.group) # What category each sample is assigned to

#Correction for multiple tests
metac$data$diff_table$wilcox_p_value <- p.adjust(metac$data$diff_table$wilcox_p_value,
                                                 method = "fdr")

range(metac$data$diff_table$log2_median_ratio, finite = TRUE) 

#Set non-significant values to zero
metac$data$diff_table$log2_median_ratio[metac$data$diff_table$wilcox_p_value > 0.05] <- 0

# plot heat tree matrix
heat_tree_matrix(metac,
                 data="diff_table",
                 seed=1,
                 node_label = taxon_names,
                 node_size = n_obs,
                 node_color = log2_median_ratio, 
                 node_color_interval = c(-6, 6),
                 node_color_range = diverging_palette(),
                 node_size_axis_label = "ASV count",
                 node_color_axis_label = "Log2 ratio of median counts",
                 layout = "davidson-harel", 
                 initial_layout = "reingold-tilford") 



#Diversidad alfa
plot_richness(qiimed, x="treatment.group", measures=c("Observed","Chao1","Shannon"), 
              color="treatment.group")+ theme(legend.position="none")+
  geom_boxplot()+ geom_point(color="black") + theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
        strip.background = element_blank(),axis.text.x = element_text(angle = 90,hjust=1))

plot_richness(heces, x="treatment.group", color="samples")+ theme(legend.position="none")+geom_boxplot()

# Calcular la diversidad alfa y hacer test de wilcoxon para ver si existen diferencias
rich = estimate_richness(qiimed)
write.table(rich, "alphadiversities.csv")

alphatest=pairwise.wilcox.test(rich$Observed, sample_data(qiimed)$treatment.group)
write.table(alphatest$p.value,file="alphatest.csv")

alphatest=pairwise.wilcox.test(rich$Shannon, sample_data(qiimed)$treatment.group)
write.table(alphatest$p.value,file="alphatest2.csv")

# Diversidad beta
wunifrac_dist = phyloseq::distance(qiimed, method="unifrac", weighted=T)
uordination = ordinate(qiimed, method="PCoA", distance=wunifrac_dist)
plot_ordination(qiimed, uordination, color="treatment.group") +theme_bw()+ theme(aspect.ratio=1)+ 
  stat_ellipse(aes(group=treatment.group))+  stat_ellipse(aes(group=treatment.group))

#plot_ordination(qiimed, uordination, color="sample.type") + theme(aspect.ratio=1)+ 
#  stat_ellipse(aes(group=sample.type))

cwunifrac_dist = phyloseq::distance(ciego, method="unifrac", weighted=T)
cuordination = ordinate(ciego, method="PCoA", distance=cwunifrac_dist)
plot_ordination(ciego, cuordination, color="treatment.group") + theme(aspect.ratio=1)+ 
  
  
  bc_dist = phyloseq::distance(ciego, method="bray")
bcordination = ordinate(ciego, method="PCoA", distance=bc_dist)
plot_ordination(ciego, bcordination, color="treatment.group") + theme(aspect.ratio=1)+  
  stat_ellipse(aes(group=treatment.group))

c_dist = phyloseq::distance(ciego, method="canberra")
cordination = ordinate(ciego, method="PCoA", distance=c_dist)
plot_ordination(ciego, cordination, color="treatment.group") + theme(aspect.ratio=1)+ 
  stat_ellipse(aes(group=treatment.group))

library(vegan)
adonis=adonis2(wunifrac_dist ~ sample_data(qiimed)$treatment.group,data = data.frame(sample_data(qiimed)))
write.table(adonis,"adonis.csv")

# Comparaciones

library(DESeq2)
sample_data(qiimed)$treatment.group <- as.factor(sample_data(qiimed)$treatment.group)
ds = phyloseq_to_deseq2(qiimed,~treatment.group)
ds = DESeq(ds)

#FRACASO
picrust=read.csv2(file="path_abun_unstrat_descrip.tsv",header=T,sep="\t")
rownames(picrust)=unlist(picrust[2])
picrust=picrust[-c(1,2,3)]
pic=as.matrix(as.data.frame(lapply(picrust,as.numeric)))
rownames(pic)=rownames(picrust)

picc=as.data.frame(cbind(as.character(qiimed@sam_data$treatment.group),t(picrust)))

picc2=as.data.frame(lapply(picc[-1],as.numeric))
rownames(picc2)=rownames(picc)
picc3=t(cbind(treatment.group=as.character(qiimed@sam_data$treatment.group),picc2))

picc4=cbind(picc3,picc3$treatment.group)

# Pathway abundance statistical tests with ALDEx2

# Making a model matrix

A=as.character(qiimed@sam_data$treatment.group)
B=rownames(qiimed@sam_data)
M=matrix(0,ncol=length(unique(A)),nrow=length(B))
colnames(M)=unique(A)
rownames(M)=B

for (i in 1:length(B)) {
  M[i,A[i]]=1
}

MM=model.matrix(~., as.data.frame(M))
colnames(MM)[2:7]=colnames(M)

# Tests estadísticos
library(ALDEx2)
aldexm = round(pic)
aldex_test = aldex.clr(aldexm, MM, verbose = TRUE)

aldex_glm=aldex.glm(aldex_test, data=as.data.frame(MM),verbose = T)

aldex_glm_effect=aldex.glm.effect(aldex_test)

aldex_test2 = aldex.clr(aldexm, as.character(qiimed@sam_data$treatment.group), verbose = TRUE)
aldex_kw=aldex.kw(aldex_test2)

effect=aldex.effect(aldex_test,glm.conds = as.character(qiimed@sam_data$treatment.group))

aldex_test2

# Plotting

phyloseq::ordinate(aldex_test, "RDA")

aldex_all=data.frame(aldex_kw,aldex_glm_effect)

aldex.plot(aldex_glm, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")

par(mfrow=c(1,1))
aldex.plot(aldex_glm_effect$`Feces WT`, test="effect", cutoff=2,main="Feces WT")
aldex.plot(aldex_glm_effect$`HET-non pregnant`, test="effect", cutoff=2)
aldex.plot(aldex_glm_effect$`Feces HET`, test="effect", cutoff=2)
aldex.plot(aldex_glm_effect$`HET-pregnant chiro`, test="effect", cutoff=2)
aldex.plot(aldex_glm_effect$`HET-pregnant`, test="effect", cutoff=2)
aldex.plot(aldex_glm_effect$`HET-pregnant folic`, test="effect", cutoff=2)

points(aldex_glm_effect$`Feces WT`$diff.win[sig],
       aldex_glm_effect$`Feces WT`$diff.btw[sig], col="blue")
sig <- aldex_glm[,20]<0.2
points(aldex_glm_effect$`Feces WT`$diff.win[sig],
       aldex_glm_effect$`Feces WT`$diff.btw[sig], col="blue")

aldex.plot(aldex_glm_effect$`Feces WT`, type="MA", test="effect")
aldex.plot(x.all, type="MW", test="glm")

aldex.plot(aldex_glm_effect$`Feces WT`, type="MA", test="glm", cutoff=2)
points(aldex_glm_effect$`Feces WT`$diff.win[sig],
       aldex_glm_effect$`Feces WT`$diff.btw[sig], col="blue")
sig <- aldex_glm[,20]<0.2
points(aldex_glm_effect$`Feces WT`$diff.win[sig],
       aldex_glm_effect$`Feces WT`$diff.btw[sig], col="blue")

plot(aldex_kw$kw.eBH, log="y", cex=0.7, col=rgb(0,0,1,0.2),
     pch=19, xlab="Difference", ylab="P value", main="Volcano plot")
points(x.all$diff.btw, x.all$we.eBH, cex=0.7, col=rgb(1,0,0,0.2),
       pch=19)
abline(h=0.05, lty=2, col="grey")



BHpVals_df = as.data.frame(aldex_kw$kw.eBH)
rownames(BHpVals_df)=rownames(aldex_kw)



ggplot(BHpVals_df) 

#importing to lefse
phyloseq2lefse(ps=qiimed,covars="treatment.group")

#???? riffomonas project


