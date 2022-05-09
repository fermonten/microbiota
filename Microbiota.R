library(phyloseq)
library(ggplot2)

# MOTHUR
mothlist  = system.file("extdata", "esophagus.fn.list.gz", package="phyloseq")
mothgroup = system.file("extdata", "esophagus.good.groups.gz", package="phyloseq")
mothtree  = system.file("extdata", "esophagus.tree.gz", package="phyloseq")
show_mothur_cutoffs(mothlist)


# Assign variables for imported data
sharedfile = "final.opti_mcc.shared"
taxfile = "final.opti_mcc.0.03.cons.taxonomy"
treefile <- "final.opti_mcc.jclass.0.03.tre" 
mapfile = "metadata.tsv"

mothur_data=import_biom("final.opti_mcc.0.03.biom",treefilename="final.opti_mcc.jclass.0.03.tre" )

# Import mothur data
mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile,
                             mothur_tree_file = treefile)

#QIIME
library(qiime2R)
qiime=qza_to_phyloseq("table.qza", tree="rep_seqs_aligned_masked_tree_rooted.qza", metadata="metadata0.tsv", taxonomy="classification.qza")

qiimed=prune_samples(sample_names(qiime)!="CN-16S",qiime)
ciego=subset_samples(qiime, sample.type=="ciego")
heces=subset_samples(qiime, sample.type=="heces")

heces0=merge_samples(heces, "treatment.group")

#tutorial
plot_tree(qiimed, color="treatment.group")

plot_richness(qiimed, x="treatment.group",color="samples")

#Observando todas las muestras
plot_bar(qiimed, fill = "Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")+facet_wrap(~treatment.group,scales = "free")

plot_bar(qiimed, fill = "Family")+
  facet_wrap(~treatment.group,scales = "free",nrow=1)


#Diversidad alfa
plot_richness(qiimed, x="treatment.group", color="samples")+ theme(legend.position="none")
plot_richness(qiimed, x="treatment.group", color="samples")+ theme(legend.position="none")+geom_boxplot()

# Calcular la diversidad alfa y hacer test de wilcoxon para ver si existen diferencias
rich = estimate_richness(qiimed)
pairwise.wilcox.test(rich$Observed, sample_data(qiimed)$treatment.group)

# Diversidad beta
wunifrac_dist = phyloseq::distance(qiimed, method="unifrac", weighted=T)
uordination = ordinate(qiimed, method="PCoA", distance=wunifrac_dist)
plot_ordination(qiimed, uordination, color="treatment.group") + theme(aspect.ratio=1)+ 
  stat_ellipse(aes(group=treatment.group))

library(vegan)
adonis(wunifrac_dist ~ sample_data(qiimed)$treatment.group)


bc_dist = phyloseq::distance(qiimed, method="bray")
bcordination = ordinate(qiimed, method="PCoA", distance=bc_dist)
plot_ordination(qiimed, bcordination, color="treatment.group") + theme(aspect.ratio=1)

c_dist = phyloseq::distance(qiimed, method="canberra")
cordination = ordinate(qiimed, method="PCoA", distance=c_dist)
plot_ordination(qiimed, cordination, color="treatment.group") + theme(aspect.ratio=1)+ 
  stat_ellipse(aes(group=treatment.group))

library(vegan)
adonis2(c_dist ~ sample_data(qiimed)$treatment.group,data = data.frame(sample_data(qiimed)))

qiimedm=qiimed$treatment.group[-which(sample_data(qiimed)$treatment.group=="Feces WT")]
a=a[-which(a=="Feces HET")]
adonis2(c_dist ~ a,data = data.frame(a))
# Comparaciones

library(DESeq2)
sample_data(qiimed)$treatment.group <- as.factor(sample_data(qiimed)$treatment.group)
ds = phyloseq_to_deseq2(qiimed,~treatment.group)
ds = DESeq(ds)

# Heces NO embarazadas: HET vs WT
res1= results(ds, contrast=c("treatment.group", "Feces HET", "Feces WT"), alpha=0.01)
res1 = res1[order(res1$padj, na.last=NA), ]
res1_sig = res1[(res1$padj < 0.01),]

res1_sig = cbind(as(res1_sig, "data.frame"), as(tax_table(qiimed)[rownames(res1_sig), ], "matrix"))
ggplot(res1_sig, aes(x=Phylum, y=log2FoldChange, color=Genus)) +  
  geom_jitter(size=3, width = 0.2) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  ggtitle("Feces HET vs Feces WT")

# Ciego vs heces de NO embarazadas HET 
res2= results(ds, contrast=c("treatment.group", "Feces HET", "HET-non pregnant"), alpha=0.01)
res2 = res2[order(res2$padj, na.last=NA), ]
res2_sig = res2[(res2$padj < 0.01),]

res2_sig = cbind(as(res2_sig, "data.frame"), as(tax_table(qiimed)[rownames(res2_sig), ], "matrix"))
ggplot(res2_sig, aes(x=Phylum, y=log2FoldChange, color=Genus)) +
  geom_jitter(size=3, width = 0.2) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  ggtitle("Feces HET vs HET-non pregnant")

# Ciego NO embarazadas HET vs embarazadas HET
res3= results(ds, contrast=c("treatment.group", "HET-non pregnant", "HET-pregnant"), alpha=0.01)
res3 = res3[order(res3$padj, na.last=NA), ]
res3_sig = res3[(res3$padj < 0.01),]

res3_sig = cbind(as(res3_sig, "data.frame"), as(tax_table(qiimed)[rownames(res3_sig), ], "matrix"))
ggplot(res3_sig, aes(x=Phylum, y=log2FoldChange, color=Genus)) +
  geom_jitter(size=3, width = 0.2) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  ggtitle("HET-non pregnant vs HET-pregnant")

# Ciego embarazadas HET dieta normal vs D-chiro-inositol
res4= results(ds, contrast=c("treatment.group", "HET-pregnant chiro", "HET-pregnant"), alpha=0.01)
res4 = res4[order(res4$padj, na.last=NA), ]
res4_sig = res4[(res4$padj < 0.01),]

res4_sig = cbind(as(res4_sig, "data.frame"), as(tax_table(qiimed)[rownames(res4_sig), ], "matrix"))
ggplot(res4_sig, aes(x=Phylum, y=log2FoldChange, color=Genus)) +
  geom_jitter(size=3, width = 0.2) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  ggtitle("HET-pregnant chiro vs HET-pregnant")

# Ciego embarazadas HET dieta normal vs Ácido fólico
res5= results(ds, contrast=c("treatment.group", "HET-pregnant folic", "HET-pregnant"), alpha=0.01)
res5 = res5[order(res5$padj, na.last=NA), ]
res5_sig = res5[(res5$padj < 0.01),]

res5_sig = cbind(as(res5_sig, "data.frame"), as(tax_table(qiimed)[rownames(res5_sig), ], "matrix"))
ggplot(res5_sig, aes(x=Phylum, y=log2FoldChange, color=Genus)) +
  geom_jitter(size=3, width = 0.2) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  ggtitle("HET-pregnant folic vs HET-pregnant")

# Agrupando en filo
phylumGlommed = tax_glom(heces0, "Phylum")
plot_bar(phylumGlommed, x="treatment.group", fill ="Phylum")

familyGlommed = tax_glom(heces0, "Family")
plot_bar(familyGlommed, x="treatment.group", fill ="Family")

plot_bar(heces0,fill="Family",x="treatment.group") + facet_wrap(~Family) + theme(legend.position="none")

plot_heatmap(qiimed,taxa.label = "phylum")+ facet_wrap(~Family) + theme(legend.position="none")
