#!/bin/bash

conda activate mothur@1.44.11

# UNIÓN DE SECUENCIAS PAREADAS Y FILTRADO

# Creación de la tabla con los pares de secuencias
#mothur '#make.file(inputdir=../raw_data, outputdir=., type=gz, prefix=stability)'

# Se unen las secuencias pareadas (esto también pasa un filtro de calidad)
mothur '#make.contigs(file=stability.files, inputdir=../raw_data, processors=50)'


# Visualizar la tabla
mothur '#summary.seqs(fasta=stability.trim.contigs.fasta, groups=stability.contigs.groups)'

# Se le puede indicar filtros más restrictivos
#mothur '#make.contigs(file=stability.files, maxambig=0, maxlength=275, maxhomop=8)'

# Alternativa: usar directamente el comando screen para filtrar
mothur '#screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, maxambig=0, minlength=400, maxlength=480, maxhomop=8)'

mothur '#summary.seqs(fasta=stability.trim.contigs.good.fasta, groups=stability.contigs.good.groups)'

# PROCESAMIENTO DE SECUENCIAS

# Se genera stability.trim.contigs.good.names y stability.trim.contigs.good.unique.fasta

mothur '#unique.seqs(fasta=stability.trim.contigs.good.fasta)'

# Se genera stability.trim.contigs.good.count_table

mothur '#count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups)'

# Resumen
mothur '#summary.seqs(fasta=stability.trim.contigs.good.unique.fasta, count=stability.trim.contigs.good.count_table)'

# Reducción de la base de datos
mothur '#pcr.seqs(fasta=../db/silva.nr_v138_1.align, oligos=oligos.oligos, keepdots=F)'
# En este caso los oligos son CCTACGGGNGGCWGCAG y GACTACHVGGGTATCTAATCC

# Alineamineto de secuencias
mothur '#align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=../db/silva.nr_v138_1.pcr.align)'

mothur '#filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=.)'

mothur '#unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.count_table)'

# Precluster y eliminacion de quimeras (????????)
mothur '#pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table, diffs=5)'

mothur '#chimera.vsearch(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)'

mothur '#remove.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)'

mothur '#summary.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)'

#48224??? NO se si ha eliminado demasiadas en pre.cluster

mothur '#classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=../db/silva.nr_v138_1.align, taxonomy=../db/silva.nr_v138_1.tax, processors=100)'

mothur '#remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v138_1.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)'

mothur '#summary.tax(taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v138_1.wang.pick.taxonomy, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)'

################################HASTA AQUI

# Seleccionar un grupo (cambiar "Mock" por el nombre del grupo)
#mothur '#get.groups(count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, groups=Mock)'

# Errores en el grupo
#mothur '#seq.error(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.fasta, count=stability.trim.contigs.goodunique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, reference=HMP_MOCK.v35.fasta, aligned=F)'

#PASOS CLUSTERING
mothur '#dist.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.03)'
mothur '#cluster(column=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)'
mothur '#make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)'
mothur '#rarefaction.single(shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared)'

#Preparando analisis

# Quitar grupo
#mothur '#remove.groups(count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, taxonomy=stability.trim.contigs.unique.good.filter.unique.precluster.pds.wang.pick.taxonomy, groups=Mock)'

mothur '#rename.file(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v138_1.wang.taxonomy, prefix=final)'

#OTU CLUSTERING

#Matriz de distancias
mothur '#dist.seqs(fasta=final.fasta, cutoff=0.03)'
mothur '#cluster(column=final.dist, count=final.count_table)'

#Alternativa a las dos lineas anteriores
#mothur '#cluster.split(fasta=final.fasta, count=final.count_table, taxonomy=final.taxonomy, taxlevel=4, cutoff=0.03)'

#Obtener informacion de secuencias por OTU por grupo
mothur '#make.shared(list=final.opti_mcc.list, count=final.count_table, label=0.03)'

#CLasificación taxonómica
mothur '#classify.otu(list=final.opti_mcc.list, count=final.count_table, taxonomy=final.taxonomy, label=0.03)'

# ALTERNATIVA: ASV
#mothur '#make.shared(count=final.count_table)'

#mothur '#classify.otu(list=final.asv.list, count=final.count_table, taxonomy=final.taxonomy, label=ASV)'

# Phylotypes

mothur '#phylotype(taxonomy=final.taxonomy)'

mothur '#make.shared(list=final.tx.list, count=final.count_table, label=1)'

mothur '#classify.otu(list=final.tx.list, count=final.count_table, taxonomy=final.taxonomy, label=1)'

# Phylogenetic

mothur '#dist.seqs(fasta=final.fasta, output=lt)'
mothur '#clearcut(phylip=final.phylip.dist)'

# ANALISIS ESTADISTICO

mothur '#count.groups(shared=final.opti_mcc.shared)'

mothur '#sub.sample(shared=final.opti_mcc.shared, size=2403)'

# OTUS

# Diversidad alfa
mothur '#rarefaction.single(shared=final.opti_mcc.shared, calc=sobs, freq=100)'

# Diversidad beta
mothur '#dist.shared(shared=final.opti_mcc.shared, calc=thetayc-jclass, subsample=t)'

mothur '#pcoa(phylip=final.opti_mcc.thetayc.0.03.lt.ave.dist)'
mothur '#nmds(phylip=final.opti_mcc.thetayc.0.03.lt.ave.dist)'

mothur '#nmds(phylip=final.opti_mcc.thetayc.0.03.lt.ave.dist, mindim=3, maxdim=3)'

########################################

mothur '#amova(phylip=final.opti_mcc.thetayc.0.03.lt.ave.dist, design=mouse.time.design)'

mothur '#homova(phylip=stability.opti_mcc.thetayc.0.03.lt.ave.dist, design=mouse.time.design)'

mothur '#corr.axes(axes=stability.opti_mcc.thetayc.0.03.lt.ave.pcoa.axes, shared=stability.opti_mcc.0.03.subsample.shared, method=spearman, numaxes=3)'

mothur '#corr.axes(axes=stability.opti_mcc.thetayc.0.03.lt.ave.pcoa.axes, metadata=mouse.dpw.metadata, method=spearman, numaxes=3)'

mothur '#get.communitytype(shared=final.opti_mcc.0.03.subsample.shared)'

mothur '#metastats(shared=final.opti_mcc.0.03.subsample.shared, design=mouse.time.design)'

mothur '#lefse(shared=final.opti_mcc.0.03.subsample.shared, design=mouse.time.design)'

# Analisis con ASV
mothur '#phylo.diversity(tree=final.phylip.tre, count=final.count_table, rarefy=T)'

mothur '#unifrac.unweighted(tree=final.phylip.tre, count=final.count_table, distance=lt,random=F, subsample=t)'
mothur '#unifrac.weighted(tree=final.phylip.tre, count=final.count_table, distance=lt, random=F, subsample=t)'
