#!/bin/bash

conda activate qiime2-2022.2

#Importando los fastq
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path raw_data \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path qiime/demux-paired-end.qza
  
#Generando gráficas para el control de calidad
qiime demux summarize \
  --i-data qiime/demux-paired-end.qza \
  --o-visualization qiime/demux-paired-end.qzv
  
#denoising con dada2  
qiime dada2 denoise-paired --verbose\
  --p-n-threads 0 \
  --i-demultiplexed-seqs qiime/demux-paired-end.qza \
  --o-table qiime/table \
  --o-representative-sequences qiime/rep-seqs \
  --o-denoising-stats qiime/denoising \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 300 \
  --p-trunc-len-r 250

#generando tablas  
qiime feature-table summarize \
  --i-table qiime/table.qza \
  --o-visualization qiime/table.qzv \
  --m-sample-metadata-file metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data qiime/rep-seqs.qza \
  --o-visualization qiime/rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file qiime/denoising.qza \
  --o-visualization qiime/denoising.qzv
  
#Construcción de un árbol filogenético
qiime alignment mafft \
--i-sequences qiime/rep-seqs.qza \
--o-alignment qiime/rep_seqs_aligned.qza

qiime alignment mask \
--i-alignment qiime/rep_seqs_aligned.qza \
--o-masked-alignment qiime/rep_seqs_aligned_masked.qza

qiime phylogeny fasttree \
--i-alignment qiime/rep_seqs_aligned_masked.qza \
--o-tree qiime/rep_seqs_aligned_masked_tree

qiime phylogeny midpoint-root \
  --i-tree qiime/rep_seqs_aligned_masked_tree.qza \
  --o-rooted-tree qiime/rep_seqs_aligned_masked_tree_rooted.qza
  
#Clasificación taxonómica
qiime feature-classifier classify-sklearn \
  --p-n-jobs -1 --verbose\
  --i-reads qiime/rep-seqs.qza \
  --i-classifier db/silva-138-99-nb-classifier.qza \
  --o-classification qiime/classification.qza
  
qiime metadata tabulate \
  --m-input-file qiime/classification.qza \
  --o-visualization qiime/classification.qzv
  
qiime taxa barplot \
  --i-table qiime/table.qza \
  --i-taxonomy qiime/classification.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization qiime/taxa_barplot.qzv  
  
qiime feature-table group \
  --i-table qiime/table.qza \
  --p-axis sample --p-mode sum \
  --m-metadata-file metadata.tsv \
  --m-metadata-column treatment-group \
  --o-grouped-table qiime/table_treatment_group.qza  

qiime taxa barplot \
  --i-table qiime/table_treatment_group.qza \
  --i-taxonomy qiime/classification.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization taxa_barplot_category.qzv  

#Eliminar secuencias de mitocondrias y cloroplastos (creo que no hace falta aqui)
qiime taxa filter-table \
  --i-table qiime/table.qza \
  --i-taxonomy qiimme/classification.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table qiime/table-no-contam.qza
  
qiime taxa filter-seqs \
  --i-sequences qiime/rep-seqs.qza \
  --i-taxonomy qiime/classification.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-sequences qiime/seq-no-contam.qza  
  
#Diversidad  
qiime diversity core-metrics-phylogenetic \
  --i-table qiime/table.qza \
  --i-phylogeny qiime/rep_seqs_aligned_masked_tree_rooted.qza \
  --p-sampling-depth 200 \
  --m-metadata-file metadata.tsv \
  --output-dir qiime/diversity
  
qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime/diversity/shannon_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization qiime/diversity/shannon_compare_groups.qzv  

qiime diversity beta-group-significance \
  --i-distance-matrix qiime/diversity/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column treatment-group \
  --p-pairwise \
  --o-visualization qiime/diversity/unweighted-unifrac-treatment-group-significance.qzv
 