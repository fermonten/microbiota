# Extraer los archivos fasta y biom a partir de los qza de qiime2
qiime tools export --input-path ../qiime/rep-seqs.qza  --output-path .
qiime tools export --input-path ../qiime/table.qza --output-path .

place_seqs.py -s dna-sequences.fasta -o out.tre -p 1 \
              --intermediate intermediate/place_seqs



hsp.py -i 16S -t out.tre -o marker_predicted_and_nsti.tsv.gz -p 1 -n

hsp.py -i EC -t out.tre -o EC_predicted.tsv.gz -p 1



#zless -S marker_predicted_and_nsti.tsv.gz

metagenome_pipeline.py -i feature-table.biom -m marker_predicted_and_nsti.tsv.gz -f EC_predicted.tsv.gz \
                       -o EC_metagenome_out --strat_out
			   
					   
#zless -S EC_metagenome_out/pred_metagenome_unstrat.tsv.gz

convert_table.py EC_metagenome_out/pred_metagenome_contrib.tsv.gz \
                 -c contrib_to_legacy \
                 -o EC_metagenome_out/pred_metagenome_contrib.legacy.tsv.gz
				 
pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_contrib.tsv.gz \
                    -o pathways_out -p 1
					

add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                    -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o pathways_out/path_abun_unstrat_descrip.tsv.gz
					
#KEGG ORTOLOGY
					
hsp.py -i KO -t out.tre -o KO_predicted.tsv.gz -p 1
										   
metagenome_pipeline.py -i feature-table.biom -m marker_predicted_and_nsti.tsv.gz -f KO_predicted.tsv.gz \
                       -o KO_metagenome_out --strat_out		
					   
convert_table.py KO_metagenome_out/pred_metagenome_contrib.tsv.gz \
                 -c contrib_to_legacy \
                 -o KO_metagenome_out/pred_metagenome_contrib.legacy.tsv.gz
				 
pathway_pipeline.py -i KO_metagenome_out/pred_metagenome_contrib.tsv.gz \
                    -o pathways_out -p 1
					
add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
                    -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m KO \
                    -o pathways_out/path_abun_unstrat_descrip.tsv.gz
					