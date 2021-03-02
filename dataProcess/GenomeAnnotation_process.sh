#!/bin/sh
### human genome annotation were downloaded from Ensembl

awk 'BEGIN{FS="gene_version"}{$2="";print $0}' Homo_sapiens.GRCh37.87.gtf > genome1

awk 'BEGIN{FS="\t"}{if($3=="exon")print $0}' genome1 > genome2

awk 'BEGIN{FS="\t"}{print "chr"$0}' genome2 > genome3

awk 'BEGIN{FS="\t";OFS="\t"}{print $NF,$4,$5}' genome3 |sort -k1,2n > genome4

### data used for bedtools merge(genome4)
#gene_id "ENSG00000000003";  	99883667	99884983
#gene_id "ENSG00000000003";  	99885756	99885863
#gene_id "ENSG00000000003";  	99887482	99887565
#gene_id "ENSG00000000003";  	99887538	99887565
