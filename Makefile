RUM_INDEX_GENE_INFO=arabidopsis_genes.txt
BAM=sorted.bam
SAM=RUM.sam

exons : $(RUM_INDEX_GENE_INFO)
	python prepade/extractexons.py --rum-gene-info $< > $@

exon_quant_old : exons $(SAM)
	perl bin/quantify_exons.pl $^ $@
