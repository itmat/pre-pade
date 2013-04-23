RUM_INDEX_GENE_INFO=arabidopsis_genes.txt
BAM=sorted.bam
SAM=RUM.sam

exons : $(RUM_INDEX_GENE_INFO)
	python prepade/extractexons.py --rum-gene-info $< > $@

exon_quant_old : exons $(SAM)
	perl bin/quantify_exons.pl $^ $@

exon_quant_new : $(RUM_INDEX_GENE_INFO) $(BAM)
	python prepade/quantifyexons.py --rum-gene-info $(RUM_INDEX_GENE_INFO) $(BAM) > $@

exon_quant_old_filtered : exon_quant_old
	python prepade/csvfilter.py -w 'x["min"] > 0' $< -o $@ --int min
