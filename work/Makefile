RUM_INDEX_GENE_INFO=arabidopsis_genes.txt
BAM=sorted.bam
SAM=RUM.sam
RUM=~/src/rum

all : exon_quant_new exon_quant_old_filtered

exons : $(RUM_INDEX_GENE_INFO)
	python -m prepade.extractexons --rum-gene-info $< > $@

exon_quant_old : exons $(SAM)
	perl $(RUM)/bin/quantify_exons.pl $^ $@

exon_quant_new : $(RUM_INDEX_GENE_INFO) $(BAM)
	python -m prepade.quantifyexons --rum-gene-info $^ -o $@

exon_quant_old_filtered : exon_quant_old
	python -m prepade.csvfilter -w 'x["min"] > 0' $< -o $@ --int min

clean :
	rm -f exons exon_quant_old exon_quant_new exon_quant_old_filtered
