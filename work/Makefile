# Started 13:31:46
#RUM_INDEX_GENE_INFO=mm10_refseq_ucsc_gene_info.txt
RUM_INDEX_GENE_INFO=arabidopsis_genes.txt
BAM=rum_sorted.bam
SAM=RUM.sam
RUM=~/src/rum

rum.bam : RUM.sam
	samtools view -bS $< > $@

rum_sorted.bam : rum.bam
	samtools sort $< rum_sorted

rum_sorted.bam.bai : rum_sorted.bam
	samtools index $< $@

all : exon_quant_new exon_quant_old_filtered

exons : $(RUM_INDEX_GENE_INFO)
	python -m prepade.extractexons --rum-gene-info $< > $@

exon_quant_old : exons $(SAM)
	perl $(RUM)/bin/quantify_exons.pl $^ $@

exon_quant_new : exons $(BAM) $(BAM).bai
#	python -m prepade.quantifyexons --exon-index exon_index --rum-gene-info $^ -o $@
	python -m prepade.quantifyexons --exon-index exon_index --exons exons $(BAM) -o $@

exon_quant_old_filtered : exon_quant_old
	python -m prepade.csvfilter -w 'x["min"] > 0' $< -o $@ --int min

joined : exon_quant_old_filtered exon_quant_new
	python -m prepade.csvjoin -j feature $^ -o $@

comparison : exon_quant_old_filtered exon_quant_new
	python -m comparequants $^ -o $@ --diffs diffs

clean :
	rm -f exons exon_quant_old_filtered exon_quant_new exon_quant_old_filtered diffs

exon_quant.diff : exon_quant_old_filtered exon_quant_new
	diff -w $^ >$@
