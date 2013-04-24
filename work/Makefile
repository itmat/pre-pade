
RUM_INDEX_GENE_INFO=mm10_refseq_ucsc_gene_info.txt
BAM=rum.sorted.bam
SAM=RUM.sam
RUM=~/src/rum

rum.bam : RUM.sam
	samtools view -bS $< > $@

rum.sorted.bam : rum.bam
	samtools sort $< rum.sorted

rum.sorted.bai : rum.sorted.bam
	samtools index $< $@

all : exon_quant_new exon_quant_old_filtered

exons : $(RUM_INDEX_GENE_INFO)
	python -m prepade.extractexons --rum-gene-info $< > $@

exon_quant_old : exons $(SAM)
	perl $(RUM)/bin/quantify_exons.pl $^ $@

exon_quant_new : $(RUM_INDEX_GENE_INFO) $(BAM)
	python -m prepade.quantifyexons --exon-index exon_index --rum-gene-info $^ -o $@

exon_quant_old_filtered : exon_quant_old
	python -m prepade.csvfilter -w 'x["min"] > 0' $< -o $@ --int min

joined : exon_quant_old_filtered exon_quant_new
	python -m prepade.csvjoin -j feature $^ -o $@

comparison : exon_quant_old_filtered exon_quant_new
	python -m comparequants $^ -o $@ --diffs diffs

clean :
	rm -f exons exon_quant_old_filtered exon_quant_new exon_quant_old_filtered

exon_quant.diff : exon_quant_old_filtered exon_quant_new
	diff -w $^ >$@
