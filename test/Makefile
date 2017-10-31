# USAGE

# Test Purge Haplotigs:
	# make all

# Clean up:
	# make clean


THREADS = 4



readhist: aligned.bam.histogram.png aligned.bam.genecov

contigcov: coverage_stats.csv

purge: curated.artefacts.fasta curated.contig_associations.log curated.fasta curated.haplotigs.fasta curated.reassignments.tsv



test: clean readhist contigcov purge

clean:
	rm -rf dotplots_* tmp_purge_haplotigs curated.* aligned.bam.genecov aligned.bam.histogram.png aligned.bam.bai contigs.fa.fai coverage_stats.csv 



aligned.bam.histogram.png aligned.bam.genecov: aligned.bam
	purge_haplotigs readhist $<

coverage_stats.csv: aligned.bam.genecov
	purge_haplotigs contigcov -i $< -l 3 -m 20 -h 50

curated.fasta curated.artefacts.fasta curated.contig_associations.log curated.haplotigs.fasta curated.reassignments.tsv: contigs.fa coverage_stats.csv aligned.bam
	purge_haplotigs purge -g $< -c coverage_stats.csv -b aligned.bam -t $(THREADS)
	@cmp curated.fasta validate.fa && echo "Purge Haplotigs successfully validated"


