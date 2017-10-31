exit; 
# If you're testing your installation of Purge Haplotigs then use the makefile
make test # to validate
make clean # to clean up

# Validation and test dataset
# The test dataset is from Redundans: https://github.com/lpryszcz/
# This isn't the best set to use for Purge Haplotigs but it will do for now


# The reads were mapped like so (intermediate files omitted here):
bwa index contigs.fa
bwa mem -x ont2d -t 16 contigs.fa nanopore.fa.gz | samtools view -bs - > nanopore.bam
bwa mem -x pacbio -t 16 contigs.fa pacbio.fq.gz | samtools view -bs - > pacbio.bam
cat nanopore.bam pacbio.bam | samtools sort - > aligned.bam


# TEST PURGE HAPLOTIGS
purge_haplotigs readhist aligned.bam

purge_haplotigs contigcov -i aligned.bam.genecov -l 3 -m 20 -h 50

purge_haplotigs purge -g contigs.fa -c coverage_stats.csv -b aligned.bam -t 16

