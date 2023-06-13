# Metagenomics analysis from ILLUMINA paired end fastq files

All fastq.gz files must be on the same directory specified on taxonomy.sh.

The taxonomy.sh script is the parent script to run. It gets the unique names of the paired end files on an array, removes adapters with AdapterRemoval, converts merged fastq to fasta with seqtk and blasts the fasta with blastN.

The taxonomy_analysis.R script reads in the blast output and produces a csv with the organisms for the alligned reads, their frequency and their percentage.

The plot_piechart.R script reads in the blast output and the taxon csv and produces a piechart with the most frequent organisms of that sample.
