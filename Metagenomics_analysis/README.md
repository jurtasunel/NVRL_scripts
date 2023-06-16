# Metagenomics analysis from ILLUMINA paired end fastq files

- This pipeline will analyze Illumina paired-end fastq.gz files and produce a taxonomy report with the most representative taxa present on each sample. The pipeline will make a Results folder for every sample that contains a taxon_piechart.pdf plot, taxon_results.csv file, paired end fastq merged pefq_merged.fasta and the blast output blastn_output.tab file.

- The taxonomy.sh script is the parent script to run. It gets the unique names of the paired end files on an array, removes adapters with AdapterRemoval, converts merged fastq to fasta with seqtk and blasts the fasta with blastN. ALL FASTQ.GZ FILES MUST BE ON THE SAME DIRECTORY SPECIFIED ON THE taxonomy.sh SCRIPT.

- The taxonomy_analysis.R script reads in the blast output and produces a csv with the organisms for the alligned reads, their frequency and their percentage.

- The plot_piechart.R script reads in the blast output and the taxon csv and produces a piechart with the most frequent organisms of that sample.

- UP TO THIS POINT, THE PIPELINE RUNS AUTOMATICALLY. From the taxon Result and taxon piechart, the most relevant organisms can be chosen for the next step.

- Do a simple alignment of the fastq files of any desired sample with Bowtie2 and get the depth.txt file using Samtools.

- The depth.txt file, the organism reference .fasta file and the annotation .gff3 file are the inputs to the plot_depth.R script. THE PLOTTING WILL NEED TO BE RESCALATED MANUALLY DEPENDING ON THE ORGANISM'S GENOME LENGTH.



