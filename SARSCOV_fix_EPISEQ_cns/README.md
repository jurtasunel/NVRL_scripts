# SARS-CoV-2 fix Episeq consensus sequence.

- This pipeline makes a consensus sequence of SARS-CoV-2 using the artic pipeline, compares it with another consensus (normally from EPISEQ) and produces a csv report with the nucleotides that are different between each consensus and a depth plot of the sample.

- fast_analysis_pipeline.sh is the parent script to run. It calls the bash script bowtie_depth.sh that depends on the Rscript plot_depth.R to produce a depth plot, then makes a consensus sequence using medaka from the artic pipeline and finally calls the Rscript find_nt_diff.R to produce a multiple sequence alignment and a csv report with the different nucleotides. 

- The Data should be on a directory with the EPISEQ consensus named {BARCODE}assembly.fa and a folder containing all the fastqs/fastq.gz files from ONT with the barcode as name of the folder.

- The pipeline will produce a Results folder containing the new generated consensus.fasta, a depth.pdf plot, a multiple sequence alignemnt mafft.fasta file with the EPISEQ consensus and the NTchanges.csv report.

- IN ORDER TO WORK, THE ARTIC PIPELINE MUST BE ACTIVATED WITH "source activate artic-ncov2019".
