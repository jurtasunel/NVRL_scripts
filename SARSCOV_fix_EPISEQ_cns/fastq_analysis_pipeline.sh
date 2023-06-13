#!bin/bash
### This script runs the basic artic pipeline steps to produce a medaka consensus. Then it calls the bowtie_depth.sh to procude a depth plot and the find_nt_diff.R to generate the NTchanges.csv report.
### It requires a directory with an alternative consensus file (normally EPISEQ) named {BARCODE}assembly.fa and a folder containing all the fastqs/fastq.gz files from ONT. The folder must have the barcode as name.
 
### Documentation:
# Artic instalation: https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
# Primer shcemes: https://github.com/phac-nml/primer-schemes
# Primer schemes documentation: https://psy-fer.github.io/interARTIC/primers/
# Artic version of muscle: https://ubuntu.pkgs.org/20.04/ubuntu-universe-amd64/muscle_3.8.1551-2build1_amd64.deb.html
# Muscle conda error: https://github.com/artic-network/artic-ncov2019/issues/93
# Artic common issues: https://github.com/artic-network/artic-ncov2019/issues

#### This script requires the artic environment to be activated with "source activate artic-ncov2019".

# Get the required variables.
data_path="/home/josemari/Desktop/OOS_fastqs/Data/"
barcodes=("N230270082" "N230270084") 
prefix="25052023"

# Loop through the barcodes
for i in ${barcodes[@]}; do

# Quality check for one barcode with artic and concatenate all fasqs in one single file.
echo "Running artic guppyplex for ${data_path}${i}"
artic guppyplex --min-length 300 --max-length 1400 --directory ${data_path}${i} --prefix ${prefix}

# Call the bash script to generate the depth plot from the concatenated fastq.
echo "Calling bowtie_depth.sh" 
bash bowtie_depth.sh ${prefix}_${i}.fastq

# Run medaka with artic primersto produce consensus sequence.
echo "Running artic medaka for ${data_path}${i}"
artic minion --medaka --medaka-model r941_min_high_g360 --normalise 200 --threads 8 --scheme-directory ~/artic-ncov2019/primer_schemes/midnight/ --read-file ${prefix}_${i}.fastq nCoV-2019/V1 ${prefix}_${i}

# Call the Rscript to generate the nucleotide changes csv report.
echo "Calling find_nt_diff.R" 
Rscript find_nt_diff.R `pwd`/${prefix}_${i}.consensus.fasta ${data_path}${barcodes}assembly.fa

# Rename the resulting csv report and the mafft alignment to match the other names.
mv Alignment_changes.csv ${i}_NTchanges.csv
mv mafft_aligned.fasta ${i}_mafft.fasta

# Move relevant files to a results folder.
mkdir ${prefix}_${i}_result
mv ${i}.fastq_depth.pdf ${i}_mafft.fasta ${i}_NTchanges.csv ${prefix}_${i}.consensus.fasta ${i}.fastq_low_depth_positions.csv ${prefix}_${i}_result
# Remove intermediate files with same prefix_barcode name but skip the result folder. ${prefix}_${i}.align*
rm *.bam *.bai *.txt *.vcf *.hdf *.depths *.tbi *.er *.in.fasta *.out.fasta *preconsensus.fasta *.vcf.gz *.fastq fasta_to_align.fasta 

done

