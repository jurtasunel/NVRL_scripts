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
barcode_path="/home/josemari/Desktop/Jose/fix_fastqs/Aisling_fastqs/15_05_2023/"
barcode=("N230570024" "N230570032" "N230570034") 
prefix="Aisling"
Rscripts="/home/josemari/Desktop/Jose/General_scripts/R"
Bashscripts="/home/josemari/Desktop/Jose/General_scripts/Bash"

# Loop through the barcodes
for i in ${barcode[@]}; do

# Quality check for one barcode with artic and concatenate all fasqs in one single file.
echo "Running artic guppyplex for ${barcode_path}${i}"
artic guppyplex --min-length 300 --max-length 1400 --directory ${barcode_path}${i} --prefix ${prefix}

# Call the bash script to generate the depth plot from the concatenated fastq.
echo "Calling bowtie_depth.sh" 
bash ${Bashscripts}/bowtie_depth.sh ${prefix}_${i}.fastq

# Rename the resulting plot to match the other names.
mv ${prefix}_depth.pdf ${prefix}_${i}_depth.pdf

# Run medaka with artic primersto produce consensus sequence.
echo "Running artic medaka for ${barcode_path}${i}"
artic minion --medaka --medaka-model r941_min_high_g360 --normalise 200 --threads 8 --scheme-directory ~/artic-ncov2019/primer_schemes/midnight/ --read-file ${prefix}_${i}.fastq nCoV-2019/V1 ${prefix}_${i}

# Call the Rscript to generate the nucleotide changes csv report.
echo "Calling find_nt_diff.R" 
Rscript ${Rscripts}/find_nt_diff.R `pwd`/${prefix}_${i}.consensus.fasta `pwd`/${barcode}assembly.fa

# Rename the resulting csv report and the mafft alignment to match the other names.
mv Alignment_changes.csv ${prefix}_${i}_NTchanges.csv
mv mafft_aligned.fasta ${prefix}_${i}_mafft.fasta

# Move relevant files to a results folder.
mkdir ${prefix}_${i}_result
mv ${prefix}_${i}_depth.pdf ${prefix}_${i}_mafft.fasta ${prefix}_${i}_NTchanges.csv ${prefix}_${i}.consensus.fasta ${prefix}_${i}_result
# Remove intermediate files with same prefix_barcode name but skip the result folder. ${prefix}_${i}.align*
rm *.bam *.bai *.txt *.vcf *.hdf *.depths *.tbi *.er *.in.fasta *.out.fasta *preconsensus.fasta *.vcf.gz *.fastq fasta_to_align.fasta ${prefix}_low_depth_positions.csv

done

