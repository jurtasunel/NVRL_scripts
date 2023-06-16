#!/bin/bash

### This script reads in paired end fastq files and performs a taxonomy analysis.
### Documentation for adapterRemoval: https://adapterremoval.readthedocs.io/en/stable/examples.html

### Define variables:
# Export blastn database for viral genomes and get the database prefix_ID.
export BLASTDB=/home/josemari/Desktop/Jose/Reference_sequences/Reference_virus_DB/ref_viruses_rep_genomes
blastN_db="ref_viruses_rep_genomes"
# Get the path to the paired end fastq files.
data_path="/home/josemari/Desktop/Metagenomics_analysis/Data"
# Define the illumina suffix for forward and reverse reads.
Fr_suffix="_R1_001.fastq.gz"
Rv_suffix="_R2_001.fastq.gz"
# Add the names of all the files to an array by:
# Creating an associative array that only stores unique values (so paired end names on two files don't repeat).
declare -A files_hash
# Loop through the files in the directory.
for file in "$data_path"/*; do
  # Remove the illumina suffixes.
  if [[ $file == *"$Fr_suffix"* ]]; then
    filename="${file/$Fr_suffix/}"
    files_hash["$(basename "${filename}")"]=1
  elif [[ $file == *"$Rv_suffix"* ]]; then
    filename="${file/$Rv_suffix/}"
    files_hash["$(basename "${filename}")"]=1
  fi
done
# Create an array to store the unique extracted strings stored on the associative array.
files=("${!files_hash[@]}")
echo -e "Fastq tags: ${files[@]}\n"

# Loop through the fastq tags.
for i in "${files[@]}"; do

  # DO TO MERGE RAW PEFQ WITH NO ADAPTERS REMOVAL USING FLASH:
  #echo -e "Merging ${i}${Fr_suffix} ${i}${Rv_suffix} with flash...\n"
  #flash ${data_path}/${i}${Fr_suffix} ${data_path}/${i}${Rv_suffix}
  #echo -e "Calling merge_fastq.R for ${i} to produce ${i}.fasta...\n"
  #Rscript merge_fastq.R `pwd`/out.extendedFrags.fastq
  
  # DO TO MERGE RAW PEFQ INTO SINGLE FASQ AFTER ADAPTERS REMOVAL USING AdapterRemoval:
  echo -e "Removing adapters and merging ${i}${Fr_suffix} ${i}${Rv_suffix} with AdapterRemoval...\n"
  AdapterRemoval --file1 ${data_path}/${i}${Fr_suffix} --file2 ${data_path}/${i}${Rv_suffix} --basename output_paired --trimns --trimqualities --collapse
  
  # Converting the merged fastq to fasta file.
  echo -e "Converting the ${i} merged fastq to fasta with seqtk...\n" 
  seqtk seq -a output_paired.collapsed > pefq_merged.fasta
  
  # Nucleotide Blast the fasta file. The outfmt string specifies which columns to get on the output. The -evalue 0.01 is a conservative value significant match reads.
  echo -e "Blasting ${i} with blastn...\n"
  blastn -db ${blastN_db} -query pefq_merged.fasta -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore slen" -max_target_seqs 10 -evalue 0.01 -out blastn_output.tab
  
  # Perform the taxonomy analysis.
  echo -e "Calling taxonomy_analysis.R...\n"
  Rscript taxonomy_analysis.R `pwd`/blastn_output.tab
  
  # Plot the taxonomy analysis.
  echo -e "Calling plot_piechart.R...\n"
  Rscript plot_piechart.R `pwd`/taxon_results.csv `pwd`/blastn_output.tab
  
  # Make a result directory for the sample, move relevant result files and remove intermediates.
  mkdir ${i}_result
  mv taxon_piechart.pdf taxon_results.csv pefq_merged.fasta output_paired.collapsed blastn_output.tab ${i}_result
  rm out*
  
  echo -e "SAMPLE ${i} analysis finished\n"
  
done
