# Produce a csv report named "Lineage_change_DATE.csv" with the lineages that have changed on the last update of pangolin.

- Pango_Lineages.sh is the parent script to run. It calls pangolin_update.sh to update pangolin first, then initializes conda and finally calls Pango_Lineages.R to do the rest of the analysis. It then finishes moving the resulting report "Lineage_change_DATE.csv" to a Reports folder.

- Pango_Lineages.R uses the constants that exist on the Pango_Lineages_source.R and sends an email with the "Lineage_change_DATE.csv" attached.

- In order to find the input files, Pango_Lineages_source.R requires the input files in specific format specified on the header of the script.

- This scripts work by updating a mysql database with only one table. The structure of that database is on the Pango_Lineages_DB.txt file.

 -THE PATH TO THE REPORT FOLDER MUST BE SPECIFIED ON THE LAST LINE OF Pango_Lineages.sh script.
 
 -THE PATH TO THE DATA FOLDER MUST BE SPECIFIED ON THE Pango_Lineages_source.R script.
 
 -THE REQUIRED DATA ARE FOUR FILES DOWNLOADED DIRECTLY FROM GSAID. SELECT THE DESIRED SAMPLES AND DOWNLOAD THE FOLLOWING (Nucleotide Sequences (FASTA), Patient status metadata, Dates and Location and Input for the Augur pipeline).
 
 -INPUT FOR THE AUGUR PIPELINE IS A .tar FILE, EXTRACT IT AND FROM THE TWO FILES THERE, TAKE THE .metadata.tsv FILE AND PLACE IT WIHT THE OTHER THREE ON THE DATA FOLDER.
 
 -RENAME Patient status metadata AND Dates and Location TO END IN _A.tsv and _B.tsv RESPECTIVELY.
  

