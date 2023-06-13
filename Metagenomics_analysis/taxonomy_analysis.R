### This script reads in an fmt6 blast output file and produces a taxonomy report file that consists on:
### A taxon_results.csv with the NCBI, taxonID and number of read matches for each organism.
### The script uses tazonomizr to do the analysis.

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(taxonomizr) # Manage ncbi taxonomy files as data bases for fast output.
# Documentation: https://cran.r-project.org/web/packages/taxonomizr/readme/README.html
# taxonomizr only needs to be set up once, then it functions can be used.
# Prepare the taxonomizr database.
prepareDatabase('/home/josemari/Desktop/Jose/Projects/MetaMIX/Scripts/accessionTaxa.sql')
accessionTaxa_sql = "/home/josemari/Desktop/Jose/Projects/MetaMIX/Scripts/accessionTaxa.sql"

# Allow argument usage.
args = commandArgs(trailingOnly = TRUE)
# Print required input file if typed help.
if (args[1] == "-h" || args[1] == "help"){
  print("Syntax: Rscript.R blastOutput.tab")
  q()
  N
}

# Get input file with depth from command line and print it.
input_file = args[1]
# Read in the blast output.
blastOut.default<- read.table(input_file, sep="\t", stringsAsFactors = FALSE)

print("Filtering blast output for read length < 50bp...")
# Filter the blast output to remove low length reads below 50 base pairs long.
blastOut.default <- blastOut.default[blastOut.default$V4 > 50, ]
# Get frequencies of each match and reorder the data frame.
acs_freq <- as.data.frame(table(blastOut.default$V2))
acs_freq <- acs_freq[order(-acs_freq$Freq),]

print("Doing taxonomy analysis with taxonomizr...")
# Get the tax_ID of the accession numbers.
taxaId <- accessionToTaxa(acs_freq$Var1,accessionTaxa_sql)
# Get the names of the organism on a variable.
taxonomy <- getTaxonomy(taxaId,accessionTaxa_sql)
organism <- as.character(taxonomy[, "species"])
print("1")
# Get the taxonomy result in a dataframe and change the columns.
tax_result <- cbind(acs_freq, taxaId, organism)
colnames(tax_result) <- c("NCBI_acs", "Frequency", "taxaID", "Organism")

print("Summarizing results to produce taxon_result.csv...")
# Make a plot df to summarize the tax results. Start storing the first line of the tax results.
filtered_res <- tax_result[1,]
# Loop through the lines of the data frame.
for (i in 2:nrow(tax_result)){
  current_row <- tax_result[i,]
  # If the organism on the current row already exist on tplot df:
  if (current_row$Organism %in% filtered_res$Organism){
    # Add the current frequency to the df plot frequency.
    filtered_res[filtered_res$Organism == current_row$Organism, "Frequency"] <- filtered_res[filtered_res$Organism == current_row$Organism, "Frequency"] + current_row$Frequency
  } else{
    # Add the new row to the df plot.
    filtered_res <- rbind(filtered_res, current_row)
  }
}
# Add a column with the percentage of each frequency.
filtered_res <- filtered_res %>% mutate(Percentage = 100 * Frequency/sum(Frequency))
# Write out csv with the filtered taxonomy results.
write.csv(filtered_res, file = "taxon_results.csv", row.names = FALSE)

