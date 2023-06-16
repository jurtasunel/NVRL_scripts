### This script reads in taxonomy csv result and blast output result from the Illumina Taxonomy to produce a piechart.

# Libraries:
library(ggplot2)
library(RColorBrewer)

# Allow argument usage.
args = commandArgs(trailingOnly = TRUE)
# Print required input file if typed help.
if (args[1] == "-h" || args[1] == "help"){
  print("Syntax: Rscript.R taxon_results.csv blastOutput.tab")
  q()
  N
}

# Get input files.
input_file = args[1]
input_file2 = args[2]
# Read in the blast output.
data_df <- read.csv(input_file, stringsAsFactors = FALSE)
blast_out<- read.table(input_file2, sep="\t", stringsAsFactors = FALSE)
aln_reads <- length(unique(blast_out$V1)) # Get the number of aligned reads. One read can align to multiple orgs, that's why the unique is required.

print("Summarizing data for plotting...")
# Summarize the data for plotting. All organisms that represent less than 1% are grouped in "others".
summary_data <- data_df[data_df$Percentage > 1,]
summary_data <- summary_data[order(-summary_data$Percentage),]
others_data <- data_df[data_df$Percentage < 1,]
others <- c("", sum(others_data$Frequency), "", "others", sum(others_data$Percentage))
summary_data <- rbind(summary_data, others)
summary_data

# Get the desired number of colours.
number_of_colors <- nrow(summary_data)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(number_of_colors)

print("Plotting taxonomy piechart...")
# Create the pie chart:
taxon_pie <- ggplot(summary_data, aes(x="", y=Percentage, fill=Organism)) + geom_bar(stat="identity", width=1) + # basic bar plot.
  coord_polar("y", start=0) + geom_text(aes(x = 1.1, label = paste0(round(as.numeric(Percentage), 2), "%")),  # x defines the text distance from the center. 
                                        position = position_stack(vjust = 0.5),
                                        check_overlap = TRUE) + # If labels overlap, don't plot them.
  theme_void() +
  labs(x = NULL, y = NULL, fill = NULL) +
  scale_fill_manual(values = mycolors, breaks = as.character(summary_data$Organism)) +  # Rearrange legend on decreasing freq.
  labs(caption = (paste0("*There are ", nrow(others_data), " organisms with less than 1% representation grouped in 'Others' \n \n *There are ", aln_reads, " aligned reads."))) +
  theme(plot.margin = unit(c(10,5,10,5), "mm"))

ggsave("taxon_piechart.pdf", taxon_pie, width = 10, height = 10, dpi = 10, limitsize = FALSE)
