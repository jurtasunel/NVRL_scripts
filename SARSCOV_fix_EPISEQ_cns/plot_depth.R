### This script creates a pdf with bar plot of a depth txt file output from samtools depth.
library(ggplot2)

# Allow argument usage.
args = commandArgs(trailingOnly = TRUE)
# Print required input file if typed help.
if (args[1] == "-h" || args[1] == "help"){
  print("Syntax: Rscript.R depth_file.txt")
  q()
  N
}

# Get input file with depth from command line and print it.
input_file = args[1]
# Separate the file name from the path.
ID = tail(unlist(strsplit(input_file, "/")), n = 1)
# Get the fastq file name and print it.
ID = unlist(strsplit(ID, "_", fixed = TRUE))[1]
print(paste0("Plotting ", ID, " depth..."))
# Read the file as table and change column names.
depth_table <- read.table(input_file,  sep = '\t', header = FALSE)
colnames(depth_table) <- c("Chrom", "Position", "Depth")

# Fill a vector with a depth score for each position, and save bad score position.
Depth_score <- c()
BC_positions <- c()
for (i in 1:nrow(depth_table)){
  
  if (depth_table$Depth[i] < 100){
    Depth_score <- c(Depth_score, "C (<100)")
    BC_positions <- c(BC_positions, depth_table$Position[i])
    
  } else if (depth_table$Depth[i] < 200){
    Depth_score <- c(Depth_score, "B (100-200)")
    BC_positions <- c(BC_positions, depth_table$Position[i])
    
  } else {Depth_score <- c(Depth_score, "A (>200)")}
  #print(tail(Depth_score, n=1))
}

# Add a column to the depth table with the depth scores vector.
depth_table <- cbind(depth_table, Depth_score)

# Get ranges of values of B an C positions.
BC_range <- split(BC_positions, cumsum(c(1, diff(BC_positions) != 1)))
# Get first and last value of each range.
Low_depth_positions <- c()
for (i in 1:length(BC_range)){
  current_range <- as.character(unlist(BC_range[i]))
  if (length(current_range) == 1){
    Low_depth_positions <- c(Low_depth_positions, current_range)
  } else{
    Low_depth_positions <- c(Low_depth_positions, paste0(head(current_range, 1), "-", tail(current_range, 1)))
  }
}
# Write the ranges of values for low depth positions to a csv file.
Low_depth_positions <- data.frame(Low_depth_positions)
write.csv(Low_depth_positions, file = paste0(ID,"_low_depth_positions.csv"), row.names = FALSE)

# Plot depth table.
p <- ggplot(depth_table, aes(x = Position, y = Depth, fill = Depth_score)) +
  geom_bar(stat = "identity") + # Fill identity for ggplot to accept the y axis data on geom bar.
  theme_minimal() + # Remove background grid.
  scale_fill_manual(values = c("#339966", "#3333FF", "#CC0033")) +
  theme(legend.position = "top") +
  ggtitle(ID) +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

# Create pdf to save the plot.
pdf(paste0(ID,"_depth.pdf"))
print(p) # Save plot on first page
#print(p2)
#print(p3).. for saving on next pages of pdf.
  
