# Set working directory to the scripts location.
setwd("path/to/scripts")
### Libraries.
source("Pango_Lineages_source.R")
### Functions
# First data upload to empty database.
first_upload <- function (){
  for (i in 1:length(fasta_file)){
    # Get the relevant information from the fasta and metadata files.
    GSAID_ID <- names(fasta_file)[i]
    seq <- fasta_file[[i]]
    last_lin <- metadata_file[metadata_file$Accession.ID == GSAID_ID, "Lineage"]
    GSAID_colldate <- metadata_file[metadata_file$Accession.ID == GSAID_ID, "Collection.date"]
    GSAID_subdate <- subdate_file[subdate_file$Accession.ID == GSAID_ID, "Submission.date"]
    submitting_lab <- sublab_file[sublab_file$gisaid_epi_isl == GSAID_ID, "submitting_lab"]
    # Make the sql query.
    sql_query = paste0('INSERT INTO Pango_Lineages (sequence, last_lineage, GSAID_ID, GSAID_colldate, GSAID_subdate, submitting_lab) VALUES ("', seq, '", "', last_lin, '", "', GSAID_ID, '", "', GSAID_colldate, '", "', GSAID_subdate, '", "', submitting_lab, '")')
    # Connect to MySQL, execute the query and clear the result.
    rs <- dbSendQuery(con_sql, sql_query)
    dbClearResult(rs);
  }
}
# Upload sequences whose ID is not already on the database.
upload_newseqs <- function(){
  # Get all the IDs existing on the database.
  sql_query = "SELECT GSAID_ID FROM Pango_Lineages"
  rs <- dbSendQuery(con_sql, sql_query)
  # Store the IDs on a temporal variable.
  existing_IDs <- dbFetch(rs)
  dbClearResult(rs)
  print("new IDs")
  # The new IDs are the downloaded IDs from GSAID that don't exist on the database.
  new_IDs <- metadata_file$Accession.ID[!(metadata_file$Accession.ID %in% existing_IDs$GSAID_ID)]
  # Loop through the new IDs.
  for (i in new_IDs){
    # Extract the relevant information.
    GSAID_ID <- i
    seq <- fasta_file[[GSAID_ID]]
    # Print a warning if the sequence is shorter than 100 base pairs.
    if(nchar(seq) < 100){
      print(paste0("There is a problem with the sequece of ", GSAID_ID))
    }
    last_lin <- metadata_file[metadata_file$Accession.ID == i, "Lineage"]
    GSAID_colldate <- metadata_file[metadata_file$Accession.ID == i, "Collection.date"]
    GSAID_subdate <- subdate_file[subdate_file$Accession.ID == i, "Submission.date"]
    submitting_lab <- sublab_file[sublab_file$gisaid_epi_isl == i, "submitting_lab"]
    # Make the sql query.
    sql_query = paste0('INSERT INTO Pango_Lineages (sequence, last_lineage, GSAID_ID, GSAID_colldate, GSAID_subdate, submitting_lab) VALUES ("', seq, '", "', last_lin, '", "', GSAID_ID, '", "', GSAID_colldate, '", "', GSAID_subdate, '", "', submitting_lab, '")')
    # Connect to MySQL, execute the query and clear the result.
    rs <- dbSendQuery(con_sql, sql_query)
    dbClearResult(rs);
 }
}
# Generate fasta from database with last month sequences.
month_seqs_fasta <- function(){
  # Create the sql query to get sequences on a time period of 31 days.
  sql_query = "SELECT GSAID_ID, sequence FROM Pango_Lineages WHERE GSAID_colldate >= CURRENT_DATE - INTERVAL 31 DAY"
  # Send the query and get the results into a variable.
  rs <- dbSendQuery(con_sql, sql_query)
  result <- dbFetch(rs)
  dbClearResult(rs)
  # Get the sequences on a list to match the format required by the write fasta format.
  seqs <- list()
  for (i in 1:length(result$sequence)) {seqs[[i]] <- result$sequence[i]}
  # Write the fasta file.
  write.fasta(sequences = seqs, names =  result$GSAID_ID, file.out = "pango.fasta")
}
# Add the new pangolin lineage to all sequences on the database.
fill_newlin <- function() {
  # Loop through the pangolin report.
  for (i in 1:nrow(pangolin_report)){
    # Get the ID and pangolin lineage of each entry of the fasta file.
    ID <- pangolin_report$taxon[i]
    lineage <- pangolin_report$lineage[i]
    # Make the sql query to update the new lineage column with the pango lineage for each ID.
    sql_query = paste0("UPDATE Pango_Lineages SET new_lineage = '", lineage, "' where (GSAID_ID = '", ID, "')")
    rs <- dbSendQuery(con_sql, sql_query)
    dbClearResult(rs)
  }
}
# Make a function that returns a data frame with the lineage changes.
report_lin <- function(){
  # Make query to get the old lineage, new lineage and ID of all entries with different lineages.
  sql_query = "SELECT GSAID_ID, last_lineage, new_lineage, submitting_lab FROM Pango_Lineages WHERE last_lineage <> new_lineage"
  rs <- dbSendQuery(con_sql, sql_query)
  result <- dbFetch(rs)
  dbClearResult(rs)
  # Make a data frame with the IDs, old and new lineages and submitting laboratory.
  lin_report <- data.frame(result$GSAID_ID, result$last_lineage, result$new_lineage, result$submitting_lab)
  # Subset the report to avoid updating lineages to "unassigned".
  unassigned_removed <- lin_report[lin_report$result.new_lineage != "Unassigned",]
  # Update the old lineage with the new one for the different ones on the database.
  for (i in 1:nrow(unassigned_removed)){
    sql_query = paste0("UPDATE Pango_Lineages SET last_lineage = '", unassigned_removed$result.new_lineage[i], "' where (GSAID_ID = '", unassigned_removed$result.GSAID_ID[i], "')")
    rs <- dbSendQuery(con_sql, sql_query)
    dbClearResult(rs)
  }
  # Return the lineage report.
  return(lin_report)
}

### Fix missing dates:
# Loop through the metadata files.
for (i in 1:nrow(metadata_file)){
  
  # Get the collection dates.
  i.metadata_date <- metadata_file$Collection.date[i]

  # Add day if the dates have less than 10 characters (normal date format has 10 yyyy-mm-dd)
  if (nchar(i.metadata_date) < 10){
    fix_date <- paste0(i.metadata_date, "-01")
    metadata_file$Collection.date <- fix_date
    subdate_file$Collection.date <- fix_date
  }
  
}

### Workflow
# Upload new sequences to database and produce fasta file with last month sequences.
db_connect()
print("Uploading new sequences to mysql database.")
upload_newseqs()
print("Generating fasta file of sequences within one month.")
month_seqs_fasta() # Produces pango.fasta as output.

# Run pangolin to get new lineages of fasta file sequences.
# Documentation: https://cov-lineages.org/resources/pangolin.html
print("Running pangoling, this might take several minutes...")
system("pangolin pango.fasta") # Produces lineage_report.csv as output.

# Read in the pangolin report and get the pangolin version and pango data version.
pangolin_report <- read.csv("lineage_report.csv", stringsAsFactors = FALSE)
pangolin_version <- pangolin_report$pangolin_version[1]
pangolearn_vstring <- pangolin_report$version[1]
# Use regexp to match the regular expression for the version, and the regmatches to extract the match from the original string.
pangolearn_version <- regmatches(pangolearn_vstring, regexpr("\\d+\\.\\d+", pangolearn_vstring, ignore.case = TRUE))
# Fill the new lineages to the database.
print("Adding new pangolin lineage to mysql database.")
fill_newlin()
# Store the output of the report into a variable.
print("Generating lineage change report and updating old lineage.")
lin_report <- report_lin()
db_disconnect()
colnames(lin_report) <- c("GSAID_ID", "previous_lineage", "updated_lineage", "submitting_lab")

# Initialize vectors to store the originating lab.
originating_lab <- c()
#Loop through the report rows ONLY if there are more than 0.
if (nrow(lin_report) > 0){
  for (i in 1:nrow(lin_report)){
    # Get the current GSAID_ID.
    current_ID <- lin_report$GSAID_ID[i]
    # Get the originating lab.
    originating_lab <- c(originating_lab, sublab_file[sublab_file$gisaid_epi_isl == current_ID, "originating_lab"])
  }
}
# Bind the originating lab vector to the lineage report data frame.
lin_report <- cbind(lin_report, originating_lab)

# Write the final csv report with the lineage changes.
write.csv(lin_report, paste0("Lineage_change_", Sys.Date(), ".csv"), row.names = FALSE)

# Configure the gmail credentials.
gm_auth_configure(path = "path/to/GmailCredentials.json")

# Create update email
update_email <- gm_mime() %>%
  # Use line below for testing.
  #gm_to("jose.urtasunelizari@ucd.ie") %>%
  # Use line below for OOS mailing list.
  gm_to("jose.urtasunelizari@ucd.ie") %>%
  gm_from("jose.urtasunelizari.reporting@gmail.com") %>%
  gm_subject(paste0("Pangolin ", pangolin_version, " (PangoLEARN ", pangolearn_version, ") Lineage update")) %>%
  gm_text_body("Hi,
  
Please find attached a csv file with the lineages that have been updated on the current version of pangoling.

Disclaimer: This is an automated message, if you have any
questions please contact me at jose.urtasunelizari@ucd.ie

Best wishes,
Josemari") %>%
  # Use line below for testing and first time run of the script.
  gm_attach_file(paste0("path/to/Scripts/Lineage_change_", Sys.Date(), ".csv"))
  # Use line below for second time run after the bash script has already moved the file to the report folder.
  #gm_attach_file(paste0("/home/gabriel/Desktop/Jose/Projects/Lineage_changes/Reports/Lineage_change_", Sys.Date(), ".csv"))

# Send email.
gm_send_message(update_email)

# Email list.
# nicola.mellotte@hpsc.ie, aislingrichmond@beaumont.ie, Wendy.Brennan@hse.ie, joannem.king@hse.ie, Aoife.Warde@hse.ie, PGRIER@stjames.ie, jonathan.dean@ucd.ie, eric.brown@ucdconnect.ie, Laura.MurilloRuiz@hse.ie, catherine.dempsey@hse.ie, sarah.mcgarry@hpsc.ie, carina.brehony@hpsc.ie, niamh.murphy@hpsc.ie, Celine.Higgins@hse.ie, Adele.Habington@olchc.ie, Karen.johnston@olchc.ie, soracha.mckinley@hpsc.ie, daniel.hare@ucd.ie, Patrick.Stapleton@hse.ie, gabo.gonzalez@ucd.ie






