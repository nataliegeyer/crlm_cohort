# Globals
base_dir <-"/home/bibu/Workspace/crlm_cohort/output"
parsed_csv_dir <- paste(base_dir, "/parsed_annotations", sep="")
combined_fn <- paste(base_dir, "/combined_annotations.csv", sep="")

# Remove previous combined file
if (file.exists(combined_fn)) {
  file.remove(combined_fn)
}

# Combined dataframe
df <- data.frame()
# List of previously parsed cvs files with annotation data, one for each slide
file_list <- list.files(path = parsed_csv_dir, pattern = "*.csv")
# Iterate csv files and combine in one large dataframe
for (fn in file_list){
    temp_df <-read.csv2(paste(parsed_csv_dir, "/",fn, sep=""), row.names=NULL)
    
    df <- rbind(df, temp_df)
}
# Write combined dataframe as csv
write.csv(df, combined_fn, row.names=FALSE)