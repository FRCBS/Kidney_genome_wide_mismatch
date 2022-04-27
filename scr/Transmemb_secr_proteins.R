###############################################################################
# Create a list of transmembrane and secretory protein ENST-IDs
###############################################################################
# Importing dataset
Transmembrane_secretory_separated <- read_excel("~/data/Transmembrane_secretory_separated.xlsx", 
                                                +     col_names = FALSE)

# Deleting extra columns
Transmem_secr_proteins <- subset(Transmembrane_secretory_separated, select = -c(...1, ...3, ...4, ...5, ...6, ...7))

# Gathering ensemble-id columns into one column
# As a result: 
# the first column is gene name, second is the column number from the previous data 
# and third is ENST-id
Transcript_list <- gather(Transmem_secr_proteins, "Column", "...8", -...2, na.rm = TRUE)

# Rename two columns
Transcript_list <- Transcript_list %>% rename(Protein = ...2, Ensemble = ...8)

# Write out the table
write.table(Transcript_list, file = "/data/Transcript_list", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

###############################################################################
# Creating a variant list from secretory and transmembrane proteins
# Import the missense variant list you created with VEP
ALL_MISSENSE_VARIANTS_cleaned <- read_delim("~/data/ALL_MISSENSE_VARIANTS_cleaned.txt", 
                                            +     "\t", escape_double = FALSE, col_names = FALSE, 
                                            +     col_types = cols(X2 = col_character()), 
                                            +     trim_ws = TRUE)

# Merge transcript-list from secretory and transmembrane proteins, and missense-variant list according to
# ENST-id column
missense_secr_transmemb <- merge(ALL_MISSENSE_VARIANTS_cleaned, Transcript_list, by.x = "X5", by.y = "Ensemble")

# Remove duplicates according to position
missense_secr_transmemb_no_duplicates <- missense_secr_transmemb[!duplicated(missense_secr_transmemb$X2), ]

# Creating a list with only positions of the variants
# First deleting extra columns
# Deleting extra columns
missense_secr_transmemb_positions <- subset(missense_secr_transmemb_no_duplicates, select = 
                                              -c(X5, X1, X3, X4, X6, X7, X8, X9, X10, X11, X12, X13, X14, Protein, Column))

# Separating chr and position into two separate columns
missense_secr_transmemb_positions <- cbind(missense_secr_transmemb_positions, 
                                                     read.table(text = as.character(missense_secr_transmemb_positions$X2), sep = ":"))

# Deleting two extra columns, leaving only positions of the variants
secr_transm_pos_FINAL <- subset(missense_secr_transmemb_positions, select = -c(X2, V1))

# Importing the bim-file of my dataset
KIDNEY_b38_no_relatives_FINAL <- read_delim("~/data/KIDNEY_b38_no_relatives_FINAL.bim", 
                                            +     "\t", escape_double = FALSE, col_names = FALSE, 
                                            +     trim_ws = TRUE)

# Merge variant table and bim-file according to position
List_of_secr_variants <- merge(KIDNEY_b38_no_relatives_FINAL, secr_transm_pos_FINAL, by.x = "X4", by.y = "V2")

# Reorder columns in created missense_variants dataframe
List_of_secr_variants <- List_of_secr_variants[, c(2, 3, 4, 1, 5, 6)]

# Remove duplicate variants according to position
List_of_secr_variants_FINAL <- List_of_secr_variants[!duplicated(List_of_secr_variants$X4), ]

# Write out the table
write.table(List_of_secr_variants_FINAL, file = "/data/List_of_secr_variants_FINAL", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
