###############################################################################
# Creating a transmembrane protein list
###############################################################################

library(tidyr)
library(dplyr)

# Import datatable
transmembrane <- read_delim("~/data/uniprot-(annotation (type transmem)+OR+keyword Transmembrane+[K--.tab", 
                            +     "\t", escape_double = FALSE, col_names = FALSE, 
                            +     trim_ws = TRUE)

# Remove extra columns
transmembrane_uniprot <- select(transmembrane, -X1, -X3, -X4, -X5, -X6, -X7)

# Write out the table
write.table(transmembrane_uniprot, file = "/data/uniprot_transmembrane.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# First remove the 'header line' from the uniprot_transmembrane.txt file,
# then paste it to excel and edit the file by separating each ENST-id to separate columns
# The list of separated ENST-ids is in /data/Transmembrane_separated.xlsx

# Importing the separated list
Transmembrane_separated <- read_excel("~/Kidney_analyses/Transmembrane_separated.xlsx", 
                                      +     col_names = FALSE)

# Gathering ensemble-id columns into one column
# As a result: first column is gene name, second is the column number from the previous data 
# and third is ENST-id
Transmembrane_list <- gather(Transmembrane_separated, "Column", "...2", -...1, na.rm = TRUE)

# Rename two columns
Transmembrane_list <- Transmembrane_list %>% rename(Protein = ...1, Ensemble = ...2)

# Write out the table
write.table(Transmembrane_list, file = "/data/Transmembrane_list", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

###############################################################################
# Creating a variant list from only transmembrane proteins

# Import the missense variant list you created with vep
ALL_MISSENSE_VARIANTS_cleaned <- read_delim("~/test_vep/ALL_MISSENSE_VARIANTS_cleaned.txt", 
                                            +     "\t", escape_double = FALSE, col_names = FALSE, 
                                            +     col_types = cols(X2 = col_character()), 
                                            +     trim_ws = TRUE)

# Merge transcript-list from transmembrane proteins, and missense-variant list according to
# ENST-id column
missense_transmemb <- merge(ALL_MISSENSE_VARIANTS_cleaned, Transmembrane_list, by.x = "X5", by.y = "Ensemble")

# Remove duplicates according to position
missense_transmemb_no_duplicates <- missense_transmemb[!duplicated(missense_transmemb$X2), ]

# Creating a list with only positions of the variants
# First deleting extra columns

# Deleting extra columns
missense_transmemb_positions <- subset(missense_transmemb_no_duplicates, select = 
                                              -c(X5, X1, X3, X4, X6, X7, X8, X9, X10, X11, X12, X13, X14, Protein, Column))

# Separating chr and position into two separate columns
missense_transmemb_positions <- cbind(missense_transmemb_positions, 
                                           read.table(text = as.character(missense_transmemb_positions$X2), sep = ":"))

# Deleting two extra columns, leaving only positions of the variants
transm_pos_FINAL <- subset(missense_transmemb_positions, select = -c(X2, V1))

# Import bim-file of my dataset
KIDNEY_b38_no_relatives_FINAL <- read_delim("~/data/KIDNEY_b38_no_relatives_FINAL.bim", 
                                            +     "\t", escape_double = FALSE, col_names = FALSE, 
                                            +     trim_ws = TRUE)

# Merge variant table and bim-file according to position
List_of_transm_variants <- merge(KIDNEY_b38_no_relatives_FINAL, transm_pos_FINAL, by.x = "X4", by.y = "V2")

# Reorder columns in created missense_variants dataframe
List_of_transm_variants <- List_of_transm_variants[, c(2, 3, 4, 1, 5, 6)]

# Remove duplicate variants according to position
List_of_transm_variants_FINAL <- List_of_transm_variants[!duplicated(List_of_transm_variants$X4), ]

# Write out the table
write.table(List_of_transm_variants_FINAL, file = "/data/List_of_transm_variants_FINAL", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
