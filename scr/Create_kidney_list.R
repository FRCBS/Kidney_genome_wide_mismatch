###############################################################################
# Creating a list of kidney-expressed proteins
###############################################################################

library(tidyr)
library(dplyr)

# Import datatable
kidney <- read_delim("~/data/uniprot-annotation (type tissue+specificity +Kidney)+reviewed yes--.tab", 
                     +     "\t", escape_double = FALSE, col_names = FALSE, 
                     +     trim_ws = TRUE)

# Remove extra columns
kidney <- select(kidney, -X1, -X3, -X4, -X5, -X6, -X7)

# Write out the table
write.table(kidney, file = "/data/uniprot_kidney.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# First remove the 'header line' from the uniprot_kidney.txt file,
# then paste it to excel and edit the file by separating each ENST-id to separate columns
# The list of separated ENST-ids is in /data/Kidney_separated.xlsx

# Importing the separated list
Kidney_separated <- read_excel("~/data/Kidney_separated.xlsx", 
                                      +     col_names = FALSE)

# Gathering ensemble-id columns into one column
# As a result: first column is gene name, second is the column number from the previous data 
# and third is ENST-id
Kidney_list <- gather(Kidney_separated, "Column", "...2", -...1, na.rm = TRUE)

# Rename two columns
Kidney_list <- Kidney_list %>% rename(Protein = ...1, Ensemble = ...2)

# Write out the table
write.table(Kidney_list, file = "/data/Kidney_list", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

###############################################################################
# Creating a variant list from only kidney-expressed proteins

# Import the missense variant list you got/created from vep
ALL_MISSENSE_VARIANTS_cleaned <- read_delim("~/data/ALL_MISSENSE_VARIANTS_cleaned.txt", 
                                            +     "\t", escape_double = FALSE, col_names = FALSE, 
                                            +     col_types = cols(X2 = col_character()), 
                                            +     trim_ws = TRUE)

# Merge transcript-list from transmembrane proteins, and missense-variant list according to
# ENST-id column
missense_kidney <- merge(ALL_MISSENSE_VARIANTS_cleaned, Kidney_list, by.x = "X5", by.y = "Ensemble")

# Remove duplicates according to position
missense_kidney_no_duplicates <- missense_kidney[!duplicated(missense_kidney$X2), ]

# Creating a list with only positions of the variants
# First deleting extra columns

# Deleting extra columns
missense_kidney_positions <- subset(missense_kidney_no_duplicates, select = 
                                         -c(X5, X1, X3, X4, X6, X7, X8, X9, X10, X11, X12, X13, X14, Protein, Column))

# Separating chr and position into two separate columns
missense_kidney_positions <- cbind(missense_kidney_positions, 
                                      read.table(text = as.character(missense_kidney_positions$X2), sep = ":"))

# Deleting two extra columns, leaving only positions of the variants
kidney_pos_FINAL <- subset(missense_kidney_positions, select = -c(X2, V1))

# Import bim-file of my dataset
KIDNEY_b38_no_relatives_FINAL <- read_delim("~/data/KIDNEY_b38_no_relatives_FINAL.bim", 
                                            +     "\t", escape_double = FALSE, col_names = FALSE, 
                                            +     trim_ws = TRUE)

# Merge variant table and bim-file according to position
List_of_kidney_variants <- merge(KIDNEY_b38_no_relatives_FINAL, kidney_pos_FINAL, by.x = "X4", by.y = "V2")

# Reorder columns in created missense_variants dataframe
List_of_kidney_variants <- List_of_kidney_variants[, c(2, 3, 4, 1, 5, 6)]

# Remove duplicate variants according to position
List_of_kidney_variants_FINAL <- List_of_kidney_variants[!duplicated(List_of_kidney_variants$X4), ]

# Write out the table
write.table(List_of_kidney_variants_FINAL, file = "/data/List_of_kidney_variants_FINAL", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
