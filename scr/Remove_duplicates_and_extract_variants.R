###############################################################################
# Removing duplicates from SNP missense list and extracting variants from 
# bim-file
###############################################################################
# Import missense-variant data
missense_pos_FINAL <- read_csv("~/data/missense_pos_FINAL.txt", 
                               +     col_names = FALSE)

# Remove duplicate variants
missense_no_duplicates <- missense_pos_FINAL[!duplicated(missense_pos_FINAL$X1), ]

# Write out the table
write.table(missense_no_duplicates, file = "/data/missense_no_duplicates.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Import bim-file 
KIDNEY_b38_no_relatives_FINAL <- read_table2("~/data/KIDNEY_b38_no_relatives_FINAL.bim", 
                                            +     col_names = FALSE)

# Merge variant table and bim-file
missense_variants <- merge(KIDNEY_b38_no_relatives_FINAL, missense_no_duplicates, by.x = "X4", by.y = "X1")

# Reorder columns in created missense_variants dataframe
missense_variants_reordered <- missense_variants[, c(2, 3, 4, 1, 5, 6)]

# Write out the table
write.table(missense_variants_reordered, file = "/data/missense_bim.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
