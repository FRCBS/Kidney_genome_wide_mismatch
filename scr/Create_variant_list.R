###############################################################################
# 1
# Creating a variant list from bim-file
# chr 1-22
###############################################################################
# Importing the data
KIDNEY_b38_FINAL_DATA <- read_table2("~/data/KIDNEY_b38_FINAL_DATA.bim", 
                                     +     col_names = FALSE)

# Adding a column to data table
x <- add_column(KIDNEY_b38_FINAL_DATA, Y = '.', .before = "X5")

# Deleting two columns from data table
y <- subset(x, select = -c(X2, X3))

# Writing out the table
write.table(y, file = "/data/variants.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

###############################################################################
# 2
# chr X
###############################################################################
# Extract all variants in chr 23
library(dplyr)
z <- slice(y, 8837301:9085413)

# Adding a column to data table
testi <- add_column(z, C = 'X', .before = "X4")

# Deleting one column from data table
chrX <- subset(testi, select = -c(X1))
write.table(chrX, file = "/data/X_variants.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

###############################################################################
# 3
# Editing the merged variant list
###############################################################################
# Importing data table
ALL_MISSENSE_VARIANTS_cleaned <- read_table2("~/data/ALL_MISSENSE_VARIANTS_cleaned.txt", 
                                             +     col_names = FALSE, col_types = cols(X2 = col_character()))

# Deleting extra columns
missense <- subset(ALL_MISSENSE_VARIANTS_cleaned, select = -c(X1, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13, X14))

# Separating chr and position into two separate columns
missense_pos <- cbind(missense, read.table(text = as.character(missense$X2), sep = ":"))

# Deleting two extra columns, leaving only positions of the variants
missense_pos_FINAL <- subset(missense_pos, select = -c(X2, V1))
# Writing out the table
write.table(missense_pos_FINAL, file = "/data/missense_pos_FINAL.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)