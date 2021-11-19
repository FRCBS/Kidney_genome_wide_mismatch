###############################################################################
# 210428 
# Creating files for both recipients and donors in which 
# the value 1 for each SNP stands for 'deletion' and value 0 for each SNP
# stands for 'no deletion'

# After creating the files, we can calculate the genomic collision between recipient
# and donor, and also calculate the sum of SNPs in both datasets 
###############################################################################

# CONTINUING WITH THE COVARIATE FILES FROM PREVIOUS survival analysis for all deletions

# Importing the covariate files 
R_dos_pheno_dels <- read_table2("~/Kidney_analyses/Results_new/Mm_and_deletion_analyses/R_covariates_deletions.txt")
D_dos_pheno_dels <- read_table2("~/Kidney_analyses/Results_new/Mm_and_deletion_analyses/D_covariates_deletions.txt")

# Recoding the SNPs accordingly:
# 0 = homozygous for major allele or heterozygous (no deletions)
# 1 = homozygous for minor allele (deletion)

###############################################################################
# RECIPIENTS FIRST

# First recoding the SNPs:
# chr7_101357735_A_G_G (rs6943474) and chr11_48548630_A_G_A (rs4882017) with different minor allele
# accordingly:
# 2 (homozygous for non-deletion allele) into value 1
# Then I have 1 = homozygous for non-deletion allele and heterozygous and
# 0 = homotzygous for deletion allele
R_dos_pheno_dels$rs6943474_0_1 <- recode(R_dos_pheno_dels$rs6943474, '2' = 1)
R_dos_pheno_dels$rs4882017_0_1 <- recode(R_dos_pheno_dels$rs4882017, '2' = 1)

# Then recoding the values 0 into value 2 (for deletion allele)
R_dos_pheno_dels$rs6943474_0_1 <- recode(R_dos_pheno_dels$rs6943474, '0' = 2)
R_dos_pheno_dels$rs4882017_0_1 <- recode(R_dos_pheno_dels$rs4882017, '0' = 2)

# Then again, recoding the values into 0 (no deletion) and 1 (deletion) according to:
# 1 -> 0 (no deletion)
R_dos_pheno_dels$rs6943474_0_1 <- recode(R_dos_pheno_dels$rs6943474_0_1, '1' = 0)
R_dos_pheno_dels$rs4882017_0_1 <- recode(R_dos_pheno_dels$rs4882017_0_1, '1' = 0)

# 2 -> 1 (deletion)
R_dos_pheno_dels$rs6943474_0_1 <- recode(R_dos_pheno_dels$rs6943474_0_1, '2' = 1)
R_dos_pheno_dels$rs4882017_0_1 <- recode(R_dos_pheno_dels$rs4882017_0_1, '2' = 1)

# Then the rest of the SNPs
# Recode 1 (heterozygous) into value 0
# (the other values remain the same (0=0 and 2=2))
# for other SNPs too
R_dos_pheno_dels$rs10927864_0_1 <- recode(R_dos_pheno_dels$rs10927864 , '1' = 0)
R_dos_pheno_dels$rs11209948_0_1 <- recode(R_dos_pheno_dels$rs11209948, '1' = 0)
R_dos_pheno_dels$rs11249248_0_1 <- recode(R_dos_pheno_dels$rs11249248, '1' = 0)
R_dos_pheno_dels$rs11587012_0_1 <- recode(R_dos_pheno_dels$rs11587012, '1' = 0)
R_dos_pheno_dels$rs158736_0_1 <- recode(R_dos_pheno_dels$rs158736, '1' = 0)
R_dos_pheno_dels$rs6693105_0_1 <- recode(R_dos_pheno_dels$rs6693105, '1' = 0)
R_dos_pheno_dels$rs7542235_0_1 <- recode(R_dos_pheno_dels$rs7542235, '1' = 0) # CFH
R_dos_pheno_dels$rs7419565_0_1 <- recode(R_dos_pheno_dels$rs7419565, '1' = 0) 
R_dos_pheno_dels$rs893403_0_1 <- recode(R_dos_pheno_dels$rs893403, '1' = 0) # LIMS1
R_dos_pheno_dels$rs10053292_0_1 <- recode(R_dos_pheno_dels$rs10053292, '1' = 0)
R_dos_pheno_dels$rs2387715_0_1 <- recode(R_dos_pheno_dels$rs2387715, '1' = 0)
R_dos_pheno_dels$rs7703761_0_1 <- recode(R_dos_pheno_dels$rs7703761, '1' = 0)
R_dos_pheno_dels$rs17654108_0_1 <- recode(R_dos_pheno_dels$rs17654108, '1' = 0)
R_dos_pheno_dels$rs2160195_0_1 <- recode(R_dos_pheno_dels$rs2160195, '1' = 0)
R_dos_pheno_dels$rs4621754_0_1 <- recode(R_dos_pheno_dels$rs4621754, '1' = 0)
R_dos_pheno_dels$rs4729606_0_1 <- recode(R_dos_pheno_dels$rs4729606, '1' = 0)
R_dos_pheno_dels$rs11985201_0_1 <- recode(R_dos_pheno_dels$rs11985201, '1' = 0)
R_dos_pheno_dels$rs4543566_0_1 <- recode(R_dos_pheno_dels$rs4543566, '1' = 0)
R_dos_pheno_dels$rs1523688_0_1 <- recode(R_dos_pheno_dels$rs1523688, '1' = 0)
R_dos_pheno_dels$rs2174926_0_1 <- recode(R_dos_pheno_dels$rs2174926, '1' = 0)
R_dos_pheno_dels$rs10885336_0_1 <- recode(R_dos_pheno_dels$rs10885336, '1' = 0)
R_dos_pheno_dels$rs2342606_0_1 <- recode(R_dos_pheno_dels$rs2342606, '1' = 0)
R_dos_pheno_dels$rs3793917_0_1 <- recode(R_dos_pheno_dels$rs3793917, '1' = 0)
R_dos_pheno_dels$rs11228868_0_1 <- recode(R_dos_pheno_dels$rs11228868, '1' = 0)
R_dos_pheno_dels$rs1944862_0_1 <- recode(R_dos_pheno_dels$rs1944862, '1' = 0)
R_dos_pheno_dels$rs1478309_0_1 <- recode(R_dos_pheno_dels$rs1478309, '1' = 0)
R_dos_pheno_dels$rs9318648_0_1 <- recode(R_dos_pheno_dels$rs9318648, '1' = 0)
R_dos_pheno_dels$rs11156875_0_1 <- recode(R_dos_pheno_dels$rs11156875, '1' = 0)
R_dos_pheno_dels$rs8007442_0_1 <- recode(R_dos_pheno_dels$rs8007442, '1' = 0)
R_dos_pheno_dels$rs8022070_0_1 <- recode(R_dos_pheno_dels$rs8022070, '1' = 0)
R_dos_pheno_dels$rs10521145_0_1 <- recode(R_dos_pheno_dels$rs10521145, '1' = 0)
R_dos_pheno_dels$rs2244613_0_1 <- recode(R_dos_pheno_dels$rs2244613, '1' = 0)
R_dos_pheno_dels$rs16966699_0_1 <- recode(R_dos_pheno_dels$rs16966699, '1' = 0)
R_dos_pheno_dels$rs8064493_0_1 <- recode(R_dos_pheno_dels$rs8064493, '1' = 0)
R_dos_pheno_dels$rs103294_0_1 <- recode(R_dos_pheno_dels$rs103294, '1' = 0)
R_dos_pheno_dels$rs324121_0_1 <- recode(R_dos_pheno_dels$rs324121, '1' = 0)
R_dos_pheno_dels$rs3810336_0_1 <- recode(R_dos_pheno_dels$rs3810336, '1' = 0)
R_dos_pheno_dels$rs4806152_0_1 <- recode(R_dos_pheno_dels$rs4806152, '1' = 0)

# Recode 2 (heterozygous) into value 1,
# after which the values are 0 = no deletion, 1 = deletion
R_dos_pheno_dels$rs10927864_0_1 <- recode(R_dos_pheno_dels$rs10927864_0_1 , '2' = 1)
R_dos_pheno_dels$rs11209948_0_1 <- recode(R_dos_pheno_dels$rs11209948_0_1, '2' = 1)
R_dos_pheno_dels$rs11249248_0_1 <- recode(R_dos_pheno_dels$rs11249248_0_1, '2' = 1)
R_dos_pheno_dels$rs11587012_0_1 <- recode(R_dos_pheno_dels$rs11587012_0_1, '2' = 1)
R_dos_pheno_dels$rs158736_0_1 <- recode(R_dos_pheno_dels$rs158736_0_1, '2' = 1)
R_dos_pheno_dels$rs6693105_0_1 <- recode(R_dos_pheno_dels$rs6693105_0_1, '2' = 1)
R_dos_pheno_dels$rs7542235_0_1 <- recode(R_dos_pheno_dels$rs7542235_0_1, '2' = 1) # CFH
R_dos_pheno_dels$rs7419565_0_1 <- recode(R_dos_pheno_dels$rs7419565_0_1, '2' = 1)
R_dos_pheno_dels$rs893403_0_1 <- recode(R_dos_pheno_dels$rs893403_0_1, '2' = 1) # LIMS1
R_dos_pheno_dels$rs10053292_0_1 <- recode(R_dos_pheno_dels$rs10053292_0_1, '2' = 1)
R_dos_pheno_dels$rs2387715_0_1 <- recode(R_dos_pheno_dels$rs2387715_0_1, '2' = 1)
R_dos_pheno_dels$rs7703761_0_1 <- recode(R_dos_pheno_dels$rs7703761_0_1, '2' = 1)
R_dos_pheno_dels$rs17654108_0_1 <- recode(R_dos_pheno_dels$rs17654108_0_1, '2' = 1)
R_dos_pheno_dels$rs2160195_0_1 <- recode(R_dos_pheno_dels$rs2160195_0_1, '2' = 1)
R_dos_pheno_dels$rs4621754_0_1 <- recode(R_dos_pheno_dels$rs4621754_0_1, '2' = 1)
R_dos_pheno_dels$rs4729606_0_1 <- recode(R_dos_pheno_dels$rs4729606_0_1, '2' = 1)
R_dos_pheno_dels$rs11985201_0_1 <- recode(R_dos_pheno_dels$rs11985201_0_1, '2' = 1)
R_dos_pheno_dels$rs4543566_0_1 <- recode(R_dos_pheno_dels$rs4543566_0_1, '2' = 1)
R_dos_pheno_dels$rs1523688_0_1 <- recode(R_dos_pheno_dels$rs1523688_0_1, '2' = 1)
R_dos_pheno_dels$rs2174926_0_1 <- recode(R_dos_pheno_dels$rs2174926_0_1, '2' = 1)
R_dos_pheno_dels$rs10885336_0_1 <- recode(R_dos_pheno_dels$rs10885336_0_1, '2' = 1)
R_dos_pheno_dels$rs2342606_0_1 <- recode(R_dos_pheno_dels$rs2342606_0_1, '2' = 1)
R_dos_pheno_dels$rs3793917_0_1 <- recode(R_dos_pheno_dels$rs3793917_0_1, '2' = 1)
R_dos_pheno_dels$rs11228868_0_1 <- recode(R_dos_pheno_dels$rs11228868_0_1, '2' = 1)
R_dos_pheno_dels$rs1944862_0_1 <- recode(R_dos_pheno_dels$rs1944862_0_1, '2' = 1)
R_dos_pheno_dels$rs1478309_0_1 <- recode(R_dos_pheno_dels$rs1478309_0_1, '2' = 1)
R_dos_pheno_dels$rs9318648_0_1 <- recode(R_dos_pheno_dels$rs9318648_0_1, '2' = 1)
R_dos_pheno_dels$rs11156875_0_1 <- recode(R_dos_pheno_dels$rs11156875_0_1, '2' = 1)
R_dos_pheno_dels$rs8007442_0_1 <- recode(R_dos_pheno_dels$rs8007442_0_1, '2' = 1)
R_dos_pheno_dels$rs8022070_0_1 <- recode(R_dos_pheno_dels$rs8022070_0_1, '2' = 1)
R_dos_pheno_dels$rs10521145_0_1 <- recode(R_dos_pheno_dels$rs10521145_0_1, '2' = 1)
R_dos_pheno_dels$rs2244613_0_1 <- recode(R_dos_pheno_dels$rs2244613_0_1, '2' = 1)
R_dos_pheno_dels$rs16966699_0_1 <- recode(R_dos_pheno_dels$rs16966699_0_1, '2' = 1)
R_dos_pheno_dels$rs8064493_0_1 <- recode(R_dos_pheno_dels$rs8064493_0_1, '2' = 1)
R_dos_pheno_dels$rs103294_0_1 <- recode(R_dos_pheno_dels$rs103294_0_1, '2' = 1)
R_dos_pheno_dels$rs324121_0_1 <- recode(R_dos_pheno_dels$rs324121_0_1, '2' = 1)
R_dos_pheno_dels$rs3810336_0_1 <- recode(R_dos_pheno_dels$rs3810336_0_1, '2' = 1)
R_dos_pheno_dels$rs4806152_0_1 <- recode(R_dos_pheno_dels$rs4806152_0_1, '2' = 1)

# Write out the table
write.table(R_dos_pheno_dels, file = "/home/markkinens/Kidney_analyses/Results_new/Mm_and_deletion_analyses/R_covariates_deletions_0_1.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

###############################################################################
# DONORS NEXT

# First recoding the SNPs:
# chr7_101357735_A_G_G (rs6943474) and chr11_48548630_A_G_A (rs4882017) with different minor allele
# accordingly:
# 2 (homozygous for non-deletion allele) into value 1
# Then I have 1 = homozygous for non-deletion allele and heterozygous and
# 0 = homotzygous for deletion allele
D_dos_pheno_dels$rs6943474_0_1 <- recode(D_dos_pheno_dels$rs6943474, '2' = 1)
D_dos_pheno_dels$rs4882017_0_1 <- recode(D_dos_pheno_dels$rs4882017, '2' = 1)

# Then recoding the values 0 into value 2 (for deletion allele)
D_dos_pheno_dels$rs6943474_0_1 <- recode(D_dos_pheno_dels$rs6943474, '0' = 2)
D_dos_pheno_dels$rs4882017_0_1 <- recode(D_dos_pheno_dels$rs4882017, '0' = 2)

# Then again, recoding the values into 0 (no deletion) and 1 (deletion) according to:
# 1 -> 0 (no deletion)
D_dos_pheno_dels$rs6943474_0_1 <- recode(D_dos_pheno_dels$rs6943474_0_1, '1' = 0)
D_dos_pheno_dels$rs4882017_0_1 <- recode(D_dos_pheno_dels$rs4882017_0_1, '1' = 0)

# 2 -> 1 (deletion)
D_dos_pheno_dels$rs6943474_0_1 <- recode(D_dos_pheno_dels$rs6943474_0_1, '2' = 1)
D_dos_pheno_dels$rs4882017_0_1 <- recode(D_dos_pheno_dels$rs4882017_0_1, '2' = 1)

# Then the rest of the SNPs
# Recode 1 (heterozygous) into value 0
# (the other values remain the same (0=0 and 2=2))
# for other SNPs too
D_dos_pheno_dels$rs10927864_0_1 <- recode(D_dos_pheno_dels$rs10927864 , '1' = 0)
D_dos_pheno_dels$rs11209948_0_1 <- recode(D_dos_pheno_dels$rs11209948, '1' = 0)
D_dos_pheno_dels$rs11249248_0_1 <- recode(D_dos_pheno_dels$rs11249248, '1' = 0)
D_dos_pheno_dels$rs11587012_0_1 <- recode(D_dos_pheno_dels$rs11587012, '1' = 0)
D_dos_pheno_dels$rs158736_0_1 <- recode(D_dos_pheno_dels$rs158736, '1' = 0)
D_dos_pheno_dels$rs6693105_0_1 <- recode(D_dos_pheno_dels$rs6693105, '1' = 0)
D_dos_pheno_dels$rs7542235_0_1 <- recode(D_dos_pheno_dels$rs7542235, '1' = 0) # CFH
D_dos_pheno_dels$rs7419565_0_1 <- recode(D_dos_pheno_dels$rs7419565, '1' = 0) 
D_dos_pheno_dels$rs893403_0_1 <- recode(D_dos_pheno_dels$rs893403, '1' = 0) # LIMS1
D_dos_pheno_dels$rs10053292_0_1 <- recode(D_dos_pheno_dels$rs10053292, '1' = 0)
D_dos_pheno_dels$rs2387715_0_1 <- recode(D_dos_pheno_dels$rs2387715, '1' = 0)
D_dos_pheno_dels$rs7703761_0_1 <- recode(D_dos_pheno_dels$rs7703761, '1' = 0)
D_dos_pheno_dels$rs17654108_0_1 <- recode(D_dos_pheno_dels$rs17654108, '1' = 0)
D_dos_pheno_dels$rs2160195_0_1 <- recode(D_dos_pheno_dels$rs2160195, '1' = 0)
D_dos_pheno_dels$rs4621754_0_1 <- recode(D_dos_pheno_dels$rs4621754, '1' = 0)
D_dos_pheno_dels$rs4729606_0_1 <- recode(D_dos_pheno_dels$rs4729606, '1' = 0)
D_dos_pheno_dels$rs11985201_0_1 <- recode(D_dos_pheno_dels$rs11985201, '1' = 0)
D_dos_pheno_dels$rs4543566_0_1 <- recode(D_dos_pheno_dels$rs4543566, '1' = 0)
D_dos_pheno_dels$rs1523688_0_1 <- recode(D_dos_pheno_dels$rs1523688, '1' = 0)
D_dos_pheno_dels$rs2174926_0_1 <- recode(D_dos_pheno_dels$rs2174926, '1' = 0)
D_dos_pheno_dels$rs10885336_0_1 <- recode(D_dos_pheno_dels$rs10885336, '1' = 0)
D_dos_pheno_dels$rs2342606_0_1 <- recode(D_dos_pheno_dels$rs2342606, '1' = 0)
D_dos_pheno_dels$rs3793917_0_1 <- recode(D_dos_pheno_dels$rs3793917, '1' = 0)
D_dos_pheno_dels$rs11228868_0_1 <- recode(D_dos_pheno_dels$rs11228868, '1' = 0)
D_dos_pheno_dels$rs1944862_0_1 <- recode(D_dos_pheno_dels$rs1944862, '1' = 0)
D_dos_pheno_dels$rs1478309_0_1 <- recode(D_dos_pheno_dels$rs1478309, '1' = 0)
D_dos_pheno_dels$rs9318648_0_1 <- recode(D_dos_pheno_dels$rs9318648, '1' = 0)
D_dos_pheno_dels$rs11156875_0_1 <- recode(D_dos_pheno_dels$rs11156875, '1' = 0)
D_dos_pheno_dels$rs8007442_0_1 <- recode(D_dos_pheno_dels$rs8007442, '1' = 0)
D_dos_pheno_dels$rs8022070_0_1 <- recode(D_dos_pheno_dels$rs8022070, '1' = 0)
D_dos_pheno_dels$rs10521145_0_1 <- recode(D_dos_pheno_dels$rs10521145, '1' = 0)
D_dos_pheno_dels$rs2244613_0_1 <- recode(D_dos_pheno_dels$rs2244613, '1' = 0)
D_dos_pheno_dels$rs16966699_0_1 <- recode(D_dos_pheno_dels$rs16966699, '1' = 0)
D_dos_pheno_dels$rs8064493_0_1 <- recode(D_dos_pheno_dels$rs8064493, '1' = 0)
D_dos_pheno_dels$rs103294_0_1 <- recode(D_dos_pheno_dels$rs103294, '1' = 0)
D_dos_pheno_dels$rs324121_0_1 <- recode(D_dos_pheno_dels$rs324121, '1' = 0)
D_dos_pheno_dels$rs3810336_0_1 <- recode(D_dos_pheno_dels$rs3810336, '1' = 0)
D_dos_pheno_dels$rs4806152_0_1 <- recode(D_dos_pheno_dels$rs4806152, '1' = 0)

# Recode 2 (heterozygous) into value 1,
# after which the values are 0 = no deletion, 1 = deletion
D_dos_pheno_dels$rs10927864_0_1 <- recode(D_dos_pheno_dels$rs10927864_0_1 , '2' = 1)
D_dos_pheno_dels$rs11209948_0_1 <- recode(D_dos_pheno_dels$rs11209948_0_1, '2' = 1)
D_dos_pheno_dels$rs11249248_0_1 <- recode(D_dos_pheno_dels$rs11249248_0_1, '2' = 1)
D_dos_pheno_dels$rs11587012_0_1 <- recode(D_dos_pheno_dels$rs11587012_0_1, '2' = 1)
D_dos_pheno_dels$rs158736_0_1 <- recode(D_dos_pheno_dels$rs158736_0_1, '2' = 1)
D_dos_pheno_dels$rs6693105_0_1 <- recode(D_dos_pheno_dels$rs6693105_0_1, '2' = 1)
D_dos_pheno_dels$rs7542235_0_1 <- recode(D_dos_pheno_dels$rs7542235_0_1, '2' = 1) # CFH
D_dos_pheno_dels$rs7419565_0_1 <- recode(D_dos_pheno_dels$rs7419565_0_1, '2' = 1)
D_dos_pheno_dels$rs893403_0_1 <- recode(D_dos_pheno_dels$rs893403_0_1, '2' = 1) # LIMS1
D_dos_pheno_dels$rs10053292_0_1 <- recode(D_dos_pheno_dels$rs10053292_0_1, '2' = 1)
D_dos_pheno_dels$rs2387715_0_1 <- recode(D_dos_pheno_dels$rs2387715_0_1, '2' = 1)
D_dos_pheno_dels$rs7703761_0_1 <- recode(D_dos_pheno_dels$rs7703761_0_1, '2' = 1)
D_dos_pheno_dels$rs17654108_0_1 <- recode(D_dos_pheno_dels$rs17654108_0_1, '2' = 1)
D_dos_pheno_dels$rs2160195_0_1 <- recode(D_dos_pheno_dels$rs2160195_0_1, '2' = 1)
D_dos_pheno_dels$rs4621754_0_1 <- recode(D_dos_pheno_dels$rs4621754_0_1, '2' = 1)
D_dos_pheno_dels$rs4729606_0_1 <- recode(D_dos_pheno_dels$rs4729606_0_1, '2' = 1)
D_dos_pheno_dels$rs11985201_0_1 <- recode(D_dos_pheno_dels$rs11985201_0_1, '2' = 1)
D_dos_pheno_dels$rs4543566_0_1 <- recode(D_dos_pheno_dels$rs4543566_0_1, '2' = 1)
D_dos_pheno_dels$rs1523688_0_1 <- recode(D_dos_pheno_dels$rs1523688_0_1, '2' = 1)
D_dos_pheno_dels$rs2174926_0_1 <- recode(D_dos_pheno_dels$rs2174926_0_1, '2' = 1)
D_dos_pheno_dels$rs10885336_0_1 <- recode(D_dos_pheno_dels$rs10885336_0_1, '2' = 1)
D_dos_pheno_dels$rs2342606_0_1 <- recode(D_dos_pheno_dels$rs2342606_0_1, '2' = 1)
D_dos_pheno_dels$rs3793917_0_1 <- recode(D_dos_pheno_dels$rs3793917_0_1, '2' = 1)
D_dos_pheno_dels$rs11228868_0_1 <- recode(D_dos_pheno_dels$rs11228868_0_1, '2' = 1)
D_dos_pheno_dels$rs1944862_0_1 <- recode(D_dos_pheno_dels$rs1944862_0_1, '2' = 1)
D_dos_pheno_dels$rs1478309_0_1 <- recode(D_dos_pheno_dels$rs1478309_0_1, '2' = 1)
D_dos_pheno_dels$rs9318648_0_1 <- recode(D_dos_pheno_dels$rs9318648_0_1, '2' = 1)
D_dos_pheno_dels$rs11156875_0_1 <- recode(D_dos_pheno_dels$rs11156875_0_1, '2' = 1)
D_dos_pheno_dels$rs8007442_0_1 <- recode(D_dos_pheno_dels$rs8007442_0_1, '2' = 1)
D_dos_pheno_dels$rs8022070_0_1 <- recode(D_dos_pheno_dels$rs8022070_0_1, '2' = 1)
D_dos_pheno_dels$rs10521145_0_1 <- recode(D_dos_pheno_dels$rs10521145_0_1, '2' = 1)
D_dos_pheno_dels$rs2244613_0_1 <- recode(D_dos_pheno_dels$rs2244613_0_1, '2' = 1)
D_dos_pheno_dels$rs16966699_0_1 <- recode(D_dos_pheno_dels$rs16966699_0_1, '2' = 1)
D_dos_pheno_dels$rs8064493_0_1 <- recode(D_dos_pheno_dels$rs8064493_0_1, '2' = 1)
D_dos_pheno_dels$rs103294_0_1 <- recode(D_dos_pheno_dels$rs103294_0_1, '2' = 1)
D_dos_pheno_dels$rs324121_0_1 <- recode(D_dos_pheno_dels$rs324121_0_1, '2' = 1)
D_dos_pheno_dels$rs3810336_0_1 <- recode(D_dos_pheno_dels$rs3810336_0_1, '2' = 1)
D_dos_pheno_dels$rs4806152_0_1 <- recode(D_dos_pheno_dels$rs4806152_0_1, '2' = 1)

# Write out the table
write.table(D_dos_pheno_dels, file = "/home/markkinens/Kidney_analyses/Results_new/Mm_and_deletion_analyses/D_covariates_deletions_0_1.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

