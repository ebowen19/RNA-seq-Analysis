#Convert Gene IDs
# library(biomaRt)
# 
# ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")  # Mouse dataset
# gene_symbols <- c("Mybbp1a", "Niban1", "Alyref", "SREBF2", "GPx3", "Ifngr1", "Cltc", "Atp1a1")
# 
# # Convert gene symbols to Ensembl Gene IDs
# gene_ids <- getBM(attributes = c("ensembl_gene_id"), filters = "external_gene_name", values = gene_symbols, mart = ensembl)
# 
# print(gene_ids)


#for normalized together:
setwd("/Users/elizabeth 1/Desktop/Wu Lab/RNA Seq practice/data")
data <- read.csv("human&renca_normalizedTogether.tsv", header = TRUE, sep = "\t")

#rearrange data
# #change first column into row names
# rownames(data) <- data[, 1]
# # Remove the first column
# data <- data[, -1]
data <- data[, c("RC1","RC2","RVN1","RVN2","A","F","h1","h2","h12","X786O","WT","L169P")]
# View the modified dataframe
View(data)


# Extract a specific row by row name
Mybbp1a <- data["ENSMUSG00000040463", ] 
Niban1 <- data["ENSMUSG00000026483", ] #correct
Alyref <- data["ENSMUSG00000025134", ] #
SREBF2 <- data["ENSMUSG00000022463", ] #
GPx3 <- data["ENSMUSG00000018339", ]  #
Ifngr1 <- data["ENSMUSG00000020009", ] #
Cltc <- data["ENSMUSG00000047126", ] #correct
Atp1a1 <- data["ENSMUSG00000033161", ] #

#display - original 8 genes
Mybbp1a
Niban1
Alyref
SREBF2
GPx3
Ifngr1
Cltc
Atp1a1

#HMGCR pathways
ACAT1 <- data["ENSMUSG00000032047", ] 
ACAT2 <- data["ENSMUSG00000023832", ] 
HMGCS1 <- data["ENSMUSG00000093930", ] 
HMGCS2 <- data["ENSMUSG00000027875", ] 
HMGCR <- data["ENSMUSG00000021670", ] 
Pmvk <- data["ENSMUSG00000027952", ] #PKM?
MVD <- data["ENSMUSG00000006517", ] 
FDPS1 <- data["ENSMUSG00000059743", ] #Fdps?
GGPS1 <- data["ENSMUSG00000021302", ] 
#display
ACAT1 
ACAT2
HMGCS1
HMGCS2
HMGCR
Pmvk
MVD
FDPS1
GGPS1

#Cholesterol uptake / scavenge pathway
LDLR <- data["ENSMUSG00000032193", ] 
SREBF1 <- data["ENSMUSG00000020538", ] 
ACACA <- data["ENSMUSG00000020532", ] 
FASN <- data["ENSMUSG00000025153", ] 
SCD <- data["ENSMUSG00000037071", ] 
Spring1 <- data["ENSMUSG00000032840", ] 
#Display
LDLR
SREBF1
ACACA
FASN
SCD
Spring1

#Hypoxia:
VHL <- data["ENSMUSG00000033933", ] 
HIF1a <- data["ENSMUSG00000021109", ] #Hif1a
Hif1b <- data["ENSMUSG00000015522", ] #*ARNT
HIF2a <- data["ENSMUSG00000024140", ] #EPAS1*
Hif2b <- data["ENSG00000172379", ] #ARNT2*
VEGFa <- data["ENSMUSG00000023951", ] #Vegfa
Vegfb <- data["ENSMUSG00000024962", ]
Vegfc <- data["ENSMUSG00000031520", ] #*
Vegfd <- data["ENSMUSG00000031380", ] #*
EPo <- data["ENSMUSG00000029711", ] 
Glut1 <- data["ENSMUSG00000028645", ] #called Slc2a1
#Display
VHL
HIF1a
Hif1b
HIF2a
Hif2b
VEGFa
Vegfb
Vegfc #
Vegfd #
EPo
Glut1

#EMT
E_cadherin <- data["ENSMUSG00000000303", ] #cdh1
N_cadherin <- data["ENSMUSG00000024304", ] #cdh2
alpha_SMA <- data["ENSMUSG00000035783", ] #Acta2
SNAI1 <- data["ENSMUSG00000042821", ] 
SNAI2 <- data["ENSMUSG00000022676", ] 
ZEB1 <- data["ENSMUSG00000024238", ] 
ZEB2 <- data["ENSMUSG00000026872", ] 
MMP9 <- data["ENSMUSG00000017737", ] 
Ntn1 <- data["ENSMUSG00000020902", ] 
#Display
E_cadherin
N_cadherin
alpha_SMA
SNAI1
SNAI2
ZEB1
ZEB2
MMP9
Ntn1

#Other Interesting
MYC <- data["ENSMUSG00000022346", ] #c-MYC
MYCn <- data["ENSMUSG00000037169", ]
MYCl <- data["ENSG00000116990", ]
KRAS <- data["ENSMUSG00000030265", ] 
NRAS <- data["ENSMUSG00000027852", ] #?? Couldn't really find in ensembl website
Ki67 <- data["ENSMUSG00000022676", ] 
Caspase3 <- data["ENSMUSG00000031628", ] #casp3
UCHL1 <- data["ENSMUSG00000029223", ] 
#Display
MYC
MYCn
MYCl
KRAS
NRAS
Ki67
Caspase3
UCHL1

#Other 2
RAB_Agfg1 <- data["ENSMUSG00000026159", ] #Agfg1
RAB_Rab11a <- data["ENSMUSG00000004771", ]
RAP_Rptor <- data["ENSMUSG00000025583", ] 
RAP_Lrap1 <- data["ENSMUSG00000029103", ] 
Rac1 <- data["ENSMUSG00000001847", ] #Rac1
Rac2 <- data["ENSMUSG00000033220", ] 
Rho <- data["ENSMUSG00000030324", ]
Rock1 <- data["ENSMUSG00000024290", ]
Rock2 <- data["ENSMUSG00000020580", ]
#Display
RAB_Agfg1
RAB_Rab11a
RAP_Rptor
RAP_Lrap1
Rac1
Rac2
Rho
Rock1
Rock2

#Other 3
Txnrd1 <- data["ENSMUSG00000020250", ] 
Txnrd2 <- data["ENSMUSG00000075704", ] 
MDM2 <- data["ENSMUSG00000020184", ] 
P14_arf <- data["ENSMUSG00000044303", ] #CDKN2A
E2F <- data["ENSMUSG00000027490", ] 
p53 <- data["ENSMUSG00000059552", ] #Trp53
PPARa <- data["ENSMUSG00000022383", ] 
PPARg <- data["ENSMUSG00000000440", ] 
PPARd <- data["ENSMUSG00000002250", ] 
PRDX5 <- data["ENSMUSG00000024953", ] 
mTOR <- data["ENSMUSG00000028991", ] 
#Display
Txnrd1
Txnrd2
MDM2
P14_arf
E2F
p53
PPARa
PPARg
PPARd
PRDX5
mTOR

#Other #4 + new genes from Eliz GSEA
MDM4 <- data["ENSMUSG00000054387", ] 
usp7 <- data["ENSMUSG00000022710", ] 
arf1 <- data["ENSMUSG00000048076", ]
calr <- data["ENSMUSG00000003814", ]
trib3 <- data["ENSMUSG00000032715", ]
ddit3 <- data["ENSMUSG00000025408", ]
mapk1 <- data["ENSMUSG00000063358", ]
sqstm1 <- data["ENSMUSG00000015837", ]
p4ha1 <- data["ENSMUSG00000019916", ]
#display
usp7
MDM4
arf1
calr
trib3
ddit3
mapk1
sqstm1
p4ha1

# # Set working directory
# setwd("/Users/elizabeth 1/Desktop/Wu Lab/RNA Seq practice/data/gene plots/renca")
# 
# gene_list <- c("Mybbp1a", "Niban1", "Alyref", "SREBF2", "GPx3", "Ifngr1", "Cltc", "Atp1a1")
# 
# samples <- c("RC1", "RC2", "vhl-ko 1", "vhl-ko 2", "Atorv", "Fluv", "vhl+hif1-ko", "vhl+hif2-ko", "vhl+hif1&2-ko")
# 
# for (gene in gene_list) {
#   # Retrieve the dataframe for the current gene & rearrange columns
#   gene_df <- get(gene)
#   gene_df <- gene_df[, c("RC1", "RC2", "RVN1", "RVN2", "F", "A", "h1", "h2", "h12")]
#     
#   bar_data <- as.numeric(gene_df)
#   
#   # Start the graphics device for saving the image
#   png(paste0(gene, ".png"))
#   
#   # Create the bar plot
#   bp <- barplot(bar_data, main = gene, xlab = "Samples", ylab = "Values", col = "blue", beside = TRUE, names.arg = NA)
#   axis(1, at=bp, labels=samples, cex.axis=0.3, las=2, srt=45)
#   
#   # Close the graphics device
#   dev.off() 
# }
# 
# 
# for (gene in gene_list) {
#   # Retrieve the dataframe for the current gene & rearrange columns
#   gene_df <- get(gene)
#   gene_df <- gene_df[, c("RC1", "RC2", "RVN1", "RVN2", "F", "A", "h1", "h2", "h12")]
#   print(gene)
#   print(gene_df)
# }
