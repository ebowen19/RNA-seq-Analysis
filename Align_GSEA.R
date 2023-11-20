#load in data
setwd("/Users/elizabeth 1/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/gsea")
reactome_sensiDown <- read.csv("Reactome GSEA sensi down.csv")
reactome_sensiDown <- reactome_sensiDown[c(1,2,3,4)] #remove NA columns
reactome_sensiUp <- read.csv("Reactome GSEA sensi up.csv")
reactome_sensiUp <- reactome_sensiUp[c(1,2,3,4)]

df <- reactome_sensiUp #****REPLACE with whichever dataframe you want to analyze & replace csv name!
name <- "reactomeSensiUp"
# Create an empty dataframe with the same number of rows as df and only NA values
df_aligned <- data.frame(matrix(NA, nrow = nrow(df), ncol = ncol(df)))
names(df_aligned) <- names(df) # Assign the column names of the original dataframe to the new one

# Fill in the reference column
df_aligned$RVN.up.vs.RC <- df$RVN.up.vs.RC

#remove empty rows
df_aligned <- df_aligned[df_aligned$RVN.up.vs.RC != "", ]
# Get the pathways from the first column
reference_pathways <- df_aligned[[1]] #returns the column as a vector rather than a dataframe
reference_pathways <- reference_pathways[reference_pathways != ""] #take out the empty entries

# Loop through each of the other columns
for (i in 2:ncol(df)) {
  # Filter and order based on reference column
  filtered_values <- df[[i]][df[[i]] %in% reference_pathways] 
  order <- match(filtered_values, reference_pathways)
  sorted_indices <- order(order) #order() returns indices of vector's elements in ascending order. 
  # Reorder filtered_values according to sorted_indices
  ordered_values <- filtered_values[sorted_indices]
  
  # Get the match indices--> indices of ordered_values (vals of subsequent column) that should appear in each entry of that column
  match_indices <- match(df_aligned$RVN.up.vs.RC, ordered_values)

  # Use the match indices to place the ordered values into the correct positions where match exists
  df_aligned[!is.na(match_indices), i] <- ordered_values[match_indices[!is.na(match_indices)]]
}

write.csv(df_aligned, paste(name, "_Aligned.csv"), row.names = FALSE) #write csv of data aligned to rvn vs rfl column
common_indices <- complete.cases(df_aligned) # boolean vector indicating rows with all complete cases
common_pathways <- df_aligned[common_indices,] # subset the data frame to only complete cases
common_pathways <- common_pathways[-c(2,3,4)]
names(common_pathways) <- paste(name, "Common Pathways")
write.csv(common_pathways, paste(name, "_Pathways.csv"), row.names = TRUE)

#with explanations...---------------------------------------------------------------------------------------------------------------
# #need dataframe to be columns of same length -- result of  lapply() 
# #call would be a list of vectors of different lengths. 
# #When you try to assign these back to the data frame with df[-1] <- ..., 
# #get an error because the replacement vectors do not all have the length that 
# #R expects for the number of rows in df.
# df <- reactome_sensiDown
# # Create an empty dataframe with the same number of rows as df and only NA values
# df_aligned <- data.frame(matrix(NA, nrow = nrow(df), ncol = ncol(df)))
# names(df_aligned) <- names(df) # Assign the column names of the original dataframe to the new one
# 
# # Fill in the reference column
# df_aligned$RVN.down.vs.RFL <- df$RVN.down.vs.RFL
# 
# #remove empty rows
# df_aligned <- df_aligned[df_aligned$RVN.down.vs.RFL != "", ]
# #subsetting dataframe syntax: [rows, columns]. condition/set of indexes before the comma
# #R interprets as row-wise operation.
# 
# # Get the pathways from the first column
# reference_pathways <- df_aligned[[1]] #returns the column as a vector rather than a dataframe
# reference_pathways <- reference_pathways[reference_pathways != ""] #take out the empty entries
# 
# # Loop through each of the other columns
# for (i in 2:ncol(df)) {
#   # Filter and order based on reference column
#   filtered_values <- df[[i]][df[[i]] %in% reference_pathways] #The outer df[[i]][] uses the logical
#   #vector to subset the original df[[i]] vector again. It selects only the elements where the 
#   #logical vector is TRUE (i.e., only those elements of the i-th column that are 
#   #also found in reference).
#   order <- match(filtered_values, reference_pathways)
#   sorted_indices <- order(order) #order() returns indices of vector's elements in ascending order. 
#   # tells positions to rearrange the elements to sort the vector. - the elements will al be in range of filtered_values
#   #bc it's coming from the # of elements (from the indices) of order, not the values themselves of order
#   # Reorder filtered_values according to sorted_indices
#   ordered_values <- filtered_values[sorted_indices]
#   
#   # Assuming 'i' is the index for the current column being processed
#   # Get the match indices--> indices of ordered_values (vals of subsequent column) that should appear in each entry of that column
#   # where column has same length as 1st column (we're corresponding/matching the vals to the vals in 1st column)
#   match_indices <- match(df_aligned$RVN.down.vs.RFL, ordered_values)
#   #match(a, b) function in R goes in the order of the elements in a. For each element in a, 
#   #match looks for its first occurrence in b and returns the position of that match in b. 
#   #If there is no match found, NA is returned for that element.
#   
#   # Use the match indices to place the ordered values into the correct positions, 
#   # but only where match_indices is not NA (i.e., there is a match)
#   df_aligned[!is.na(match_indices), i] <- ordered_values[match_indices[!is.na(match_indices)]]
#   
#   #match gives positions in the df_aligned$V1 column where the values of ordered_values are found.
#   #assigning ordered_values to the i-th column of df_aligned at the positions that 
#   #align with the df_aligned$V1 column.
# }
# 
