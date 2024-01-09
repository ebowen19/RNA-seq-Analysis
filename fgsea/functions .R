#PART 1: Dif Exp Analysis


# Function for differential expression analysis
# performDEA Function
#
# Description:
#   This function performs differential expression analysis (DEA) on gene expression data.
#   It assigns each sample to predefined groups, creates a design matrix, fits a linear model,
#   and computes contrasts to identify differentially expressed genes between the groups.
#   It saves the results of these comparisons and generates log fold change density plots.
#
# Parameters:
#   data - A data frame where each row represents a gene and each column a sample.
#          The first two columns should be gene identifiers and are excluded from analysis.
#
#   comparisonGroups - A named list of character vectors. Each name represents a group
#                      (e.g., "Sensitive", "Insensitive"), and each vector contains
#                      patterns to match column names in 'data' to this group.
#
#   resultPath - A string specifying the path where the results (CSV files and plots) 
#                will be saved.
#
#   contrastName - A string used as a base for naming output files, reflecting the
#                  nature of the comparison being performed.
#
# Usage:
#   performDEA(data, comparisonGroups, resultPath, contrastName)
#
# Example:
#   comparisonGroups <- list(
#     Sensitive = c("RVN", "H2"),
#     Insensitive = c("H1", "H12", "RC")
#   )
#   performDEA(data, comparisonGroups, "path/to/results", "Sensitive_vs_Insensitive")
#
# Output:
#   The function outputs CSV files containing the results of the differential expression analysis
#   for each contrast and PNG files of the log fold change density plots. These files are saved
#   in the specified 'resultPath' directory.

# Function for differential expression analysis with voom
performDEA <- function(data, comparisonGroups, resultPath, contrastName) {
  # Create a DGEList object from the normalized count data
  dge <- DGEList(counts = data)

  # Apply the voom transformation
  v <- voom(dge, plot = TRUE)  # plot=TRUE generates a mean-variance plot

  # Assign each column in 'v' to the respective comparison group
  groups_vector <- sapply(colnames(v$E), function(col) {
    group_assignment <- NA  # Default assignment
    for (group in names(comparisonGroups)) {
      # Create a regex pattern that matches the group name followed by a dot and a number
      patterns <- paste0("^", comparisonGroups[[group]], "\\.[0-9]+$")
      if (any(sapply(patterns, grepl, col))) {
        group_assignment <- group
        break  # Exit the loop once a match is found
      }
    }
    return(group_assignment)
  })

  # Debugging: print out the groups_vector
  print("Group assignments for each column:")
  print(groups_vector)

  # Exclude columns that do not match any group
  valid_columns <- !is.na(groups_vector)
  v <- v$E[, valid_columns]  # Modify as needed based on the structure of 'v' post-voom
  groups_vector <- groups_vector[valid_columns]
  groups <- factor(groups_vector, levels = names(comparisonGroups))

  # Debugging: print out the factor levels of groups
  print("Factor levels of groups:")
  print(levels(groups))

  # Create a design matrix with no intercept
  design <- model.matrix(~ 0 + groups)

  # Debugging: print out the design matrix
  print("Design matrix:")
  print(design)

  # Fit the model using voom output
  fit <- lmFit(v, design)

  # Construct and evaluate the contrast within makeContrasts
  # Assuming 'groupsExperimental' and 'groupsBaseline' are defined correctly
  cont.matrix <- makeContrasts(Contrasts = groupsExperimental - groupsBaseline, levels = design)

  # Debugging: print out the contrast matrix
  print("Contrast matrix:")
  print(cont.matrix)

  fit2 <- eBayes(contrasts.fit(fit, cont.matrix))

  # Get the results for the contrast
  results <- topTable(fit2, coef = "Contrasts", n = Inf)

  # Save results
  csv_filename <- paste0(resultPath, "/", contrastName, ".csv")
  write.csv(results, file = csv_filename)

  # Uncomment to save plot
  # png_filename <- paste0(resultPath, "/", contrastName, "_LogFC_Density.png")
  # png(file = png_filename)
  # plot(density(results$logFC), main = paste("Density of Log Fold Change -", contrastName))
  # dev.off()
}


convertToSymbols <- function(df, species) {
  if (species == "Mus musculus") {
    library(org.Mm.eg.db)
    symbols <- mapIds(org.Mm.eg.db, keys=df$X, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  } else if (species == "Homo sapiens") {
    library(org.Hs.eg.db)
    symbols <- mapIds(org.Hs.eg.db, keys=df$X, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  } else {
    stop("Invalid species specified: ", species)
  }
  df$X <- symbols
  return(df)
}

generateVolcanoPlotsAndSaveData <- function(data_frames) {
  for(name in names(data_frames)) {
    volcano_data <- data_frames[[name]]
    
    # Determine significance and direction (up or down regulated)
    volcano_data$significant <- ifelse(volcano_data$adj.P.Val < 0.05 & abs(volcano_data$logFC) > 3.5, "yes", "no")
    volcano_data$direction <- ifelse(volcano_data$logFC > 0, "up", "down")
    
    # Plot
    p <- ggplot(volcano_data, aes(x = logFC, y = -log10(P.Value))) +
      geom_point(aes(color = ifelse(significant == "yes" & direction == "up", "Significant Up", 
                                    ifelse(significant == "yes" & direction == "down", "Significant Down", "Not Significant"))), 
                 alpha = 0.5) +
      theme_minimal() +
      labs(title = paste(name, "volcano plot"),
           x = "Log2 Fold Change",
           y = "-Log10 p-value",
           color = "Gene Regulation") +
      scale_color_manual(values = c("Not Significant" = "gray", 
                                    "Significant Up" = "red", 
                                    "Significant Down" = "blue")) +
      theme(plot.title = element_text(hjust = 0.5))
    
    # Print and save the plot
    print(p)
    ggsave(filename = paste0(name, " volcano plot.jpg"), plot = last_plot())
    
    # Filter the data
    filtered_data <- volcano_data %>%
      filter(significant == "yes")
    
    # Count the number of significant points that show up on the graph
    num_significant_points <- nrow(filtered_data)
    print(paste(num_significant_points, "significant points for", name))
    write.csv(filtered_data, paste0(name, " significant points.csv"), row.names = TRUE)
    
    # Save list of upregulated & downregulated genes
    upregulated <- filtered_data %>%
      filter(direction == "up")
    print(paste(nrow(upregulated), "upregulated"))
    downregulated <- filtered_data %>%
      filter(direction == "down")
    print(paste(nrow(downregulated), "downregulated"))
    
    # make & save data frames for the up/downregulated genes with both the gene id & symbols
    upregulated <- data.frame(gene_id = upregulated[1], symbol = convertToSymbols(upregulated[1], "Mus musculus"))
    colnames(upregulated) <- c("gene id", "symbol")
    downregulated <- data.frame(gene_id = downregulated[1], symbol = convertToSymbols(downregulated[1], "Mus musculus"))
    colnames(downregulated) <- c("gene id", "symbol")
    write.csv(upregulated, paste(name, "upregulated genes.csv"))
    write.csv(downregulated, paste(name, "downregulated genes.csv"))
  }
}

# PART 2: GSEA

# Function: loadAndOrganizeGeneSets
#
# Description:
#   This function loads and organizes gene sets for specified species and categories.
#   It utilizes the msigdbr package to fetch gene sets for Homo sapiens and Mus musculus 
#   across different categories like GO, REACTOME, KEGG, and HALLMARK. The function 
#   organizes these gene sets into a structured list, making it convenient for further analysis.
#
# Returns:
#   A list containing organized gene sets. The list is structured by species, each containing 
#   a nested list of categories. Each category list contains gene sets split by gene symbol and 
#   gene set name.
#
# Example of use:
#   geneSets <- loadAndOrganizeGeneSets()

loadAndOrganizeGeneSets <- function() {
  # Load gene sets
  hs_go_all <- msigdbr(species = "Homo sapiens", category = "C5")
  mm_go_all <- msigdbr(species = "Mus musculus", category = "C5")
  hs_reactome_all <- msigdbr(species = "Homo sapiens",  subcategory = "CP:REACTOME")
  mm_reactome_all <- msigdbr(species = "Mus musculus", subcategory = "CP:REACTOME")
  hs_kegg_all <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
  mm_kegg_all <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
  hs_hallmark_all <- msigdbr(species = "Homo sapiens", category = "H")
  mm_hallmark_all <- msigdbr(species = "Mus musculus", category = "H")
  
  # Organize gene sets
  hs_hallmark_sets_list <- split(hs_hallmark_all$gene_symbol, hs_hallmark_all$gs_name)
  hs_reactome_sets_list <- split(hs_reactome_all$gene_symbol, hs_reactome_all$gs_name)
  hs_go_sets_list <- split(hs_go_all$gene_symbol, hs_go_all$gs_name)
  hs_kegg_sets_list <- split(hs_kegg_all$gene_symbol, hs_kegg_all$gs_name)
  
  mm_hallmark_sets_list <- split(mm_hallmark_all$gene_symbol, mm_hallmark_all$gs_name)
  mm_reactome_sets_list <- split(mm_reactome_all$gene_symbol, mm_reactome_all$gs_name)
  mm_go_sets_list <- split(mm_go_all$gene_symbol, mm_go_all$gs_name)
  mm_kegg_sets_list <- split(mm_kegg_all$gene_symbol, mm_kegg_all$gs_name)
  
  # Combine them into a single list
  geneSets <- list(
    "Homo sapiens" = list("HALLMARK" = hs_hallmark_sets_list, "REACTOME" = hs_reactome_sets_list, "GO" = hs_go_sets_list, "KEGG" = hs_kegg_sets_list),
    "Mus musculus" = list("HALLMARK" = mm_hallmark_sets_list, "REACTOME" = mm_reactome_sets_list, "GO" = mm_go_sets_list, "KEGG" = mm_kegg_sets_list)
  )
  
  return(geneSets)
}


# Conversion function from Ensembl ID to Gene Symbol
#
# Arguments:
#   df: A dataframe containing gene expression data. The dataframe is expected
#       to have a column 'X' which contains the Ensembl IDs.
#   species: A string specifying the species of the genes. It can be either
#            'Mus musculus' or 'Homo sapiens'.
#
# Return:
#   The function returns a dataframe similar to the input but with the 'X' column
#   now containing gene symbols instead of Ensembl IDs.

convertToSymbols <- function(df, species) {
  if (species == "Mus musculus") {
    library(org.Mm.eg.db)
    symbols <- mapIds(org.Mm.eg.db, keys=df$X, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  } else if (species == "Homo sapiens") {
    library(org.Hs.eg.db)
    symbols <- mapIds(org.Hs.eg.db, keys=df$X, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  } else {
    stop("Invalid species specified: ", species)
  }
  df$X <- symbols
  return(df)
}


# The performGSEA function performs Gene Set Enrichment Analysis using the fgsea package.
#
# Arguments:
#   dataframe: A data frame containing gene expression data. Expected to have columns 'logFC' for 
#              log fold change and 'X' for gene names.
#   geneSetList: A list of gene sets from a certain database where each gene set is a character vector of gene symbols.
#   comparisonName: A string representing the name of the comparison being analyzed.
#   outputPath: A string representing the path where the GSEA output files will be saved.
#   species: specify species you're working with
#   capping: option to cap extreme/biologically implausible logFC values
#
# The function sorts genes based on their log fold change, uses that ranking to perform GSEA, and saves the results
# as a CSV file in the specified output path.
performGSEA <- function(dataframe, geneSetList, comparisonName, outputPath, species, capping = TRUE) {
  # Check if the first gene identifier looks like an Ensembl ID
  if (startsWith(dataframe$X[1], "ENS")) {
    # Convert Ensembl IDs to Gene Symbols
    dataframe <- convertToSymbols(dataframe, species)
  }
  
  # # Debugging: Print summary of the input data
  # print(paste("Analysis for:", comparisonName))
  # print(summary(dataframe))
  
  geneList <- dataframe$logFC
  names(geneList) <- dataframe$X
  rankedGenes <- sort(geneList, decreasing = TRUE)
  
  # Capping logFC values in the rankedGenes vector at +/- 20
  if (capping == TRUE) {
    rankedGenes <- sapply(rankedGenes, function(gene) {
      if (gene > 20) {
        return(20)
      } else if (gene < -20) {
        return(-20)
      } else {
        return(gene)
      }
    })
  } else {
    break
  }
  
  
  # # Debugging -- closer look @ ranking process
  # print("Top ranked genes:")
  # print(head(rankedGenes))
  # 
  # print("Bottom ranked genes:")
  # print(tail(rankedGenes))
  
  # Save a plot of the log fold change distribution
  plotFilename <- paste0(outputPath, "/", comparisonName, "_LogFC_Density.png")
  png(filename = plotFilename)
  plot(density(rankedGenes), main=paste("Density of Log Fold Change -", comparisonName))
  dev.off() # Close the plotting device
  
  
  # run the fgsea
  fgseaResults <- fgsea(
    pathways = geneSetList,
    stats = rankedGenes, # using the capped vals
    minSize = 15,
    maxSize = 500
  )
  
  if (nrow(fgseaResults) == 0) {
    print(paste("No significant pathways found for ", comparisonName))
  } else {
    fgseaResults <- as.data.frame(fgseaResults[order(fgseaResults$NES, decreasing = TRUE), ])
    fgseaResults$leadingEdge <- sapply(fgseaResults$leadingEdge, function(x) paste(x, collapse = ";"))
    #   print(head(fgseaResults))
    write.csv(fgseaResults, file = paste0(outputPath, "/fgseaResults_", comparisonName, ".csv"))
    # output # of up/downregulated pathways 
    downregulated <- sum(fgseaResults$NES < 0)
    upregulated <- sum(fgseaResults$NES > 0)
    print(paste0("# of significant pathways for ", comparisonName, ": ", nrow(fgseaResults)))
    print(paste(upregulated, "upregulated &", downregulated, "downregulated"))
    
  }
}

# The runAnalysis function calls the performGSEA function to perform Gene Set Enrichment Analysis (GSEA) 
# for a set of comparisons with the specified gene set type (e.g., GO, KEGG, Hallmark) for each comparison.
#
# Arguments:
#   comparisons: A list of comparisons. Each comparison is a list containing:
#                - dataPath: A string representing the file path to the CSV file for the comparison.
#                - species: A string indicating the species ("Homo sapiens" or "Mus musculus").
#   geneSetType: A string representing the type of gene set to be used for GSEA (e.g., "GO", "KEGG").
#   geneSets: A nested list of gene sets. The outer list is keyed by species, and each inner list
#             contains the actual gene sets keyed by type (e.g., "GO", "KEGG").
#   outputPath: A string representing the path where the output files will be saved.
#
# Returns:
#   resultFiles: A list of file paths to the saved result files. Each element in the list corresponds
#                to a comparison, with the list names being the comparison names and the values being
#                the paths to the respective result files. This list can be used for further processing
#                or for generating reports.
#
# The function reads the data for each comparison, performs GSEA using the specified gene set type,
# and saves the results to the specified output path.
runAnalysis <- function(comparisons, geneSetType, geneSets, outputPath) {
  resultFiles <- list()
  
  # iterate through the comparisons u want to perform gsea on
  for (comparisonName in names(comparisons)) {
    comparisonDetails <- comparisons[[comparisonName]]
    species <- comparisonDetails$species
    
    # Check if species exists in geneSets & gene set type exists  for the  species
    if (!species %in% names(geneSets)) {
      stop(paste("Species not found in geneSets:", species)) #error line
    }
    if (!geneSetType %in% names(geneSets[[species]])) {
      stop(paste("Gene set type not found for species", species, ":", geneSetType)) #error line
    }
    # load in data
    comparisonData <- read.csv(comparisonDetails$dataPath)
    geneSetList <- geneSets[[species]][[geneSetType]]
    
    # Perform GSEA
    performGSEA(comparisonData, geneSetList, comparisonName, outputPath, species)
    
    # Add the path to the result file to the list
    resultFilePath <- paste0(outputPath, "/fgseaResults_", comparisonName, ".csv")
    resultFiles[[comparisonName]] <- resultFilePath
  }
  
  # Return the list of result file paths
  return(resultFiles)
}

# The extractCommonPathways function processes the results of GSEA for multiple comparisons,
# aligns the dataframes, finds common pathways, and saves the results.
#
# Arguments:
#   dataframes: A named list of data frames, where each data frame corresponds to fgsea output
#               for a specific comparison. Each data frame should contain gene pathway data with
#               columns such as 'pathway', 'NES', etc. The names of the list elements should
#               correspond to the names of the comparisons.
#   outputPath: A string representing the path where the aligned output files will be saved.
#
# Returns:
#   pathways_list: A list where each element corresponds to a specific comparison. Each element is a 
#                  vector containing the common pathways identified in that comparison. This list can
#                  be used for further analysis or reporting on the common pathways across different 
#                  comparisons.
#
# The function aligns each comparison's pathways based on a reference, identifies common pathways
# across comparisons, and writes the aligned data and common pathways to CSV files.

extractCommonPathways <- function(dataframes, referenceComparison, outputPath) {
  # Check for NULL dataframes and remove them
  dataframes <- dataframes[sapply(dataframes, function(x) !is.null(x))]
  
  # Check if referenceComparison exists in dataframes
  if (!(referenceComparison %in% names(dataframes))) {
    stop(paste("Reference comparison not found:", referenceComparison))
  }
  
  # Separate into upregulated & downregulated pathways
  up_down_all <- lapply(dataframes, function(df) {
    list(
      up = df %>% dplyr::filter(NES > 0) %>% dplyr::arrange(desc(NES)) %>% dplyr::select(pathway),
      down = df %>% dplyr::filter(NES < 0) %>% dplyr::arrange(NES) %>% dplyr::select(pathway)
    )
  })
  
  # Define function to align pathways based on reference
  alignPathways <- function(pathways, referencePathways) {
    # Create aligned column
    aligned_column <- rep(NA, length(referencePathways))
    
    # Loop through each reference pathway
    for (i in seq_along(referencePathways)) {
      if (referencePathways[i] %in% pathways) {
        aligned_column[i] <- referencePathways[i]
      }
    }
    return(aligned_column)
  }
  
  # Initialize lists to store aligned pathways
  aligned_up_pathways <- list()
  aligned_down_pathways <- list()
  
  # Reference pathways
  reference_up_pathways <- up_down_all[[referenceComparison]]$up$pathway
  reference_down_pathways <- up_down_all[[referenceComparison]]$down$pathway
  
  # Align pathways
  for (comparison in names(up_down_all)) {
    aligned_up_pathways[[comparison]] <- alignPathways(up_down_all[[comparison]]$up$pathway, reference_up_pathways)
    aligned_down_pathways[[comparison]] <- alignPathways(up_down_all[[comparison]]$down$pathway, reference_down_pathways)
  }
  
  # Convert to dataframes
  aligned_up_df <- data.frame(aligned_up_pathways)
  aligned_down_df <- data.frame(aligned_down_pathways)
  
  # Set 1st column to be the reference comparison
  aligned_up_df <- aligned_up_df[c(referenceComparison, setdiff(names(aligned_up_df), referenceComparison))]
  aligned_down_df <- aligned_down_df[c(referenceComparison, setdiff(names(aligned_down_df), referenceComparison))]
  
  # Save the aligned dataframes
  write.csv(aligned_up_df, paste0(outputPath, "/Aligned_Upregulated_Pathways.csv"), row.names = FALSE)
  write.csv(aligned_down_df, paste0(outputPath, "/Aligned_Downregulated_Pathways.csv"), row.names = FALSE)
  
  # IDENTIFYING COMMON PATHWAYS after aligned pathway dataframes are generated
  # Define Function to find pathways present in all but one comparison
  findAlmostCommonPathways <- function(df) {
    n <- ncol(df)
    almost_common <- apply(df, 1, function(x) {
      missing_index <- which(is.na(x))
      if (length(missing_index) == 1 && sum(!is.na(x)) == (n - 1)) {
        return(c(pathway = x[1], missing_from = names(df)[missing_index]))
      }
      return(NULL)
    })
    
    # Filter out NULL entries and convert to dataframe
    if (length(almost_common) == 0) { # error message for case when there's no common or almost common pathways
      stop(paste("No common or almost common pathways found for", deparse(substitute(df)),"."))
    }
    else {
      almost_common <- do.call(rbind, almost_common)
    }
    return(almost_common)
  }
  
  # Initialize the return object
  pathways_list <- list(up = NULL, down = NULL)
  
  # Process for upregulated pathways, different process based on if there's any common pathways 
  # in all comparisons or not.
  if (any(complete.cases(aligned_up_df))) {
    common_up_pathways <- data.frame(aligned_up_df[complete.cases(aligned_up_df), 1])
    colnames(common_up_pathways)[1] <- "upregulated pathways"
    pathways_list$up <- common_up_pathways
    write.csv(common_up_pathways, paste0(outputPath, "/Common_Upregulated_Pathways.csv"), row.names = FALSE)
    print(paste(nrow(common_up_pathways), "common upregulated"))
  } else {
    almost_common_up_pathways <- findAlmostCommonPathways(aligned_up_df)
    colnames(almost_common_up_pathways)[1] <- "upregulated pathways"
    pathways_list$up <- almost_common_up_pathways
    write.csv(almost_common_up_pathways, paste0(outputPath, "/Almost_Common_Upregulated_Pathways.csv"), row.names = FALSE)
    print(paste(nrow(almost_common_up_pathways), "almost common upregulated"))
  }
  
  # Same Process for downregulated pathways
  if (any(complete.cases(aligned_down_df))) {
    common_down_pathways <- as.data.frame(aligned_down_df[complete.cases(aligned_down_df), 1])
    colnames(common_down_pathways)[1] <- "downregulated pathways"
    pathways_list$down <- common_down_pathways
    write.csv(common_down_pathways, paste0(outputPath, "/Common_Downregulated_Pathways.csv"), row.names = FALSE)
    print(paste(nrow(common_down_pathways), "common downregulated"))
  } else {
    almost_common_down_pathways <- findAlmostCommonPathways(aligned_down_df)
    colnames(almost_common_down_pathways)[1] <- "downregulated pathways"
    pathways_list$down <- almost_common_down_pathways
    write.csv(almost_common_down_pathways, paste0(outputPath, "/Almost_Common_Downregulated_Pathways.csv"), row.names = FALSE)
    print(paste(nrow(almost_common_down_pathways), "almost common downregulated"))
  }
  
  # Return the pathways list along with aligned dataframes
  return(list(pathways_list = pathways_list, aligned_up_df = aligned_up_df, aligned_down_df = aligned_down_df))
}

# Usage example:
# alignedPathways <- extractCommonPathways(gseaResults, "RVN_vs_RC", path_results)

# Function to read in (create dataframes) & group files by 2 last subfolders 
# apply to the files that have "pattern" in their file name, or if there's no pattern argument do all the files
# rowNames specifies if u wanna include row names when reading in the CSVs
group_files_by_subfolder <- function(main_folder, pattern = "") {
  # Recursively list all files in the directory
  files <- list.files(main_folder, full.names = TRUE, recursive = TRUE)
  
  # Initialize an empty list to store the list of dataframes grouped by subfolder
  grouped_list <- list()
  
  # Loop through the files
  for (file_path in files) {
    # Extract the folder and subfolder parts
    path_parts <- str_split(dirname(file_path), "/")[[1]]
    folder_key <- path_parts[length(path_parts) - 1]  # Second last element is the folder
    subfolder_key <- path_parts[length(path_parts)]  # Last element is the subfolder
    
    # Check if the file matches the pattern or if pattern is empty
    if (pattern == "" || grepl(pattern, file_path)) {
      # Regular expression to capture the part between the last "/" and the last "."
      df_key <- sub("^.*/(.*?)\\.[^./]*$", "\\1", file_path)
      
      # Extract the contents of the file into a dataframe
      df <- read.csv(file_path, row.names = NULL)
      
      # Check if the folder_key already exists, if not add it as an element key
      if (!folder_key %in% names(grouped_list)) {
        grouped_list[[folder_key]] <- list()
      }
      
      # Check if the subfolder_key already exists under the folder_key, if not add it
      if (!subfolder_key %in% names(grouped_list[[folder_key]])) {
        grouped_list[[folder_key]][[subfolder_key]] <- list()
      }
      
      # Add the dataframe to the list under the appropriate folder and subfolder
      grouped_list[[folder_key]][[subfolder_key]][[df_key]] <- df
    }
  }
  
  # Return the grouped list of dataframes
  return(grouped_list)
}


# Example usage
# grouped_dataframes <- group_files_by_subfolder('path_to_main_folder', 'fgseaResults')


#PART 3: ANALYZE COMMON GENES
# extract_leading_edge Function
# Description:
#   Extracts the leading edge genes from an fgsea results dataframe for a given pathway.
#   Returns NULL if the dataframe is NULL. Converts all gene names to lowercase if any 
#   are found to be in uppercase.
#
# Args:
#   fgsea_results_df: A dataframe containing fgsea results, expected to have a 'pathway' column
#                     with pathway names and a 'leadingEdge' column with semicolon-separated gene names.
#   pathway: A string representing the name of the pathway for which leading edge genes are to be extracted.
#
# Returns:
#   A character vector containing the leading edge genes for the specified pathway in lowercase, 
#   if the dataframe is not NULL and the pathway is found; otherwise, NULL.

extract_leading_edge <- function(fgsea_results_df, pathway) {
  if (!is.null(fgsea_results_df)) {
    leading_edge_genes <- fgsea_results_df[fgsea_results_df$pathway == pathway, "leadingEdge"]
    
    # If there are no genes, return NULL
    if (is.null(leading_edge_genes) || length(leading_edge_genes) == 0) {
      return(NULL)
    }
    
    leading_edge_genes <- unlist(strsplit(leading_edge_genes, ";")) %>% tolower()
    
    return(leading_edge_genes)
  }
  return(NULL)
}

#' Extract Common Leading Edge Genes
#'
#' This function iterates over fgsea results and extracts leading edge genes for pathways
#' categorized as common or almost common, and as upregulated or downregulated. It handles
#' the exclusion of specific comparisons for almost common pathways.
#'
#' @param fgseaResults A nested list of fgsea result dataframes, structured by comparison categories and subcategories.
#' @param commonPathways_master A nested list containing dataframes of pathway names,
#'        structured by comparison categories and subcategories. Dataframes include a column
#'        for pathways and, for almost common pathways, an additional column indicating the
#'        comparison from which the pathway is missing.
#' @return A nested list mirroring the structure of `commonPathways_master`, where each
#'         pathway name is replaced with a vector of common leading edge genes across the valid comparisons.
#'         The list is further categorized into upregulated and downregulated, as well as common and almost common.
#'
#' @examples
#' # fgseaResults_list and common_pathways_master should be predefined
#' results <- extract_common_leading_edges(fgseaResults_list, common_pathways_master)
#'
extract_common_leading_edges <- function(fgseaResults, commonPathways_master, exclude_comparison=NULL) {
  final_results <- list()
  # browser()
  # Loop over categories and subcategories
  for (category in names(commonPathways_master)) {
    final_results[[category]] <- list()
    
    for (subcategory in names(commonPathways_master[[category]])) {
      final_results[[category]][[subcategory]] <- list()
      
      # Loop through all dataframes in the subcategory
      for (df_name in names(commonPathways_master[[category]][[subcategory]])) {  
        df <- commonPathways_master[[category]][[subcategory]][[df_name]] # df of common or almost common pathways down or up
        
        #check & store names for common/almost commom & upregulated/downregulated
        is_almost_common <- grepl("Almost_Common", df_name)
        is_upregulated <- grepl("Upregulated", df_name)
        
        pathway_type <- paste0(
          ifelse(is_upregulated, "upregulated_", "downregulated_"),
          ifelse(is_almost_common, "almost_common", "common")
        )
        
        for (pathway in df[[1]]) { # make sure you're just referring to the 1st column (if there's 2 columns for almost common)
          # browser()
          exclude_comparison <- if (pathway_type == "Almost_Common") df[df[[1]] == pathway, 2] else NULL
          subcategory_comps <- if (!is.null(exclude_comparison)) {
            # if there's a comparison to exclude
            names(fgseaResults[[category]][[subcategory]])[names(fgseaResults[[category]][[subcategory]]) != exclude_comparison] 
          } else {
            names(fgseaResults[[category]][[subcategory]]) # if there's no comparison to exclude
          }
          # get the names of the subcategory comparisons we want to consider when finding the common genes (ignore the exclude_comparison)
          
          leading_edges <- setNames(
            lapply(seq_along(subcategory_comps), function(i) {
              fgsea_df <- fgseaResults[[category]][[subcategory]][[i]]
              comp <- subcategory_comps[i]
              if (!is.null(fgsea_df) && (is.null(exclude_comparison) || comp != exclude_comparison)) {
                extract_leading_edge(fgsea_df, pathway)
              }
            }), # set content of leading_edges list
            subcategory_comps  # Assign the names of comparisons as the names of the list elements
          )
          
          # Find the intersection of leading edge genes across the valid comparisons
          common_genes <- Reduce(intersect, leading_edges)
          
          final_results[[category]][[subcategory]][[pathway_type]][[pathway]] <- common_genes
        }
      }
    }
  }
  
  return(final_results)
}

# save_final_results Function
# Description:
#   Saves the final results of common leading edge gene analysis into a structured folder hierarchy
#   with CSV files. The folder structure is created based on categories, subcategories, and pathway types,
#   mirroring the structure of the 'final_results' list.
#
# Args:
#   final_results: A nested list containing common leading edge genes, structured by categories, 
#                  subcategories, and pathway types.
#   base_path: The base directory path where the folder structure and CSV files will be created and saved.
#
# Example usage:
#   base_path <- "~/Desktop/YourFolder"  # Replace with your actual folder path
#   save_final_results(final_results, base_path)
save_final_results <- function(final_results, base_path, file_name_template = "{pathway_type}_genes.csv") {
  for (category in names(final_results)) {
    cat_folder <- file.path(base_path, category)
    if (!dir.exists(cat_folder)) {
      dir.create(cat_folder, recursive = TRUE)
    }
    for (subcategory in names(final_results[[category]])) {
      subcat_folder <- file.path(cat_folder, subcategory)
      if (!dir.exists(subcat_folder)) {
        dir.create(subcat_folder)
      }
      for (pathway_type in names(final_results[[category]][[subcategory]])) {
        pathway_dfs <- list()
        for (pathway in names(final_results[[category]][[subcategory]][[pathway_type]])) {
          genes <- final_results[[category]][[subcategory]][[pathway_type]][[pathway]]
          pathway_df <- data.frame(
            Pathway = pathway,
            Genes = paste(genes, collapse = ", ")
          )
          pathway_dfs[[pathway]] <- pathway_df
        }
        combined_df <- do.call(rbind, pathway_dfs)
        rownames(combined_df) <- NULL
        final_results[[category]][[subcategory]][[pathway_type]] <- combined_df
        
        # Construct the filename using the template and replacing the placeholder
        file_name <- gsub("\\{pathway_type\\}", pathway_type, file_name_template)
        write.csv(combined_df, file.path(subcat_folder, file_name), row.names = FALSE)
      }
    }
  }
  return(final_results)  # Return the modified final_results with the new dataframes
}


# TALLY UP GENES INTO HTML FILE

# Function to tally up genes in a dataframe and create a gene count table
# returns table with col of gene & their counts
tally_genes <- function(df) {
  # browser()
  gene_counts <- df %>%
    separate_rows(Genes, sep = ",\\s*") %>% # separate vals into multiple rows
    filter(Genes != "") %>% # filter empty rows
    count(Genes, name = "Count") # count occurences of each unique gene & assign to a new col called count 
  
  return(gene_counts)
}

# Function to process a list of dataframes and create tallied gene count tables
# based on "downregulated" or "upregulated" in the dataframe names
process_list <- function(list_df) {
  # browser()
  # initialize a list called "results" with down & upregulated elements (both initially empty/null)
  results <- list(Downregulated = NULL, Upregulated = NULL)
  
  # loop over ea dataframe using name
  for (df_name in names(list_df)) {
    if (grepl("downregulated", df_name, ignore.case = TRUE)) {  # not case-sensitive
      # grepl returns a bool vector based on if a pattern exists in each element of a char vector
      results$Downregulated <- bind_rows(results$Downregulated, tally_genes(list_df[[df_name]]))
      # results$Downregulated: This is the target data frame to which rows are being added or extended. 
      # It is assumed to be a data frame or tibble that is already initialized
      # tally_genes(list_df[[df_name]]): This is the data frame generated by calling the tally_genes function on the 
      # contains the gene counts for a specific category, and it is being added to the results$Downregulated data frame
      
    } else if (grepl("upregulated", df_name, ignore.case = TRUE)) {
      results$Upregulated <- bind_rows(results$Upregulated, tally_genes(list_df[[df_name]]))
    }
  }
  
  return(results)
}

# Compile results from different categories into a single *list of tables*
# apply the process_list function to each category in the input list, 
# resulting in a list of tallied gene count tables for each category.
compile_results <- function(data_list) {
  all_results <- lapply(data_list, process_list)
  return(all_results)
}

# process_and_combine_gene_data Function
#
# Summary:
#   Consolidates gene count data into separate downregulated and upregulated data frames
#   with a total count per gene across all categories.
#
# Description:
#   Takes a nested list of gene count data frames and processes them to create two wide-format data frames: 
#   one for downregulated genes and one for upregulated genes. Each data frame includes genes as rows, 
#   categories as columns, and adds a column for the total count per gene.
#
# Parameters:
#   gene_data_list - A nested list where each element is a list with two data frames named "Downregulated" 
#                    and "Upregulated", containing 'Genes' and 'Count' columns.
#
# Returns:
#   A list with two elements:
#     - Downregulated: A data frame of genes vs. categories with total downregulated counts.
#     - Upregulated: A data frame of genes vs. categories with total upregulated counts.
#
# Example:
#   combined_gene_data <- process_and_combine_gene_data(gene_counts_list)
process_and_combine_gene_data <- function(gene_data_list) {
  
  # Function:
  # Processes a list of gene count data frames for either downregulated or upregulated genes,
  # consolidating counts by gene and adding a total count column for the specified regulation type.
  #
  # Inputs:
  #   data_list: A list where each element is a data frame with 'Genes' and 'Count' columns.
  #   regulation_type: String indicating the regulation type ('Downregulated_Total' or 'Upregulated_Total').
  #
  # Output:
  #   A wide-format data frame with genes as rows, categories as columns, and a total count column.
  process_category <- function(data_list, regulation_type) {
    combined_df <- do.call(rbind, lapply(names(data_list), function(category_name) {
      df <- data_list[[category_name]]
      if (is.null(df)) {
        warning(paste("The data frame for", category_name, regulation_type, "is NULL. Skipping this category."))
        return(NULL)  # Return NULL to be filtered out later
      }
      # Add the category name as a new column in the dataframe
      df %>%
        mutate(Category = category_name) %>%
        select(Genes, Count, Category)
    })) %>%
      na.omit()  # Remove rows with NULL values resulting from the above check
    
    # Pivot the dataframe from long to wide format and calculate totals
    wide_df <- combined_df %>%
      # make the category column vals into the colnames of the wide df. fill it w/ vals from count col & turn na vals to 0
      pivot_wider(names_from = Category, values_from = Count, values_fill = list(Count = 0)) %>%
      # add col called Total tallying up counts for each gene
      mutate(Total = rowSums(select(., -Genes), na.rm = TRUE)) %>%
      rename(!!regulation_type := Total)  # !! is used to unquote the variable holding the column name for dynamic naming
    
    return(wide_df)
  }
  # Process downregulated and upregulated data frames
  # '[[' operator extracts elements by name -->  retrieve data frames associated with  "Downregulated" category from each list element
  downregulated_df <- process_category(lapply(gene_data_list, `[[`, "Downregulated"), "Downregulated_Total") %>% arrange(desc(Downregulated_Total))
  upregulated_df <- process_category(lapply(gene_data_list, `[[`, "Upregulated"), "Upregulated_Total") %>% arrange(desc(Upregulated_Total))
  
  # Return a list containing both the sorted data frames
  return(
    list(
      Downregulated = downregulated_df, 
      Upregulated = upregulated_df
    )
  )
}

# Generates HTML tables for downregulated and upregulated gene count data frames.
#
# Parameters:
#   df: list w/ up & downregulated dataframes of gene counts 
#   comparison_type: String prefix for naming the output HTML file.
#
# Returns:
#   The file path of the generated HTML file containing the tables.
#
# Example:
#   html_file_og <- generate_html_output(down_df, up_df, "Original_Comparisons")

generate_html_geneTables <- function(comparison_results, comparison_type) {
  # Define custom CSS for table styling and title appearance
  # CSS for the gene regulation tables
  custom_css <- "
<style>
  body { font-family: 'Arial', sans-serif; }
  .table { font-size: 14px; }
  .table td, .table th {
    border: 1px solid #ddd;
    text-align: center;
    padding: 8px;
  }
  .table { border-collapse: collapse; }
  .table-striped > tbody > tr:nth-of-type(odd) { background-color: #f9f9f9; }
  .table-hover > tbody > tr:hover { background-color: #f5f5f5; }
  caption {
    font-family: 'Times New Roman', Times, serif;
    color: #007bff; /* Adjust this blue to the shade you want */
    font-size: 24px; /* Keep font size large */
    font-weight: bold; /* Keep font bold */
    text-align: center;
    margin-top: 20px; /* Add space above the caption */
  }
</style>"
  
  
  
  # Function to create an HTML table from a dataframe
  create_html_table <- function(df, file_name, comparison_type) {
    # Replace underscores with spaces and capitalize for title
    title <- paste(comparison_type, file_name, "Genes") %>%
      gsub("_", " ", .) %>%
      tools::toTitleCase()
    kable(df, caption = title, format = "html", align = "c") %>%
      kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>%
      column_spec(1:ncol(df), width = "2in") %>%
      as.character()
  }
  
  # Extract downregulated and upregulated data frames
  downregulated_df <- comparison_results$Downregulated
  upregulated_df <- comparison_results$Upregulated
  
  # Create HTML tables for both downregulated and upregulated data frames
  down_html <- create_html_table(downregulated_df, "Downregulated", comparison_type)
  up_html <- create_html_table(upregulated_df, "Upregulated", comparison_type)
  
  # Combine custom CSS with HTML tables into one HTML content string
  html_content_down <- paste(custom_css, down_html)
  html_content_up <- paste(custom_css, up_html)
  
  # Save to HTML files
  html_file_down <- paste0(gsub(" ", "", comparison_type), "_Down_GeneTable.html")
  html_file_up <- paste0(gsub(" ", "", comparison_type), "_Up_GeneTable.html")
  
  writeLines(html_content_down, html_file_down)
  writeLines(html_content_up, html_file_up)
  
  return(list(Downregulated = html_file_down, Upregulated = html_file_up))
}

# Generate an index HTML file for a collection of Gene Regulation Reports.
# This function creates an index page with links to different categories of reports.
#
# Args:
#   main_folder: The main folder where the index.html file will be created.
#
# Returns:
#   The path to the generated index.html file.
#
# Side Effects:
#   - Writes the HTML content to the index.html file.
#   - Checks if the linked HTML files exist in the specified directory.
#
# Example Usage:
#   main_folder <- "/path/to/reports"
#   index_html_file <- generate_index_html(main_folder)
#   # Creates index.html with links to report categories.

generate_index_html <- function(main_folder) {
  # Define custom CSS for the index page styling
  custom_css <- "
  <style>
    body { font-family: 'Times New Roman', Times, serif; margin: 20px; }
    h1 { color: #333; }
    h2 { font-family: 'Times New Roman', Times, serif; font-weight: normal; }
    .report-link { margin-left: 20px; font-size: 18px; }
    .category { margin-top: 20px; }
    a:link, a:visited { color: #007bff; }
  </style>"
  
  # Start building the HTML content
  index_html_content <- paste0("
  <html>
  <head>
    <title>Gene Regulation Reports</title>", custom_css, "
  </head>
  <body>
    <h1>Gene Regulation Reports</h1>
    <div class='category'>
      <h2>Original Comparisons</h2>
      <div class='report-link'><a href='Original Comparisons/Original_Comparisons_Down_GeneTable.html'>Downregulated Genes</a></div>
      <div class='report-link'><a href='Original Comparisons/Original_Comparisons_Up_GeneTable.html'>Upregulated Genes</a></div>
    </div>
    <div class='category'>
      <h2>New Comparisons</h2>
      <div class='report-link'><a href='New Comparisons/New_Comparisons_Down_GeneTable.html'>Downregulated Genes</a></div>
      <div class='report-link'><a href='New Comparisons/New_Comparisons_Up_GeneTable.html'>Upregulated Genes</a></div>
    </div>
  </body>
  </html>")
  
  # Define the output file path
  output_file <- file.path(main_folder, "index.html")
  
  # Write the HTML content to the file
  writeLines(index_html_content, output_file)
  
  # Debug: Check if the file exists after writing
  if (file.exists(output_file)) {
    message("Index file created at: ", output_file)
  } else {
    stop("Failed to create the index file.")
  }
  
  # Debug: Check if the linked files exist
  message("Checking if linked files exist:")
  linked_files <- c(
    'Original Comparisons/OriginalComparisons_Down_GeneTable.html',
    'Original Comparisons/OriginalComparisons_Up_GeneTable.html',
    'New Comparisons/NewComparisons_Down_GeneTable.html',
    'New Comparisons/NewComparisons_Up_GeneTable.html'
  )
  
  for (file in linked_files) {
    file_path <- file.path(main_folder, file)
    if (!file.exists(file_path)) {
      warning("Linked file does not exist: ", file_path)
    } else {
      message("Linked file exists: ", file_path)
    }
  }
  
  return(output_file)
}

# Function to convert CSV to HTML with RStudio-like table styling
convert_csv_to_html <- function(csv_file_path, output_html_file_path) {
  # Read the CSV file
  data <- read.csv(csv_file_path)
  
  # Remove the first column if its name is 'X'
  if ("X" %in% names(data)) {
    data <- data[, -1]
  }
  
  # Remove rows where the 'Genes' column is an empty string
  data <- subset(data, Genes != "")
  
  # Generate a title based on the file path segments
  path_parts <- unlist(strsplit(csv_file_path, split = "/"))
  comparison_type <- path_parts[length(path_parts) - 2] # Second last segment for comparison type
  analysis_type <- path_parts[length(path_parts) - 1] # Last segment for analysis type
  analysis_type <- sub("\\.csv$", "", analysis_type) # Remove file extension
  analysis_type <- gsub("_(.*)", " \\1", analysis_type) # Replace underscores with spaces
  title <- paste(comparison_type, analysis_type, sep = ": ") # Combine for title
  
  # Define CSS styles to mimic RStudio's data frame appearance
  table_css <- "
  <style>
  body { font-family: Arial, sans-serif; }
  table {
    border-collapse: collapse;
    width: 100%;
    margin-top: 20px;
  }
  th, td {
    border: 1px solid #dddddd;
    text-align: left;
    padding: 8px;
  }
  th {
    background-color: #f2f2f2;
  }
  tr:nth-child(even) {
    background-color: #f9f9f9;
  }
  .title {
    font-size: 24px;
    font-weight: bold;
    margin-top: 20px;
  }
  </style>"
  
  # Create the title for the HTML table
  title_html <- paste0("<div class='title'>", title, "</div>")
  
  # Convert to HTML using kable from knitr
  html_table <- kable(data, format = "html", table.attr = "class='table' style='border-collapse: collapse;'")
  
  # Combine the CSS, title, and the HTML table
  html_content <- paste0(table_css, title_html, html_table)
  
  # Write the combined HTML content to file
  writeLines(html_content, output_html_file_path)
}

# Function to find CSV files and convert them
convert_csv_files_to_html <- function(folder_path) {
  # List all CSV files in the folder and subfolders
  csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  
  # Loop over each CSV file and convert it to HTML
  for (csv_file in csv_files) {
    # Create the HTML file path by replacing the .csv extension with .html
    html_file <- sub("\\.csv$", ".html", csv_file)
    
    # Convert the CSV to HTML
    convert_csv_to_html(csv_file, html_file)
    
    # Print out the path to the HTML file
    message("Created HTML file: ", html_file)
  }
}


