source('~/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/scripts/fgsea/functions.R')
library(knitr)
library(kableExtra)
library(rvest)
library(xml2)

main_folder <- "/Users/elizabeth 1/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/LW15 analysis/gsea/Common Genes"
# create & save html files for each gene table (containing total gene counts across gene sets)
filePath_og <- paste0(main_folder, "/Original Comparisons")
filePath_new <- paste0(main_folder, "/New Comparisons") 
setwd(filePath_og)
html_files_original <- generate_html_geneTables(final_results$original_comparisons, "Original Comparisons")
setwd(filePath_new)
html_files_new <- generate_html_geneTables(final_results$new_comparisons, "New Comparisons")

# generate index
setwd(main_folder)
index_html_file <- generate_index_html(main_folder)

# next link the colnames to the appropriate gene set data telling which genes 
# come from each pathway

# create html files for the pathway-gene data
# convert_csv_files_to_html(main_folder)

# Function to link column headers to HTML files with debugging
link_headers_to_files <- function(report_html_path, folder_path) {
  report_html <- read_html(report_html_path)
  cat("Loaded HTML content from:", report_html_path, "\n")
  
  # Get all the column header nodes
  header_nodes <- xml2::xml_find_all(report_html, "//th")
  
  for (header_node in header_nodes) {
    # Extract the text of the node
    header_text <- xml2::xml_text(header_node)
    cat("Processing header:", header_text, "\n")
    
    # Check if header text matches expected values and skip if it does
    if (header_text %in% c("Genes", "Downregulated_Total", "Upregulated_Total")) {
      cat("Skipping header:", header_text, "\n")
      next
    }
    
    # Clean up the header text to match file paths
    clean_header_text <- gsub(" ", "", header_text)
    
    # Determine if it is upregulated or downregulated from the file name
    reg_type <- ifelse(grepl("Up", report_html_path), "upregulated", "downregulated")
    comparison_type <- ifelse(grepl("New", report_html_path), "New Comparisons", "Original Comparisons")
    
    # Construct file name based on pattern observed
    file_patterns <- c("almost_common", "common")
    html_file_found <- FALSE
    
    for (pattern in file_patterns) {
      html_file_name <- paste0(reg_type, "_", pattern, " genes.html")
      html_file_path <- file.path(folder_path, comparison_type, clean_header_text, html_file_name)
      print(paste("file path:", html_file_path))
      if (file.exists(html_file_path)) {
        cat("Found HTML file for header:", clean_header_text, "\n")
        
        # Create the absolute link path
        absolute_link_path <- normalizePath(html_file_path, winslash = "/")
        cat("Absolute link path:", absolute_link_path, "\n")
        
        # Replace existing header text or link with the new hyperlink
        xml2::xml_set_text(header_node, "") # Clear the text content of the header node
        new_a_node <- xml2::xml_add_child(header_node, "a", header_text)
        xml2::xml_set_attr(new_a_node, "href", absolute_link_path)
        
        html_file_found <- TRUE
        break
      }
    }
    
    if (!html_file_found) {
      cat("No HTML file found for header:", clean_header_text, "\n")
    }
  }
  
  # Save the entire HTML document back to the file
  xml2::write_html(report_html, report_html_path)
  cat("Saved modified HTML content to:", report_html_path, "\n")
}



# Folder path where 'New Comparisons' and 'Original Comparisons' directories are located
folder_path <- main_folder

# List all report HTML files that need to be linked (you'll need to fill in these paths)
report_html_files <- list.files(path = folder_path, pattern = "Comparisons", full.names = TRUE, recursive = TRUE)

# Apply the function to each report
for (file in report_html_files) {
  link_headers_to_files(file, folder_path)
}
