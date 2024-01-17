import os
import pandas as pd
from bs4 import BeautifulSoup

#%%

root_folder = ('/Users/elizabeth 1/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/LW15 '
               'analysis/LW15-Target-Genes/Common Genes/Volc & PCA Data')
os.chdir(root_folder)
csv_files = {}

# Loop through all files in the directory
# Loop through all files in the directory
for filename in os.listdir(root_folder):
    # Check if the file is a CSV file
    if filename.endswith('.csv'):
        # Construct the full path to the CSV file
        csv_file_path = os.path.join(root_folder, filename)
        # Read the CSV file
        df = pd.read_csv(csv_file_path)

        # Convert the CSV file name to HTML file name
        html_filename = filename.replace('.csv', '.html')
        # Construct the full path to the HTML file
        html_file_path = os.path.join(root_folder, html_filename)

        # Extract the base filename without the .csv extension to use as a title
        file_title = os.path.splitext(filename)[0]

        # Convert the DataFrame to HTML
        html_string = df.to_html(index=False)

        # Add a heading to the HTML with the file title
        html_string_with_title = f"<h1>{file_title}</h1>\n{html_string}"

        # Write the HTML content to a file
        with open(html_file_path, 'w') as f:
            f.write(html_string_with_title)

        print(f'Converted {filename} to HTML and saved to {html_file_path}')


#%%
filename = 'sensi_vs_non_total_genes.csv'
csv_file_path = os.path.join(root_folder, filename)
# Read the CSV file
df = pd.read_csv(csv_file_path)

# Convert the CSV file name to HTML file name
html_filename = filename.replace('.csv', '.html')
# Construct the full path to the HTML file
html_file_path = os.path.join(root_folder, html_filename)

# Extract the base filename without the .csv extension to use as a title
file_title = os.path.splitext(filename)[0]

# Convert the DataFrame to HTML
html_string = df.to_html(index=False)

# Add a heading to the HTML with the file title
html_string_with_title = f"<h1>{file_title}</h1>\n{html_string}"

# Write the HTML content to a file
with open(html_file_path, 'w') as f:
    f.write(html_string_with_title)
#%%
