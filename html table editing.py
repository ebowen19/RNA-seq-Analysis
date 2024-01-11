import pandas as pd
import os

main_path = '/Users/elizabeth 1/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/LW15 analysis/LW15-Target-Genes/Common Genes'

#%%

# Load HTML files
tables = {}

os.chdir(main_path+'/Original Comparisons')
files_og = ['OriginalComparisons_Up_GeneTable.html', 'OriginalComparisons_Down_GeneTable.html']
for file in files_og:
    tables[file] = pd.read_html(file)

os.chdir(main_path+'/New Comparisons')
files_new = ['NewComparisons_Up_GeneTable.html', 'NewComparisons_Down_GeneTable.html']
for file in files_new:
    tables[file] = pd.read_html(file)

os.chdir(main_path+'/New Comparisons/Without trop2')
files_trop = ['Up_GeneTable.html', 'Down_GeneTable.html']
for file in files_trop:
    tables[file] = pd.read_html(file)

#%%

# add column "Gene Sets Found In" that counts up non-zero values from 2nd to 5th columns for each row
for name,table in tables.items():
    df = table[0]
    columns_to_check = df.columns[1:5]
    df['Gene Sets Found In'] = df[columns_to_check].astype(bool).sum(axis=1)
    print(name)
    print(list(df['Gene Sets Found In']))
    print("")

#%%
for table in tables:
    print(table)

#%%

# Step 3: Write the modified table to an HTML string
modified_table_html = df.to_html(escape=False)  # escape=False to keep links intact

# Step 4: Preserve non-table HTML elements
file_path = 'OriginalComparisons_Up_GeneTable.html'
with open(file_path, 'r') as file:
    original_html = file.read()

# Assuming the table is the first element after the title
# Find the start of the original table and split the HTML
split_html = original_html.split('<table', 1)
before_table = split_html[0]

# Combine everything
final_html = before_table + modified_table_html

# Write the final HTML to a new file
with open('modified_file.html', 'w') as file:
    file.write(final_html)
#%%

main_path = '/Users/elizabeth 1/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/LW15 analysis/LW15-Target-Genes/Common Genes'

# Function to modify the HTML table
# Function to modify the HTML table
def modify_html_table(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        original_html = file.read()

    # Parse the HTML
    soup = BeautifulSoup(original_html, 'html.parser')

    # Find the table
    table = soup.find('table')

    # Add a new header cell for the new column
    header_row = table.find('tr')
    new_header_cell = soup.new_tag('th')
    new_header_cell.string = "Gene Sets Found In"
    header_row.append(new_header_cell)

    # Iterate through each row of the table, excluding the header
    for row in table.find_all('tr')[1:]:
        cells = row.find_all('td')
        if len(cells) > 4:  # Ensure there are enough cells
            # Count non-zero values from 2nd to 5th columns
            non_zero_count = sum(bool(int(cell.text)) for cell in cells[1:5])
            # Create a new cell with this count
            new_cell = soup.new_tag('td')
            new_cell.string = str(non_zero_count)
            row.append(new_cell)

    # Write the modified HTML back to file
    with open(file_path, 'w', encoding='utf-8') as file:
        file.write(str(soup))

#%%
# List of directories and their respective files to process
directories_files = {
    main_path + '/Original Comparisons': ['OriginalComparisons_Up_GeneTable.html', 'OriginalComparisons_Down_GeneTable.html'],
    main_path + '/New Comparisons': ['NewComparisons_Up_GeneTable.html', 'NewComparisons_Down_GeneTable.html'],
    main_path + '/New Comparisons/Without trop2': ['Up_GeneTable.html', 'Down_GeneTable.html']
}

# Iterate over directories and files
for directory, files in directories_files.items():
    os.chdir(directory)
    for file in files:
        modify_html_table(file)
