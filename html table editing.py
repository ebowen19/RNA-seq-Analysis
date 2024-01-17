import pandas as pd
import os
from bs4 import BeautifulSoup


# %%
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


# %%
# Function to create Interactive Data tables
def add_datatables(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        html_content = file.read()

    soup = BeautifulSoup(html_content, 'html.parser')

    # Add an ID to the table
    table = soup.find('table')
    table['id'] = 'dataTable'

    # Find the index of the "Upregulated_Total" or "Downregulated_Total" column
    headers = table.find('thead').find_all('th')
    upregulated_index = next((i for i, th in enumerate(headers) if "Upregulated_Total" in th.text), -1)
    downregulated_index = next((i for i, th in enumerate(headers) if "Downregulated_Total" in th.text), -1)
    sort_index = upregulated_index if upregulated_index != -1 else downregulated_index

    # Add DataTables CSS and JS in the <head>
    head = soup.head
    jquery_script = soup.new_tag("script", src="https://code.jquery.com/jquery-3.5.1.js")
    head.append(jquery_script)
    datatables_css = soup.new_tag("link", rel="stylesheet", type="text/css",
                                  href="https://cdn.datatables.net/1.10.21/css/jquery.dataTables.css")
    head.append(datatables_css)
    datatables_script = soup.new_tag("script", src="https://cdn.datatables.net/1.10.21/js/jquery.dataTables.js",
                                     charset="utf8")
    head.append(datatables_script)

    # Initialize DataTables with dynamic sorting
    init_script = soup.new_tag("script")
    order_string = f"\"order\": [[{sort_index}, 'desc']], " if sort_index != -1 else ""
    init_script.string = "$(document).ready(function() { $('#dataTable').DataTable({ " + order_string + "\"lengthMenu\": [[-1], ['All']] }); });"
    soup.body.append(init_script)

    # Write the modified content back to a new file
    new_file_path = file_path.replace('.html', '_interactive.html')
    with open(new_file_path, 'w', encoding='utf-8') as file:
        file.write(str(soup))


# %%
def add_datatables_search(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        html_content = file.read()

    soup = BeautifulSoup(html_content, 'html.parser')

    # Ensure <html>, <head>, and <body> tags are present
    if not soup.html:
        html_tag = soup.new_tag('html')
        soup.append(html_tag)
    if not soup.head:
        head_tag = soup.new_tag('head')
        soup.html.insert(0, head_tag)
    if not soup.body:
        body_tag = soup.new_tag('body')
        soup.html.append(body_tag)

    head = soup.head

    # Add an ID to the table if not already present
    table = soup.find('table')
    if table and not table.has_attr('id'):
        table['id'] = 'dataTable'

    # Add DataTables CSS and JS in the <head>
    if not head.find("script", {"src": "https://code.jquery.com/jquery-3.5.1.js"}):
        jquery_script = soup.new_tag("script", src="https://code.jquery.com/jquery-3.5.1.js")
        head.append(jquery_script)
    if not head.find("link", {"href": "https://cdn.datatables.net/1.10.21/css/jquery.dataTables.css"}):
        datatables_css = soup.new_tag("link", rel="stylesheet", type="text/css",
                                      href="https://cdn.datatables.net/1.10.21/css/jquery.dataTables.css")
        head.append(datatables_css)
    if not head.find("script", {"src": "https://cdn.datatables.net/1.10.21/js/jquery.dataTables.js"}):
        datatables_script = soup.new_tag("script", src="https://cdn.datatables.net/1.10.21/js/jquery.dataTables.js",
                                         charset="utf8")
        head.append(datatables_script)

    # Initialize DataTables with search functionality and showing all entries
    init_script = soup.new_tag("script")
    init_script.string = "$(document).ready(function() { $('#dataTable').DataTable({ 'paging': false, 'searching': true, 'info': false }); });"
    soup.body.append(init_script)

    # Write the modified content back to a new file
    with open(file_path, 'w', encoding='utf-8') as file:
        file.write(str(soup))


# %%
main_path = ('/Users/elizabeth 1/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/LW15 '
             'analysis/LW15-Target-Genes/Common Genes')

# List of directories and their respective files to process
directories_files = {
    main_path + '/Original Comparisons': ['OriginalComparisons_Up_GeneTable.html',
                                          'OriginalComparisons_Down_GeneTable.html'],
    main_path + '/New Comparisons': ['NewComparisons_Up_GeneTable.html', 'NewComparisons_Down_GeneTable.html'],
    main_path + '/New Comparisons/Without trop2': ['Up_GeneTable.html', 'Down_GeneTable.html']
}
# %%
# Iterate over directories and files
for directory, files in directories_files.items():
    os.chdir(directory)
    for file in files:
        # modify_html_table(file)
        add_datatables(file)

# %%
# List Pathway Files (Which we want to add search ability to)

# Define the directory path where you want to start searching
start_directory = ('/Users/elizabeth 1/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/LW15 '
                   'analysis/LW15-Target-Genes/Common Genes')

# Create an empty list to store the HTML file paths
pathway_files = []

# Use os.walk() to traverse through the directory and its subfolders
for root, _, files in os.walk(start_directory):
    for file in files:
        # Check if the file has a .html extension
        if file.endswith('genes.html'):
            # Create the full file path and add it to the list
            file_path = os.path.join(root, file)
            pathway_files.append(file_path)

# Now, html_files contains a list of all HTML files in the specified directory and its subfolders
print(pathway_files)

for file_path in pathway_files:
    add_datatables_search(file_path)


# %%

def replace_string_in_files(file_paths, old_string, new_string):
    for file_path in file_paths:
        # Read the content of the file
        with open(file_path, 'r', encoding='utf-8') as file:
            content = file.read()

        # Replace the old string with the new string
        modified_content = content.replace(old_string, new_string)

        # Write the modified content back to the file
        with open(file_path, 'w', encoding='utf-8') as file:
            file.write(modified_content)


# List of file paths to modify
file_paths = [
    '/path/to/your/firstfile.html',
    '/path/to/your/secondfile.html',
    # Add more file paths as needed
]

# Strings to be replaced
old_string = 'The string to be replaced'
new_string = 'The new string'

# Call the function
replace_string_in_files(file_paths, old_string, new_string)
# Say 'Key Genes' & 'Enriched Pathways' in title
