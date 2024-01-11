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
def add_datatables(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        html_content = file.read()

    soup = BeautifulSoup(html_content, 'html.parser')

    # Add an ID to the table
    table = soup.find('table')
    table['id'] = 'dataTable'

    # Find the index of the "Upregulated_Total" column
    headers = table.find('thead').find_all('th')
    upregulated_total_index = next((i for i, th in enumerate(headers) if "Upregulated_Total" in th.text), -1)

    # Add DataTables CSS and JS in the <head>
    head = soup.head
    jquery_script = soup.new_tag("script", src="https://code.jquery.com/jquery-3.5.1.js")
    head.append(jquery_script)
    datatables_css = soup.new_tag("link", rel="stylesheet", type="text/css", href="https://cdn.datatables.net/1.10.21/css/jquery.dataTables.css")
    head.append(datatables_css)
    datatables_script = soup.new_tag("script", src="https://cdn.datatables.net/1.10.21/js/jquery.dataTables.js", charset="utf8")
    head.append(datatables_script)

    # Initialize DataTables with dynamic sorting
    init_script = soup.new_tag("script")
    init_script.string = f"""
    $(document).ready(function() {{
        $('#dataTable').DataTable({{
            "order": [[{upregulated_total_index}, "desc"]],
            "lengthMenu": [[-1], ["All"]]
        }});
    }});
    """
    soup.body.append(init_script)

    # Write the modified content back to a new file
    new_file_path = file_path.replace('.html', '_interactive.html')
    with open(new_file_path, 'w', encoding='utf-8') as file:
        file.write(str(soup))


# %%
main_path = ('/Users/elizabeth 1/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/LW15 '
             'analysis/LW15-Target-Genes/Common Genes')

# List of directories and their respective files to process
directories_files = {
    main_path + '/Original Comparisons': ['OriginalComparisons_Up_GeneTable.html'],
    main_path + '/New Comparisons': ['NewComparisons_Up_GeneTable.html', 'NewComparisons_Down_GeneTable.html'],
    main_path + '/New Comparisons/Without trop2': ['Up_GeneTable.html']
}
#%%
# Iterate over directories and files
for directory, files in directories_files.items():
    os.chdir(directory)
    for file in files:
        # modify_html_table(file)
        add_datatables(file)

