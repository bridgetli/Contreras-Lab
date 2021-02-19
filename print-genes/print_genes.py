''' 
    print_genes.py
    1/22/2020
'''

# Use Anaconda, a distribution of Python for data science
# Run on Spyder

# In order to run, make sure Excel file is closed
# Otherwise, it will return an error saying "Permission Denied"

# Import Pandas package, used for data analysis (Excel)
import pandas as pd

# Read in Excel file given file name
file = pd.ExcelFile('data_sets_GEO_specs.xlsx')

# Loops through sheets in Excel file
for sheet_name in file.sheet_names:
    
    # Creates a DataFrame for each sheet, assuming there is no header row
    df = file.parse(sheet_name, header=None)
    
    # If there is a header row, updates the DataFrame to use the headers
    if not isinstance(df.loc[0][1], str):
        df = file.parse(sheet_name)
    
    # Creates a new text file for each sheet, using first word of the sheet as the file name
    txtFile = open(sheet_name.split(' ', 1)[0] + ".txt", "w")
    
    # Loop through rows in the DataFrame and append gene names to the text file
    for i, row in df.iterrows():
        txtFile.write("./sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump " + df.loc[i][1] + "\n")
    
    txtFile.close()
