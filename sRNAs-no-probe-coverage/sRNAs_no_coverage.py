''' 
    sRNAs_no_coverage.py
    1/24/2020
    This program finds regions of sRNAs that don't have asRNA coverage. This will help us
    design probes for those uncovered regions.
'''

# Use Anaconda, a distribution of Python for data science.
# Run on Spyder.

# In order to run, make sure Excel file is closed.
# Otherwise, it will return an error saying "Permission Denied".

# Import Pandas package, used for data analysis (Excel).
import pandas as pd
from openpyxl import load_workbook

# Load workbook so you can add new sheets w/o overwriting existing ones.
# Use the specific name of your Excel file.
book = load_workbook('myco_length_for_BL.xlsx')
writer = pd.ExcelWriter('myco_length_for_BL.xlsx', engine='openpyxl') 
writer.book = book
writer.sheets = dict((ws.title, ws) for ws in book.worksheets)

# Read in Excel file given file name.
mycoLengthFile = pd.ExcelFile('myco_length_for_BL.xlsx')
mycoLengthFile.sheet_names
# Parse Excel sheet into DataFrame.
mycoLengthDF = mycoLengthFile.parse('myco_length_for_BL') 
# Parse csv file directly into DataFrame.
mycoProbesDF = pd.read_csv('passing_myco_probes_for_BL.csv')
regionsList = []

# User can input 'n' - sequence +/- n nucleotides.
# This is because probes can overlap a little with the covered regions.
n = int(input("Number of nucleotides outside uncovered region: "))

# Loop through every row in myco_length DataFrame.
for i, row in mycoLengthDF.iterrows():
    
    # Set of all ints up to sRNA length:
    length = mycoLengthDF.loc[i]['sRNA len'] + 1
    sRNAPos = set(range(0, length))
    # New column for regions without coverage (output):
    notCoveredCol = 'regions without coverage'
    # String to append new regions without coverage:
    notCovered = ""

    # Loop through every row in myco_probes DataFrame.
    for j, row in mycoProbesDF.iterrows():
        
        # If sRNA names are the same, create set of ints between probe start and end positions.
        # Remove ints in probe position set from sRNA position set,
        # so only ints not covered by probe remain in sRNA position set.
        if mycoLengthDF.loc[i]['sRNA'] == mycoProbesDF.loc[j]['sRNA']:
            probePos = set(range(int(mycoProbesDF.loc[j]['start nt']), int(mycoProbesDF.loc[j]['end nt'] + 1)))
            sRNAPos = sRNAPos - probePos
    
    # Convert sRNA position set into list and sort list, so it can be easily iterated over.
    posList = list(sRNAPos)
    posList.sort()
    
    # Create new list for 'bounds'- starts and ends of regions not covered.
    gaps = [[s, e] for s, e in zip(posList, posList[1:]) if s+1 < e]
    edges = iter(posList[:1] + sum(gaps, []) + posList[-1:])
    bounds = list(zip(edges, edges)) # Each set of bounds is stored as a tuple.

    # Iterate through even indices (start of region) in bounds list.
    for j in range(len(bounds)):
        tup = bounds[j]
        # Only consider region if it is more than 4 nucleotides long.
        if tup[1]-tup[0] > 4:
            
            # Get entire sRNA sequence from Excel sheet.
            sequence = mycoLengthDF.loc[i]['sequence']
            # Get sequence of uncovered region +/- n nucleotides.
            # If region is near the beginning, sequence starts at the beginning.
            if tup[0] <= n:
                if length - tup[1] > n:
                    sequence = sequence[:tup[1] + n+1]
            # If region is near the end, sequence goes to the end.
            elif length - tup[1] <= n:
                sequence = sequence[tup[0]-n-1:]
            # For all other regions in the middle:
            else:
                sequence = sequence[tup[0] - n:tup[1] + n+1]
            
            # Append sRNA name, bounds of uncovered region, and sequence of uncovered region
            # to new region.
            region = []
            region.append(mycoLengthDF.loc[i]['sRNA'])
            region.append(tup[0])
            region.append(tup[1])
            region.append(sequence)
        
            # Add new region to list of all uncovered regions for all sRNAs.
            regionsList.append(region)

# Create a new DataFrame from regionsList, so it can be exported to Excel.
regionsDF = pd.DataFrame(regionsList, columns = ['sRNA','start','end','sequence'])
    
# Create new sheet for regions DataFrame.
regionsDF.to_excel(writer, 'regions_without_coverage', index=False)
writer.save()
