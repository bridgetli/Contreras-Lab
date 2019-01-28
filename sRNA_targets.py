''' 
    sRNAHiddenTargets.py
    1/12/2019
    This program finds matches in the sRNA positions of "hot regions" and the predicted targets
    to determine which targets to test experimentally for interactions
'''

# Use Anaconda, a distribution of Python for data science
# Run on Spyder

# In order to run, make sure Excel file is closed
# Otherwise, it will return an error saying "Permission Denied"

# Import Pandas package, used for data analysis (Excel)
import pandas as pd
from openpyxl import load_workbook
from xlrd import XLRDError

# Read in Excel file given file name
file = pd.ExcelFile('example.intarna.for.overlap.code.xlsx')
file.sheet_names

# Load workbook so you can add new sheets w/o overwriting existing ones
# Use the specific name of your Excel file
book = load_workbook('example.intarna.for.overlap.code.xlsx')
writer = pd.ExcelWriter('example.intarna.for.overlap.code.xlsx', engine='openpyxl') 
writer.book = book
writer.sheets = dict((ws.title, ws) for ws in book.worksheets)

# Parse Excel sheets into separate DataFrames for predictions and hotRegions 
# Use the specific names of the sheets in your file
predictions = file.parse('example.intarna.for.overlap.cod')
hotRegions = file.parse('hot.region.matrix')

# Loop through every row in hotRegions DataFrame
for i, row in hotRegions.iterrows():
    
    # Column for the number of matches in positions of prediction and hotRegion
    matchCol = hotRegions.loc[i]['sRNA'] + ' [' + str(hotRegions.loc[i]['region start']) + ',' + str(hotRegions.loc[i]['region end']) + ']'
    
    # Column for percent: number of matches / length of hotRegion
    percentCol = '% [' + str(hotRegions.loc[i]['region start']) + ',' + str(hotRegions.loc[i]['region end']) + ']'
        
    # Set of all ints between hotRegion start and end positions
    hotRegionsPos = set(range(hotRegions.loc[i]['region start'], hotRegions.loc[i]['region end'] + 1))

    # try except else statement determines whether to make new sheet for sRNA or add to an existing sRNA sheet
    try:

        # Check if there is an existing sheet with sRNA name
        existing = pd.read_excel(writer, sheet_name = hotRegions.loc[i]['sRNA'])
    
    # Exception if there is no existing sheet for the sRNA --> create new sRNA sheet
    # Sheet belonging to each sRNA will only contain information relevant to it
    except XLRDError:

        # Copy predictions into new DataFrame so each sRNA can be modified separately
        each = predictions.copy(deep=False)
    
        # Loop through every row in predictions DataFrame to find matches in sRNA and position
        for j, row in each.iterrows():
            
            # Make both lowercase and removes whitespaces to check if equal
            if hotRegions.loc[i]['sRNA'].lower().strip() == each.loc[j]['sRNA_name'].lower().strip():
                
                # Using sets, find matches between the sRNA start and end and the hotRegion positions
                predictionsPos = set(range(predictions.loc[j]['sRNA_start_pos'], predictions.loc[j]['sRNA_end_pos'] + 1))
                match = hotRegionsPos & predictionsPos
    
                # Fill matchCol with the number of matches
                # More matches means higher likelihood of interaction with predicted targets
                each.loc[j, matchCol] = len(match)
                
                # Add percentage matches (# of matches/total length of hotRegion) to percentCol
                each.loc[j, percentCol] = float(len(match))/len(hotRegionsPos)
            
            else:
                
                # If sRNA names not matching, delete row
                each.drop(j, inplace=True)
        
        # Reset index because rows were deleted
        each.reset_index(drop=True)
        
        # Rank the sRNA sequences by energy and adds to new 'Rank' column
        # Lowest energy is the highest rank because more favorable
        each['Rank'] = each['energy'].rank()
            
        # Reorder columns so 'Rank' before other info
        cols = list(each.columns.values) # Make a list of columns in the DataFrame
        cols.pop(cols.index(matchCol)) # Remove columns from list
        cols.pop(cols.index(percentCol))
        cols.pop(cols.index('Rank'))
        each = each[cols+['Rank', matchCol, percentCol]] # Add columns back in desired order
            
        # Sort sequences first by matches, then by energy
        # At the top, the most likely targets will be shown
        each = each.sort_values(by = [matchCol, 'Rank'], ascending=[0,1])
            
        # Add new sheet for sRNA
        each.to_excel(writer, hotRegions.loc[i]['sRNA'])
    
    # If there is already a sheet for that sRNA, the process is easier
    # Just add two new columns to the existing sheet containing new information for the hotRegion
    else:
        
        # Loop through every row in existing DataFrame to find matches in sRNA and position
        for j, row in each.iterrows():
            
            # Find matches between the sRNA start and end and the hotRegion positions
            predictionsPos = set(range(existing.loc[j]['sRNA_start_pos'], existing.loc[j]['sRNA_end_pos'] + 1))
            match = hotRegionsPos & predictionsPos
            
            # Fill matchCol with the number of matches
            existing.loc[j, matchCol] = len(match)
            
            # Add percentage matches to percentCol
            existing.loc[j, percentCol] = float(len(match))/len(hotRegionsPos)
        
        # Add columns for that hotRegion to the existing DataFrame from the sRNA's sheet
        existing = existing.sort_values(by = matchCol, ascending = 0)
        
        # Rewrite existing sheet with the DataFrame with the new hotRegion added
        existing.to_excel(writer, hotRegions.loc[i]['sRNA'])

    writer.save()
