''' 
    sRNA_targets.py
    1/12/2019
    This program finds matches in the sRNA positions of "hot regions" and the predicted targets
    to determine which targets to test experimentally for interactions
'''

# Use Anaconda, a distribution of Python for data science
# Run on Spyder

# In order to run, make sure Excel file is closed
# Otherwise, it will return an error saying "Permission Denied"

# Import Pandas package, used for data analysis (Excel).
import pandas as pd
from openpyxl import load_workbook
from xlrd import XLRDError

# Read in Excel file given file name.
file_name = 'example.intarna.for.overlap.code.xlsx'
file = pd.ExcelFile(file_name)

# Load workbook so you can add new sheets w/o overwriting existing ones.
# Use the specific name of your Excel file.
book = load_workbook(file_name)
writer = pd.ExcelWriter(file_name, engine='openpyxl') 
writer.book = book
writer.sheets = dict((ws.title, ws) for ws in book.worksheets)


# Parse Excel sheets into separate DataFrames for predictions and hotRegions.
# Use the specific names of the sheets in your file.
predictions = file.parse('example.intarna.for.overlap.cod')
hotRegions = file.parse('hot.region.matrix')


# Read in Excel sheets that contain GO terms of interest for each sRNA and GO terms
# for each mRNA.
sRNA_GO = pd.read_excel('sRNA_GO.xlsx')
mRNA_GO = pd.read_excel('Ecoli_K12_MG1655_GO_for_BL.xlsx')

            
# Function for determining if mRNA has sRNA GO term of interest:
def has_GO_term(j, mRNA_GO, GO_terms, GO_row, each):
    
    empty = True
    for r in range(4,7):
        
        if pd.isnull(GO_row.iloc[0,r]) == False:
            empty = False
            
            for g in GO_terms:
                if g in GO_row.iloc[0, r]:
                    # If mRNA has GO term of interest, put 1 in GO term column.
                    each.loc[j, 'GO term'] = 1
                    return each
                
    # If mRNA has no GO terms, leave GO term column empty.
    if empty == True:
        each.loc[j, 'GO term'] = ''
    # If mRNA has GO terms but not GO term of interest, put 0 in GO term column.
    else:
        each.loc[j, 'GO term'] = 0
    return each

                   
# Loop through every row in hotRegions DataFrame.
for i, row in hotRegions.iterrows():
    
    # Column for the number of matches in positions of prediction and hotRegion:
    matchCol = hotRegions.loc[i]['sRNA'] + ' [' + str(hotRegions.loc[i]['region start']) + ',' + str(hotRegions.loc[i]['region end']) + ']'
    
    # Get all GO terms of interest for each sRNA in hotRegions DataFrame.
    GO_terms = []
    for s, row in sRNA_GO.iterrows():
        if hotRegions.loc[i,'sRNA'] == sRNA_GO.loc[s,'sRNA name']:
            GO_terms.append(str(sRNA_GO.loc[s]['GO_of_interest']))
        
        
    # Column for percent: number of matches / length of hotRegion:
    percentCol = '% [' + str(hotRegions.loc[i]['region start']) + ',' + str(hotRegions.loc[i]['region end']) + ']'
        
    # Set of all ints between hotRegion start and end positions:
    hotRegionsPos = set(range(hotRegions.loc[i]['region start'], \
                              hotRegions.loc[i]['region end'] + 1))


    # try except else statement determines whether to make new sheet for sRNA or add
    # to an existing sRNA sheet.
    try:

        # Check if there is an existing sheet with sRNA name.
        existing = pd.read_excel(writer, sheet_name = hotRegions.loc[i]['sRNA'])
    
    
    # Exception if there is no existing sheet for the sRNA --> create new sRNA sheet.
    # Sheet belonging to each sRNA will only contain information relevant to it.
    except XLRDError:
        
        # Copy predictions into new DataFrame so each sRNA can be modified separately.
        each = predictions.copy()
        
        # Make both lowercase and removes whitespaces to check if sRNA names are equal.
        # where() function changes rows that don't satisfy condition to NaN.
        each.where(each['sRNA_name'].str.lower().str.strip() == \
                   hotRegions.loc[i]['sRNA'].lower().strip(), inplace=True)
        # Drop NaN rows so only rows with matching sRNA names remain.
        each.dropna(subset=['target_name'], inplace=True)
        
        # Drop rows in original predictions DataFrames that are in each DataFrame.
        # This saves time because further iterations won't loop through same mRNAs again.
        predictions.drop(each.index, axis=0, inplace=True)
        
        # Reset indices because rows were deleted.
        each.reset_index(drop=True, inplace=True)
        
        
        # Add columns for match, match percentage, and GO term.
        each = pd.concat([each, pd.DataFrame(columns=[matchCol, percentCol, 'GO term'])], sort=False)
        
        # Loop through every row in predictions DataFrame to find matches in sRNA
        # and position.
        for j, row in each.iterrows():

            # Using sets, find matches between the sRNA start and end and the hotRegion
            # positions.
            predictionsPos = set(range(int(each.loc[j]['sRNA_start_pos']), \
                                       int(each.loc[j]['sRNA_end_pos'] + 1)))
            match = hotRegionsPos & predictionsPos

            # Fill matchCol with the number of matches. More matches means higher
            # likelihood of interaction with predicted targets.
            each.loc[j, matchCol] = len(match)
            
            # Add percentage matches (# of matches/total length of hotRegion) to
            # percentCol.
            each.loc[j, percentCol] = float(len(match))/len(hotRegionsPos)
            
            
            # Get row for current mRNA in GO term database, returns a 1-row DataFrame.
            GO_row = mRNA_GO.loc[mRNA_GO['CommonGene Name'].str.lower().str.strip() \
                                  == each.loc[j]['target_name'].lower().strip()]

            # If mRNA not found in GO term database, leave GO term column empty.
            if GO_row.empty:
                each.loc[j, 'GO term'] = ''
            # If mRNA is found in GO term database, update each DataFrame with whether
            # or not mRNA has GO term of interest using has_GO_term() method.
            else:
                each = has_GO_term(j, mRNA_GO, GO_terms, GO_row, each) 


        # Rank the sRNA sequences by energy and adds to new 'Rank' column
        # Lowest energy is the highest rank because more favorable
        each['Rank'] = each['energy'].rank()
            
        # Reorder columns so 'Rank' before other info
        # Make a list of columns in the DataFrame
        cols = list(each.columns.values)
        # Remove columns from list
        cols.pop(cols.index(matchCol))
        cols.pop(cols.index(percentCol))
        cols.pop(cols.index('GO term'))
        cols.pop(cols.index('Rank'))
        # Add columns back in desired order
        each = each[cols+['Rank', 'GO term', matchCol, percentCol]] 
        # Sort sequences first by matches, then by energy
        # At the top, the most likely targets will be shown
        each = each.sort_values(by = [matchCol, 'Rank'], ascending=[0,1])
        
        # Add new sheet for sRNA
        each.to_excel(writer, hotRegions.loc[i]['sRNA'], index=False)
    
    
    # If there is already a sheet for that sRNA, the process is easier. Just add
    # two new columns to the existing sheet containing new information for the hotRegion.
    else:
        
        # Loop through every row in existing DataFrame to find matches in sRNA and
        # position.
        for j, row in existing.iterrows():
            
            # Find matches between the sRNA start and end and the hotRegion positions.
            predictionsPos = set(range(existing.loc[j]['sRNA_start_pos'], \
                                       existing.loc[j]['sRNA_end_pos'] + 1))
            match = hotRegionsPos & predictionsPos
            
            # Fill matchCol with the number of matches.
            existing.loc[j, matchCol] = len(match)
            
            # Add percentage matches to percentCol
            existing.loc[j, percentCol] = float(len(match))/len(hotRegionsPos)
        

        # Add columns for that hotRegion to the existing DataFrame from the sRNA's
        # sheet.
        existing = existing.sort_values(by = [matchCol, 'Rank'], ascending=[0,1])
        
        # Rewrite existing sheet with the DataFrame with the new hotRegion added.
        existing.to_excel(writer, hotRegions.loc[i]['sRNA'], index=False)

    writer.save()
