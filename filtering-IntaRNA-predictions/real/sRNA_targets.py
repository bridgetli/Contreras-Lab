''' 
    sRNA_targets.py
    8/06/2020
    This program finds matches in the sRNA positions of "hot regions" and the predicted targets
    to determine which targets to test experimentally for interactions.
'''

# Use Anaconda, a distribution of Python for data science
# Run on Spyder

# In order to run, make sure Excel file is closed
# Otherwise, it will return an error saying "Permission Denied"

# Import Pandas package, used for data analysis (Excel).
import pandas as pd
from openpyxl import load_workbook
import glob


# Read in hot regions Excel file given file name.
hot_regions_file = 'hot_region_matrix.xlsx'
file = pd.ExcelFile(hot_regions_file)

# Load workbook so you can add new sheets w/o overwriting existing ones.
# Use the specific name of your Excel file.
book = load_workbook(hot_regions_file)
writer = pd.ExcelWriter(hot_regions_file, engine='openpyxl') 
writer.book = book
writer.sheets = dict((ws.title, ws) for ws in book.worksheets)

hotRegions = pd.read_excel(hot_regions_file)


# Read in Excel sheets that contain GO terms of interest for each sRNA and GO terms
# for each mRNA.
sRNA_GO = pd.read_excel('sRNA_GO.xlsx')
mRNA_GO = pd.read_excel('Ecoli_K12_MG1655_GO_for_BL.xlsx')


# Get files for individual sRNAs from inta_outputs folder.
path = 'inta_outputs'
filenames = glob.glob(path + '/*.csv')

# New DataFrame for top 5 predictions based on energy for each sRNA:
top5 = pd.DataFrame(columns = ['id1','start1','end1','id2','region','start2','end2',
                               'E','Rank','GO term','# overlap','% overlap'])


# Function for determining if mRNA has sRNA GO term of interest:
def has_GO_term(j, mRNA_GO, GO_terms, GO_row, df):
    
    empty = True
    for r in range(4,7):
        
        if pd.isnull(GO_row.iloc[0,r]) == False:
            empty = False
            
            for g in GO_terms:
                if g in GO_row.iloc[0, r]:
                    # If mRNA has GO term of interest, put 1 in GO term column.
                    df.loc[j, 'GO term'] = 1
                    return df
                
    # If mRNA has no GO terms, leave GO term column empty.
    if empty == True:
        df.loc[j, 'GO term'] = ''
    # If mRNA has GO terms but not GO term of interest, put 0 in GO term column.
    else:
        df.loc[j, 'GO term'] = 0
    return df


# Loop through each row of hotRegions dataframe to find matching name.
for i, row in hotRegions.iterrows():
    
    # Loop through all files in inta_outputs folder.
    for filename in filenames:
    
        # Parse Excel sheets into DataFrames for predictions.
        predictions = pd.read_csv(filename)
        
        # Sort sequences first by matches, then by energy
        # At the top, the most likely targets will be shown
        predictions = predictions.sort_values(by=['E'])
        predictions.reset_index(inplace=True)
        
        
        if hotRegions.loc[i, 'sRNA'].lower() == predictions.loc[0, 'id2'].lower():
            
            # Add column for original IntaRNA rank.
            cols = list(predictions.columns.values)
            if cols[-1] == 'E':
                rank_list = list(range(1, len(predictions)+1))
                predictions['Rank'] = rank_list
            else:
                cols = list(predictions.columns.values)
                predictions.columns = cols[:-1] + ['Rank']
        
            # Store hot region coordinates.
            start = str(hotRegions.loc[i]['accepted_region_start'])
            end = str(hotRegions.loc[i]['accepted_region_end'])
            
            # Column for the number of matches in positions of prediction and hotRegion:
            matchCol = '[' + start + ',' + end + ']'
            
            # Get all GO terms of interest for each sRNA in hotRegions DataFrame.
            GO_terms = []
            for s, row in sRNA_GO.iterrows():
                if hotRegions.loc[i,'sRNA'] == sRNA_GO.loc[s,'sRNA name']:
                    GO_terms.append(str(sRNA_GO.loc[s]['GO_of_interest']))
                
            # Column for percent: number of matches / length of hotRegion:
            percentCol = '% [' + start + ',' + end + ']'
                
            # Set of all ints between hotRegion start and end positions:
            hotRegionsPos = set(range(int(start), int(end) + 1))
            
            # Create column containing hot region coordinates.
            predictions['region'] = '[' + start + ',' + end + ']'
            
            # Add columns for match, match percentage, and GO term.
            predictions = pd.concat([predictions, pd.DataFrame(columns=[matchCol, percentCol, 'GO term'])], sort=False)
    
    
            # Loop through every row in sRNA DataFrame to find matches in sRNA and position.
            j = 0
            while j < len(predictions):
            
                # Find matches between the sRNA start and end and the hotRegion positions.
                predictionsPos = set(range(int(predictions.loc[j,'start2']), \
                                           int(predictions.loc[j,'end2']) + 1))
                match = hotRegionsPos & predictionsPos
                
                # Calculate percent overlap between hot region and sRNA region.
                percent_overlap = float(len(match))/len(hotRegionsPos)
                
                # If percent overlap is less than 0.7, filter out the prediction.
                if percent_overlap < 0.8:
                    predictions.drop(predictions.index[j], inplace=True)
                    predictions.reset_index(drop=True, inplace=True)
                    
                else:
                            
                    # Add percentage matches to percentCol.
                    predictions.loc[j, percentCol] = percent_overlap
                    
                    # Add number of matches to matchCol.
                    predictions.loc[j, matchCol] = len(match)
            
                    # Get row for current mRNA in GO term database, returns a 1-row DataFrame.
                    GO_row = mRNA_GO.loc[mRNA_GO['CommonGene Name'].str.lower().str.strip() \
                                         == predictions.loc[j]['id1'].lower().strip()]
                
                    # If mRNA not found in GO term database, leave GO term column empty.
                    if GO_row.empty:
                        predictions.loc[j, 'GO term'] = ''
                    # If mRNA is found in GO term database, update each DataFrame with whether
                    # or not mRNA has GO term of interest using has_GO_term() method.
                    else:
                        predictions = has_GO_term(j, mRNA_GO, GO_terms, GO_row, predictions)
                
                    j += 1
            
            
            # Reorder columns
            predictions = predictions[['id1','start1','end1','id2','region','start2','end2',
                                       'E','Rank','GO term',matchCol,percentCol]]

            # Export predictions DataFrame for current sRNA to Excel.
            predictions.to_excel(writer, sheet_name=filename[13:-4]+ ' (' + start + ',' + end + ')',
                                 index=False)
            
            # Add first 5 rows from current sRNA DataFrame to top5 DataFrame.
            predictions = predictions.rename(columns={matchCol:"# overlap", percentCol:"% overlap"})
            top5 = pd.concat([top5, predictions.head()])
        
            writer.save()
            
            break

top5.to_excel(writer, sheet_name='top 5 merged', index=False)
writer.save()


# =============================================================================
# Find sRNA:mRNA pairs in RIL-seq datasets that match top INTERFACE predictions.
# =============================================================================

# Read in Excel file with RIL-seq datasets.
datasets = pd.ExcelFile('info_for_RILseq_support.xlsx')
file = pd.ExcelFile(hot_regions_file)
top5 = file.parse('top 5 merged')


# Function for formatting datasets_found list as a string for entry into top5 DataFrame:
def format_datasets(datasets_found):
    # Input: list of datasets where sRNA:mRNA pair was found
    # Output: string of datasets with commas and spaces in between
    string = datasets_found[0]
    
    i = 1
    while i < len(datasets_found):
        string += ', ' + datasets_found[i]
        i += 1
        
    return string
        
# Loop through each row in top5 DataFrame.
for i, row in top5.iterrows():
    
    # List to store which datasets sRNA:mRNA pair is found in
    datasets_found = []
    
    # Boolean to distinguish whether or not sRNA exists in datasets
    sRNA_found = False
    
    # Get sRNA:mRNA pair from top5 DataFrame
    sRNA = top5.loc[i,'id2'].lower()
    mRNA = top5.loc[i,'id1'].lower().split('_')[0]
    
    
    # Loop through each RIL-seq dataset to look for matching pair.
    for sheet in datasets.sheet_names:
        
        dataset = datasets.parse(sheet)
        
        # Loop through each row in RIL-seq dataset.
        for j, row in dataset.iterrows():
            
            # If top5 sRNA is in RNA1 and top5 mRNA is in RNA2, there is a match.
            # Append name of dataset to datasets_found list.
            if sRNA in dataset.loc[j,'RNA1 name'].lower():
                sRNA_found = True
                if mRNA in dataset.loc[j,'RNA2 name'].lower():
                    datasets_found.append(sheet)
            
            # The order of sRNA and mRNA may be switched, so also need to check
            # if sRNA is in RNA2 and mRNA is in RNA1 to identify all matches.
            elif sRNA in dataset.loc[j,'RNA2 name'].lower():
                sRNA_found = True
                if mRNA in dataset.loc[j,'RNA1 name'].lower():
                    datasets_found.append(sheet)
    
    
    # If pair is not found in any dataset, put 'N/A' in 'RIL-seq' column of top5
    # if sRNA found. Put 'not found' if sRNA not found in any dataset.
    if not datasets_found:
        if sRNA_found:
            top5.loc[i,'RIL-seq'] = 'N/A'
        else:
            top5.loc[i,'RIL-seq'] = 'sRNA not found'
            
    # If pair is found, put names of datasets where it was found in
    # 'RIL-seq' column.
    else:
        top5.loc[i,'RIL-seq'] = format_datasets(datasets_found)


# Export top5 DataFrame to Excel.
top5.to_excel(writer, sheet_name='top 5 merged', index=False)
writer.save()
