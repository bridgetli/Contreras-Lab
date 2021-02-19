''' 
    sRNA_CsrA_targets.py
    1/15/2019
    
    This program finds matches of CsrA and sRNA targets in large datasets of RNA:RNA
    interactions to identify potential sRNA:CsrA interactions to test experimentally
'''

# Use Anaconda, a distribution of Python for data science

# In order to run, make sure Excel file is closed
# Otherwise, it will return an error saying "Permission Denied"

# Import Pandas package, used for data analysis (Excel)
import pandas as pd
from openpyxl import load_workbook

file = pd.ExcelFile('info_for_identifying_shared_CsrA_sRNA_targets.xlsx')
file.sheet_names

# Load a workbook so you can add new sheets w/o overwriting existing ones
# Use the specific name of your Excel file
book = load_workbook('info_for_identifying_shared_CsrA_sRNA_targets.xlsx')
writer = pd.ExcelWriter('info_for_identifying_shared_CsrA_sRNA_targets.xlsx', engine='openpyxl') 
writer.book = book
writer.sheets = dict((ws.title, ws) for ws in book.worksheets)

# Parse sheets with sRNA and CsrA targets into separate DataFrames
CsrA_df = file.parse('CsrA_targets')
sRNA_df = file.parse('sRNAs_of_interest')

# Create a list of all CsrA targets (the rest of the DataFrame is not used)
CsrA_list = CsrA_df['mRNA'].tolist()

# Input [RNA dataset name, index of RNA_1 column, index of RNA_2 column] for all datasets
# Make sure you input the column indexes (starting at 0)
# You can add additional datasets and change current ones
input_list = [['Hfq_RIL_seq_Margalit', 0, 1], ['RNase_CLASH_Tree', 1, 6], ['Hfq_CLASH_Granneman', 1, 2]]

# Loop through each dataset inputted
for i in input_list:
    
    # Create a list of sRNA specific to that dataset using column from sRNA_df
    # Drop empty Strings from sRNA_list (sRNA not in the dataset) to save time looping
    sRNA_list = sRNA_df[i[0]].dropna().unique().tolist()
    # Make a DataFrame for that dataset by parsing its specific sheet
    dataset_df = file.parse(i[0])
    # Make a new DataFrame to store all sRNA:CsrA matches for that dataset
    results_df = pd.DataFrame(columns = list(dataset_df))
    
    # Loop through each sRNA in sRNA_list
    for sRNA in sRNA_list:
        
        # Most "N/A" show up as empty Strings in the list and have already been removed
        # However, if there are "N/A" remaining, continues (goes to next iteration of for loop)
        # to save time looping through the dataset
        if sRNA == "N/A":
            continue
            
        # Loop through each row in the RNA dataset
        for j, row in dataset_df.iterrows():
            
            # If the sRNA matches the RNA in the RNA_1 column, and if the RNA in the RNA_2 column
            # is contained in CsrA_list, there is a match
            # Append rows with matches to results_df

            if sRNA == dataset_df.iloc[j][i[1]]:
                if dataset_df.iloc[j][i[2]] in CsrA_list:
                    results_df = results_df.append(row, ignore_index=True)
                
            # The order of sRNA and CsrA may be switched, so also need to check if sRNA is in RNA_2
            # and CsrA is in RNA_1 to identify all matches
            # Append rows with matches to results_df
            elif sRNA == dataset_df.iloc[j][i[2]]:
                
                # For Hfq_RIL_seq_Margalit dataset, there may be extra info in the CsrA name
                # So to check if it matches a target, split the CsrA name by periods,
                # and check if any individual element is contained in CsrA_list
                if '.' in dataset_df.iloc[j][i[1]]:
                    split = dataset_df.iloc[j][i[1]].split(".")
                    if not set(split).isdisjoint(CsrA_list):
                        results_df = results_df.append(row, ignore_index=True)
                
                # Check normally if no periods in the CsrA name       
                elif dataset_df.iloc[j][i[1]] in CsrA_list:
                    results_df = results_df.append(row, ignore_index=True)

    # Add new sheet to Excel file for results_df
    results_df.to_excel(writer, "results_" + i[0])

writer.save()
