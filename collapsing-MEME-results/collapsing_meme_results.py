''' 
    collapsing_meme_results.py
    7/16/2020
    This program organizes MEME results from three different motif databases (PRODORIC2,
    DPInteract, and SWISS). Redundant motifs are collapsed together.
'''

# Use Anaconda, a distribution of Python for data science
# Run on Spyder

# Import Pandas package, used for data analysis (Excel)
import pandas as pd
from openpyxl import load_workbook

file_name = 'FNR_combined_FNR_1.xlsx'
file = pd.ExcelFile(file_name)

# Load workbook so you can add new sheets w/o overwriting existing ones.
# Use the specific name of your Excel file.
book = load_workbook(file_name)
writer = pd.ExcelWriter(file_name, engine='openpyxl') 
writer.book = book
writer.sheets = dict((ws.title, ws) for ws in book.worksheets)

# Parse Excel sheets into data DataFrame. Use the specific name of the input sheet
# in your file.
MEME_results_df = file.parse('FNR_combined_FNR_1')
# Sort DataFrame by 'sequence name' so rows for the same sRNA will be adjacent.
MEME_results_df.sort_values(by='sequence name')

merged_df = pd.read_excel('merged_peaks_SNR_1_min.xlsx')


# Create new DataFrames for outputs.
output_1_df = pd.DataFrame(columns = ['sRNA_name','unique_motif','peak_name', 'min_FNR',
                                      'max_FNR', 'start_gen_coord', 'end_gen_coord', 
                                      '# databases'])
output_2_df = pd.DataFrame(columns = ['sRNA_name','proposed_TFs'])


# Determine which database the motif is from.
def which_database(i):
    
    motif = MEME_results_df.loc[i,'#pattern name']
    # If there is something in the 2nd column, the motif is from PRODORIC2.
    if pd.isnull(MEME_results_df.iloc[i,1]) == False:
        return 'PRODORIC2'
    # Else if the motif's first letter is uppercase, it is from SWISS.
    elif motif[0].isupper():
        return 'SWISS'
    # Otherwise the motif is from DPInteract.
    else:
        return 'DPInteract'


# Format motif names so they can be compared for equality and outputted correctly.
def format_motif(name):
    
    # Drop extra info at the end and make name lowercase.
    if '_' in name:
        name = name.rsplit('_',1)[0]
    elif ' ' in name:
        name = name[:name.index(' ')]
    name = name.lower()
    
    # Format based on length of motif name.
    if len(name) == 3:
        return name.capitalize()
    elif len(name) == 4:
        return name[:-1].capitalize() + name[-1].upper()
    elif len(name) > 4:
        return name[:3].capitalize() + name[3].upper()
    else:
        return name.capitalize()


# Format sRNA names so they can be compared for equality and outputted correctly.
def format_sRNA(name):
    
    name = name.rsplit('_',1)[0]
    if '_long' in name:
        return name[0].upper() + name[1:-6] + name[-6].upper() + '_long'
    else:
        return name[0].upper() + name[1:-1] + name[-1].upper()


# Retrieve information from merged_df given a peak from MEME_results_df.
def retrieve_from_merged_peaks(sRNA):
    
    # Find row with matching peak.
    i = merged_df.index[merged_df['sRNA peak'] == sRNA].tolist()[0]
    
    # List will contain all SNRs for the given peak.
    snr_list = []
    
    # Loop through each condition.
    dfNum = 5
    while dfNum < len(merged_df.columns)-6:
        
        # Append SNR for condition to SNR list.
        snr_list.append(merged_df.iloc[i, dfNum])
        dfNum += 3

    # Find the minimum and maximum SNRs in SNR list. If one of the SNRs is 'na',
    # minimum SNR is 'na'.
    if 'na' in snr_list:
        min_SNR = 'na'
        snr_list = list(filter(lambda a: a != 'na', snr_list))
    else:
        min_SNR = float(min(snr_list))
    max_SNR = float(max(snr_list))
    
    # Get the start and end genome coordinates for the given peak.
    start = merged_df.loc[i, 'Merged L']
    end = merged_df.loc[i, 'Merged R']

    return(min_SNR, max_SNR, start, end)


# Update output DataFrames with info for new peak.
def update_for_prev(motifs, sRNA, output_1_df, output_2_df):

    proposed_TFs = ''
    # Loop through all motifs for that peak.
    for motif in motifs:
        
        # Retrieve peak's minimum SNR, maximum SNR, start coord, and end coord from
        # merged_df
        min_SNR, max_SNR, start, end = retrieve_from_merged_peaks(sRNA)
        
        # Add a new row to output_1_df containing the sRNA name, motif, peak name,
        # min FNR, max FNR, start gen coord, end gen coord, and # of databases.
        new_row = {'sRNA_name':format_sRNA(sRNA), 'unique_motif':motif, 'peak_name':sRNA,
                   'min_FNR':min_SNR, 'max_FNR':max_SNR, 'start_gen_coord':start,
                   'end_gen_coord':end, '# databases':len(motifs[motif])}
        output_1_df = output_1_df.append(new_row, ignore_index=True)
        # Keep track of proposed transcription factors for output_2_df
        proposed_TFs = proposed_TFs + motif + ', '
    
    # Add a new row to output_2_df containing the peak and its proposed TFs.
    output_2_df = output_2_df.append({'sRNA_name':sRNA,
                        'proposed_TFs':proposed_TFs[:-2]}, ignore_index=True)
    
    # Return output_1_df and output_2_df so they can be modified outside the function.
    return output_1_df, output_2_df


# Reset data (motifs and their databases) for new sRNA.
def start_new(i, MEME_results_df):
    
    databases = [which_database(i)]
    # Dictionary that pairs each motif (key) to its list of databases (value):
    return {format_motif(MEME_results_df.loc[i,'#pattern name']):databases}


# Keep track of motifs and databases in a dictionary.
motifs = start_new(0, MEME_results_df)
# curr_sRNA is the sRNA peak in the row you are currently on, starting on row 0.
curr_sRNA = MEME_results_df.loc[0, 'sequence name']

i = 1
# Loop through all rows in DataFrame.
while i < len(MEME_results_df.index):
    
    # Set prev_sRNA to old curr_sRNA. prev_sRNA is the sRNA peak in the row before.
    prev_sRNA = curr_sRNA
    # Update curr_sRNA to sRNA peak in current row.
    curr_sRNA = MEME_results_df.loc[i, 'sequence name']
    
    # If current name is the same as the previous name, you are on the same peak.
    if curr_sRNA == prev_sRNA:
        
        # curr_motif is the motif in the row you are currently on.
        curr_motif = format_motif(MEME_results_df.loc[i,'#pattern name'])
        # Get database that the entry is from using the which_database() function.
        curr_database = which_database(i)
        
        # If the motif is already in the motifs dictionary:
        if curr_motif in motifs:
            
            # If the database is not already in that motif's value list, append the
            # database to the list. If the database is already in the list, no
            # changes needed.
            if curr_database not in motifs[curr_motif]:
                motifs[curr_motif].append(curr_database)
                
        # If curr_motif is not already in motifs, add an entry (motif:[database])
        else:
            motifs[curr_motif] = [curr_database]
        
    # If current name is not the same as the previous name, we are on a new peak.
    else:
        
        # Update outputs for previous peak, since we're finished with it.
        output_1_df, output_2_df = update_for_prev(motifs, prev_sRNA, output_1_df,
                                                   output_2_df)

        # Start new dictionary for current peak.
        motifs = start_new(i, MEME_results_df)

    i += 1

# Update outputs for final motif, since the loop cannot do that.
output_1_df, output_2_df = update_for_prev(motifs, curr_sRNA, output_1_df,
                                           output_2_df)


# Export output DataFrames to Excel.
output_1_df.to_excel(writer, 'output_1', index=False)
output_2_df.to_excel(writer, 'output_2', index=False)
writer.save()
