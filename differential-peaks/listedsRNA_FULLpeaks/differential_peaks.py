''' 
    differential_peaks.py
    5/24/2020
    This program finds which sRNA peaks are expressed under certain conditions.
    This will help determine which transcription factors are involved for each sRNA.
'''

# Use Anaconda, a distribution of Python for data science.
# Run on Spyder.

# In order to run, make sure Excel file is closed.
# Otherwise, it will return an error saying "Permission Denied".


# =============================================================================
# Import packages.
# =============================================================================

# Import Pandas package, used for data analysis (Excel).
import pandas as pd
# Import numpy and math for dealing with empty DataFrame cells.
import numpy as np
import math
import itertools
import re


# =============================================================================
# List of files for each condition.
# To add new conditions, add name of file to end of array.
# =============================================================================

files = ['listedsRNA_FULLpeaks_rdmWT.csv',
         'listedsRNA_FULLpeaks_rdmStatWT.csv',
         'listedsRNA_FULLpeaks_minWT.csv']


# =============================================================================
# User input
# =============================================================================

# User can input 'n' - minimum number of nucleotides overlap in sRNA peaks.
# Default to 50 nt
n = 50
user_input = input("Minimum number of nucleotides overlap: ")
# User has to input an integer to be valid.
try:
    n = int(user_input)
except ValueError:
    print("Input not valid. Defaulting to 50 nucleotides.")

# User can input minimum SNR to consider a peak.
# Default to 0.5
min_SNR = 0.5
user_input = input("Minimum SNR to consider a peak: ")
# User must input a float to be valid.
try:
    min_SNR = float(user_input)
except ValueError:
    print("Input not valid. Defaulting to minimum SNR of 0.5")


# =============================================================================
# Various functions
# =============================================================================

# Function for adding new columns for each new condition/file:
def add_columns(df, dfNum):
    # Column for left coordinate:
    left = files[dfNum][21:-4] + ' L'
    # Column for right coordinate:
    right = files[dfNum][21:-4] + ' R'
    # Column for SNR:
    snr = files[dfNum][21:-4] + ' SNR'
    # Return new DataFrame with columns added.
    return pd.concat([df, pd.DataFrame(columns=[left, right, snr])], sort=False)


# DataFrame for merged data:
merged = pd.DataFrame(columns=['sRNA', 'sRNA peak'])
# Add columns for rdmWT condition.
merged = add_columns(merged, 0)


# Function for adding sRNA data to cells, works for all sRNAs and conditions:
def add_sRNA(merged, i, j, dfNum, df):
    # Add left coordinate.
    merged.iloc[i, dfNum*3+2] = int(df.iloc[j, 2])
    # Add right coordinate.
    merged.iloc[i, dfNum*3+3] = int(df.iloc[j, 3])
    # Add SNR.
    if str(df.iloc[j, 4])[5:-3] != '':
        snr_temp = re.findall('\d+\.\d+', df.iloc[j, 4])
        merged.iloc[i, dfNum*3+4] = float(snr_temp[0])
        
    return merged

# Function for inserting an empty row into merged DataFrame for new peak:
def insert_row(row_number, df):
    # Split DataFrame in two so new row can be added in between.
    df_1 = df.iloc[0:row_number]
    df_2 = df.iloc[row_number:]
    # Create new empty row.
    empty_row = pd.DataFrame([[np.nan] * len(df.columns)], columns=df.columns)
    # Concatenate DataFrames with empty row in between.
    newDF = pd.concat([df_1, empty_row, df_2], ignore_index = True)
    newDF.reset_index()
    return newDF

# Function for inserting row with sRNA data added:
def insert_row_final(merged, newDF, i, j, k):
    # Add empty row at the end.
    merged = insert_row(i+1, merged)
    # Add rdmStat peak data to new empty row.
    merged = add_sRNA(merged, i+1, j, k, newDF)
    # Add sRNA name.
    merged.loc[i+1, 'sRNA'] = newDF.iloc[j, 0]
    return merged


# =============================================================================
# Add rdmWT data to merged DataFrame. 
# =============================================================================

# Create DataFrame for rdmWT condition.
rdmWT_DF = pd.read_csv(files[0])

# Loop through rdmWT DataFrame and add data to merged DataFrame.
for i, row in rdmWT_DF.iterrows():

    # Add sRNA names.
    merged.loc[i, 'sRNA'] = rdmWT_DF.iloc[i, 0]
    # Add remaining data (left coordinate, right coordinate, SNR).
    merged = add_sRNA(merged, i, i, 0, rdmWT_DF)


# =============================================================================
# Loop through remaining files for different conditions. Add data to merged
# DataFrame. Peaks are compared to rdmWT peak.
# =============================================================================
    
k = 1
while k < len(files):
    
    # Create DataFrame for given condition.
    condition_DF = pd.read_csv(files[k])
    # Add columns for given condition.
    merged = add_columns(merged, k)
    
    # Loop through peaks and add data to merged DataFrame.
    for j, row in condition_DF.iterrows():
        
        # Create set of positions for new peak.
        newPos = set(range(condition_DF.iloc[j, 2], condition_DF.iloc[j, 3] + 1))
        
        # Loop through existing peaks to compare.
        i = 0
        while i < len(merged.index):
            
            # If sRNAs have the same name, compare peak positions.
            if merged.loc[i, 'sRNA'] == condition_DF.iloc[j, 0]:
                
                # Boolean to keep track of whether the peak was added or not:
                # Default to False because peak not added yet.
                added = False
                
                # If there is a rdmWT peak:
                if math.isnan(merged.iloc[i, 2]) == False:
                    # Create set of positions for rdmWT peak.
                    WTPos = set(range(int(merged.iloc[i, 2]), int(merged.iloc[i, 3]) + 1))
                    # List of matching nt between positions for new and rdmWT peaks.
                    match = newPos & WTPos
                    
                    # If the number of matching nucleotides is greater than n, the peak is the same.
                    # Add new peak data to existing row.
                    if len(match) >= n:
                        # Since peak was added, added = True.
                        merged = add_sRNA(merged, i, j, k, condition_DF)
                        added = True
                        # Break out of while loop since peak already added.
                        break
                    
                # If on the last existing peak and new peak not added yet
                if i == len(merged.index) - 1 and added == False:
                    merged = insert_row_final(merged, condition_DF, i, j, k)
                
                    # Increase index by 1 so while loop doesn't go over new row.
                    i += 1
                    
                # If next rdmWT peak is a different sRNA and peak not added yet.
                elif merged.loc[i+1, 'sRNA'] != condition_DF.iloc[j, 0] and added == False:
                    merged = insert_row_final(merged, condition_DF, i, j, k)
                        
                    # Increase index by 1 so while loop doesn't go over new row.
                    i += 1
                
            i += 1
  
    k += 1        
    
    
# =============================================================================
# Add leftmost coordinate, rightmost coordinate, and maximum SNR for peaks.
# =============================================================================

# Add new columns for merged data.
# Merged L is minimum left coordinate of the peak across different conditions.
# Merged R is maximum right coordinate of the peak across different conditions.
merged = pd.concat([merged, pd.DataFrame(columns=['Merged L', 'Merged R'])], sort=False)

# Loop through all peaks in merged DataFrame.
for i, row in merged.iterrows():
    
    # Lists to keep track of each peak's data under different condition.
    left = []
    right = []
    snr = []
    
    # Loop through each condition.
    dfNum = 0
    while dfNum < len(files):
        # If there is a peak for that condition (not empty):
        if math.isnan(merged.iloc[i, dfNum*3+2]) == False:
            # Add left coord, right coord, & SNR to respective lists.
            left.append(merged.iloc[i, dfNum*3+2])
            right.append(merged.iloc[i, dfNum*3+3])
            if math.isnan(merged.iloc[i, dfNum*3+4]) == False:
                snr.append(merged.iloc[i, dfNum*3+4])
        dfNum += 1

    # If there are SNRs, find maximum SNR. Add to merged DataFrame.
    if snr:
        merged.loc[i, 'Max SNR'] = float(max(snr))
    # Find leftmost coordinate. Add to merged DataFrame.
    merged.iloc[i, -3] = int(min(left))
    # Find rightmost coordinate. Add to merged DataFrame.
    merged.iloc[i, -2] = int(max(right))


# =============================================================================
# Sort peaks by sRNA and position.
# =============================================================================
    
# Order by sRNA name, keep current order.
order = 0
# New column to keep track of name order.
merged.loc[0, 'Order by name'] = order
i = 1
# Loop through all peaks in merged DataFrame.
while i < len(merged.index):
    
    # If new sRNA, order increases by 1.
    if merged.loc[i, 'sRNA'] != merged.loc[i-1, 'sRNA']:
        order += 1
    merged.loc[i, 'Order by name'] = order
    i += 1

# Loop for order by position.
# Loop through all peaks in merged DataFrame.
for i, row in merged.iterrows():
    
    # Loop through each condition.
    dfNum = 0
    while dfNum < len(files):
        # If there is a peak for that condition (not empty):
        if math.isnan(merged.iloc[i, dfNum*3+2]) == False:
            # Add left coordinate to new column 'Order by position'.
            merged.loc[i, 'Order by position'] = merged.iloc[i, dfNum*3+2]
            break
        dfNum += 1

# Sort first by sRNA name, keep current order.
# Sort second by position of left coordinate.
merged = merged.sort_values(by=['Order by name', 'Order by position'])
# Delete additional columns created.
merged.drop(['Order by name', 'Order by position'], axis=1, inplace=True)
# Reset indexes so they aren't confusing.
# drop = True so new 'index' column is not created.
# inplace = True to modify existing DataFrame.f
merged.reset_index(drop = True, inplace = True)


# =============================================================================
# Number peaks for each sRNA.
# =============================================================================

# Number peaks for each sRNA.
num = 1
# Loop through all peaks in merged DataFrame.
for i, row in merged.iterrows():
    # Start the first peak at 1.
    if i != 0:
        # If new sRNA, reset numbering.
        if merged.loc[i, 'sRNA'] != merged.loc[i-1, 'sRNA']:
            num = 1
    merged.loc[i, 'sRNA peak'] = merged.loc[i, 'sRNA'] + "_" + str(num)
    # Add 1 for next peak.
    num += 1


# Replace all numpy NaN values with string "na".
merged = merged.fillna("na")


# =============================================================================
# Exclude peaks below minimum SNR from analysis.
# =============================================================================

# Copy of merged DataFrame that excludes peaks below minimum SNR from analysis:
temp = merged.copy()

# Loop through all peaks in temp DataFrame.
for i, row in temp.iterrows():
    
    # Loop through each condition.
    dfNum = 0
    while dfNum < len(files):
        
        # If there is a peak for that condition (not empty).
        if temp.iloc[i, dfNum*3+4] != "na":
            
            # If the SNR is below the minimum SNR, count peak as non-existent.
            if temp.iloc[i, dfNum*3+4] < min_SNR:
                temp.iloc[i, dfNum*3+2] = "na"
                temp.iloc[i, dfNum*3+3] = "na"
                temp.iloc[i, dfNum*3+4] = "na"
                
        dfNum += 1
    

# If peak is nonexistent in any condition, peak is differential.
# Unless peak is nonexistent for all conditions (below minimum SNR).
for i, row in temp.iterrows():
    if temp.iloc[i].str.contains("na").any() and temp.iloc[i].str.contains("na").sum() < len(files)*3:
        
        # 1 means differential peak. 0 means nondifferential peak.
        merged.loc[i, 'Differential peak?'] = 1

# Fill all other peaks with 0.
merged = merged.fillna(0)


# =============================================================================
# Find and compare maximum SNR ratios and differences to determine if peak is
# differential.
# =============================================================================

# Create list of condition combinations to compare each 2 conditions separately.
file_nums = list(range(len(files)))
file_subsets = list(itertools.combinations(file_nums, 2))

# Loop through all peaks in temp DataFrame.
for i, row in temp.iterrows():

    # Lists for all SNR ratios and differences to find maximum later:
    snr_ratios = []
    snr_diffs = []
    
    # For each combination, find the maximum SNR ratio and SNR difference.
    # A big enough SNR ratio means the peak is differential.
    # The SNR ratio is the absolute value of log (snr1/snr2). Doesn't matter which is being divided.
    for subset in file_subsets:

        # If both SNRs are not empty, add SNR ratio and difference to lists.
        if temp.iloc[i, subset[0]*3+4] != "na" and temp.iloc[i, subset[1]*3+4] != "na":
            snr1 = temp.iloc[i, subset[0]*3+4]
            snr2 = temp.iloc[i, subset[1]*3+4]
            snr_diffs.append((abs(snr1 - snr2)))
            if merged.iloc[i, subset[0]*3+4] != 0 and merged.iloc[i, subset[1]*3+4] != 0:
                snr_ratios.append((abs(math.log10(snr1/snr2))))

    # If SNR ratio is greater than log(2), peak counts as differential.
    # This means one SNR is at least twice as large as the other.
    if any(y >= math.log10(2) for y in snr_ratios):
        merged.loc[i, 'Differential peak?'] = 1
    
    # If peak is nonexistent under any condition, max SNR difference is max SNR.
    # Max SNR ratio is empty.
    if merged.iloc[i].str.contains("na").any():
        merged.loc[i, 'Max SNR difference'] = merged.loc[i, 'Max SNR']
    
    # If peak is existent for all conditions:
    else:
        # Find maximum SNR ratio and add to merged DataFrame.
        if snr_ratios:
            merged.loc[i, 'Max SNR ratio'] = max(snr_ratios)
            
        # Find maximum SNR difference and add to merged DataFrame.
        if snr_diffs:
            merged.loc[i, 'Max SNR difference'] = max(snr_diffs)
    
# Fill all other cells with "na".
merged = merged.fillna("na")


# =============================================================================
# Add directions of all peaks (forward or reverse).
# =============================================================================

# Read file with directions of all sRNAs
sRNAs_DF = pd.read_excel('sRNAs_list_from-lib_050919.xlsx')

# Loop through each peak in merged DataFrame.
for i, row in merged.iterrows():
    
    # For each peak, loop through sRNAs DataFrame to find match
    for j, row in sRNAs_DF.iterrows():
        
        # Add direction for each sRNA from sRNAs DataFrame to merged DataFrame
        if sRNAs_DF.loc[j, 'sRNA'].lower().strip() in merged.loc[i, 'sRNA'].lower().strip():
            merged.loc[i, 'Direction'] = sRNAs_DF.loc[j, 'Direction']
            break


# =============================================================================
# Extract DNA sequences of peak regions.
# =============================================================================

# Open file that contains the K-12 genome for extracting DNA sequences.
with open("GCF_000005845.2_ASM584v2_genomic (1).fna") as f:
# Create list of all lines in file.
    content = f.readlines()
content = [x.strip() for x in content]

# Loop through strings in content list and add to one long string.
# This long string is the entire genome.
# Start at i = 1 since genome starts on second line.
i = 1
# Genome starts with "0" so indices will match nucleotide coordinates.
genome = "0"
while i < len(content):
    genome += content[i]
    i += 1

# Loop through each peak in merged DataFrame.
for i, row in merged.iterrows():
        
    # Get leftmost and rightmost coordinates.
    left_coord = merged.loc[i, 'Merged L']
    right_coord = merged.loc[i, 'Merged R']
    
    # Extract forward sequence from genome string using the left and right coordinates.
    seq_F = genome[left_coord:right_coord+1]
    
    # If the direction is forward, add forward sequence to merged DataFrame.
    if merged.loc[i, 'Direction'] == 'F':
        merged.loc[i, 'Sequence'] = seq_F
    
    # If the direction is reverse, take reverse complement of forward sequence.
    else:
        seq_R = ""
        for nt in seq_F:
            if nt == 'A':
                seq_R = 'T' + seq_R
            elif nt == 'T':
                seq_R = 'A' + seq_R
            elif nt == 'G':
                seq_R = 'C' + seq_R
            elif nt == 'C':
                seq_R = 'G' + seq_R
        # Add reverse sequence to merged DataFrame.
        merged.loc[i, 'Sequence'] = seq_R


# These columns aren't needed for final output.
merged.drop(['Max SNR', 'Direction'], axis=1, inplace=True)


# Export merged DataFrame to new Excel file
# Or rewrite existing Excel file with new merged DataFrame
merged.to_excel("merged_peaks.xlsx", sheet_name='Merged Peaks', index=False)


# =============================================================================
# Create FASTA file containing all differential peaks.
# =============================================================================

# Create new FASTA file for writing (or overwrite existing file).
fasta_file = open("differential_peaks.fasta", "w")

# Loop through each peak in merged DataFrame.
for i, row in merged.iterrows():
    
    # If the peak is differential:
    if merged.loc[i, 'Differential peak?'] == 1:
        # Write the region's name on first line.
        fasta_file.write(">" + merged.loc[i, 'sRNA peak'] + "\n")
        # Write the sequence on the second line.
        if i != len(merged.index) - 1:
            fasta_file.write(merged.loc[i, 'Sequence'] + "\n")
        else:
            fasta_file.write(merged.loc[i, 'Sequence'])
        
fasta_file.close()


# =============================================================================
# Create files containing commands to run MEME (fimo) on each differential peak.
# =============================================================================

# Create new FASTA file for writing (or overwrite existing file).
txtFile = open("run_fimo_dp_prod.txt", "w")
txtFile_no_redund = open("run_fimo_dp_prod_no_redund.txt", "w")

# Loop through each peak in merged DataFrame.
for i, row in merged.iterrows():
    
    # If the peak is differential:
    if merged.loc[i, 'Differential peak?'] == 1:
        
        if i != len(merged.index) - 1:
            
            # Write fimo command to 'all' text file.
            txtFile.write("fimo --oc . --verbosity 1 --thresh 1.0E-4 " +
                          "dpinteract_prod2.meme " + merged.loc[i, 'sRNA peak'] +
                          ".txt\n")
            # Write fimo command to 'not redundant' text file.
            txtFile_no_redund.write("fimo --oc . --verbosity 1 --thresh 1.0E-4 " +
                                    "dpinteract_prod2_no_redund.meme " +
                                    merged.loc[i, 'sRNA peak'] + ".txt\n")
            
        # Don't add blank line for last row.
        else:
            txtFile.write("fimo --oc . --verbosity 1 --thresh 1.0E-4 " +
                          "dpinteract_prod2.meme " + merged.loc[i, 'sRNA peak'] +
                          ".txt")
            txtFile_no_redund.write("fimo --oc . --verbosity 1 --thresh 1.0E-4 " +
                                    "dpinteract_prod2_no_redund.meme " +
                                    merged.loc[i, 'sRNA peak'] + ".txt")

txtFile.close()
txtFile_no_redund.close()
