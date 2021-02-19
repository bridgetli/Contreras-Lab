"""
update_sequences.py
    6/23/2020
    This program extracts sequences based on updated sRNA coordinates and makes a
    FASTA file for each sRNA. It also updates the starts of ends of each IPOD region,
    then creates IntaRNA commands for each region.
"""

# Use Anaconda, a distribution of Python for data science
# Run on Spyder

# Import Pandas package, used for data analysis (Excel).
import pandas as pd
from openpyxl import load_workbook
# Import os package for making directories and adding text files to them.
import os


# =============================================================================
# Extract sRNA sequences using updated coordinates and directions. Make a FASTA
# file for each sRNA.
# =============================================================================

# Read in Excel file given file name.
file = pd.ExcelFile('final_INTERFACE_analysis.xlsx')

# Load workbook so you can add new sheets w/o overwriting existing ones
# Use the specific name of your Excel file
book = load_workbook('final_INTERFACE_analysis.xlsx')
writer = pd.ExcelWriter('final_INTERFACE_analysis.xlsx', engine='openpyxl') 
writer.book = book
writer.sheets = dict((ws.title, ws) for ws in book.worksheets)

# Parse IPOD sheet into DataFrame.
ipod_df = file.parse('IPOD')


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


# Create directory for FASTA files if it doesn't already exist.
if not os.path.exists('FASTA Files'):
    os.mkdir('FASTA Files')


# Loop through each sRNA in IPOD DataFrame.
for i, row in ipod_df.iterrows():
    
    # Create a new text file in directory for each sRNA.
    txtFile = open('FASTA Files\\' + ipod_df.loc[i, 'sRNA'] + '.fasta', "w")
    
    # Write the sRNA's name on first line.
    txtFile.write('>' + ipod_df.loc[i, 'sRNA'] + "\n")
        
    # Get left and right coordinates of text file.
    left_coord = ipod_df.loc[i, 'L_Gen_Coord']
    right_coord = ipod_df.loc[i, 'R_Gen_Coord']
    
    # Extract forward sequence from genome string using the left and right coordinates.
    seq_F = genome[left_coord:right_coord+1]
    
    # If the direction is forward, write forward sequence on  the second line of
    # the text file.
    if ipod_df.loc[i, 'Direction'] == 'F':
        txtFile.write(seq_F)
    
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
                
        # Write reverse sequence on the second line of the text file.
        txtFile.write(seq_R)
    
    # Close text file so it can't be written anymore.
    txtFile.close()
    
    
# =============================================================================
# Update start and end coordinates of seed regions.
# =============================================================================
    
# Parse seed region sheet into DataFrame.
seed_regions_df = file.parse('pos_ctrl_seeds')

# Function for updating starts and ends of seed regions.
def update_seed_regions(column_num):
    
    # If 0 mismatch between old and new coordinates, seed region coordinates stay
    # the same

    if ipod_df.iloc[i, column_num] == 0:
        
        seed_regions_df.loc[j, 'accepted_region_start'] = \
        seed_regions_df.loc[j, 'region_start']
        
        seed_regions_df.loc[j, 'accepted_region_end'] = \
        seed_regions_df.loc[j, 'region_end']
    
    # If there is a mismatch, add length of mismatch to old seed region coordinates
    # to get new coordinates. Length of region stays the same.
    else:
        
        seed_regions_df.loc[j, 'accepted_region_start'] = \
        seed_regions_df.loc[j, 'region_start'] + \
        ipod_df.iloc[i, column_num]
        
        seed_regions_df.loc[j, 'accepted_region_end'] = \
        seed_regions_df.loc[j, 'region_end'] + \
        ipod_df.iloc[i, column_num]
        
    return(seed_regions_df)


# Loop through each seed region in seed region DataFrame.
for j, row in seed_regions_df.iterrows():
    
    # For each seed region, loop through IPOD DataFrame to find match.
    for i, row in ipod_df.iterrows():

        # If sRNA names match:        
        if seed_regions_df.loc[j, 'IPOD_name'].lower() == ipod_df.loc[i, 'sRNA'].lower():
            
            # Only left mismatch matters for forward sequences. Use left mismatch
            # to update seed region coordinates.
            if ipod_df.loc[i, 'Direction'] == 'F':
                seed_regions_df = update_seed_regions(-2)
                
            # Only right mismatch matters for forward sequences. Use right mismatch
            # to update seed region coordinates.
            else:
                seed_regions_df = update_seed_regions(-1)
            
            break


# Rewrite existing sheet with the DataFrame with updated seed region coordinates.
# Save sheet.
seed_regions_df.to_excel(writer, 'pos_ctrl_seeds', index=False)
writer.save()


# =============================================================================
# Create IntaRNA commands for each sRNA.
# =============================================================================

# Create a new text file that contains commands.
txtFile = open('IntaRNA_commands.txt', 'w')

# Function for formatting and concatenating IntaRNA command:
def write_IntaRNA_command(sRNA_name):
    
    # Make first letter of sRNA name lowercase so it matches FASTA file name.
    sRNA_name = sRNA_name[0].lower() + sRNA_name[1:]
    
    # Write IntaRNA command to text file.
    txtFile.write("IntaRNA -q subdirectory/" + sRNA_name + ".fasta -t " + 
                  "annotated_5_utr_plus_100.fasta --outMode=C --outCsvCols " +
                  "'id1,start1,end1,id2,start2,end2,E' > $SCRATCH/" + sRNA_name +
                  ".csv \nIntaRNA -q subdirectory/" + sRNA_name + ".fasta -t " +
                  "annotated_5_utr_plus_100.fasta " + 
                  " --outMode=C --outCsvCols 'id1,start1,end1,id2,start2,end2,E' > $SCRATCH/" + 
                  sRNA_name + ".csv")


# Sort DataFrame by IPOD_name first. This is so seed regions of the sRNA will be
# adjacent.
# Sort DataFrame by coordinates second. This is so end and start coordinates can
# be compared to determine whether two regions overlap.
seed_regions_df = seed_regions_df.sort_values(by=['IPOD_name', 'accepted_region_start', 'accepted_region_end'])

j = 1
# Loop through all peaks in seed regions DataFrame.
while j < len(seed_regions_df.index):
    
    current_name = seed_regions_df.loc[j, 'IPOD_name']
    previous_name = seed_regions_df.loc[j-1, 'IPOD_name']
    
    # If current name is the same as the previous name, we are on the same sRNA.
    if current_name != previous_name:
        
        # Write IntaRNA command for previous sRNA, since we're finished with it.
        write_IntaRNA_command(previous_name)
        txtFile.write("\n")
        
    j+=1

# Write IntaRNA command for final sRNA, since the for loop doesn't do that.
write_IntaRNA_command(seed_regions_df.loc[j-1, 'IPOD_name'])

txtFile.close()
