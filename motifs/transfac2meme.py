"""
transfac2meme.py
    6/21/2020
    This program creates text files that contain transfac2meme commands for MG1655
    motifs from the Prodoric2 database. These commands will convert the TRANSFAC
    files of each motif into MEME motif files suitable for use with the MEME Suite.
"""

# Use Anaconda, a distribution of Python for data science
# Run on Spyder

# Import Pandas package, used for data analysis (Excel)
import pandas as pd

# Read in CSV Excel file containing Prodoric2 motifs
prodoric2_df = pd.read_csv('prodoric2.csv')

# Writes transfac2meme command for given motif
def write_command(i, txt_file):
        # The TRANSFAC file follows this format:
        # prodoric_<accession number>.txt
        file_name = 'prodoric_' + prodoric2_df.iloc[i, 0].lower() + '.txt'
        
        # Write transfac2meme command to text file.
        txt_file.write('transfac2meme [options] ' + file_name + '\n')

# =============================================================================
# ALL
# All MG1655 motifs in Prodoric2 database
# =============================================================================

# Create a new text file that contains commands for ALL Prodoric2 motifs.
txt_all = open("transfac2meme_all.txt", "w")

# Loop through all motifs in Prodoric2 DataFrame.
# Create a transfac2meme command for each motif.
for i, row in prodoric2_df.iterrows():
    
    # Write transfac2meme command for given motif.
    write_command(i, txt_all)

txt_all.close()


# =============================================================================
# NOT REDUNDANT
# MG1655 motifs found in Prodoric2 database but not in DPInteract database
# =============================================================================

# Read in XLSX Excel file containing DPInteract motifs.
dpinteract_df = pd.read_excel('dpinteract.xlsx')

# Create a new text file that contains commands.
txt_not_redundant = open("transfac2meme_not_redundant.txt", "w")

# Loop through all motifs in Prodoric2 DataFrame.
for i, row in prodoric2_df.iterrows():
    redundant = False
    
    # Loop through all motifs in DPInteract DataFrame
    for j, row in dpinteract_df.iterrows():
        
        # If Prodoric2 motif contains DPInteract motif, Prodoric2 motif is redundant.
        # Make all strings in list lowercase, so they will be compared correctly.
        if dpinteract_df.loc[j, 'motif'].lower() in prodoric2_df.iloc[i, 1].lower():
            redundant = True
            break
    
    # If Prodoric2 motif doesn't match any DPInteract motifs, write transfac2meme
    # command for given motif.
    if redundant == False:
        write_command(i, txt_not_redundant)

txt_not_redundant.close()
