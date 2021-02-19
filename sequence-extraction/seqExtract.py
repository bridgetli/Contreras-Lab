''' Bridget Li
    10/10/18
    seqExtract.py
    This program extracts the desired DNA sequence given an entire genome to determine
    where the sRNA interacts with the DNA
'''

# Imports Pandas package, used for data analysis (Excel)
import pandas as pd

# Reads in Excel file given file name
file = pd.ExcelFile('example_extraction.xlsx')
file.sheet_names

# Parses sheet into DataFrame
extract = file.parse('Input')

# Creates new seqColumn for the sequences we are extracting
seqColumn = 'Seq'
extract[seqColumn] = 0

# Loops through each row in extract DataFrame
for i in range(extract.shape[0]):    #for each line
    
    # Reads in genome from file
    genomeFile = open("ecoligenome_MKM.txt","r")

    # Sets orientation as 'fwd' or 'rev' based on what is in the 'Orientation' column
    orientation = extract.loc[i]['Orientation']
    # Sets start and length of sequence
    UTR_start = extract.loc[i]['UTR_Start']
    length = extract.loc[i]['Transcript_Length']
        
    genome = ""
    position = 0
    location = 0
    linecount = 0
    seq = ""
    revcomp = ""
    revseq = ""
        
    for raw_line in genomeFile:
        line = raw_line.rstrip("\r\n")
    
        if linecount == 0:
            name = line.split(' ',1)[0] # Save the name of the genome
            for char in line:
                # If position > (len(name)+1): 
                if position > 13:
                    genome = genome + char
                    location +=1
                else:
                    genome = genome
                position+=1
        linecount += 1
    
        if linecount > 0:
            for char in line:
                genome = genome + char
                location +=1
        linecount += 1
    file.close()
    
    # Complementary bases for reverse sequences
    tcomp = "A"
    ccomp = "G"
    acomp = "T"
    gcomp = "C"

    # If reverse orientation, flips sequence and replaces with complementary base
    if orientation == "rev":
        seq = genome[(UTR_start-length):UTR_start]
        revseq = seq[::-1]
        for nuc in revseq:
            if nuc == "A":
                revcomp = revcomp + acomp
            if nuc == "T":
                revcomp= revcomp + tcomp
            if nuc == "C":
                revcomp = revcomp + ccomp
            if nuc == "G":
                revcomp= revcomp + gcomp

        # Adds revcomp to 'Seq' column
        extract.loc[i, seqColumn] = revcomp
            
    else:
        # Adds fwd sequence to 'Seq' column
        extract.loc[i, seqColumn] = genome[UTR_start:(UTR_start+length)]
                    
# Exports DataFrame back to original Excel file, adding a new 'Output' sheet
from openpyxl import load_workbook

book = load_workbook('example_extraction.xlsx')
writer = pd.ExcelWriter('example_extraction.xlsx', engine='openpyxl') 
writer.book = book
writer.sheets = dict((ws.title, ws) for ws in book.worksheets)

extract.to_excel(writer, "Output")

writer.save()
