#!/usr/bin/env python3

#A function to clean up DNA files
def clean_seq(input_seq):
  clean = input_seq.upper()
  clean = clean.replace("N", "")
  return clean

#Function to calculate nucleotide percentage.
def nuc_freq(sequence, base, sig_figs=2):
#Length of the sequence.
  length = len(sequence)
#Count the number of nucleotides in sequence.
  nuc_count = sequence.count(base)
  nuc_freq = nuc_count/length
#Return the frequency and length.
  return (length, round(nuc_freq, sig_figs))

#Load the system module.
import sys

print('This is the name of the script: ', sys.argv[0])
print('Number of arguments: ', len(sys.argv))
print('The arguments are: ', str(sys.argv))

#Set the name of the input file to be opened.
fsa_file = sys.argv[1];
gff_file = sys.argv[2];

#Open the file and read it in.
gff_in = open(gff_file, 'r')
fsa_in = open(fsa_file, 'r')

#declare variable that will hold the genome sequence
genome = ''

#Count the number of lines and make that a variable.
first_line = 0

#Read the file line by line to stdout.
for line in fsa_in:
#    print(str(first_line) + ": " + line)
    line = line.rstrip("\n")
    if first_line > 0:
        genome = genome + line
#counts the lines incrementally
    first_line += 1
genome_len = len(genome)

#Close the file.
fsa_in.close()

cds = ''
tRNA = ''
rRNA = ''
intron = ''
repeat = ''
misc = ''

for line in gff_in:
    types = line.split('type ')
    fields = line.split('\t')
    type = fields[2]
    start = int(fields[3])
    stop = int(fields[4])
    fragment = genome[start-1:stop]
    fragment = clean_seq(fragment)

    if type == 'CDS':
        cds += fragment
    elif type == 'intron':
        intron += fragment
    elif type == 'tRNA':
        tRNA += fragment
    elif type == 'rRNA':
        rRNA += fragment
    elif type == 'misc_feature':
        misc += fragment
    elif type == 'repeat_region':
        repeat += fragment

#Close the gff file.
gff_in.close()

#Use our shiny new function instead of all that math we worked so hard on.
types=["CDS", "rRNA", "tRNA", "Intron", "Misc","Repeats"]
i=0

for feature_type in [cds, rRNA, tRNA, intron, misc, repeat]:
  for nucleotide in ['A','T','G','C']:
    (feature_length, feature_comp) = nuc_freq(sequence=feature_type, base=nucleotide)
    print(types[i] + "\t" + "length= " + str(feature_length) + "\t" + "and the percentage of " + nucleotide + "'s in this feature is " + str(feature_comp*100) + "%.")
  i = i + 1
