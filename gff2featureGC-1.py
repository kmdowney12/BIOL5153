#!/usr/bin/env python3

#A function to clean up DNA files
def clean_seq(input_seq):
  clean = input_seq.upper()
  clean = clean.replace("N", "")
  return clean

#Function to calculate nucleotide percentage.
def nuc_freq(sequence, b1, b2):
#Length of the sequence.
  length = len(sequence)
#Count the number of nucleotides in sequence.
  b1_count = sequence.count(b1)
  b2_count = sequence.count(b2)
  total_bases = b1_count + b2_count
  b_perc = round(((total_bases/length)) *100, 2)
#Return the frequency and length.
  return (b_perc)

#Function for creating the complementary strands.
def complement(seq):
  complement = (seq.replace("A", "t"))
  complement = (complement.replace("C", "g"))
  complement = (complement.replace("T", "a"))
  complement = (complement.replace("G", "c"))
  return complement.upper()

#Function to isolate just the name of the gene we want.
def snip_exon(gff_line):
	break_1 = line.split('Gene')
	break_2 = break_1[1].split(';')
	exon = break_2[0].strip()
	return exon

# Our key is the feature type which includes all the sequences for that particular feature.
feature_seq = {}
# This is just the exon sequences for each feature (*cough* and maybe a little tRNA *cough*).
exon_seq = {}
# All the exons in the gene.
gene_seq = {}

#Load the system module.
import sys
import collections

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
	# gather full_seq on feature_types
	full_seq = line.split("	")
	feature = full_seq[2].strip()
	start_here = int(full_seq[3]) - 1
	strand = clean_seq(genome[start_here:int(full_seq[4])])
	if feature in feature_seq:
		feature_seq[feature] = feature_seq[feature] + strand
	else:
		feature_seq[feature] = strand
	# create dictionary of exons and their sequences
	full_seq = line.split("\t")
	if full_seq[2].strip() == 'CDS':
		exon = snip_exon(line)
		start_here = int(full_seq[3]) - 1
		strand = clean_seq(genome[start_here:int(full_seq[4])])
		if full_seq[6].strip() == '-':
			complement_strand = complement(strand)
			exon_seq[exon] = complement_strand
		else:
			exon_seq[exon] = strand

#Exon dictionary
exon_dict = collections.OrderedDict(sorted(exon_seq.items()))

#Build the exon
for exon, sequence in exon_dict.items():
	exon_bits = exon.split(' ')
	if '-' in exon_bits[0]:
		exon_bits = exon_bits[0].split('-')
	if exon_bits[0] in gene_seq:
		gene_seq[exon_bits[0]] += sequence
	else:
		gene_seq[exon_bits[0]] = sequence

# For-loop to get the features.
for type, sequence in feature_seq.items():
	GC_perc = nuc_freq(sequence, "C", "G")
	type_perc = (len(sequence)/genome_len) * 100
	output = "(" + str(round(type_perc, 1)) + "%)"
	info_junk = "{0:15}{1:7} {2:7}{3:15}"
	print(info_junk.format(type, len(sequence), output, GC_perc))

print ("\n")

#Make the sequences.
for gene, sequence in gene_seq.items():
	print(">" + gene + "\n" + sequence + "\n")
