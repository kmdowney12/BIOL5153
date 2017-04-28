#!/usr/bin/env python3

#This script reads and parses collection data for a new ctenophore species.

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

# All the GC junk for CDS.
CDS_len = len(cds)
cds_g = cds.count('G')
cds_c = cds.count('C')
gc_cds = round((cds_g + cds_c) / CDS_len *100,2)
cds_perc = round((CDS_len / genome_len) *100,1)

# All the GC junk for introns.
intron_len = len(intron)
intron_g = intron.count('G')
intron_c = intron.count('C')
gc_intron = round((intron_g + intron_c) / intron_len *100,2)
intron_perc = round((intron_len / genome_len) *100,1)

# All the GC junk for tRNA.
tRNA_len = len(tRNA)
tRNA_g = tRNA.count('G')
tRNA_c = tRNA.count('C')
gc_tRNA = round((tRNA_g + tRNA_c) / tRNA_len *100,2)
tRNA_perc = round((tRNA_len / genome_len) *100,1)

# All the GC junk for rRNA.
rRNA_len = len(rRNA)
rRNA_g = rRNA.count('G')
rRNA_c = rRNA.count('C')
gc_rRNA = round((rRNA_g + rRNA_c) / rRNA_len *100,2)
rRNA_perc = round((rRNA_len / genome_len) *100,1)

# All the GC junk for repeat.
repeat_len = len(repeat)
repeat_g = repeat.count('G')
repeat_c = repeat.count('C')
gc_repeat = round((repeat_g + repeat_c) / repeat_len *100,2)
repeat_perc = round((repeat_len / genome_len) *100,1)

# All the GC junk for misc.
misc_len = len(misc)
misc_g = misc.count('G')
misc_c = misc.count('C')
gc_misc = round((misc_g + misc_c) / misc_len *100,2)
misc_perc = round((misc_len / genome_len) *100,1)

#Close the gff file.
gff_in.close()

print("exon              ", CDS_len, "(" + str(cds_perc) + "%)",'\t', gc_cds)
print("intron            ", intron_len, "(" + str(intron_perc) + "%)",'\t', gc_intron)
<<<<<<< HEAD
print("misc_feature      ", misc_len, "(" + str(misc_perc) + "%)",'\t', gc_misc)
=======
print("misc_feature       ", misc_len, "(" + str(misc_perc) + "%)",'\t', gc_misc)
>>>>>>> 090cab1986a7c4414c0d84f20a20542d13397bb5
print("repeat_region     ", repeat_len, "(" + str(repeat_perc) + "%)",'\t', gc_repeat)
print("rRNA              ", rRNA_len," " "(" + str(rRNA_perc) + "%)",'\t', gc_rRNA)
print("tRNA              ", tRNA_len," " "(" + str(tRNA_perc) + "%)",'\t', gc_tRNA)
