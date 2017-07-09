#! /usr/bin/python
# generate_KASPAR_markers.py
# description:	
#	1) Import IUPAC sequences

import commands

import optparse
from optparse import OptionParser 

import re

import sets

import string

# global variables
cmp_DNA = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'M':'K', 'R':'Y', 'S':'S', 'W':'W', 'Y':'R', 'K':'M', 'V':'B', 'H':'D', 'D':'H', 'B':'V', 'X':'X', 'N':'N'}
IUPAC = {'AG':'R', 'GA':'R', 'CT':'Y', 'TC':'Y', 'CA':'M', 'AC':'M', 'TG':'K', 'GT':'K', 'TA':'W', 'AT':'W', 'CG':'S', 'GC':'S'}
DNA = ['A', 'C', 'G', 'T']
IUPAC_dinucleotides = ['R', 'Y', 'M', 'K', 'W', 'S']
IUPAC_SNP = {'R':['A', 'G'], 'Y':['C', 'T'], 'M':['A', 'C'], 'K':['G', 'T'], 'W':['A', 'T'], 'S':['C', 'G']}
VIC = 'GAAGGTCGGAGTCAACGGATT'
FAM = 'GAAGGTGACCAAGTTCATGCT'

# FUNCTIONS
# reverse
# input : DNA sequence
# output : reverse of said DNA sequence
def reverse(orig_DNA):
	rev_DNA = ''
	for index in range(len(orig_DNA)):
		rev_DNA += orig_DNA[len(orig_DNA)-index-1]
	return rev_DNA

# reverse complement
# input : DNA sequence
# output : reverse complement of said DNA sequence
def reverse_complement(orig_DNA):
	rev_comp_DNA = ''
	for index in range(len(orig_DNA)):
		rev_comp_DNA += cmp_DNA[orig_DNA[len(orig_DNA)-index-1]]
	return rev_comp_DNA

# primer3 output
def primer3_output(output_stream, markerID, template, right_or_force):
	output_stream.write('SEQUENCE_ID=' + markerID + '\n')
	output_stream.write('SEQUENCE_TEMPLATE=' + template + '\n')
	
	if right_or_force == 'right':
		output_stream.write('PRIMER_PICK_LEFT_PRIMER=0' + '\n')
	elif right_or_force == 'force':
		output_stream.write('PRIMER_PICK_LEFT_PRIMER=1' + '\n')
		output_stream.write('SEQUENCE_FORCE_LEFT_START=0' + '\n')
		output_stream.write('SEQUENCE_FORCE_LEFT_END=19' + '\n')

	output_stream.write('PRIMER_PICK_INTERNAL_OLIGO=0' + '\n')
	output_stream.write('PRIMER_PICK_RIGHT_PRIMER=1' + '\n')
	output_stream.write('PRIMER_OPT_SIZE=18' + '\n')
	output_stream.write('PRIMER_MIN_SIZE=15' + '\n')
	output_stream.write('PRIMER_MAX_SIZE=21' + '\n')
	output_stream.write('PRIMER_MAX_NS_ACCEPTED=0' + '\n')
	output_stream.write('PRIMER_LIBERAL_BASE=1' + '\n')

	if right_or_force == 'right':
		output_stream.write('PRIMER_PRODUCT_SIZE_RANGE=21-100' + '\n')
	elif right_or_force == 'force':
		output_stream.write('PRIMER_PRODUCT_SIZE_RANGE=41-100' + '\n')
	
	output_stream.write('P3_FILE_FLAG=0' + '\n')
	output_stream.write('PRIMER_EXPLAIN_FLAG=0' + '\n')
	output_stream.write('=' + '\n')

	return 1

# import arguments and options
parser = OptionParser()
(options, args) = parser.parse_args()

# dictionaries
marker_sequence = {}		# marker -> sequence
marker_SNP_positions = {}	# marker -> [SNP positions, ...]
SNP_marker_information = {}	# SNP_ID -> [marker, SNP position, {F:[s1, s2, r1], R:[s1, s2, r1]}]

# STEP 1. Import sequence assembly
gene_sequence = open(args[0], 'r')

line = gene_sequence.readline()

while line:
	sline = string.split(line)
	
	if len(line) > 0:
		if line[0] == '>':
			marker = sline[0][1:]
			marker_sequence[marker] = ''
		else:
			marker_sequence[marker] += sline[0]

	line = gene_sequence.readline()

gene_sequence.close()

# STEP 2. Import SNP positions and allelic states
SNP_file = open(args[1], 'r')

for line in SNP_file.readlines():
	sline = string.split(line)

	if sline[0] not in marker_SNP_positions.keys():
		marker_SNP_positions[sline[0]] = []

	marker_SNP_positions[sline[0]].append([sline[1], sline[2], sline[3]])

SNP_file.close()

# STEP 3. Output
polymaker_gene_sequence = open(string.split(args[1], '.')[0] + '_polymaker.txt', 'w')
formatted_gene_sequence = open(string.split(args[1], '.')[0] + '_IUPAC.fa', 'w')

for marker in marker_SNP_positions.keys():
	gene_IUPAC = marker_sequence[marker]

	for SNP in marker_SNP_positions[marker]:
		polymaker_gene_sequence.write(marker + ',' + 'XX' + ',' + marker_sequence[marker][(int(SNP[0]) - 101):(int(SNP[0]) - 1)] + '[' + SNP[1] + '/' + SNP[2] + ']' + marker_sequence[marker][int(SNP[0]):(int(SNP[0]) + 100)] + '\n')

		# confirm quality of all SNPs (reference matches expected reference)
		if gene_IUPAC[int(SNP[0]) - 1] == SNP[1]:
			print marker, 'OK', SNP[0], gene_IUPAC[int(SNP[0]) - 1], SNP[1]
		else:
			print marker, 'BAD', SNP[0], gene_IUPAC[int(SNP[0]) - 1], SNP[1]

		gene_IUPAC = gene_IUPAC[:(int(SNP[0]) - 1)] + IUPAC[SNP[1] + SNP[2]] + gene_IUPAC[int(SNP[0]):]

	formatted_gene_sequence.write('>' + marker + '\n')
	formatted_gene_sequence.write(gene_IUPAC + '\n')

polymaker_gene_sequence.close()
formatted_gene_sequence.close()
