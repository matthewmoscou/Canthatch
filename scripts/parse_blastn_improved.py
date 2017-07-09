#! /usr/bin/python

# parse_blastn.py
# description:	input a blastn output file and output a variety of information
# input:	input blastn output

import optparse
from optparse import OptionParser 

import pickle

import string

# import arguments and options
parser = OptionParser()
(options, args) = parser.parse_args()

blastn_file = open(args[0], 'r')
output_file = open(args[1], 'w')

current_id = ''
subject = ''
id = 0

while True:
	line = blastn_file.readline()
	
	if len(line) == 0:
		break

	sline = string.split(line)

	if string.find(line, 'BLASTN') >= 0:
		id += 1
		first_truth = 0
		identity_truth = 0
	
	if len(line) >= 6 and len(sline) > 0:
		if sline[0] == 'Query=':
			# no white space (and no suffix characters)
			current_id = sline[1]

			# white space wanted (and suffix characters)
#			current_id = line[7:len(line)-1]
	
		if sline[0] == 'Score' and first_truth == 1:
			score = sline[2]
			e_value = sline[7]

		if sline[0] == 'Strand' and first_truth == 1:
			if sline[2] == sline[4]:
				strand = '+'
			else:
				strand = '-'
			

	if len(line) >= 0 and len(sline) > 0:
		if line[0] == '>' and first_truth == 0:
			first_truth += 1

		if line[0] == '>':
			subject = sline[0][1:]

		if sline[0] == 'Query:' and first_truth == 1:
			query_start = sline[1]

		if sline[0] == 'Sbjct:' and first_truth == 1:
			first_truth += 1
			subject_start = sline[1]

			if len(subject) > 0:
				output_file.write(current_id + '\t' + subject + '\t' + score + '\t' + e_value + '\t' + positives[0] + '\t' + positives[1] + '\t' + query_start + '\t' + subject_start + '\t' + strand + '\t' + '\n')
			
			subject = ''

	
	if len(line) >= 11 and len(sline) > 0:
		if sline[0] == 'Identities':
			if len(sline) > 3:
				positives = string.split(sline[2], '/')
			else:
				print current_id, subject
	
#	if len(line) > 2:
#		if line[0] == '>':
#			print current_id, first_truth, subject
#
#			raw_input('...')

blastn_file.close()
output_file.close()
