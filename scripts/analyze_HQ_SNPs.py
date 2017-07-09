#! /usr/bin/python

import string

DNA = ['A', 'C', 'G', 'T']

directory = 'analysis_95_5/'

HQ_SNP_file = open(directory + 'Canthatch_mutant_flow_sorting_Can.7DL_HQ_SNPs.txt', 'r')

truth = False

contig_data = {}
contig_pdata = {}

for line in HQ_SNP_file.readlines():
	sline = string.split(line)

	if truth:
		if sline[0] not in contig_data.keys():
			contig_data[sline[0]] = []
			contig_pdata[sline[0]] = []
	
		contig_data[sline[0]].append(line)
		contig_pdata[sline[0]].append(sline)
	
	truth = True

HQ_SNP_file.close()

multiple_SNP_contigs = open(directory + 'Canthatch_mutant_flow_sorting_Can.7DL_HQ_contigs_no_selection.txt', 'w')

for contig in contig_data.keys():
	if len(contig_data[contig]) > 1:
		for line in contig_data[contig]:
			multiple_SNP_contigs.write(line)

multiple_SNP_contigs.close()

multiple_SNP_contigs_selection = open(directory + 'Canthatch_mutant_flow_sorting_Can.7DL_HQ_contigs_selection.txt', 'w')
nonref_SNP_contigs_selection = open(directory + 'Canthatch_mutant_flow_sorting_Can.7DL_HQ_contigs_nonreference_SNP.txt', 'w')

all_SNPs = []
selected_contigs = []
both_SNP_SNPs = []
state_transition = {}

for baseX in DNA:
	state_transition[baseX] = {}

	for baseY in DNA:
		state_transition[baseX][baseY] = [0, 0]

for contig in contig_pdata.keys():
	both_genotypes = [False, False]

	for SNP in contig_pdata[contig]:
		all_SNPs.append(contig + '_' + str(SNP))
		both_SNPs = False

		if SNP[2] != SNP[4]:
			both_genotypes[0] = True

			if SNP[2] != SNP[6]:
				both_SNPs = True
			else:
				state_transition[SNP[2]][SNP[4]][0] += 1

		if SNP[2] != SNP[6]:
			both_genotypes[1] = True
			
			if SNP[2] != SNP[4]:
				both_SNPs = True
			else:
				state_transition[SNP[2]][SNP[6]][1] += 1
	
		if both_SNPs:
			both_SNP_SNPs.append(SNP)
			both_genotypes = [False, False]

	if both_genotypes[0] and both_genotypes[1]:
		selected_contigs.append(contig)

		for line in contig_data[contig]:
			multiple_SNP_contigs_selection.write(line)
	
	if both_genotypes[0] or both_genotypes[1]:
		for line in contig_data[contig]:
			nonref_SNP_contigs_selection.write(line)
	
multiple_SNP_contigs_selection.close()
nonref_SNP_contigs_selection.close()

print 'All SNPs:', len(all_SNPs)
print 'Bad SNPs:', len(both_SNP_SNPs)
print 'Selected contigs:', len(selected_contigs)

for baseX in DNA:
	for baseY in DNA:
		if ''.join([baseX, baseY]) in ['CT', 'GA']:
			print baseX + ' -> ' + baseY + '\t' + str(state_transition[baseX][baseY][0]) + '\t' + str(state_transition[baseX][baseY][1]) + '\t' + str(sum(state_transition[baseX][baseY])) + '\t' + 'EMS'
		else:
			print baseX + ' -> ' + baseY + '\t' + str(state_transition[baseX][baseY][0]) + '\t' + str(state_transition[baseX][baseY][1]) + '\t' + str(sum(state_transition[baseX][baseY]))

print sum(state_transition['C']['T']) + sum(state_transition['G']['A'])
print len(all_SNPs) - len(both_SNP_SNPs)

data_file = open('Can.7DL_edena_clean_v1.fasta.masked', 'r')
output_file = open(directory + 'Can.7DL_edena_clean_v1.fasta.masked_Can.7DL_HQ_contigs_selection.fa', 'w')

export_flag = False

for line in data_file.readlines():
	sline = string.split(line)

	if sline[0][0] == '>':
		if sline[0][1:] in selected_contigs:
			export_flag = True
		else:
			export_flag = False

	if export_flag:
		output_file.write(line)

data_file.close()
output_file.close()

data_file = open('Can.7DL_edena_clean_v1.fasta', 'r')
output_file = open(directory + 'Can.7DL_edena_clean_v1_Can.7DL_HQ_contigs_selection.fa', 'w')

export_flag = False

for line in data_file.readlines():
	sline = string.split(line)

	if sline[0][0] == '>':
		if sline[0][1:] in selected_contigs:
			export_flag = True
		else:
			export_flag = False

	if export_flag:
		output_file.write(line)

data_file.close()
output_file.close()
