#! /usr/bin/python

# 0_Canthatch_mutant_flow_sorting_analysis_MM2.py
# description:	Relaxed parameters for SNP identification, improve EMS rate, estimates over variable thresholds

import sets
import string

# import_SNPs
# reads through a VarScan output file and only retains SNPs that are above the read count and percent threshold
# contig -> position of SNP -> [reference allele, alternate allele, percent variant frequency]
def import_SNPs(varscan, read_threshold, percent_threshold):
	contig_position_allele = {}

	varscan_file = open(varscan, 'r')

	truth = False

	for line in varscan_file.readlines():
		sline = string.split(line)

		if truth:
			if (int(sline[4]) + int(sline[5])) >= read_threshold:
				if float(string.replace(sline[6], '%', '')) >= percent_threshold:
					if sline[0] not in contig_position_allele.keys():
						contig_position_allele[sline[0]] = {}

					contig_position_allele[sline[0]][int(sline[1])] = [sline[2], sline[18], string.replace(sline[6], '%', '')]

		truth = True

	varscan_file.close()

	return contig_position_allele

# function could be converted to only look in contig_SNPs once between changes in contig ID
def complete_SNPs(reference_genotype, genotype, contig_SNPs):
	contig_position_coverage = {}

	genomecov_file = open(reference_genotype + '_' + genotype + '.sorted.genomecov.threshold', 'r')

	line = genomecov_file.readline()
	current_ID = ''
	truth = False

	while line:
		sline = string.split(line)

		if sline[0] != current_ID:
			if sline[0] in contig_SNPs.keys():
				truth = True
			else:
				truth = False
			current_ID = sline[0]

		if truth:
			if sline[0] not in contig_position_coverage.keys():
				contig_position_coverage[sline[0]] = {}

			if int(sline[1]) in contig_SNPs[sline[0]]:
				contig_position_coverage[sline[0]][int(sline[1])] = float(sline[2])
			
		line = genomecov_file.readline()

	genomecov_file.close()

	return contig_position_coverage

genotypes = ['Can.7DL', 'NS1M', 'NS2N']
mutant_genotypes = ['NS1M', 'NS2N']
genotype_contig_position_allele = {}
genotype_contig_position_coverage = {}
contig_SNPs = {}

reference_genotype = 'Can.7DL'
global_coverage = 10
wildtype_percent_threshold = 5.0
mutant_percent_threshold = 95.0
directory = 'analysis_95_5/'

for genotype in mutant_genotypes:
	print 'reads variants from', genotype
	genotype_contig_position_allele[genotype] = import_SNPs(genotype + '.sorted.rmdup.pileup2snp.txt', global_coverage, mutant_percent_threshold)

	for contig in genotype_contig_position_allele[genotype].keys():
		if contig not in contig_SNPs.keys():
			contig_SNPs[contig] = {}

		for position in genotype_contig_position_allele[genotype][contig].keys():
			contig_SNPs[contig][position] = genotype_contig_position_allele[genotype][contig][position][0]

# evaluate reference accession for variants based on relaxed parameters
print 'reads variants from reference'
reference_contig_position_allele = import_SNPs(reference_genotype + '.sorted.rmdup.pileup2snp.txt', global_coverage, wildtype_percent_threshold)

# for every contig and position, remove all SNPs that are different based on self-masking
for contig in reference_contig_position_allele.keys():
	if contig in contig_SNPs.keys():
		for position in reference_contig_position_allele[contig].keys():
			if position in contig_SNPs[contig].keys():
				del contig_SNPs[contig][position]

for contig in reference_contig_position_allele.keys():
	if contig in contig_SNPs.keys():
		if len(contig_SNPs[contig]) == 0:
			del contig_SNPs[contig]

#raw_input('...')

#print len(contig_SNPs.keys())
#print contig_SNPs.keys()[0]
#print contig_SNPs['out_33682']

#raw_input('...')

HQ_SNPs = open(directory + 'Canthatch_mutant_flow_sorting_' + reference_genotype + '_HQ_SNPs.txt', 'w')

HQ_SNPs.write('contig' + '\t' + 'position')

for genotype in genotypes:
	HQ_SNPs.write('\t' + genotype)

HQ_SNPs.write('\n')

# allele should not be N in reference
for contig in contig_SNPs.keys():
	#print contig

	for SNP in contig_SNPs[contig]:
		truth = True

# contig -> position of SNP -> [reference allele, alternate allele, percent variant frequency]
		for genotype in mutant_genotypes:
			if contig in genotype_contig_position_allele[genotype].keys():
				if SNP in genotype_contig_position_allele[genotype][contig].keys():
					if genotype_contig_position_allele[genotype][contig][SNP][0] == 'N':
						truth = False
						
		if contig in reference_contig_position_allele.keys():
			if SNP in reference_contig_position_allele[contig].keys():
				if reference_contig_position_allele[contig][SNP][0] == 'N':
					truth = False
					
		if truth:
			HQ_SNPs.write(contig + '\t' + str(SNP))

			for genotype in genotypes:
				if genotype == reference_genotype:
					if contig in reference_contig_position_allele.keys():
						if SNP in reference_contig_position_allele[contig].keys():
							HQ_SNPs.write('\t' + reference_contig_position_allele[contig][SNP][1] + '\t' + reference_contig_position_allele[contig][SNP][2])
						else:
							HQ_SNPs.write('\t' + contig_SNPs[contig][SNP] + '\t' + 'NA')
					else:
						HQ_SNPs.write('\t' + contig_SNPs[contig][SNP] + '\t' + 'NA')
				else:
					if contig in genotype_contig_position_allele[genotype].keys():
						if SNP in genotype_contig_position_allele[genotype][contig].keys():
							HQ_SNPs.write('\t' + genotype_contig_position_allele[genotype][contig][SNP][1] + '\t' + genotype_contig_position_allele[genotype][contig][SNP][2])
						else:
							HQ_SNPs.write('\t' + contig_SNPs[contig][SNP] + '\t' + 'NA')
					else:
						HQ_SNPs.write('\t' + contig_SNPs[contig][SNP] + '\t' + 'NA')

			HQ_SNPs.write('\n')
		
HQ_SNPs.close()

# need to create a script that merges genomecov and pileup2snp to create a matrix of conserved SNPs relative to the reference genome
# as a special case, make morex the default whenever inputing a SNP for the first time... or store both? this is easiest...
# makes sense to start with SNPs, as these are a minority of the observations
# read all SNPs, then quality assess against coverage in other transcriptomes
