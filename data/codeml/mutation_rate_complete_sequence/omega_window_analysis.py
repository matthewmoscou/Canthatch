#! /usr/bin/python
# omega_window_analysis.py

import commands
import string

window_size = 420
step_size = 30

gene_sequence = {}

phylip_file = open('Med15_phylogeny_Os_outgroup_complete.phy', 'r')

line = phylip_file.readline()
sequences, proteinlen = string.split(line)
line = phylip_file.readline()

while line:
	sline = string.split(line)
	gene_sequence[sline[0]] = sline[1]
	line = phylip_file.readline()

phylip_file.close()

window_result_file = open('window_' + str(window_size) + '_stepsize_' + str(step_size) + '_results.txt', 'w')

window_result_file.write('index' + '\t' + 'H0_lnL' + '\t' + 'H0_omega' + '\t' + 'H2_lnL' + '\t' + 'H2_omega1' + '\t' + 'H2_omega2' + '\n')

index = 0

lnL_H0 = []
lnL_H2 = []

omega_H0 = []
omega_H2 = []

while index < (int(proteinlen) - window_size):
	print index

	window_result_file.write(str(index) + '\t')

	# generate sequence file
	phylip_file = open('Med15_phylogeny_Os_outgroup_complete_' + str(index) + '.phy', 'w')
	phylip_file.write(' ' + sequences + ' ' + str(window_size) + '\n')

	for gene in gene_sequence.keys():
		phylip_file.write(gene + ' ' * (12 - len(gene)) + gene_sequence[gene][index:(index + window_size)] + '\n')
	
	phylip_file.close()

	# generate control file H0
	control_file = open('codeml.ctl', 'w')

	control_file.write('      seqfile = Med15_phylogeny_Os_outgroup_complete_' + str(index) + '.phy' + '\n')
	control_file.write('     treefile = RAxML_bestTree.Med15_phylogeny_Os_outgroup_complete_H0' + '\n')
	control_file.write('      outfile = results.' + str(index) + '.H0.txt' + '\n')
	control_file.write('        noisy = 1      * 0,1,2,3,9: how much rubbish on the screen' + '\n')
	control_file.write('      verbose = 1      * 1:detailed output' + '\n')
	control_file.write('      runmode = 0      * 0:user defined tree' + '\n')
	control_file.write('      seqtype = 1      * 1:codons' + '\n')
	control_file.write('    CodonFreq = 2      * 0:equal, 1:F1X4, 2:F3X4, 3:F61' + '\n')
	control_file.write('        model = 0      * 0:one omega ratio for all branches' + '\n')
	control_file.write('      NSsites = 0      * ' + '\n')
	control_file.write('        icode = 0      * 0:universal code' + '\n')
	control_file.write('    fix_kappa = 0      * 1:kappa fixed, 0:kappa to be estimated' + '\n')
	control_file.write('        kappa = 2      * initial or fixed kappa' + '\n')
	control_file.write('    fix_omega = 0      * 1:omega fixed, 0:omega to be estimated ' + '\n')
	control_file.write('        omega = 0.2    * initial omega' + '\n')
	control_file.write('    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?' + '\n')

	control_file.close()
	
	# run H0
	commands.getstatusoutput('codeml')
	
	# parse output, move file
	results_file = open('results.' + str(index) + '.H0.txt', 'r')

	for line in results_file.readlines():
		sline = string.split(line)

		if len(line) > 5:
			if line[:3] == 'lnL':
				lnL_H0.append(sline[4])
				window_result_file.write(sline[4] + '\t')
			elif line[:5] == 'omega':
				omega_H0.append(sline[3])
				window_result_file.write(sline[3] + '\t')

	# generate control file H2window
	control_file = open('codeml.ctl', 'w')

	control_file.write('      seqfile = Med15_phylogeny_Os_outgroup_complete_' + str(index) + '.phy' + '\n')
	control_file.write('     treefile = RAxML_bestTree.Med15_phylogeny_Os_outgroup_complete_H2window' + '\n')
	control_file.write('      outfile = results.' + str(index) + '.H2.txt' + '\n')
	control_file.write('        noisy = 1      * 0,1,2,3,9: how much rubbish on the screen' + '\n')
	control_file.write('      verbose = 1      * 1:detailed output' + '\n')
	control_file.write('      runmode = 0      * 0:user defined tree' + '\n')
	control_file.write('      seqtype = 1      * 1:codons' + '\n')
	control_file.write('    CodonFreq = 2      * 0:equal, 1:F1X4, 2:F3X4, 3:F61' + '\n')
	control_file.write('        model = 2      * 0:one omega ratio for all branches' + '\n')
	control_file.write('      NSsites = 0      * ' + '\n')
	control_file.write('        icode = 0      * 0:universal code' + '\n')
	control_file.write('    fix_kappa = 0      * 1:kappa fixed, 0:kappa to be estimated' + '\n')
	control_file.write('        kappa = 2      * initial or fixed kappa' + '\n')
	control_file.write('    fix_omega = 0      * 1:omega fixed, 0:omega to be estimated ' + '\n')
	control_file.write('        omega = 0.2    * initial omega' + '\n')
	control_file.write('    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?' + '\n')

	control_file.close()
	
	# run H2window
	commands.getstatusoutput('codeml')
	
	# parse output, move file
	results_file = open('results.' + str(index) + '.H2.txt', 'r')

	for line in results_file.readlines():
		sline = string.split(line)

		if len(line) > 9:
			if line[:3] == 'lnL':
				lnL_H0.append(sline[4])
				window_result_file.write(sline[4] + '\t')
			elif line[:9] == 'w (dN/dS)':
				omega_H2.append([sline[4], sline[5]])
				window_result_file.write(sline[4] + '\t' + sline[5] + '\n')

	index += step_size

window_result_file.close()
