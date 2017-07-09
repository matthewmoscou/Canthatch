# link_position_expression.py
# description:	link position information from several sources for selecting SNPs between wild-type Canthatch and mutants, including:
#			Position (relative to NRgene 7D)
#			Expression (expressed gene in leaf from Canthatch and mutants
#			SNPs in genes with expression with non-synonymous changes

import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

import sets

import string

# global variables
cmp_DNA = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'M':'K', 'R':'Y', 'W':'S', 'S':'W', 'Y':'R', 'K':'M', 'V':'B', 'H':'D', 'D':'H', 'B':'V', 'X':'X', 'N':'N'}
IUPAC_convert = {'AG':'R', 'GA':'R', 'CT':'Y', 'TC':'Y', 'CA':'M', 'AC':'M', 'TG':'K', 'GT':'K', 'TA':'W', 'AT':'W', 'CG':'S', 'GC':'S'}


directory = 'analysis_95_5/'

# reverse complement
# input : DNA sequence
# output : reverse complement of said DNA sequence
def reverse_complement(orig_DNA):
	rev_comp_DNA = ''
	for index in range(len(orig_DNA)):
		rev_comp_DNA += cmp_DNA[orig_DNA[len(orig_DNA)-index-1]]
	return rev_comp_DNA

# read in reference genome, determine size of all contigs
fasta = open('Can.7DL_edena_clean_v1.fasta', 'r')

ID_sequence = {}

for line in fasta.readlines():
	if len(line) > 0:
		if line[0] == '>':
			ID = string.split(line)[0][1:]
			ID_sequence[ID] = ''
		else:
			ID_sequence[ID] += string.split(line)[0]

fasta.close()

# read in BLAST information
BLAST_file = open('BLAST_NRgene/Can.7DL_edena_clean_v1.fasta_chr7D_blastn_analysis.txt', 'r')

contig_BLAST = {}

for line in BLAST_file.readlines():
	sline = string.split(line)

	contig_BLAST[sline[0]] = [float(sline[4]) / float(sline[5]), sline[6], sline[7], sline[8]]

	#print sline[0] + '\t' + str(float(sline[4]) / float(sline[5])) + '\t' + sline[5] + '\t' + str(len(ID_sequence[sline[0]]))

BLAST_file.close()

# read in cufflinks
contig_transcriptID_transcript_start_stop_strand = {}
transcriptID_exons = {}
contig_gene_clusters = {}
transcriptID_CDS_start_stop = {}

# initialize contigs
GTF_file = open('cufflinks_all/transcripts.gtf', 'r')

for line in GTF_file.readlines():
	sline = string.replace(line, '\n', '')
	sline = string.split(line, '\t')

	contig_transcriptID_transcript_start_stop_strand[sline[0]] = {}
	
GTF_file.close()

# reads information
GTF_file = open('cufflinks_all/transcripts.gtf', 'r')

for line in GTF_file.readlines():
	sline = string.replace(line, '\n', '')
	sline = string.split(line, '\t')
		
	transcriptID = string.replace(string.replace(string.replace(string.split(sline[8])[3], '"', ''), ';', ''), 'CUFF', 'ALL')

	if sline[2] == 'transcript':
		contig_transcriptID_transcript_start_stop_strand[sline[0]][transcriptID] = [int(sline[3]) - 1, int(sline[4]) - 1, sline[6]]
		transcriptID_exons[transcriptID] = []

	if sline[2] == 'exon':
		transcriptID_exons[transcriptID].append([int(sline[3]) - 1, int(sline[4]) - 1])


	transcriptID_CDS_start_stop[transcriptID] = {}

GTF_file.close()

# import SNPs
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

# identify polymorphisms between wild-type and mutant
# need to understand the data files I want out of this part of the script and adjust accordingly
# generate independent modified set within this scan
all_SNPs = []
selected_contigs = []
both_SNP_SNPs = []

NS1_ID_sequence = {}
NS2_ID_sequence = {}
IUPAC_ID_sequence = {}

NS1_contig_SNP_file = open(directory + 'Canthatch_mutant_flow_sorting_Can.7DL_HQ_SNPs_NS1_position_cufflinks.txt', 'w')
NS2_contig_SNP_file = open(directory + 'Canthatch_mutant_flow_sorting_Can.7DL_HQ_SNPs_NS2_position_cufflinks.txt', 'w')

NS1_contig_SNP_file.write('contig' + '\t' + 'position' + '\t' + 'reference_allele' + '\t' + 'mutant' + '\t' + 'NS1_allele' + '\t' + 'chr' + '\t' + 'chr_position' + '\t' + 'expressed_genes' + '\n')
NS2_contig_SNP_file.write('contig' + '\t' + 'position' + '\t' + 'reference_allele' + '\t' + 'mutant' + '\t' + 'NS2_allele' + '\t' + 'chr' + '\t' + 'chr_position' + '\t' + 'expressed_genes' + '\n')

for contig in contig_pdata.keys():
	both_genotypes = [False, False]

	for SNP in contig_pdata[contig]:
		all_SNPs.append(contig + '_' + str(SNP))
		both_SNPs = False

		# Reference is different from NS1
		if SNP[2] != SNP[4]:
			both_genotypes[0] = True

			# Reference is different from NS2
			if SNP[2] != SNP[6]:
				both_SNPs = True

		# Reference is different from NS2
		if SNP[2] != SNP[6]:
			both_genotypes[1] = True
			
			# Reference is different from NS1
			if SNP[2] != SNP[4]:
				both_SNPs = True
	
		if both_SNPs:
			both_SNP_SNPs.append(SNP)
			both_genotypes = [False, False]

		if not (both_genotypes[0] and both_genotypes[1]):
			# if a gene is expressed in Canthatch, provide names
			if contig in contig_transcriptID_transcript_start_stop_strand.keys():
				expressed_gene = '|'.join(contig_transcriptID_transcript_start_stop_strand[contig])
			else:
				expressed_gene = ''

			# if sequence maps onto chr 7D of NRGene assembly, provide position
			if contig in contig_BLAST.keys():
				contig_position = 'chr7D' + '\t' + contig_BLAST[contig][2]
			else:
				contig_position = '\t'

			if both_genotypes[0]:
				if contig not in IUPAC_ID_sequence.keys():
					IUPAC_ID_sequence[contig] = ID_sequence[contig]

				if contig not in NS1_ID_sequence.keys():
					NS1_ID_sequence[contig] = ID_sequence[contig]

				NS1_contig_SNP_file.write(contig + '\t' + SNP[1] + '\t' + SNP[2] + '\t' + 'NS1' + '\t' + SNP[4] + '\t' + contig_position + '\t' + expressed_gene + '\n')
				NS1_ID_sequence[contig] = NS1_ID_sequence[contig][:int(SNP[1]) - 1] + SNP[4] + NS1_ID_sequence[contig][int(SNP[1]):]
				IUPAC_ID_sequence[contig] = IUPAC_ID_sequence[contig][:int(SNP[1]) - 1] + IUPAC_convert[IUPAC_ID_sequence[contig][int(SNP[1]) - 1] + SNP[4]] + IUPAC_ID_sequence[contig][int(SNP[1]):]
			if both_genotypes[1]:
				if contig not in IUPAC_ID_sequence.keys():
					IUPAC_ID_sequence[contig] = ID_sequence[contig]

				if contig not in NS2_ID_sequence.keys():
					NS2_ID_sequence[contig] = ID_sequence[contig]

				NS2_contig_SNP_file.write(contig + '\t' + SNP[1] + '\t' + SNP[2] + '\t' + 'NS2' + '\t' + SNP[6] + '\t' + contig_position + '\t' + expressed_gene + '\n')
				NS2_ID_sequence[contig] = NS2_ID_sequence[contig][:int(SNP[1]) - 1] + SNP[6] + NS2_ID_sequence[contig][int(SNP[1]):]
				IUPAC_ID_sequence[contig] = IUPAC_ID_sequence[contig][:int(SNP[1]) - 1] + IUPAC_convert[IUPAC_ID_sequence[contig][int(SNP[1]) - 1] + SNP[6]] + IUPAC_ID_sequence[contig][int(SNP[1]):]

		# QC01
		#print contig + '\t' + SNP[2] + '\t' + ID_sequence[contig][int(SNP[1]) - 1]

	if both_genotypes[0] and both_genotypes[1]:
		selected_contigs.append(contig)

NS1_contig_SNP_file.close()
NS2_contig_SNP_file.close()

print len(sets.Set(ID_sequence.keys()) & sets.Set(contig_BLAST.keys()))
print len(sets.Set(contig_BLAST.keys()) - sets.Set(ID_sequence.keys()))
print len(sets.Set(ID_sequence.keys()) - sets.Set(contig_BLAST.keys()))
print list(sets.Set(ID_sequence.keys()) - sets.Set(contig_BLAST.keys()))

# export contigs that do not have hits
missing_BLAST = open(directory + 'missing_BLAST_contigs.fa', 'w')

for contig in list(sets.Set(ID_sequence.keys()) - sets.Set(contig_BLAST.keys())):
	missing_BLAST.write('>' + contig + '\n')
	missing_BLAST.write(ID_sequence[contig] + '\n')

missing_BLAST.close()

# export IUPAC FASTA files for generating KASP markers
IUPAC_KASP = open(directory + 'IUPAC_KASP_contigs.fa', 'w')

for contig in IUPAC_ID_sequence.keys():
	IUPAC_KASP.write('>' + contig + '\n')
	IUPAC_KASP.write(IUPAC_ID_sequence[contig] + '\n')

IUPAC_KASP.close()

# export transcript files for transdecoder analysis
fasta_strand_file = open(directory + 'Canthatch_cufflinks_ALL_strand.fa', 'w')
fasta_nostrand_file = open(directory + 'Canthatch_cufflinks_ALL_nostrand.fa', 'w')

transcriptID_wt_mt_sequence = {}

for contig in contig_transcriptID_transcript_start_stop_strand.keys():
	contig_length = len(ID_sequence[contig])
	contig_gene_clusters[contig] = []

	contig_transcript_space = []

	for base_index in range(contig_length):
		contig_transcript_space.append([])

	# populating transcript coverage space
	for transcriptID in contig_transcriptID_transcript_start_stop_strand[contig].keys():
		for base_index in range(contig_transcriptID_transcript_start_stop_strand[contig][transcriptID][0], contig_transcriptID_transcript_start_stop_strand[contig][transcriptID][1]):
			contig_transcript_space[base_index].append(transcriptID)


	# identify non-overlapping segments
	gene_clusters = []
	truth = False
	cluster_index = 0

	for base_index in range(contig_length):
		if len(contig_transcript_space[base_index]) > 0:
			if not truth:
				gene_clusters.append([])

			truth = True

			for transcriptID in contig_transcript_space[base_index]:
				gene_clusters[cluster_index].append(transcriptID)
		else:
			if truth:
				cluster_index += 1

			truth = False

	# export concatenated exons of gene clusters
	for gene_index in range(len(gene_clusters)):
		contig_gene_clusters[contig].append(list(sets.Set(gene_clusters[gene_index])))

		for transcriptID in list(sets.Set(gene_clusters[gene_index])):
			coding_sequence = ''
			coding_sequence_NS1 = ''
			coding_sequence_NS2 = ''

			for exon in transcriptID_exons[transcriptID]:
				coding_sequence += ID_sequence[contig][exon[0]:exon[1] + 1]
				
				if contig in NS1_ID_sequence.keys():
					coding_sequence_NS1 += NS1_ID_sequence[contig][exon[0]:exon[1] + 1]
				else:
					coding_sequence_NS1 += ID_sequence[contig][exon[0]:exon[1] + 1]

				if contig in NS2_ID_sequence.keys():
					coding_sequence_NS2 += NS2_ID_sequence[contig][exon[0]:exon[1] + 1]
				else:
					coding_sequence_NS2 += ID_sequence[contig][exon[0]:exon[1] + 1]

			if contig_transcriptID_transcript_start_stop_strand[contig][transcriptID][2] == '+':
				fasta_strand_file.write('>' + transcriptID + '\n')
				fasta_strand_file.write(coding_sequence + '\n')
				transcriptID_wt_mt_sequence[transcriptID] = [contig, coding_sequence, coding_sequence_NS1, coding_sequence_NS2]

			elif contig_transcriptID_transcript_start_stop_strand[contig][transcriptID][2] == '-':
				fasta_strand_file.write('>' + transcriptID + '\n')
				fasta_strand_file.write(reverse_complement(coding_sequence) + '\n')
				transcriptID_wt_mt_sequence[transcriptID] = [contig, reverse_complement(coding_sequence), reverse_complement(coding_sequence_NS1), reverse_complement(coding_sequence_NS2)]
			else:
				fasta_nostrand_file.write('>' + transcriptID + '\n')
				fasta_nostrand_file.write(coding_sequence + '\n')
				transcriptID_wt_mt_sequence[transcriptID] = [contig, coding_sequence, coding_sequence_NS1, coding_sequence_NS2]

fasta_strand_file.close()
fasta_nostrand_file.close()

raw_input('Please press enter if Transdecoder analysis is complete')

# import GFF file for longest ORFs
# need to pick up from here, not functional after alcohol is consumed...
# need to import the position of the CDS and determine if mutations occur within the CDS (relative to original sequence)
# create mutant form and parse sequence and translate, compare
gff_file = open('Canthatch_cufflinks_ALL_nostrand.fa.transdecoder_dir/longest_orfs.gff3', 'r')

# for each locus
# identify the longest ORF and second longest non-overlapping ORF with >50% length of longest ORF

for line in gff_file.readlines():
	sline = string.split(line)

	if len(sline) > 2:
		if sline[2] == 'gene':
			truncated_description = sline[8][sline[8].index('%3A'):]
			gene_model = truncated_description[3:truncated_description.index('%20')]
		if sline[2] == 'CDS':
			transcriptID_CDS_start_stop[sline[0]][sline[8][(sline[8].index('|') + 1):sline[8].index(';')]] = [int(sline[3]), int(sline[4]), gene_model]

gff_file.close()

gff_file = open('Canthatch_cufflinks_ALL_strand.fa.transdecoder_dir/longest_orfs.gff3', 'r')

# for each locus
# identify the longest ORF and second longest non-overlapping ORF with >50% length of longest ORF

for line in gff_file.readlines():
	sline = string.split(line)

	if len(sline) > 2:
		if sline[2] == 'gene':
			truncated_description = sline[8][sline[8].index('%3A'):]
			gene_model = truncated_description[3:truncated_description.index('%20')]
		if sline[2] == 'CDS':
			transcriptID_CDS_start_stop[sline[0]][sline[8][(sline[8].index('|') + 1):sline[8].index(';')]] = [int(sline[3]), int(sline[4]), gene_model]

gff_file.close()


# identify if SNPs lead to protein sequence changes
POI = open(directory + 'Canthatch_mutant_flow_sorting_Can.7DL_HQ_SNPs_POI.fa', 'w')

for transcriptID in transcriptID_wt_mt_sequence.keys():
	for gene_model in transcriptID_CDS_start_stop[transcriptID].keys():
		start = transcriptID_CDS_start_stop[transcriptID][gene_model][0]
		stop = transcriptID_CDS_start_stop[transcriptID][gene_model][1]

		wildtype_CDS =  Seq(transcriptID_wt_mt_sequence[transcriptID][1][(start - 1):stop], IUPAC.unambiguous_dna)
		NS1_CDS =  Seq(transcriptID_wt_mt_sequence[transcriptID][2][(start - 1):stop], IUPAC.unambiguous_dna)
		NS2_CDS =  Seq(transcriptID_wt_mt_sequence[transcriptID][3][(start - 1):stop], IUPAC.unambiguous_dna)

		wildtype_pep = wildtype_CDS.translate(to_stop=True)
		NS1_pep = NS1_CDS.translate(to_stop=True)
		NS2_pep = NS2_CDS.translate(to_stop=True)

		if len(sets.Set([str(wildtype_pep), str(NS1_pep), str(NS2_pep)])) > 1:
			print transcriptID_wt_mt_sequence[transcriptID][0], transcriptID, len(sets.Set([str(wildtype_pep), str(NS1_pep), str(NS2_pep)]))
			POI.write('>' + transcriptID_wt_mt_sequence[transcriptID][0] + '|' + transcriptID + '\n')
			POI.write(str(wildtype_pep) + '\n')

POI.close()
