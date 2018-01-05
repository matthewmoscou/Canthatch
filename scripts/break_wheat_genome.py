import string

fasta = open('161010_Chinese_Spring_v1.0_pseudomolecules.fasta', 'r')

Agenome = open('161010_Chinese_Spring_v1.0_pseudomoleculesA.fasta', 'w')
Bgenome = open('161010_Chinese_Spring_v1.0_pseudomoleculesB.fasta', 'w')
Dgenome = open('161010_Chinese_Spring_v1.0_pseudomoleculesD.fasta', 'w')

for line in fasta.readlines():
	if len(line) > 0:
		if line[0] == '>':
			genome = line[5]

		if genome == 'A':
			Agenome.write(line)
		elif genome == 'B':
			Bgenome.write(line)
		elif genome == 'D':
			Dgenome.write(line)

Agenome.close()
Bgenome.close()
Dgenome.close()
fasta.close()
