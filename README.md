# Canthatch
In 1980, Kerber and Green identified the presence of a suppressor of resistance in wheat to wheat stem rust (*Puccinia graminis* f. sp. *tritici*). Kerber and Green (1980) mapped the suppression to chromosome 7D, and Kerber (1991) established that a single locus conferred suppression on chromosome 7D. We set out to identify the suppressor gene in Canthatch by applying chromosome flow sorting and high throughput sequencing to Canthatch and two EMS-derived mutants.

# Table of Contents


## Defining the genetic interval encompassing *Srs1*
Our initial approach was to identify SNPs along the long arm of chromosome 7D in order to develop SNP markers that will be applied to the Canthatch x NS1 and Canthatch x NS2 doubled-haploid mapping populations. We used an approach that integrates multiple data sets including *de novo* assembly of flow sorted chromosomes of Canthatch, RNAseq data derived from Canthatch, NS1, and NS2, alignment-based SNP calling, and the physical assembly of chromosome 7D from the IWGSC (NRGene assembly).

### *De novo* assembly of Canthatch, NS1, and NS2 flow sorted chromosome arm 7DL reads
We sequenced flow sorted chromosome arm 7DL from Canthatch (wild-type) and mutants NS1 and NS2. We used high stringency parameters in `Trimmomatic` to identify high quality reads. `Edena` requires that all reads used for assembly have identical size, therefore a considerable number of reads were removed before assembly.

```bash
java -jar trimmomatic-0.36.jar PE -phred33 Can.7DL_DDPL00006_H32VYALXX_L6_1.clean.fq Can.7DL_DDPL00006_H32VYALXX_L6_2.clean.fq Can.7DL_gDNA_forward_paired.fq.gz Can.7DL_gDNA_forward_unpaired.fq.gz Can.7DL_gDNA_reverse_paired.fq.gz Can.7DL_gDNA_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:150 > Can.7DL.trimmomatic.run.log 2>&1 &

java -jar trimmomatic-0.36.jar PE -phred33 Can.7DL_DDPL00006_H32VYALXX_L7_1.clean.fq Can.7DL_DDPL00006_H32VYALXX_L7_2.clean.fq Can.7DL_2_gDNA_forward_paired.fq.gz Can.7DL_2_gDNA_forward_unpaired.fq.gz Can.7DL_2_gDNA_reverse_paired.fq.gz Can.7DL_2_gDNA_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:150 > Can.7DL_2.trimmomatic.run.log 2>&1 &

java -jar trimmomatic-0.36.jar PE -phred33 NS1M_DDPL00004_H32VYALXX_L6_1.clean.fq.gz NS1M_DDPL00004_H32VYALXX_L6_2.clean.fq.gz NS1M_gDNA_forward_paired.fq.gz NS1M_gDNA_forward_unpaired.fq.gz NS1M_gDNA_reverse_paired.fq.gz NS1M_gDNA_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:150 > NS1M.trimmomatic.run.log 2>&1 &

java -jar trimmomatic-0.36.jar PE -phred33 NS1M_DDPL00004_H32VYALXX_L7_1.clean.fq.gz NS1M_DDPL00004_H32VYALXX_L7_2.clean.fq.gz NS1M_2_gDNA_forward_paired.fq.gz NS1M_2_gDNA_forward_unpaired.fq.gz NS1M_2_gDNA_reverse_paired.fq.gz NS1M_2_gDNA_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:150 > NS1M_2.trimmomatic.run.log 2>&1 &

java -jar trimmomatic-0.36.jar PE -phred33 NS2N_DDPL00005_H32VYALXX_L6_1.clean.fq.gz NS2N_DDPL00005_H32VYALXX_L6_2.clean.fq.gz NS2N_gDNA_forward_paired.fq.gz NS2N_gDNA_forward_unpaired.fq.gz NS2N_gDNA_reverse_paired.fq.gz NS2N_gDNA_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:150 > NS2N.trimmomatic.run.log 2>&1 &

java -jar trimmomatic-0.36.jar PE -phred33 NS2N_DDPL00005_H32VYALXX_L7_1.clean.fq.gz NS2N_DDPL00005_H32VYALXX_L7_2.clean.fq.gz NS2N_2_gDNA_forward_paired.fq.gz NS2N_2_gDNA_forward_unpaired.fq.gz NS2N_2_gDNA_reverse_paired.fq.gz NS2N_2_gDNA_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:150 > NS2N_2.trimmomatic.run.log 2>&1 &

java -jar trimmomatic-0.36.jar PE -phred33 Can.K_DDPL00003_H32VYALXX_L6_1.clean.fq.gz Can.K_DDPL00003_H32VYALXX_L6_2.clean.fq.gz Can.K_gDNA_forward_paired.fq.gz Can.K_gDNA_forward_unpaired.fq.gz Can.K_gDNA_reverse_paired.fq.gz Can.K_gDNA_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:150 > Can.K.trimmomatic.run.log 2>&1 &

java -jar trimmomatic-0.36.jar PE -phred33 Can.K_DDPL00003_H32VYALXX_L7_1.clean.fq.gz Can.K_DDPL00003_H32VYALXX_L7_2.clean.fq.gz Can.K_2_gDNA_forward_paired.fq.gz Can.K_2_gDNA_forward_unpaired.fq.gz Can.K_2_gDNA_reverse_paired.fq.gz Can.K_2_gDNA_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:150 > Can.K_2.trimmomatic.run.log 2>&1 &

java -jar trimmomatic-0.36.jar PE -phred33 Can.K_DDPL00003_H32VYALXX_L8_1.clean.fq.gz Can.K_DDPL00003_H32VYALXX_L8_2.clean.fq.gz Can.K_3_gDNA_forward_paired.fq.gz Can.K_3_gDNA_forward_unpaired.fq.gz Can.K_3_gDNA_reverse_paired.fq.gz Can.K_3_gDNA_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:150 > Can.K_3.trimmomatic.run.log 2>&1 &

java -jar trimmomatic-0.36.jar PE -phred33 Tratcher_DDPL00007-w_H3223ALXX_L7_1.clean.fq.gz Tratcher_DDPL00007-w_H3223ALXX_L7_2.clean.fq.gz Thatcher_gDNA_forward_paired.fq.gz Thatcher_gDNA_forward_unpaired.fq.gz Thatcher_gDNA_reverse_paired.fq.gz Thatcher_gDNA_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:150 > Thatcher.trimmomatic.run.log 2>&1 &

java -jar trimmomatic-0.36.jar PE -phred33 Tratcher_DDPL00007-w_H3223ALXX_L8_1.clean.fq.gz Tratcher_DDPL00007-w_H3223ALXX_L8_2.clean.fq.gz Thatcher_2_gDNA_forward_paired.fq.gz Thatcher_2_gDNA_forward_unpaired.fq.gz Thatcher_2_gDNA_reverse_paired.fq.gz Thatcher_2_gDNA_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:150 > Thatcher_2.trimmomatic.run.log 2>&1 &

java -jar trimmomatic-0.36.jar PE -phred33 Tratcher_DDPL00007-w_H3727ALXX_L2_1.clean.fq.gz Tratcher_DDPL00007-w_H3727ALXX_L2_2.clean.fq.gz Thatcher_3_gDNA_forward_paired.fq.gz Thatcher_3_gDNA_forward_unpaired.fq.gz Thatcher_3_gDNA_reverse_paired.fq.gz Thatcher_3_gDNA_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:150 > Thatcher_3.trimmomatic.run.log 2>&1 &
```

Next, we assembled individual chromosome arms genomes using `Edena` for Canthatch, NS1, and NS2.

```bash
# Canthatch
./edena -paired Can.7DL_gDNA_forward_paired.fq Can.7DL_gDNA_reverse_paired.fq Can.7DL_2_gDNA_forward_paired.fq Can.7DL_2_gDNA_reverse_paired.fq -p Can.7DL.clean -nThreads 32 > Can.7DL.edena.clean.run.log 2>&1 &
./edena -e Can.7DL.clean.ovl -m 100

# NS1
./bin/edena -paired NS1M_gDNA_forward_paired.fq NS1M_gDNA_reverse_paired.fq NS1M_2_gDNA_forward_paired.fq NS1M_2_gDNA_reverse_paired.fq -p NS1M.clean -nThreads 32 > NS1M.edena.clean.run.log 2>&1 &
./bin/edena -e NS1M.clean.ovl -m 100

# NS2
./bin/edena -paired NS2N_gDNA_forward_paired.fq NS2N_gDNA_reverse_paired.fq NS2N_2_gDNA_forward_paired.fq NS2N_2_gDNA_reverse_paired.fq -p NS2N.clean -nThreads 32 > NS2N.edena.clean.run.log 2>&1 &
./bin/edena -e NS2N.clean.ovl -m 100
```

The Canthatch 7DL assembly was composed of approximately 129k contigs over 241 Mb.

### Identification of EMS generated SNPs in NS1 and NS2 relative to Canthatch
An alignment-based strategy was taken to identify SNPs between Canthatch and mutants NS1 and NS2, using `bwa` for alignment.

```bash
# Initialize assembly for bwa alignment
bwa index Can.7DL_edena_clean_v1.fasta.masked

# Align reads to reference sequence
bwa aln -t 16 Can.7DL_edena_clean_v1.fasta.masked Can.7DL_gDNA_forward_paired.fq > Can.7DL_1_1.aln
bwa aln -t 16 Can.7DL_edena_clean_v1.fasta.masked Can.7DL_gDNA_reverse_paired.fq > Can.7DL_1_2.aln
bwa aln -t 16 Can.7DL_edena_clean_v1.fasta.masked Can.7DL_2_gDNA_forward_paired.fq > Can.7DL_2_1.aln
bwa aln -t 16 Can.7DL_edena_clean_v1.fasta.masked Can.7DL_2_gDNA_reverse_paired.fq > Can.7DL_2_2.aln

bwa aln -t 16 Can.7DL_edena_clean_v1.fasta.masked NS1M_gDNA_forward_paired.fq > NS1M_1_1.aln
bwa aln -t 16 Can.7DL_edena_clean_v1.fasta.masked NS1M_gDNA_reverse_paired.fq > NS1M_1_2.aln
bwa aln -t 16 Can.7DL_edena_clean_v1.fasta.masked NS1M_2_gDNA_forward_paired.fq > NS1M_2_1.aln
bwa aln -t 16 Can.7DL_edena_clean_v1.fasta.masked NS1M_2_gDNA_reverse_paired.fq > NS1M_2_2.aln

bwa aln -t 16 Can.7DL_edena_clean_v1.fasta.masked NS2N_gDNA_forward_paired.fq > NS2N_1_1.aln
bwa aln -t 16 Can.7DL_edena_clean_v1.fasta.masked NS2N_gDNA_reverse_paired.fq > NS2N_1_2.aln
bwa aln -t 16 Can.7DL_edena_clean_v1.fasta.masked NS2N_2_gDNA_forward_paired.fq > NS2N_2_1.aln
bwa aln -t 16 Can.7DL_edena_clean_v1.fasta.masked NS2N_2_gDNA_reverse_paired.fq > NS2N_2_2.aln

# Merge forward and reverse reads
bwa sampe Can.7DL_edena_clean_v1.fasta.masked Can.7DL_1_1.aln Can.7DL_1_2.aln Can.7DL_gDNA_forward_paired.fq Can.7DL_gDNA_reverse_paired.fq > Can.7DL_1.sam
bwa sampe Can.7DL_edena_clean_v1.fasta.masked Can.7DL_2_1.aln Can.7DL_2_2.aln Can.7DL_2_gDNA_forward_paired.fq Can.7DL_2_gDNA_reverse_paired.fq > Can.7DL_2.sam

bwa sampe Can.7DL_edena_clean_v1.fasta.masked NS1M_1_1.aln NS1M_1_2.aln NS1M_gDNA_forward_paired.fq NS1M_gDNA_reverse_paired.fq > NS1M_1.sam
bwa sampe Can.7DL_edena_clean_v1.fasta.masked NS1M_2_1.aln NS1M_2_2.aln NS1M_2_gDNA_forward_paired.fq NS1M_2_gDNA_reverse_paired.fq > NS1M_2.sam

bwa sampe Can.7DL_edena_clean_v1.fasta.masked NS2N_1_1.aln NS2N_1_2.aln NS2N_gDNA_forward_paired.fq NS2N_gDNA_reverse_paired.fq > NS2N_1.sam
bwa sampe Can.7DL_edena_clean_v1.fasta.masked NS2N_2_1.aln NS2N_2_2.aln NS2N_2_gDNA_forward_paired.fq NS2N_2_gDNA_reverse_paired.fq > NS2N_2.sam

# Select only paired end reads, sort, and remove duplicate reads
samtools view -f2 -Shub -o Can.7DL_1.bam Can.7DL_1.sam
samtools view -f2 -Shub -o Can.7DL_2.bam Can.7DL_2.sam

samtools view -f2 -Shub -o NS1M_1.bam NS1M_1.sam
samtools view -f2 -Shub -o NS1M_2.bam NS1M_2.sam

samtools view -f2 -Shub -o NS2N_1.bam NS2N_1.sam
samtools view -f2 -Shub -o NS2N_2.bam NS2N_2.sam

samtools sort Can.7DL_1.bam Can.7DL_1.sorted
samtools sort Can.7DL_2.bam Can.7DL_2.sorted

samtools sort NS1M_1.bam NS1M_1.sorted
samtools sort NS1M_2.bam NS1M_2.sorted

samtools sort NS2N_1.bam NS2N_1.sorted
samtools sort NS2N_2.bam NS2N_2.sorted

samtools rmdup Can.7DL_1.sorted.bam Can.7DL_1.sorted.rmdup.bam
samtools rmdup Can.7DL_2.sorted.bam Can.7DL_2.sorted.rmdup.bam

samtools rmdup NS1M_1.sorted.bam NS1M_1.sorted.rmdup.bam
samtools rmdup NS1M_2.sorted.bam NS1M_2.sorted.rmdup.bam

samtools rmdup NS2N_1.sorted.bam NS2N_1.sorted.rmdup.bam
samtools rmdup NS2N_2.sorted.bam NS2N_2.sorted.rmdup.bam

samtools merge Can.7DL.sorted.rmdup.bam Can.7DL_1.sorted.rmdup.bam Can.7DL_2.sorted.rmdup.bam
samtools merge NS1M.sorted.rmdup.bam NS1M_1.sorted.rmdup.bam NS1M_2.sorted.rmdup.bam
samtools merge NS2N.sorted.rmdup.bam NS2N_1.sorted.rmdup.bam NS2N_2.sorted.rmdup.bam

samtools index Can.7DL.sorted.rmdup.bam
samtools index NS1M.sorted.rmdup.bam
samtools index NS2N.sorted.rmdup.bam

samtools faidx Can.7DL_edena_clean_v1.fasta.masked

samtools mpileup -f Can.7DL_edena_clean_v1.fasta.masked -BQ0 Can.7DL.sorted.rmdup.bam > Can.7DL.sorted.rmdup.pileup.txt
samtools mpileup -f Can.7DL_edena_clean_v1.fasta.masked -BQ0 NS1M.sorted.rmdup.bam > NS1M.sorted.rmdup.pileup.txt
samtools mpileup -f Can.7DL_edena_clean_v1.fasta.masked -BQ0 NS2N.sorted.rmdup.bam > NS2N.sorted.rmdup.pileup.txt

bedtools genomecov -d -split -ibam Can.7DL.sorted.rmdup.bam > Can.7DL.sorted.rmdup.genomecov.txt
bedtools genomecov -d -split -ibam NS1M.sorted.rmdup.bam > NS1M.sorted.rmdup.genomecov.txt
bedtools genomecov -d -split -ibam NS2N.sorted.rmdup.bam > NS2N.sorted.rmdup.genomecov.txt

java -jar VarScan.v2.3.8.jar mpileup2snp Can.7DL.sorted.rmdup.pileup.txt > Can.7DL.sorted.rmdup.mpileup2snp.txt
java -jar VarScan.v2.3.8.jar mpileup2indel Can.7DL.sorted.rmdup.pileup.txt > Can.7DL.sorted.rmdup.mpileup2indel.txt
java -jar VarScan.v2.3.8.jar mpileup2snp NS1.sorted.rmdup.pileup.txt > NS1.sorted.rmdup.mpileup2snp.txt
java -jar VarScan.v2.3.8.jar mpileup2indel NS1.sorted.rmdup.pileup.txt > NS1.sorted.rmdup.mpileup2indel.txt
java -jar VarScan.v2.3.8.jar mpileup2snp NS2.sorted.rmdup.pileup.txt > NS2.sorted.rmdup.mpileup2snp.txt
java -jar VarScan.v2.3.8.jar mpileup2indel NS2.sorted.rmdup.pileup.txt > NS2.sorted.rmdup.mpileup2indel.txt
```

### Anchoring of genomic contigs onto the NRGene assembly of chromosome 7D
Our approach was to anchor genomic contigs relative to the NRGene assembly in order to make efficient use of the physical space in selecting equidistant markers. We anchored contigs using BLAST.

```bash
blastall -p blastn -d chr7D.fa -i Can.7DL_edena_clean_v1.fasta -o Can.7DL_edena_clean_v1_chr7D_blastn.txt -a 16 -F F -v 1 -b 1
```

### Annotation of leaf expressed genes on *de novo* Canthatch assembly
RNAseq alignment was performd against the *de novo* Canthatch assembly with spliced alignment using `Tophat`. RNAseq data derived from Canthatch, NS1, and NS2 were used for gene prediction.

```bash
bowtie2-build Can.7DL_edena_clean_v1.fasta Can.7DL_edena_clean_v1
tophat --max-intron-length 20000 --num-threads 4 Can.7DL_edena_clean_v1 TA_01918_L1_1.fq.gz,TA_01918_L2_1.fq.gz,TA_01915_L1_1.fq.gz,TA_01915_L2_1.fq.gz,TA_01916_L1_1.fq.gz,TA_01916_L2_1.fq.gz,TA_01916_L3_1.fq.gz,TA_01916_L4_1.fq.gz TA_01918_L1_2.fq.gz,TA_01918_L2_2.fq.gz,TA_01915_L1_2.fq.gz,TA_01915_L2_2.fq.gz,TA_01916_L1_2.fq.gz,TA_01916_L2_2.fq.gz,TA_01916_L3_2.fq.gz,TA_01916_L4_2.fq.gz
```

Next, we identify gene models using Cufflinks.

```bash
cufflinks -p 4 -o cufflinks_all accepted_hits_Can.7DL_Can.7DL.bam
```

### Analysis to merge multiple data sources to identify non-synonymous SNPs between Canthatch and mutants
The analyses above were integrated to identify SNPs between Canthatch and mutants for marker development, as well as uncover non-synonymous mutations in candidate genes. The pipeline is described below using a set of scripts that identify high quality SNPs based on either a stringent or relaxed set of parameters (see table below).

|Parameter set|Coverage|SNP freqency in wild-type|SNP frequency in mutant|
|:-----------:|:------:|:-----------------------:|:---------------------:|
|Stringent    |   10   |          <= 5%          |        >= 95%         |   
|Relaxed      |   10   |         <= 20%          |        >= 70%         |   

```bash
# Identify mutations that meet either strigent or relaxed sets of parameters
python mutant_flow_sorting_analysis.py

# Identify the types of identified mutations and export contigs containing a high quality SNP
python analyze_HQ_SNPs.py

# Link multiple sources of information into a single tab delimited file
python link_position_expression.py
```

The Python script `link_position_expression.py` will generate two FASTA files that identify transcripts with known or unknown strand orientation. Transdecoder is used to identfy the longest open reading frames and peptide sequence. Then, `link_position_expression.py` is ran again, where it can now integrate this information with the entire data set.

```bash
TransDecoder.LongOrfs -t Canthatch_cufflinks_ALL_strand.fa -S
TransDecoder.LongOrfs -t Canthatch_cufflinks_ALL_nostrand.fa

python link_position_expression.py
```

The output includes the following files:

**Canthatch_mutant_flow_sorting_Can.7DL_HQ_SNPs_NS1_position_cufflinks.txt**
**Canthatch_mutant_flow_sorting_Can.7DL_HQ_SNPs_NS2_position_cufflinks.txt**
All critical information is linked within these text files, including contigs, reference and mutant SNP, position in contig, position on chromosome 7D, and expressed genes in contig.

**IUPAC_KASP_contigs.fa**
Input file for generating KASP markers on EMS generated SNPs.

**Canthatch_mutant_flow_sorting_Can.7DL_HQ_SNPs_POI.fa**
Protein sequence for candidate genes (where SNPs lead to a non-synonymous change in sequence).


## Establishing the gene content and EMS-derived polymorphisms within the *Srs1* physical interval
The high identify found between Chinese Spring and Canthatch in the *Srs1* region permitted the use of a high strigency mapping approach to the IWGSC reference genome (v1.0). Two approaches were used:
   * Alignment of *de novo* assembled reads from Canthatch, NS1, and NS2
   * High strigent alignment of reads to the Chinese Spring physical interval

### High strigency alignment of Canthatch flow sorted chromosome arm reads to *Srs1* region
Our first step is to convert the Chinese Spring reference sequence into the Canthatch haplotype. The [QKgenome](https://github.com/matthewmoscou/QKgenome) pipeline was used, although read mapping was performed using `BBMap` with parameters that will allow a maximum of 2 SNPs and 1 InDel between individually mapped reads and the reference genome sequence.

```bash
./bbmap/bbsplit.sh ref=chr7D_621250000_622360000.fa in1=Can_gDNA_forward_paired.fq.gz in2=Can_gDNA_reverse_paired.fq.gz minid=0.985 maxindel=1 outm=Srs_Can_1.sam
./bbmap/bbsplit.sh ref=chr7D_621250000_622360000.fa in1=Can_2_gDNA_forward_paired.fq.gz in2=Can_2_gDNA_reverse_paired.fq.gz minid=0.985 maxindel=1 outm=Srs_Can_2.sam
```

Next, we convert the SAM files into BAM files, remove all unpaired reads (*i.e.* requiring that all reads map as pairs), and remove duplicate reads generated from PCR.

```bash
samtools view -f2 -Shub -o Srs_Can_1.bam Srs_Can_1.sam
samtools sort Srs_Can_1.bam Srs_Can_1.sorted
samtools rmdup Srs_Can_1.sorted.bam Srs_Can_1.sorted.rmdup.bam

samtools view -f2 -Shub -o Srs_Can_2.bam Srs_Can_2.sam
samtools sort Srs_Can_2.bam Srs_Can_2.sorted
samtools rmdup Srs_Can_2.sorted.bam Srs_Can_2.sorted.rmdup.bam

samtools merge Srs_Can.sorted.bam Srs_Can_1.sorted.bam Srs_Can_2.sorted.bam
```

### Conversion of the *Srs1* region from the Chinese Spring haplotype to the Canthatch haplotype
`QKgenome` was used to convert the Chinese Spring haplotype into the Canthatch haplotype. Parameters include a coverage of 10 reads and 80% alternate allele frequency.

```bash
bedtools genomecov -d -split -ibam Srs_Can.sorted.rmdup.bam > Srs_Can.sorted.rmdup.genomecov.txt

samtools index Srs_Can.sorted.rmdup.bam
samtools mpileup -f chr7D_621250000_622360000.fa -BQ0 Srs_Can.sorted.rmdup.bam > Srs_Can.sorted.rmdup.mpileup.txt

java -jar VarScan.v2.3.8.jar mpileup2snp Srs_Can.sorted.rmdup.mpileup.txt > Srs_Can.sorted.rmdup.mpileup2snp.txt
java -jar VarScan.v2.3.8.jar mpileup2indel Srs_Can.sorted.rmdup.mpileup.txt > Srs_Can.sorted.rmdup.mpileup2indel.txt

python QKgenome_conversion.py 10 80.0 chr7D_621250000_622360000.fa chr7D_621250000_622360000.gff3 Srs_Can.sorted.rmdup.mpileup2snp.txt Srs_Can.sorted.rmdup.mpileup2indel.txt Srs_Can.sorted.rmdup.genomecov.txt chr7D_621250000_622360000_Canthatch
```

### Defining the repeat sequence within the *Srs1* interval
The repetitive content of the *Srs1* interval was determined using `RepeatMasker`.

```bash
RepeatMasker -species monocotyledons chr7D_621250000_622360000_Canthatch.fa
```

### Realignment of Canthatch, NS1, and NS2 flow sorted chromosome arm reads to Canthatch converted *Srs1* region
```bash
./bbmap/bbsplit.sh ref=chr7D_621250000_622360000_Canthatch.fa in1=Can.7DL_gDNA_forward_paired.fq.gz in2=Can.7DL_gDNA_reverse_paired.fq.gz minid=0.985 maxindel=1 outm=SrsQK_Can_1.sam
./bbmap/bbsplit.sh ref=chr7D_621250000_622360000_Canthatch.fa in1=Can.7DL_2_gDNA_forward_paired.fq.gz in2=Can.7DL_2_gDNA_reverse_paired.fq.gz minid=0.985 maxindel=1 outm=SrsQK_Can_2.sam

./bbmap/bbsplit.sh ref=chr7D_621250000_622360000_Canthatch.fa in1=NS1M_gDNA_forward_paired.fq.gz in2=NS1M_gDNA_reverse_paired.fq.gz minid=0.985 maxindel=1 outm=SrsQK_NS1M_1.sam
./bbmap/bbsplit.sh ref=chr7D_621250000_622360000_Canthatch.fa in1=NS1M_2_gDNA_forward_paired.fq.gz in2=NS1M_2_gDNA_reverse_paired.fq.gz minid=0.985 maxindel=1 outm=SrsQK_NS1M_2.sam

./bbmap/bbsplit.sh ref=chr7D_621250000_622360000_Canthatch.fa in1=NS2N_gDNA_forward_paired.fq.gz in2=NS2N_gDNA_reverse_paired.fq.gz minid=0.985 maxindel=1 outm=SrsQK_NS2N_1.sam
./bbmap/bbsplit.sh ref=chr7D_621250000_622360000_Canthatch.fa in1=NS2N_2_gDNA_forward_paired.fq.gz in2=NS2N_2_gDNA_reverse_paired.fq.gz minid=0.985 maxindel=1 outm=SrsQK_NS2N_2.sam
```

```bash
# Canthatch
samtools view -f2 -Shub -o SrsQK_Can_1.bam SrsQK_Can_1.sam
samtools sort SrsQK_Can_1.bam SrsQK_Can_1.sorted
samtools rmdup SrsQK_Can_1.sorted.bam SrsQK_Can_1.sorted.rmdup.bam

samtools view -f2 -Shub -o SrsQK_Can_2.bam SrsQK_Can_2.sam
samtools sort SrsQK_Can_2.bam SrsQK_Can_2.sorted
samtools rmdup SrsQK_Can_2.sorted.bam SrsQK_Can_2.sorted.rmdup.bam

samtools merge SrsQK_Can.sorted.bam SrsQK_Can_1.sorted.bam SrsQK_Can_2.sorted.bam

# NS1
samtools view -f2 -Shub -o SrsQK_NS1M_1.bam SrsQK_NS1M_1.sam
samtools sort SrsQK_NS1M_1.bam SrsQK_NS1M_1.sorted
samtools rmdup SrsQK_NS1M_1.sorted.bam SrsQK_NS1M_1.sorted.rmdup.bam

samtools view -f2 -Shub -o SrsQK_NS1M_2.bam SrsQK_NS1M_2.sam
samtools sort SrsQK_NS1M_2.bam SrsQK_NS1M_2.sorted
samtools rmdup SrsQK_NS1M_2.sorted.bam SrsQK_NS1M_2.sorted.rmdup.bam

samtools merge SrsQK_NS1M.sorted.bam SrsQK_NS1M_1.sorted.bam SrsQK_NS1M_2.sorted.bam

# NS2
samtools view -f2 -Shub -o SrsQK_NS2N_1.bam SrsQK_NS2N_1.sam
samtools sort SrsQK_NS2N_1.bam SrsQK_NS2N_1.sorted
samtools rmdup SrsQK_NS2N_1.sorted.bam SrsQK_NS2N_1.sorted.rmdup.bam

samtools view -f2 -Shub -o SrsQK_NS2N_2.bam SrsQK_NS2N_2.sam
samtools sort SrsQK_NS2N_2.bam SrsQK_NS2N_2.sorted
samtools rmdup SrsQK_NS2N_2.sorted.bam SrsQK_NS2N_2.sorted.rmdup.bam

samtools merge SrsQK_NS2N.sorted.bam SrsQK_NS2N_1.sorted.bam SrsQK_NS2N_2.sorted.bam
```

```bash
bedtools genomecov -d -split -ibam SrsQK_NS1M.sorted.bam > SrsQK_NS1M.sorted.genomecov.txt
bedtools genomecov -d -split -ibam SrsQK_NS2N.sorted.bam > SrsQK_NS2N.sorted.genomecov.txt

samtools index SrsQK_NS1M.sorted.bam
samtools index SrsQK_NS2N.sorted.bam

samtools mpileup -f chr7D_621250000_622360000_Canthatch.fa -BQ0 SrsQK_NS1M.sorted.bam > SrsQK_NS1M.sorted.mpileup.txt
samtools mpileup -f chr7D_621250000_622360000_Canthatch.fa -BQ0 SrsQK_NS2N.sorted.bam > SrsQK_NS2N.sorted.mpileup.txt

java -jar VarScan.v2.3.8.jar mpileup2snp SrsQK_NS1M.sorted.mpileup.txt > SrsQK_NS1M.sorted.mpileup2snp.txt
java -jar VarScan.v2.3.8.jar mpileup2snp SrsQK_NS2N.sorted.mpileup.txt > SrsQK_NS2N.sorted.mpileup2snp.txt

java -jar VarScan.v2.3.8.jar mpileup2indel SrsQK_NS1M.sorted.mpileup.txt > SrsQK_NS1M.sorted.mpileup2indel.txt
java -jar VarScan.v2.3.8.jar mpileup2indel SrsQK_NS2N.sorted.mpileup.txt > SrsQK_NS2N.sorted.mpileup2indel.txt

python QKgenome_conversion.py 10 80.0 chr7D_621250000_622360000_Canthatch.fa chr7D_621250000_622360000_Canthatch.gff3 SrsQK_NS1M.sorted.mpileup2snp.txt SrsQK_NS1M.sorted.mpileup2indel.txt SrsQK_NS1M.sorted.genomecov.txt chr7D_621250000_622360000_Canthatch_NS1M
python QKgenome_conversion.py 10 80.0 chr7D_621250000_622360000_Canthatch.fa chr7D_621250000_622360000_Canthatch.gff3 SrsQK_NS2N.sorted.mpileup2snp.txt SrsQK_NS2N.sorted.mpileup2indel.txt SrsQK_NS2N.sorted.genomecov.txt chr7D_621250000_622360000_Canthatch_NS2N
```

**TODO** Need to perform hisat2 alignments on SrsQK reference dataset, alternatively, I could use existing data set as it is just a featureCounts output. No, better to use converted genome as it will ensure the maximum number of SNPs will map to the reference sequence.



### Gene expression of the homeologous *Med15a* gene family
Initial mapping of RNAseq reads was performed using `HISAT2`.

```bash
./hisat2-2.1.0/hisat2-build chr7ABD_Med15.fa chr7ABD_Med15
./hisat2-2.1.0/hisat2 --max-intronlen 20000 -p 16 -x chr7ABD_Med15 -1 TA_01918_L1_1.fq.gz,TA_01918_L2_1.fq.gz -2 TA_01918_L1_2.fq.gz,TA_01918_L2_2.fq.gz -S Med15_Can_RNAseq.sam
./hisat2-2.1.0/hisat2 --max-intronlen 20000 -p 16 -x chr7ABD_Med15 -1 TA_01915_L1_1.fq.gz,TA_01915_L2_1.fq.gz -2 TA_01915_L1_2.fq.gz,TA_01915_L2_2.fq.gz -S Med15_NS1M_RNAseq.sam
./hisat2-2.1.0/hisat2 --max-intronlen 20000 -p 16 -x chr7ABD_Med15 -1 TA_01916_L1_1.fq.gz,TA_01916_L2_1.fq.gz,TA_01916_L3_1.fq.gz,TA_01916_L4_1.fq.gz -2 TA_01916_L1_2.fq.gz,TA_01916_L2_2.fq.gz,TA_01916_L3_2.fq.gz,TA_01916_L4_2.fq.gz -S Med15_NS2N_RNAseq.sam
```

The challenge with the hexaploid genome is aligning reads that specifically map to individual sub-genomes. The initial `HISAT2` run was used to identify reads mapping to the Med15 7DL gene family. Next, we extract the aligned RNAseq reads and realign them using `tophat2` using highlt stringent parameters (`-N 0`, *i.e.* no differences between the read and the reference).

```bash
samtools view -F 4 -Shub -o Med15_Can_RNAseq.bam Med15_Can_RNAseq.sam
samtools sort Med15_Can_RNAseq.bam Med15_Can_RNAseq.sorted
samtools rmdup Med15_Can_RNAseq.sorted.bam Med15_Can_RNAseq.sorted.rmdup.bam

samtools view -F 4 -Shub -o Med15_NS1M_RNAseq.bam Med15_NS1M_RNAseq.sam
samtools sort Med15_NS1M_RNAseq.bam Med15_NS1M_RNAseq.sorted
samtools rmdup Med15_NS1M_RNAseq.sorted.bam Med15_NS1M_RNAseq.sorted.rmdup.bam

samtools view -F 4 -Shub -o Med15_NS2N_RNAseq.bam Med15_NS2N_RNAseq.sam
samtools sort Med15_NS2N_RNAseq.bam Med15_NS2N_RNAseq.sorted
samtools rmdup Med15_NS2N_RNAseq.sorted.bam Med15_NS2N_RNAseq.sorted.rmdup.bam

java -jar picard.jar SamToFastq I=Med15_Can_RNAseq.sorted.rmdup.pairs.bam FASTQ=Med15_Can_RNAseq_1.fastq SECOND_END_FASTQ=Med15_Can_RNAseq_2.fastq UNPAIRED_FASTQ=temp.fastq

bowtie2-build chr7ABD_Med15.fa chr7ABD_Med15
tophat2 -N 0 -p 4 --report-secondary-alignments chr7ABD_Med15 Med15_Can_RNAseq_1.fastq Med15_Can_RNAseq_2.fastq
featureCounts -T 4 -M -O -t exon -g ID -a chr7ABD.gff3 -o Med15_Can_RNAseq.sorted.rmdup.tophat2_readCounts.txt tophat_out/accepted_hits.bam
```

After this initial run, we found that several regions within all Med15.7L genes lacked mapped reads. This was due to SNPs between Canthatch and Chinese Spring alleles. Manual curation using `HISAT2` aligned RNAseq reads was used to convert the coding sequences for *Med15.7AL* and *Med15.7BL* from Chinese Spring to Canthatch. The `QKgenome` pipeline was used to curate *Med15.7DL*. `HISAT2` was performed again to identify reads mapping to the Med15 7DL gene family. We extract the aligned RNAseq reads and realigned them using `tophat2` using highlt stringent parameters (`-N 0`, *i.e.* no differences between the read and the reference).

```bash
./hisat2-2.1.0/hisat2-build Med15_7L_Canthatch.fa Med15_7L_Canthatch
./hisat2-2.1.0/hisat2 --max-intronlen 20000 -p 16 -x Med15_7L_Canthatch -1 TA_01918_L1_1.fq.gz,TA_01918_L2_1.fq.gz -2 TA_01918_L1_2.fq.gz,TA_01918_L2_2.fq.gz -S Med15_Can_Can_RNAseq.sam
./hisat2-2.1.0/hisat2 --max-intronlen 20000 -p 16 -x Med15_7L_Canthatch -1 TA_01915_L1_1.fq.gz,TA_01915_L2_1.fq.gz -2 TA_01915_L1_2.fq.gz,TA_01915_L2_2.fq.gz -S Med15_Can_NS1M_RNAseq.sam
./hisat2-2.1.0/hisat2 --max-intronlen 20000 -p 16 -x Med15_7L_Canthatch -1 TA_01916_L1_1.fq.gz,TA_01916_L2_1.fq.gz,TA_01916_L3_1.fq.gz,TA_01916_L4_1.fq.gz -2 TA_01916_L1_2.fq.gz,TA_01916_L2_2.fq.gz,TA_01916_L3_2.fq.gz,TA_01916_L4_2.fq.gz -S Med15_Can_NS2N_RNAseq.sam

samtools view -F 4 -Shub -o Med15_Can_Can_RNAseq.bam Med15_Can_Can_RNAseq.sam
samtools sort Med15_Can_Can_RNAseq.bam Med15_Can_Can_RNAseq.sorted
samtools rmdup Med15_Can_Can_RNAseq.sorted.bam Med15_Can_Can_RNAseq.sorted.rmdup.bam

samtools view -F 4 -Shub -o Med15_Can_NS1M_RNAseq.bam Med15_Can_NS1M_RNAseq.sam
samtools sort Med15_Can_NS1M_RNAseq.bam Med15_Can_NS1M_RNAseq.sorted
samtools rmdup Med15_Can_NS1M_RNAseq.sorted.bam Med15_Can_NS1M_RNAseq.sorted.rmdup.bam

samtools view -F 4 -Shub -o Med15_Can_NS2N_RNAseq.bam Med15_Can_NS2N_RNAseq.sam
samtools sort Med15_Can_NS2N_RNAseq.bam Med15_Can_NS2N_RNAseq.sorted
samtools rmdup Med15_Can_NS2N_RNAseq.sorted.bam Med15_Can_NS2N_RNAseq.sorted.rmdup.bam

java -jar picard.jar SamToFastq I=Med15_Can_Can_RNAseq.sorted.rmdup.pairs.bam FASTQ=Med15_Can_Can_RNAseq_1.fastq SECOND_END_FASTQ=Med15_Can_Can_RNAseq_2.fastq UNPAIRED_FASTQ=temp.fastq

bowtie2-build Med15_7L_Canthatch.fa Med15_7L_Canthatch
tophat2 -N 0 -p 4 --report-secondary-alignments Med15_7L_Canthatch Med15_Can_Can_RNAseq_1.fastq Med15_Can_Can_RNAseq_2.fastq
```

### Does loss of *Srs1* lead to increased expression of *Lr34*?

```bash
cat chr7A_50000000_50100000.fa chr7D_47400000_47430000.fa > Lr34.fa

./hisat2-2.1.0/hisat2-build Lr34.fa Lr34

./hisat2-2.1.0/hisat2 --max-intronlen 20000 -p 16 -x Lr34 -1 TA_01918_L1_1.fq.gz,TA_01918_L2_1.fq.gz -2 TA_01918_L1_2.fq.gz,TA_01918_L2_2.fq.gz -S Lr34_Can_RNAseq.sam

samtools view -f2 -Shub -o Lr34_Can_RNAseq.bam Lr34_Can_RNAseq.sam
samtools sort Lr34_Can_RNAseq.bam Lr34_Can_RNAseq.sorted
samtools rmdup Lr34_Can_RNAseq.sorted.bam Lr34_Can_RNAseq.sorted.rmdup.bam

./hisat2-2.1.0/hisat2 --max-intronlen 20000 -p 16 -x Lr34 -1 TA_01915_L1_1.fq.gz,TA_01915_L2_1.fq.gz -2 TA_01915_L1_2.fq.gz,TA_01915_L2_2.fq.gz -S Lr34_NS1M_RNAseq.sam

samtools view -f2 -Shub -o Lr34_NS1M_RNAseq.bam Lr34_NS1M_RNAseq.sam
samtools sort Lr34_NS1M_RNAseq.bam Lr34_NS1M_RNAseq.sorted
samtools rmdup Lr34_NS1M_RNAseq.sorted.bam Lr34_NS1M_RNAseq.sorted.rmdup.bam

./hisat2-2.1.0/hisat2 --max-intronlen 20000 -p 16 -x Lr34 -1 TA_01916_L1_1.fq.gz,TA_01916_L2_1.fq.gz,TA_01916_L3_1.fq.gz,TA_01916_L4_1.fq.gz -2 TA_01916_L1_2.fq.gz,TA_01916_L2_2.fq.gz,TA_01916_L3_2.fq.gz,TA_01916_L4_2.fq.gz -S Lr34_NS2N_RNAseq.sam

samtools view -f2 -Shub -o Lr34_NS2N_RNAseq.bam Lr34_NS2N_RNAseq.sam
samtools sort Lr34_NS2N_RNAseq.bam Lr34_NS2N_RNAseq.sorted
samtools rmdup Lr34_NS2N_RNAseq.sorted.bam Lr34_NS2N_RNAseq.sorted.rmdup.bam
```

Do reads cross map between homeologs?


**TODO**
* Need to run BBmap on the Canthatch converted region of Chinese Spring, DONE
* Need to run RNAseq alignments on everything again, DONE
* Need to setup a Github with all code, DONE
* Need to generate a figure to describe read coverage in the *Srs1* region
* Need to incorporate read coverage information form already generated figure

Sequence flanking mutations in the *Srs1* interval identified using aligned *de novo* assemblies.

```
>NS2_158661
CGCCTCCGCCTCCGTGTACATGGTGACGCGAGGTGAGTCAAACTCAATGTCACTTGGGATGACGGCCTCGGCACCGTACACGAGGAAGAATGGAGTGAAGYCGGTTGATTTGTTCGAAGTGGTACGGAGGCTCCAGAGGACGGCCGGCAGCTCGTCGATCCAGCAGTCGGCCGACCGCTCCAGTGGCTCGACCAGTCGGGG
>NS2_390486
CGGGAGCTATTTGTTGCACATCCCAGTTATTTTTCACATTTAGTATATACATACATAATGTAAGTTAGGGCAAGTTCTTTAGGGCAGCTTAAAAATTAAGYTGCCTCTCCCCAGCTTAAAAAATAATCCGTCCTCTAAGTTTGTTGAGACTTCTAAATTAGCTATTTCATAACTAGTTTAGAAGCCCCAATGAACAAGAGA
>NS1_393242
ATCGATGTCTTCGCTGGTGGTGGTGCGGTGAAGAACGTCGGCGGGGTCGTAGTGAAGGAGGCTGGCGGTGGCATGGCCATTGTCGCCCTCGGCTTGTTCGYTGCCGCCCGGCATGCCGCCTCCATTCTCTCCCACGGCGCCGGCGGGGATGGGACACGTGCTTAGCGGTGCTCCTGGCCTGGGTGTCACCGCCGGCTGAGG
>NS1_394896
GCCACAACAATTGGGCTCTCAGGCAAACATGTCAAGTTTACAGCAGCAGCAACAAAATCAACAGCAGCAGCGGATGCATATGCTACAAATGAAAGCTCAGYAAACGCAGCAACAACAGCATGCTCAACAACCACCAATGGGTTTGATGCAACCTCAGTCCCAACACAACCAACTTCAGCAATCGCAGCAACATCTTATGTC
>NS2_397520
CTCCATGCAAGCAAATGCAAGTTCATTACAGCAGCTGAAGCAGCAACAGCAGGATCATCATATGATGCAGAGTCAGCAAATGAAGCGTCAGATGTTTCAGYAGTACCAGCAAAAGCAACAAATGCTTCAACAGCAGTTCCCAATACAGCAACAATTACAGAAACAGCAGCAAGTACAGATGCAGGTTCCACAGCTTCATGC
>NS1_483619
ACACTCCGAGTTTGCGCAGCTGGGTGAGTTTGTTGAGCTCTTTGACGATGTCCTTCCCGCCTGAAGCACCAATGTTGACAACACCGAGCGTGTGCAATGCYGTCAGTTCACCAATCCTGCTTGGCACTACAACACCAACTAGGCGACGACATCTGCGGAACTCCGGCAAGAAGCTAGTTGAGGCATCTGGCGTTGATGCTG
>NS2_569625
TATAGGTACATCCATTTTTGGACAAATGAAAGCCAAGTATTTTGGAACGGAGGGAGTAGAACTGAACTTTTTGGTTTTGTGAATGAAATTGGGGAGATATYGCTTTCTTCAACAAATAAATGATCACCCAGAACAAATACATTTTAAAATTATTTTAATATGGAGTACCGCAATATAATCAGTAGATTAGGTTGGGTTGCG
>NS1_1093857
TGCTCCTGTTGAACTAAAACTGTTTCAAATAAATTTGATAAAAATGCCACAAAATACTGTGAGACCCTCTGATGTGTATAGATCACTTTTAGGTAAAAATYCATCTTGTTCATTTGAAAAGTTTGGCCACAACATTCTACCATGGCATTACGAAAGTTTGCCACAGCATTACTAAAGTTTGCCACAACATTACTAAAGTTC
>NS1_1100289
GTAGTTGGTTAGCTGCATGTGCGGTTGTAAGCTAGTTTGATGGGTGAGTGGCATAGGAGACGGGCAGGTTGTTGCCCTGGTAGTGGTAGCGTGTGCATGCRTGTGCCATTTAAACAGCCGTGGTGCTGTGTAGAGAGCGCGGAAGCCAGAAGAAGAGATGTAGTCTCCAATGGTCGGTGCGTGTGTGTTCGGAGAGAAGCT
```


### *Lr34*

|Gene           |Chromosome|Start   |End     |
|:-------------:|:--------:|:------:|:------:|
|*Lr34* homeolog|    7A    |50000000|50100000|
|*Lr34*         |    7D    |47400000|47430000|

Read coverage of Lr34 and homeolog
Gene	Canthatch	NS1	NS2
Lr34	1599	864	792
Reads	41000991	29237025	27664945


## Phylogenetic and diversity analysis of *TaMed15* gene family

### Natural variation in *Med15* in *Aegilops tauschii*

```bash
~/src/hisat2-2.0.5/hisat2-build chr7D_Med15.fa chr7D_Med15

~/src/hisat2-2.0.5/hisat2 --max-intronlen 20000 -p 4 -x chr7D_Med15 -1 ~/temp/DRR058959_1.fastq.gz -2 ~/temp/DRR058959_2.fastq.gz -S AetMed15_AT76_RNAseq.sam
samtools view -f 2 -Shub AetMed15_AT76_RNAseq.sam > AetMed15_AT76_RNAseq.bam
samtools sort AetMed15_AT76_RNAseq.bam AetMed15_AT76_RNAseq_sorted
samtools rmdup AetMed15_AT76_RNAseq_sorted.bam AetMed15_AT76_RNAseq_sorted_rmdup.bam

~/src/hisat2-2.0.5/hisat2 --max-intronlen 20000 -p 4 -x chr7D_Med15 -1 ~/temp/DRR058960_1.fastq.gz -2 ~/temp/DRR058960_2.fastq.gz -S AetMed15_KU-2003_RNAseq.sam
samtools view -f 2 -Shub AetMed15_KU-2003_RNAseq.sam > AetMed15_KU-2003_RNAseq.bam
samtools sort AetMed15_KU-2003_RNAseq.bam AetMed15_KU-2003_RNAseq_sorted
samtools rmdup AetMed15_KU-2003_RNAseq_sorted.bam AetMed15_KU-2003_RNAseq_sorted_rmdup.bam

~/src/hisat2-2.0.5/hisat2 --max-intronlen 20000 -p 4 -x chr7D_Med15 -1 ~/temp/DRR058961_1.fastq.gz -2 ~/temp/DRR058961_2.fastq.gz -S AetMed15_KU-2025_RNAseq.sam
samtools view -f 2 -Shub AetMed15_KU-2025_RNAseq.sam > AetMed15_KU-2025_RNAseq.bam
samtools sort AetMed15_KU-2025_RNAseq.bam AetMed15_KU-2025_RNAseq_sorted
samtools rmdup AetMed15_KU-2025_RNAseq_sorted.bam AetMed15_KU-2025_RNAseq_sorted_rmdup.bam

~/src/hisat2-2.0.5/hisat2 --max-intronlen 20000 -p 4 -x chr7D_Med15 -1 ~/temp/DRR058962_1.fastq.gz -2 ~/temp/DRR058962_2.fastq.gz -S AetMed15_KU-2075_RNAseq.sam
samtools view -f 2 -Shub AetMed15_KU-2075_RNAseq.sam > AetMed15_KU-2075_RNAseq.bam
samtools sort AetMed15_KU-2075_RNAseq.bam AetMed15_KU-2075_RNAseq_sorted
samtools rmdup AetMed15_KU-2075_RNAseq_sorted.bam AetMed15_KU-2075_RNAseq_sorted_rmdup.bam

~/src/hisat2-2.0.5/hisat2 --max-intronlen 20000 -p 4 -x chr7D_Med15 -1 ~/temp/DRR058963_1.fastq.gz -2 ~/temp/DRR058963_2.fastq.gz -S AetMed15_KU-2078_RNAseq.sam
samtools view -f 2 -Shub AetMed15_KU-2078_RNAseq.sam > AetMed15_KU-2078_RNAseq.bam
samtools sort AetMed15_KU-2078_RNAseq.bam AetMed15_KU-2078_RNAseq_sorted
samtools rmdup AetMed15_KU-2078_RNAseq_sorted.bam AetMed15_KU-2078_RNAseq_sorted_rmdup.bam

~/src/hisat2-2.0.5/hisat2 --max-intronlen 20000 -p 4 -x chr7D_Med15 -1 ~/temp/DRR058964_1.fastq.gz -2 ~/temp/DRR058964_2.fastq.gz -S AetMed15_KU-2087_RNAseq.sam
samtools view -f 2 -Shub AetMed15_KU-2087_RNAseq.sam > AetMed15_KU-2087_RNAseq.bam
samtools sort AetMed15_KU-2087_RNAseq.bam AetMed15_KU-2087_RNAseq_sorted
samtools rmdup AetMed15_KU-2087_RNAseq_sorted.bam AetMed15_KU-2087_RNAseq_sorted_rmdup.bam

~/src/hisat2-2.0.5/hisat2 --max-intronlen 20000 -p 4 -x chr7D_Med15 -1 ~/temp/DRR058965_1.fastq.gz -2 ~/temp/DRR058965_2.fastq.gz -S AetMed15_KU-2093_RNAseq.sam
samtools view -f 2 -Shub AetMed15_KU-2093_RNAseq.sam > AetMed15_KU-2093_RNAseq.bam
samtools sort AetMed15_KU-2093_RNAseq.bam AetMed15_KU-2093_RNAseq_sorted
samtools rmdup AetMed15_KU-2093_RNAseq_sorted.bam AetMed15_KU-2093_RNAseq_sorted_rmdup.bam

~/src/hisat2-2.0.5/hisat2 --max-intronlen 20000 -p 4 -x chr7D_Med15 -1 ~/temp/DRR058966_1.fastq.gz -2 ~/temp/DRR058966_2.fastq.gz -S AetMed15_KU-2124_RNAseq.sam
samtools view -f 2 -Shub AetMed15_KU-2124_RNAseq.sam > AetMed15_KU-2124_RNAseq.bam
samtools sort AetMed15_KU-2124_RNAseq.bam AetMed15_KU-2124_RNAseq_sorted
samtools rmdup AetMed15_KU-2124_RNAseq_sorted.bam AetMed15_KU-2124_RNAseq_sorted_rmdup.bam

~/src/hisat2-2.0.5/hisat2 --max-intronlen 20000 -p 4 -x chr7D_Med15 -1 ~/temp/DRR058967_1.fastq.gz -2 ~/temp/DRR058967_2.fastq.gz -S AetMed15_KU-2627_RNAseq.sam
samtools view -f 2 -Shub AetMed15_KU-2627_RNAseq.sam > AetMed15_KU-2627_RNAseq.bam
samtools sort AetMed15_KU-2627_RNAseq.bam AetMed15_KU-2627_RNAseq_sorted
samtools rmdup AetMed15_KU-2627_RNAseq_sorted.bam AetMed15_KU-2627_RNAseq_sorted_rmdup.bam

~/src/hisat2-2.0.5/hisat2 --max-intronlen 20000 -p 4 -x chr7D_Med15 -1 ~/temp/DRR058968_1.fastq.gz -2 ~/temp/DRR058968_2.fastq.gz -S AetMed15_PI499262_RNAseq.sam
samtools view -f 2 -Shub AetMed15_PI499262_RNAseq.sam > AetMed15_PI499262_RNAseq.bam
samtools sort AetMed15_PI499262_RNAseq.bam AetMed15_PI499262_RNAseq_sorted
samtools rmdup AetMed15_PI499262_RNAseq_sorted.bam AetMed15_PI499262_RNAseq_sorted_rmdup.bam
```

### RNAseq coverage over *TaMed15.7DL*
bedtools genomecov -d -split -ibam Med15_Can_RNAseq.sorted.rmdup.bam > Med15_Can.7DL.sorted.rmdup.genomecov.txt
bedtools genomecov -d -split -ibam Med15_NS1M_RNAseq.sorted.rmdup.bam > Med15_NS1M.7DL.sorted.rmdup.genomecov.txt
bedtools genomecov -d -split -ibam Med15_NS2N_RNAseq.sorted.rmdup.bam > Med15_NS2N.7DL.sorted.rmdup.genomecov.txt

Note: In Trinity assembly of Canthatch, all three homeologs of collapse into a single RNAseq contig

```R
library(ggplot2)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

setwd("~/Desktop/bioinformatics/Canthatch/mutchromseq/physical_analysis_v2/Canthatch_Med15_RNAseq/")

data = read.table(file="Canthatch_Med15_RNAseq.txt", header=T)
data = data.frame(data)

ggplot(data, aes(position, normexpression, color=genotype)) + geom_area(aes(fill=genotype)) + scale_fill_manual(values=cbPalette) + scale_color_manual(values=cbPalette)

ggsave("Med15_RNAseq.eps", width=10, height=4)
```
