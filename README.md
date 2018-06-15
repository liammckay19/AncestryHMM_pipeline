# AncestryHMM_pipeline
Undergraduate research Adrian Salguero and Liam McKay under guidance by Russ Corbett-Detig PhD<br>

# What This Does
- Takes a VCF (Variant Call Format) SNP data file (NGS read pileup data/genotype data/...)
- Converts it into a file (like example.panel) for input into russcd/Ancestry_HMM

| 1 | 1105366 | . | T | C | . | PASS | AA=T;AC=4;AN=114;DP=3251 | GT:DP | ./.:0 | ./.:0 | 0/0:107 | ./.:0 | ./.:0 | 0/0:25 | 0/0:30 | 0/0:31 | 0/0:57 | 0/0:69 | 0/0:53 | 0/0:225 | ./.:0 | 0/0:6 | 0/0:79 | ./.:0 | 0/0:110 | 0/0:79 | ./.:0 | ./.:0 | ./.:0 | 0/0:43 | 1/0:54 | 0/0:7 | 0/0:89 | 0/0:87 | 0/0:98 | 0/0:83 | ./.:0 | 0/0:62 | 0/0:1 | 0/0:4 | ./.:0 | 0/0:97 | ./.:0 | 0/0:115 | ./.:0 | 0/0:77 | 0/0:8 | 0/0:63 | ./.:0 | 0/0:92 | ./.:0 | 0/0:1 | 0/0:1 | ./.:0 | ./.:0 | ./.:0 | ./.:0 | 0/0:76 | ./.:0 | ./.:0 | ./.:0 | 0/0:41 | 0/0:35 | 1/0:135 | ./.:0 | 1/0:116 | 0/0:6 | ./.:0 | 0/0:147 | ./.:0 | ./.:0 | 0/0:4 | 0/0:40 | 1/0:23 | ./.:0 | 0/0:1 | 0/0:2 | ./.:0 | 0/0:7 | 0/0:1 | 0/0:90 | 0/0:49 | ./.:0 | 0/0:6 | ./.:0 | 0/0:82 | 0/0:31 | 0/0:7 | 0/0:9 | 0/0:7 | ./.:0 | ./.:0 | ./.:0 | 0/0:176 | 0/0:3 | 0/0:81 | 0/0:67 | 0/0:156 | 
<br> Converted to: <br>
| 1 | 20916748 | 1 | 0  | 0 | 4  | 0.20916748000000002 | 4   | 0 | 0 | 0 | 6  | 1  | 0  | 0 | 1  | 0  | 1  | 0 | 4  | 3  | 0 | 0  | 0 | 0   | 1   | 0  | 3 | 0  | 0  | 0 | 1 | 44  | 3 | 11  | 21 | 53  | 
|   | 0        | 0 | 39 | 0 | 26 | 29                  | 0   | 0 | 0 | 3 | 1  | 2  | 36 | 0 | 19 | 27 | 0  | 2 | 0  | 0  | 0 | 5  | 0 | 29  | 0   |    |   |    |    |   |   |     |   |     |    |     | 
| 1 | 37098064 | 1 | 1  | 8 | 0  | 0.16181316          | 344 | 5 | 2 | 1 | 76 | 49 | 1  | 2 | 14 | 0  | 24 | 2 | 87 | 57 | 8 | 12 | 9 | 121 | 109 | 35 | 3 | 53 | 13 | 8 | 8 | 504 | 4 | 215 | 94 | 291 | 

# Dependencies
- Install Anacaonda for python 3 https://anaconda.org/anaconda/python
- (To check what version of python you have type `python --version`)

# Quickstart:
- set config.ini with these parameters:<br>
allelefreq cutoff = 0.5<br>
min locus distance = 10<br>
number of reference panels = 4<br>
recombination_rate = 1e-8<br>
minChrom = 1<br>
filename = vcfDownloadTestData/CEU.exon.2010_03.genotypes.vcf<br>
populationNames = NA121,NA122,NA123,NA124,NA125,NA127,NA128,NA10,NA11,NA120<br>
- Then type in command line: <br>
`python createAncestryHMM-Input.py`<br>
- Type a name for output at prompt
- Creates a tsv file with 136 lines 
- To see it in terminal type:<br>
`cat <filename>`<br>

# Documentation
To use this program, edit config.ini for your input file in VCF format<br>
<br>[allelefreq cutoff]<br> (Float) cutoff value for reference panel allele frequency calculation<br>
<br>[min locus distance]<br> (Integer) minimum distance between each allele locus<br> 
<br>[number of reference panels]<br> (Integer) number of reference panel species in the VCF<br> 
<br>[recombination_rate]<br> (Float) estimated recombination rate for recombination probability for Ancestry_HMM input.<br> Average recombinations per base pairs<br> 
<br>[minChrom]<br> (Integer) the minimum amount of chromosomes that must be present in the reference panel alleles to make it through the threshold<br> 
<br>[filename]<br> (String no quotes) Name of VCF file on local machine<br> 
<br>[refPanel names(requires 2 names)]<br> (at least 2 Strings no quotes) Names of reference panels in the VCF file.<br> Reference panels should be named like guanaco0 guanaco1 guanaco2 etc.<br> Example argument: guanaco<br> 
<br>[sample names]<br> (at least 2 Strings no quotes) Names of sample panels to be run in Ancestry_HMM. <br> Should be named like llama1 llama2 llama3 etc.<br> Example argument: llama<br>
