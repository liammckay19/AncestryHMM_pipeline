# AncestryHMM_pipeline
Undergraduate research Adrian Salguero and Liam McKay under guidance by Russ Corbett-Detig PhD<br>

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
