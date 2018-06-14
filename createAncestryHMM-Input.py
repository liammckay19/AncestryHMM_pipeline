import sys
import re
import io
import configparser



class VCFreader :
    ''' 
    Define objects to read VCF files.
    instantiation: 
    thisReader = VCFreader("NewChr1.snp.vcf")
    usage:
    for line in thisReader.readVCF():
        #### print (line)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
        self.doOpen()
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
    def readVCF (self):
        ''' Read an entire VCF record and return the sequence header/sequence'''
        line = ''
        with self.doOpen() as fileH:
            for line in fileH:
                while not line.startswith('##'):
                    yield line
                    line = fileH.readline()
def pretty(d, indent=0): #https://stackoverflow.com/questions/3229419/how-to-pretty-print-nested-dictionaries
   for key, value in d.items():
      print('\t' * indent + str(key))
      if isinstance(value, dict):
         pretty(value, indent+1)
      else:
         print('\t' * (indent+1) + str(value))
def file_len(fname): #https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
class FreqDistanceCalculator :
    def __init__(self):
        config = configparser.ConfigParser()
        config.read('config.ini')
        print('Your config file input:\n')
        pretty(dict(config['DEFAULT']),0)
        self.c = float(config['DEFAULT']['allelefreq cutoff'])
        self.locusLength = int(config['DEFAULT']['min locus distance'])
        self.refPanelNumber = int(config['DEFAULT']['number of reference panels'])
        self.recombinationRate = float(config['DEFAULT']['recombination_rate'])
        self.minChrom = int(config['DEFAULT']['minChrom'])
        self.fileName = str(config['DEFAULT']['filename'])
        self.populationNames = config['DEFAULT']['populationNames'].split(',')
        self.lineCount = 0

    def createAncestryHMMInputFile(self, outFileName):
        '''
        Create outfile for AncestryHMM Input 
        Instantiation:
        myOutput = createAncestryHMMInputFile('output.tsv')
        '''
        # attempt to read file for input; if unable, exit program
        try:
            self.read = VCFreader(self.fileName)
        except FileNotFoundError as err:
            print('{} \ncannot open'.format(err), self.fileName)
            exit()
        else:
            print(self.fileName, 'to be read')

        numAboveC = 0
        outFile = open(outFileName, 'w')
        currentLocus = 0

        for line in self.read.readVCF(): # read line by line of VCF
            reference1_count = 0;reference2_count = 0;reference1AlleleA = 0;
            reference2AlleleA = 0;sampleAlleleA=0;sample_count=0
            chromosomeNumber = 0
            position = 0
            recombinationFrequency = 0.0

            lineList = list(line.split('\t'))
            pastChromosomeNumber = lineList[0].split('[A-Z][a-z]')[0]

            alleleCountRefAltList = []
            readCountRefAltList = []
            alleleFreqList = []
            if lineList[0] == '': # reached EOF
                break
            if lineList[0].startswith('##'): # comment of VCF
                continue
            elif lineList[0].startswith('#'): # header of VCF
                populationColumnsList = []
                for populationName in self.populationNames:
                    populationColumnsList.append([ i for i, column in enumerate(lineList) if re.search('^'+populationName, column) ])
                
                if self.refPanelNumber > len(populationColumnsList):
                    print('Count of references in config is larger than columns in vcf: ref='+str(self.refPanelNumber)+' col='+str(len(populationColumnsList))+'\nExiting')
                    exit()
                # check if all columns have been read
                p = 0
                initPopulationCount=len(populationColumnsList[0])
                for populationColumn in populationColumnsList :
                    if populationColumn == [] :
                        print("Could not find population "+self.populationNames[p]+'\nContinuing')
                    p += 1
                # ref lists [0] [1] ... sample list [2] [3] ...
            else: # actual data of VCF
               
                for refpopulation in range(0,self.refPanelNumber-1): # compute allele count for at least two reference populations
                    reference1_count=0
                    reference1AlleleA=0
                    bothRefList = list(zip(populationColumnsList[refpopulation+1],populationColumnsList[refpopulation]))
                    # print(populationColumnsList)
                    # print((populationColumnsList[refpopulation+1],populationColumnsList[refpopulation]))
                    # print(bothRefList)

                    for column1,column2 in bothRefList:
                        # get 1st chrom - 0/1
                        if (lineList[column1][:1] != '.'): 
                            if int(lineList[column1][:1])<=1:
                                reference1AlleleA += int(lineList[column1][:1])
                                reference1_count += 1
                        # get 2nd chrom - 0/1
                        if (lineList[column1][2:3] != '.'): 
                            if int(lineList[column1][2:3])<=1:
                                reference1AlleleA += int(lineList[column1][2:3])
                                reference1_count += 1

                        # compare to next population
                        if (lineList[column2][:1] != '.'):
                            if int(lineList[column2][:1])<=1:
                                reference2AlleleA += int(lineList[column2][:1])
                                reference2_count += 1
                        
                        if (lineList[column2][2:3] != '.'):
                            if int(lineList[column2][2:3])<=1:
                                reference2AlleleA += int(lineList[column2][2:3])
                                reference2_count += 1
                                    

                    if reference1_count == 0:
                        reference1_count = 1 # so theres no dividing by zero
                    if reference2_count == 0:
                        reference2_count = 1 # so theres no dividing by zero
                    chrom_number1 = reference1_count
                    chrom_number2 = reference2_count
                    freq1=reference1AlleleA/reference1_count
                    freq2=reference2AlleleA/reference2_count
                    diff12=abs(freq1-freq2)
                    if chrom_number1 < self.minChrom or chrom_number2 < self.minChrom: # make parameter the code for min chromosome number
                        continue # returns back to the top of the for loop
                    elif diff12 < self.c:
                        continue
                    else : 
                        chromosomeNumber= lineList[0].split('[A-Z][a-z]')[0]
                        if chromosomeNumber != pastChromosomeNumber:
                            currentLocus = 0
                        distanceToNextLocus = abs(currentLocus - int(lineList[1]))
                        
                        if (distanceToNextLocus > self.locusLength): # if passed all thresholds specified from config.ini
                            # for printing
                            alleleCountRefAltList.append(abs(reference1_count-reference1AlleleA))
                            alleleCountRefAltList.append(reference1AlleleA)
                            alleleCountRefAltList.append(abs(reference2_count-reference2AlleleA))
                            alleleCountRefAltList.append(reference2AlleleA)
                            for samplePopulation in range(self.refPanelNumber,len(populationColumnsList)): # find read count for each sample 
                                for column in range(0,len(populationColumnsList[samplePopulation])):
                                    try:
                                        samPopCol = populationColumnsList[samplePopulation][column]
                                    except IndexError as err:
                                        print("Failed to load sample panel entry for (col:row)\t"+str(samplePopulation)+':'+str(column))
                                        continue
                                    # GT:AD:DP:GQ:PGT:PID:PL

                                    readCounts = lineList[samPopCol].split(':')[1].split(',')
                                    for read in readCounts :
                                        readCountRefAltList.append(read)
                    
                            currentLocus = int(lineList[1])
                            recombinationFrequency=distanceToNextLocus*self.recombinationRate 
                            # compute recombination frequency -- simple formula = (distance*recombinationRate)
                            alleleCountRefAltListString ='\t'.join(str(e) for e in alleleCountRefAltList)
                            readCountRefAltListString='\t'.join(str(e) for e in readCountRefAltList)
                            outFile.write("{}\t{}\t{}\t{}\t{}\n".format(\
                                        chromosomeNumber,currentLocus,alleleCountRefAltListString,recombinationFrequency,readCountRefAltListString))
                            self.lineCount +=1

def main():
    usage="\nUsage: python createAncestryHMM-Input.py  \n"\
                    +"Type python createAncestryHMM-Input.py -help to get descriptions of arguments\n "\
                    +"You must have edited the config.ini file to read in parameters"
    helpString = "[allelefreq cutoff]\n\t (Float) cutoff value for reference panel allele frequency calculation\n"\
                +"\n[min locus distance]\n\t (Integer) minimum distance between each allele locus\n" \
                +"\n[number of reference panels]\n\t (Integer) number of reference panel columns in the VCF.\n" \
                +"\n[recombination_rate]\n\t (Float) estimated recombination rate for recombination probability for Ancestry_HMM input.\n\t Average recombinations per base pairs\n" \
                +"\n[minChrom]\n\t (Integer) the minimum amount of chromosomes that must be present in the reference panel alleles to make it through the threshold\n" \
                +"\n[filename]\n\t (String no quotes) Name of VCF file on local machine\n" \
                +"\n[refPanel names(requires 2 names)]\n\t (at least 2 Strings no quotes) Names of reference panels in the VCF file.\n\t Reference panels should be named like guanaco0 guanaco1 guanaco2 etc.\n\t Example argument: guanaco\n" \
                +"\n[sample names]\n\t (at least 2 Strings no quotes) Names of sample panels to be run in Ancestry_HMM. \n\t Should be named like llama1 llama2 llama3 etc.\n\t Example argument: llama\n"
    numAboveC = 0
    
    
    if (len(sys.argv) == 2):
        if '-help' in sys.argv[1]:
            print(helpString)
            exit()
    vcfInputProcessor = FreqDistanceCalculator()
    fileOut= input("Name of output file:")
    vcfInputProcessor.createAncestryHMMInputFile(fileOut)    
    print('Created output file: '+fileOut+" \n(for Ancestry_HMM input)\n"+str(vcfInputProcessor.lineCount)+" Lines")

main()
