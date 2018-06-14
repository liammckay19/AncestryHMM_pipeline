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

class FreqDistanceCalculator :
    def __init__(self):
        config = configparser.ConfigParser()
        config.read('config.ini')
        print('Your config file input:\n',dict(config['DEFAULT'].items()))
        usage="\nUsage: python createAncestryHMM-Input.py  \n"\
                    +"Type python createAncestryHMM-Input.py -help to get descriptions of arguments\n "\
                    +"You must have edited the config.ini file to read in parameters"
        helpString = "[allelefreq cutoff]\n\t (Float) cutoff value for reference panel allele frequency calculation\n"\
                    +"\n[min locus distance]\n\t (Integer) minimum distance between each allele locus\n" \
                    +"\n[number of reference panels]\n\t (Integer) number of reference panel species in the VCF\n" \
                    +"\n[recombination_rate]\n\t (Float) estimated recombination rate for recombination probability for Ancestry_HMM input.\n\t Average recombinations per base pairs\n" \
                    +"\n[minChrom]\n\t (Integer) the minimum amount of chromosomes that must be present in the reference panel alleles to make it through the threshold\n" \
                    +"\n[filename]\n\t (String no quotes) Name of VCF file on local machine\n" \
                    +"\n[refPanel names(requires 2 names)]\n\t (at least 2 Strings no quotes) Names of reference panels in the VCF file.\n\t Reference panels should be named like guanaco0 guanaco1 guanaco2 etc.\n\t Example argument: guanaco\n" \
                    +"\n[sample names]\n\t (at least 2 Strings no quotes) Names of sample panels to be run in Ancestry_HMM. \n\t Should be named like llama1 llama2 llama3 etc.\n\t Example argument: llama\n"
        numAboveC = 0
        
        
        if (len(sys.argv) <0):
            if '-help' in sys.argv[1]:
                print(helpString)
                exit()
        self.c = float(config['DEFAULT']['allelefreq cutoff'])
        self.locusLength = int(config['DEFAULT']['min locus distance'])
        self.refPanelNumber = int(config['DEFAULT']['number of reference panels'])
        self.recombinationRate = float(config['DEFAULT']['recombination_rate'])
        self.minChrom = int(config['DEFAULT']['minChrom'])
        self.fileName = str(config['DEFAULT']['filename'])
        self.populationNames = config['DEFAULT']['populationNames'].split(',')
        print(self.populationNames)

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
            reference0_count = 0;reference1_count = 0;reference0AlleleA = 0;
            reference1AlleleA = 0;sampleAlleleA=0;sample_count=0
            chromosomeNumber = 0
            position = 0
            recombinationFrequency = 0.0

            lineList = list(line.split('\t'))
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
                
                # check if all columns have been read
                p = 0
                initPopulationCount=len(populationColumnsList[0])
                for populationColumn in populationColumnsList :
                    if populationColumn == [] :
                        print("Could not find population "+self.populationNames[p]+'\nExiting')
                        exit()
                    elif len(populationColumn) != initPopulationCount:
                        print("Different population sizes are not allowed "+self.populationNames[p]+'\nExiting')
                        # exit()
                    p += 1
                # ref lists [0] [1] ... sample list [2] [3] ...
                print(populationColumnsList)
            else: # actual data of VCF
               
                for refpopulation in range(0,self.refPanelNumber-1): # compute allele count for at least two reference populations
                    reference0_count=0
                    reference0AlleleA=0
                    print("refpopulation"+str(refpopulation)+'-'+str(refpopulation+1))
                    for column in range(0,len(populationColumnsList[refpopulation])):
                        refPop1 = populationColumnsList[refpopulation+1][column]
                        print(refPop1)
                        refPop2 = populationColumnsList[refpopulation][column]
                        print(refPop2)

                        if (lineList[refPop1][:1] != '.'): 
                            if int(lineList[refPop1][:1])<=1:
                                reference0AlleleA += int(lineList[refPop1][:1])
                                reference0_count += 1
                        
                        if (lineList[refPop1][2:3] != '.'): 
                            if int(lineList[refPop1][2:3])<=1:
                                reference0AlleleA += int(lineList[refPop1][2:3])
                                reference0_count += 1

                        if (lineList[refPop2][:1] != '.'):
                            if int(lineList[refPop2][:1])<=1:
                                reference1AlleleA += int(lineList[refPop2][:1])
                                reference1_count += 1
                        
                        if (lineList[refPop2][2:3] != '.'):
                            if int(lineList[refPop2][2:3])<=1:
                                reference1AlleleA += int(lineList[refPop2][2:3])
                                reference1_count += 1
                                
                    alleleCountRefAltList.append(abs(reference0_count-reference0AlleleA))
                    alleleCountRefAltList.append(reference0AlleleA)
                    alleleCountRefAltList.append(abs(reference1_count-reference1AlleleA))
                    alleleCountRefAltList.append(reference1AlleleA)

                if reference0_count == 0:
                    reference0_count = 1 # so theres no dividing by zero
                if reference1_count == 0:
                    reference1_count = 1 # so theres no dividing by zero
                chrom_number1 = reference0_count
                chrom_number2 = reference1_count
                if chrom_number1 < self.minChrom or chrom_number2 < self.minChrom: # make parameter the code for min chromosome number
                    print("low chrom")
                    continue # returns back to the top of the for loop
                elif abs(reference0AlleleA/reference0_count-reference1AlleleA/reference1_count) < self.c:
                    print('low allele freq')
                    continue
                else : # if passed all thresholds specified from config.ini
                    chromosomeNumber= lineList[0][lineList[0].find('^[A-Z]'):]
                    distanceToNextLocus = abs(currentLocus - int(lineList[1]))
                    if (distanceToNextLocus > self.locusLength):
                        print("low dist")
                        for samplePopulation in range(self.refPanelNumber,len(populationColumnsList)): # find read count for each sample 
                            for column in range(0,len(populationColumnsList[len(populationColumnsList)-1])):
                                samPopCol = populationColumnsList[samplePopulation][column]
                                # GT:AD:DP:GQ:PGT:PID:PL
                                readCounts = lineList[samPopCol].split(':')[1].split(',')
                                if len(readCounts) == 2:
                                    for read in readCounts :
                                        readCountRefAltList.append(read)
                
                        currentLocus = int(lineList[1])
                        recombinationFrequency=distanceToNextLocus*self.recombinationRate 
                        # compute recombination frequency -- simple formula = (distance*recombinationRate)
                        alleleCountRefAltListString ='\t'.join(str(e) for e in alleleCountRefAltList)
                        readCountRefAltListString='\t'.join(str(e) for e in readCountRefAltList)
                        print("{}\t{}\t{}\t{}\t{}\n".format(\
                                    chromosomeNumber,currentLocus,alleleCountRefAltListString,recombinationFrequency,readCountRefAltListString))
def main():
    testObject = FreqDistanceCalculator()
    fileOut= input("Name of output file:")
    testObject.createAncestryHMMInputFile(fileOut)    
    print('Created output file:'+fileOut+" for Ancestry_HMM input")

main()
