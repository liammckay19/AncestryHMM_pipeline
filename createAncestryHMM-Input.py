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
        self.config = configparser.ConfigParser()
        self.config.read('config.ini')
        print('Your config file input:\n')
        pretty(dict(self.config['DEFAULT']),0)
        self.c = float(self.config['DEFAULT']['allelefreq cutoff'])
        self.locusLength = int(self.config['DEFAULT']['min locus distance'])
        self.refPanelNumber = -1
        self.recombinationRate = float(self.config['DEFAULT']['recombination_rate'])
        self.minChrom = int(self.config['DEFAULT']['minChrom'])
        self.fileName = str(self.config['DEFAULT']['filename'])
        self.refPopulationNames = self.config['DEFAULT']['refPopulationNames'].split(',')
        self.samplePopulationNames = self.config['DEFAULT']['samplePopulationNames'].split(',')
        self.populationNames = self.refPopulationNames + self.samplePopulationNames
        self.refPopColI=[]
        self.samPopColI=[]  
        refColIndexString=self.config['DEFAULT']['refPopulationColumnIndices']
        samColIndexString=self.config['DEFAULT']['samplePopulationColumnIndices']
        try :
            self.parseColumnIndices(refColIndexString, self.refPopColI)
        except:
            print('Unable to parse column integer numbers from', refColIndexString)
        try :
            self.parseColumnIndices(samColIndexString, self.samPopColI)
        except:
            print('Unable to parse column integer numbers from', samColIndexString)

        if self.populationNames != [] and (self.samPopColI == [] or self.refPopColI == []):
            print('Using names of columns not column indices specified')
        elif self.populationNames == [] and (self.samPopColI == [] or self.refPopColI == []):
            print('Unable to load columns with information given. \nExiting')
            exit()

        print()
        self.lineCount = 0
   
    def parseColumnIndices (self, string, a):
        string=string.split(';')
        for tupleString in string:
            if ',' in tupleString:
                ranges = tupleString.split(',')
                rangeLow = int(ranges[0])
                rangeHigh = int(ranges[1])
                a += [list(range(rangeLow,rangeHigh+1))]
            else:
                a += [[int(tupleString)]]

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
                print()
                if self.refPopColI == [] and self.samPopColI == []:
                    for populationName in self.populationNames:
                        populationColumnsList.append([ i for i, column in enumerate(lineList) if re.search('^'+populationName, column) ])

                    # If no reference panel number of columns specified, 
                    # set to length of reference panels found with names specified
                    self.refPanelNumber = len(self.refPopulationNames)

                    if self.refPanelNumber > len(self.refPopulationNames) : print('Using more than named reference panels for filtering alleles:\t')
                    else : print('Using these reference panels for filtering alleles:\t')
                    for a in range(self.refPanelNumber):
                        print('\t'+self.populationNames[a]+' | Columns Found: '+str(populationColumnsList[a]))

                    print('Using these samples panels:\t')
                    # check if all columns have been read
                    p = 0
                    initPopulationCount=len(populationColumnsList[0])
                    for populationColumn in populationColumnsList :
                        if populationColumn == [] :
                            print("Could not find population "+self.populationNames[p]+'\nContinuing')
                        if p>self.refPanelNumber:
                            print('\t'+self.populationNames[p]+' | Columns Found: '+str(populationColumnsList[p]))
                        p += 1
                else :
                    self.refPanelNumber = len(self.refPopColI)
                    populationColumnsList = self.refPopColI + self.samPopColI
                    print('Using these reference panels for filtering alleles:\t')
                    for x in range(len(self.refPopColI)):
                        for y in range(len(self.refPopColI[x])) :
                            print(lineList[self.refPopColI[x][y]],end ="\t")

                    print('Using these sample panels:\t')
                    for x in range(len(self.samPopColI)):
                        for y in range(len(self.samPopColI[x])) :
                            print(lineList[self.samPopColI[x][y]],end ="\t")
                    print('\n')

            else: # actual data of VCF
               
                for refpopulation in range(0,self.refPanelNumber-1): # compute allele count for at least two reference populations
                    reference1_count=0
                    reference1AlleleA=0
                    bothRefList = list(zip(populationColumnsList[refpopulation+1],populationColumnsList[refpopulation]))

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
                +"\nThere must be at least a reference and panel specified each in one of two ways\n" \
                +"\n[refPopulationNames]\n\t (at least 2 Strings no quotes) Names of reference panels in the VCF file.\n\t Reference panels should be named like guanaco0 guanaco1 guanaco2 etc.\n\t Example argument: guanaco,vicugna\n" \
                +"\n[samplePopulationNames]\n\t (at least 2 Strings no quotes) Names of sample panels to be run in Ancestry_HMM. \n\t Should be named like llama1 llama2 llama3 etc.\n\t Example argument: llama,alpaca\n"\
                +"\n[refPopulationColumnIndices]\n\t (at least 2 Strings no quotes) Column index of reference panels to be run in Ancestry_HMM. \n\t [Syntax] beginning,end;beginning,end;...\n\t Example argument: 46,49;50,58;59,64\n"\
                +"\n[samplePopulationColumnIndices]\n\t (at least 2 Strings no quotes) Column index of sample panels to be run in Ancestry_HMM. \n\t [Syntax] beginning,end;beginning,end;...\n\t Example argument: 69;70,82;83,98\n"
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
