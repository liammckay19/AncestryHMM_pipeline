import sys
import re
import io


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
        usage="\nUsage: python alleleFreqCorrected.py [allelefreq cutoff] [min locus distance]\n" \
                    +"      [number of reference panels] [recombination_rate] [minChrom] [filename] \n" \
                    +"      [refPanel names(requires 2 names)] [sample names][...][...] \n"\
                    +"Type -help to get descriptions of arguments\n "
        helpString = "[allelefreq cutoff]\n\t (Float) cutoff value for reference panel allele frequency calculation\n"\
                    +"\n[min locus distance]\n\t (Integer) minimum distance between each allele locus\n" \
                    +"\n[number of reference panels]\n\t (Integer) number of reference panel species in the VCF\n" \
                    +"\n[recombination_rate]\n\t (Float) estimated recombination rate for recombination probability for Ancestry_HMM input.\n\t Average recombinations per base pairs\n" \
                    +"\n[minChrom]\n\t (Integer) the minimum amount of chromosomes that must be present in the reference panel alleles to make it through the threshold\n" \
                    +"\n[filename]\n\t (String no quotes) Name of VCF file on local machine\n" \
                    +"\n[refPanel names(requires 2 names)]\n\t (at least 2 Strings no quotes) Names of reference panels in the VCF file.\n\t Reference panels should be named like guanaco0 guanaco1 guanaco2 etc.\n\t Example argument: guanaco\n" \
                    +"\n[sample names]\n\t (at least 2 Strings no quotes) Names of sample panels to be run in Ancestry_HMM. \n\t Should be named like llama1 llama2 llama3 etc.\n\t Example argument: llama\n"
        numAboveC = 0
        
        try:
            if '-help' in sys.argv[1]:
                print(helpString)
                exit()
            self.c = float(sys.argv[1])
            self.locusLength = int(sys.argv[2])
            self.refPanelNumber = int(sys.argv[3])
            self.recombinationRate = float(sys.argv[4])
            self.minChrom = int(sys.argv[5])
            self.fileName = str(sys.argv[6])
            self.populationNames = list(sys.argv[7:])
            try:
                self.read = VCFreader(self.fileName)
            except FileNotFoundError as err:
                print('{} \ncannot open'.format(err), self.fileName)
                exit()
            else:
                print(self.fileName, 'to be read')

        except IndexError:
            if(len(sys.argv) != 8):
                print("Incorrect number of arguments."+usage)
            exit()
        except ValueError as err:
            print("Could not parse {}\n".format(err)+usage)
            exit()
        except UnboundLocalError as err:
            print("{}\nCould not read file: {}".format(err, self.fileName))
            exit()
        # except:
        #     print(usage)
        #     exit()

        

    def createPrunedFreqDistanceDataFile(self, outFileStr):
        if(len(self.populationNames)!=2):
            print(self.populationNames)
            print('cannot create pruned data only from other than two reference panel populations')
            return

        numAboveC = 0
        outFile = open(outFileStr, 'w')
        currentLocus = 0

        for line in self.read.readVCF():
            reference0_count = 0
            numer1 = 0
            reference1_count = 0
            numer2 = 0
            lineList = list(line.split('\t'))
            if lineList[0] == '':
                break
            if lineList[0].startswith('##'):
                continue
            elif lineList[0].startswith('#'):
                populationColumns1 = [ i for i, column in enumerate(lineList) if re.search('^'+self.populationNames[0], column) ]
                populationColumns2 = [ i for i, column in enumerate(lineList) if re.search('^'+self.populationNames[1], column) ]
            else:
                for column in range(0,max(len(populationColumns1),len(populationColumns2))):
                    colPop_1 = populationColumns1[column]
                    colPop_2 = populationColumns2[column]
                    if (lineList[colPop_1][2:3] != '.'): 
                        if int(lineList[colPop_1][2:3])<=1:
                            reference0_count += 1

                    if (lineList[colPop_1][:1] != '.'): 
                        if int(lineList[colPop_1][:1])<=1:
                            numer1 += int(lineList[colPop_1][:1])
                            reference0_count += 1
                    if (lineList[colPop_2][2:3] != '.'):
                        if int(lineList[colPop_2][2:3])<=1:
                            reference1_count += 1

                    if (lineList[colPop_2][:1] != '.'):
                        if int(lineList[colPop_2][:1])<=1:
                            numer2 += int(lineList[colPop_2][:1])
                            reference1_count += 1

                if reference0_count == 0:
                    reference0_count = 1 # so theres no dividing by zero
                if reference1_count == 0:
                    reference1_count = 1 # so theres no dividing by zero
                chrom_number1 = reference0_count
                chrom_number2 = reference1_count
                if chrom_number1 <= 14 or chrom_number2 <= 14: # make parameter the code for min chromosome number
                    continue
                elif abs(numer1/reference0_count-numer2/reference1_count) < self.c:
                    continue
                else :
                    if (abs(currentLocus - int(lineList[1])) > self.locusLength):
                        currentLocus = int(lineList[1])
                        outFile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tallelefreq={}\n".format(\
                                    self.populationNames[0],self.populationNames[1],lineList[1],lineList[3],lineList[4],\
                                    numer1,reference0_count,round(numer1/reference0_count,3),numer2,reference1_count,round(numer2/reference1_count,3),\
                                    round(abs(numer1/reference0_count-numer2/reference1_count),3)))
                    numAboveC += 1
        outFile.write("#"+str(numAboveC))

    def createAncestryHMMInputFile(self, outFileName):
        numAboveC = 0
        outFile = open(outFileName, 'w')
        currentLocus = 0
        for line in self.read.readVCF():
            reference0_count = 0
            numer1 = 0
            reference1_count = 0
            numer2 = 0
            reference0AlleleA = 0
            reference0Allelea = 0
            reference1AlleleA = 0
            reference1Allelea = 0

            lineList = list(line.split('\t'))
            alleleCountRefAltList = []
            readCountRefAltList = []
            chromosomeNumber = 0
            position = 0
            recombinationFrequency = 0.0
            if lineList[0] == '':
                break
            if lineList[0].startswith('##'):
                continue
            elif lineList[0].startswith('#'):
                populationColumnsList = []
                for populationName in self.populationNames:
                    populationColumnsList.append([ i for i, column in enumerate(lineList) if re.search('^'+populationName, column) ])
                # print(populationColumnsList) 
                # ref lists [0] [1] sample list [2] [3] 
            else:
                for refpopulation in range(1,self.refPanelNumber):
                    for column in range(0,len(populationColumnsList[refpopulation])):
                        refPop1 = populationColumnsList[refpopulation-1][column]
                        refPop2 = populationColumnsList[refpopulation][column]
                        if (lineList[refPop1][:1] != '.'): 
                            if int(lineList[refPop1][:1])<=1:
                                numer1 += int(lineList[refPop1][:1])
                                reference0AlleleA += int(lineList[refPop1][:1])
                                reference0_count += 1
                        if (lineList[refPop1][2:3] != '.'): 
                            if int(lineList[refPop1][2:3])<=1:
                                reference0AlleleA += int(lineList[refPop1][2:3])
                                reference0_count += 1

                        if (lineList[refPop2][:1] != '.'):
                            if int(lineList[refPop2][:1])<=1:
                                reference1AlleleA += int(lineList[refPop2][:1])
                                numer2 += int(lineList[refPop2][:1])
                                reference1_count += 1
                        if (lineList[refPop2][2:3] != '.'):
                            if int(lineList[refPop2][2:3])<=1:
                                reference1AlleleA += int(lineList[refPop2][2:3])
                                reference1_count += 1
                                
                    alleleCountRefAltList.append(reference0_count-reference0AlleleA)
                    alleleCountRefAltList.append(reference0AlleleA)
                    alleleCountRefAltList.append(reference1_count-reference1AlleleA)
                    alleleCountRefAltList.append(reference1AlleleA)

                    if reference0_count == 0:
                        reference0_count = 1 # so theres no dividing by zero
                    if reference1_count == 0:
                        reference1_count = 1 # so theres no dividing by zero
                    chrom_number1 = reference0_count
                    chrom_number2 = reference1_count
                    if chrom_number1 < self.minChrom or chrom_number2 < self.minChrom: # make parameter the code for min chromosome number
                        continue
                    elif abs(reference0AlleleA/reference0_count-reference1AlleleA/reference1_count) < self.c:
                        continue
                    else :
                        # alleleCountRefAltList.append()
                        chromosomeNumber= lineList[0][lineList[0].find('^[A-Z]'):]
                        distanceToNextLocus = abs(currentLocus - int(lineList[1]))
                        if (distanceToNextLocus > self.locusLength):
                            currentLocus = int(lineList[1])
                            recombinationFrequency=distanceToNextLocus*self.recombinationRate
                            alleleCountRefAltListString ='\t'.join(str(e) for e in alleleCountRefAltList)
                            readCountRefAltListString='\t'.join(str(e) for e in readCountRefAltList)
                            print("{}\t{}\t{}\t{}\t{}\n".format(\
                                        chromosomeNumber,currentLocus,alleleCountRefAltListString,recombinationFrequency,readCountRefAltListString))

# 1. Chromosome
# 2. Position in basepairs
# 3. Allele counts of allele A in reference panel 0
# 4. Allele counts of allele a in reference panel 0
# 5. Allele counts of allele A in reference panel 1
# 6. Allele counts of allele a in reference panel 1
# If there are additional reference panels, they are included following these columns. So, the next panel would be columns 7 and 8, the next 9 and 10, and so on. All of the following column numbers would be augmented by two for each additional reference panel included. The number of reference panels provided must match the number specified on the command line using –a.
# 7. Distance in Morgans between the previous marker and this position. For the first position on a given chromosome, this may take any value, as it will be ignored.
# Following this, an option column is the site-specific error rates for each allele. Here, the first column denotes the error rate where a read (or genotype) that is really an A is reported as a, and vice versa for the following column. If this is provided (via specifying –E on the command line), the following columns numbers should also be augmented by two as well.
# Each sample is then represented by two columns with counts corresponding to
# 8. Read counts of allele A in sample 1
# 9. Read counts of allele a in sample 1
# 10. Read counts of allele A in sample 2 
# 11. Read counts of allele a in sample 2


      

def main():
    testObject = FreqDistanceCalculator()
    # testObject.createPrunedFreqDistanceDataFile('prunedAlleleFreq-Distance.tsv')
    fileOut= input("Name of output file:")
    testObject.createAncestryHMMInputFile(fileOut)
    print('Created '+fileOut+" for Ancestry_HMM input")
    # create a commandline flag for making a pruned data set for two populations or for computing an input file for Russ's program

main()
