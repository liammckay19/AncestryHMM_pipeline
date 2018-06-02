import sys
import re
import io

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
# 10. Read counts of allele A in sample 2 11. Read counts of allele a in sample 2


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

class inputCreation :
    def __init__(self, outputFile):
        self.numAboveC = 0
        self.outputFile = outputFile
        try:
            self.c = float(sys.argv[1])
            self.locusLength = int(sys.argv[2])

            self.populationName1 = str(sys.argv[3])
            self.populationName2 = str(sys.argv[4])
            self.fileName = str(sys.argv[5])
            try:
                self.read = VCFreader(fileName)
            except FileNotFoundError:
                print('cannot open', fileName)
                return
            else:
                print(fileName, 'to be read')

        except IndexError:
            if(len(sys.argv) != 6):
                print("Incorrect number of arguments. Usage:\n \
        python alleleFreqCorrected.py [float c-value] [min locus distance] [population name 1] [population name 2] [filename]\n")
            else:
                print("Usage: python alleleFreqCorrected.py [float c-value] [population name 1] [population name 2] [filename]")
            return
        except UnboundLocalError:
            print("Could not read file: {}".format(fileName))
            return
        

    def calcAlleleFreq(self):
        outFile = open(self.outputFile, 'w')
        for line in self.read.readVCF():
            denom1 = 0
            numer1 = 0
            denom2 = 0
            numer2 = 0
            lineList = list(line.split('\t'))
            if lineList[0] == '':
                break
            if lineList[0].startswith('##'):
                continue
            elif lineList[0].startswith('#'):
                populationColumns1 = [ i for i, column in enumerate(lineList) if re.search('^'+populationName1, column) ]
                populationColumns2 = [ i for i, column in enumerate(lineList) if re.search('^'+populationName2, column) ]
            else:
                for column in range(0,max(len(populationColumns1),len(populationColumns2))):
                    colPop_1 = populationColumns1[column]
                    colPop_2 = populationColumns2[column]
                    if (lineList[colPop_1][2:3] != '.'): 
                        if int(lineList[colPop_1][2:3])<=1:
                            denom1 += 1

                    if (lineList[colPop_1][:1] != '.'): 
                        if int(lineList[colPop_1][:1])<=1:
                            numer1 += int(lineList[colPop_1][:1])
                            denom1 += 1
                    if (lineList[colPop_2][2:3] != '.'):
                        if int(lineList[colPop_2][2:3])<=1:
                            denom2 += 1

                    if (lineList[colPop_2][:1] != '.'):
                        if int(lineList[colPop_2][:1])<=1:
                            numer2 += int(lineList[colPop_2][:1])
                            denom2 += 1

                if denom1 == 0:
                    denom1 = 1 # so theres no dividing by zero
                if denom2 == 0:
                    denom2 = 1 # so theres no dividing by zero
                chrom_number1 = denom1
                chrom_number2 = denom2
                if chrom_number1 <= 14 or chrom_number2 <= 14: # make parameter the code for min chromosome number
                    continue
                elif abs(numer1/denom1-numer2/denom2) < float(c):
                    continue
                else :
                    if (abs(currentchrom - int(lineList[1])) > locusLength):
                        currentchrom = int(lineList[1])
                        outFile.write("{}\t{}\t{}\t{}\t{}\t{}/{}\t{}\t{}/{}\t{}\tallelefreq={}".format(\
                                    populationName1,populationName2,lineList[1],lineList[3],lineList[4],\
                                    numer1,denom1,round(numer1/denom1,3),numer2,denom2,round(numer2/denom2,3),\
                                    round(abs(numer1/denom1-numer2/denom2),3)))
                    numAboveC += 1
        outFile.write("#",numAboveC)

def main():
    testObject = inputCreation('prunedAlleleFreq-Distance.tsv')
    testObject.calcAlleleFreq()

main()
