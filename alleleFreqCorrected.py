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
            line = fileH.readline()
            print(line)
            for line in fileH:
                yield line
        yield line

def main():
    numAboveC = 0

    
    try:
        c = float(sys.argv[1])
        populationName1 = str(sys.argv[2])
        populationName2 = str(sys.argv[3])
        fileName = str(sys.argv[4])
        read = VCFreader(fileName)
        try:
            f = open(fileName, 'r')
        except OSError:
            print('cannot open', fileName)
            return
        else:
            print(fileName, 'has', len(f.readlines()), 'lines')
            f.close()

    except IndexError:
        if(len(sys.argv) != 5):
            print("Incorrect number of arguments. Usage:\n \
    python alleleFreqCorrected.py [float c-value] [population name 1] [population name 2] [filename]\n")
        else:
            print("Usage: python alleleFreqCorrected.py [float c-value] [population name 1] [population name 2] [filename]")
        return
    except UnboundLocalError:
        print("Could not read file: {}".format(fileName))
        return
        
    

    for line in read.readVCF():
        denom1 = 0
        numer1 = 0
        denom2 = 0
        numer2 = 0
        lineList = list(line.split('\t'))
        print(lineList)
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
            if chrom_number1 < 14 or chrom_number2 < 14: # make parameter the code for min chromosome number
                continue
            elif abs(numer1/denom1-numer2/denom2) < float(c):
                continue
            else :
                print("{}\t{}\t{}\t{}\t{}\t{}/{}\t{}\t{}/{}\t{}\tallelefreq={}".format(\
                            populationName1,populationName2,lineList[1],lineList[3],lineList[4],\
                            numer1,denom1,round(numer1/denom1,3),numer2,denom2,round(numer2/denom2,3),\
                            round(abs(numer1/denom1-numer2/denom2),3)))
                # if numAboveC % 10000==0: print(numAboveC)
                numAboveC += 1
    print("#",numAboveC)
main()

