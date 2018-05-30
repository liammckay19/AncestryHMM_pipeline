import sys
import re
import io

class VCFreader :
    ''' 
    Define objects to read FastA files.
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        #### print (head,seq)
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
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        with self.doOpen() as fileH:
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('##') :
                line = fileH.readline()
            for line in fileH:
                yield line
        yield line

def main():
    numAboveC = 0
    read = VCFreader("NewChr1.snp.vcf")
    c = sys.argv[1]
    populationName1 = str(sys.argv[2])
    populationName2 = str(sys.argv[3])
    for line in read.readVCF():
        denom1 = 0
        numer1 = 0
        denom2 = 0
        numer2 = 0
        lineList = list(line.split('\t'))
        if lineList[0].startswith('##'):
            continue
        elif lineList[0].startswith('#'):
            populationColumns1 = [ i for i, column in enumerate(lineList) if re.search('^'+populationName1, column) ]
            populationColumns2 = [ i for i, column in enumerate(lineList) if re.search('^'+populationName2, column) ]
        else:
            for column in range(0,max(len(populationColumns1),len(populationColumns2))):
                colPop_1 = populationColumns1[column]
                colPop_2 = populationColumns2[column]
                # print(populationColumns1,populationColumns2)
                # print((lineList[column][:1]),lineList[column][2:3])
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
            # print(chrom_number2,chrom_number1)
            if chrom_number1 < 14 or chrom_number2 < 14:
                continue
                # print("moving on1")
            elif abs(numer1/denom1-numer2/denom2) < float(c):
                continue
            else :
                # print("found allele: {}".format(lineList))
                # print("{},{}: {} {} {} {}/{}={}|{}/{}={}, allelefreq={}".format(\
                #             populationName1,populationName2,lineList[1],lineList[3],lineList[4],\
                #             numer1,denom1,round(numer1/denom1,3),numer2,denom2,round(numer2/denom2,3),\
                #             round(abs(numer1/denom1-numer2/denom2),3)))
                # if numAboveC % 10000==0: print(numAboveC)
                numAboveC += 1
            # if numer1+numer2 >= 14:
            #     if denom1 > len(populationColumns1)/2 and denom2 > len(populationColumns2)/2:
            #         if abs(numer1/denom1-numer2/denom2) > float(c):
                        # numAboveC+=1
                        # print(lineList[colPop_1],lineList[colPop_2])

                        # print("{},{}: {} {} {} {}/{}={}|{}/{}={}, allelefreq={}".format(\
                        #     populationName1,populationName2,lineList[1],lineList[3],lineList[4],\
                        #     numer1,denom1,round(numer1/denom1,3),numer2,denom2,round(numer2/denom2,3),\
                        #     round(abs(numer1/denom1-numer2/denom2),3)))
    print(numAboveC)
main()

