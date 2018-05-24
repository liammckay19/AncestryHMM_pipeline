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
    populationName = str(sys.argv[2])
    for line in read.readVCF():
        denom = 0
        numer = 0
        lineList = list(line.split('\t'))
        if lineList[0].startswith('##'):
            continue
        elif lineList[0].startswith('#'):
            populationColumns = [ i for i, column in enumerate(lineList) if re.search('^'+populationName, column) ]
        else:
            for column in populationColumns:
                # print(lineList[column])
                # print(int(lineList[column][:1]),lineList[column][2:3])
                if (lineList[column][2:3] != '.'):
                    denom += 1
                if (lineList[column][:1] != '.'):
                    numer += int(lineList[column][:1])
            if denom == 0:
                denom = 1 # so theres no dividing by zero
            if denom > len(populationColumns)/2:
                if numer/denom > float(c):
                    numAboveC+=1
                    # print("{}: {} {} {} {} {}/{}={}".format(populationName,lineList[0],lineList[1],lineList[3],lineList[4],numer,denom,numer/denom))
    print(numAboveC)
main()

