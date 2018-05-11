
import matplotlib.pyplot as plt
import csv
import numpy as np
import pandas as pd
import re

# run python within the file containing this tsv file
fname = 'NewChr1.snp.vcf'
print(fname)
newChr1 = pd.read_csv(fname, skiprows=range(0,43),sep='\t')
# newChr1['INFO'].to_csv('newChr1-INFO.tsv',sep='\t')


newChr1['Bactrian_Camelfreq'] = newChr1['Bactrian_Camel'].str.replace('[:].*','',0)
newChr1['Bactrian_Camelfreq'] = newChr1['Bactrian_Camelfreq'].str.replace('[.]','0',0)
newChr1['Bactrian_Camelfreq'] = newChr1['Bactrian_Camelfreq'].str.replace('[2-9]','0',0)
newChr1['Bactrian_Camelfreq1st'] = newChr1['Bactrian_Camelfreq'].str[0:1]
newChr1['Bactrian_Camelfreq2nd'] = newChr1['Bactrian_Camelfreq'].str[2:3]
firstAlleleFreqCamel = sum(list(map(int,newChr1['Bactrian_Camelfreq1st'].tolist())))
secondAlleleFreqCamel = sum(list(map(int,newChr1['Bactrian_Camelfreq2nd'].tolist())))
print(firstAlleleFreqCamel,secondAlleleFreqCamel)

# vicugna18
newChr1['vicugna18freq'] = newChr1['vicugna18'].str.replace('[:].*','',0)
newChr1['vicugna18freq'] = newChr1['vicugna18freq'].str.replace('[.]','0',0)
newChr1['vicugna18freq'] = newChr1['vicugna18freq'].str.replace('[2-9]','0',0)
newChr1['vicugna18freq1st'] = newChr1['vicugna18freq'].str[0:1]
newChr1['vicugna18freq2nd'] = newChr1['vicugna18freq'].str[2:3]
firstAlleleVicugnaEighteen = sum(list(map(int,newChr1['vicugna18freq1st'].tolist())))
secondAlleleVicugnaEighteen = sum(list(map(int,newChr1['vicugna18freq2nd'].tolist())))
print(firstAlleleVicugnaEighteen,secondAlleleVicugnaEighteen)


# vicugna20
newChr1['vicugna20freq'] = newChr1['vicugna20'].str.replace('[:].*','',0)
newChr1['vicugna20freq'] = newChr1['vicugna20freq'].str.replace('[.]','0',0)
newChr1['vicugna20freq'] = newChr1['vicugna20freq'].str.replace('[2-9]','0',0)
newChr1['vicugna20freq1st'] = newChr1['vicugna20freq'].str[0:1]
newChr1['vicugna20freq2nd'] = newChr1['vicugna20freq'].str[2:3]
firstAlleleVicugnaTwenty = sum(list(map(int,newChr1['vicugna20freq1st'].tolist())))
secondAlleleVicugnaTwenty = sum(list(map(int,newChr1['vicugna20freq2nd'].tolist())))
print(firstAlleleVicugnaTwenty,secondAlleleVicugnaTwenty)


# vicugna33
newChr1['vicugna33freq'] = newChr1['vicugna33'].str.replace('[:].*','',0)
newChr1['vicugna33freq'] = newChr1['vicugna33freq'].str.replace('[.]','0',0)
newChr1['vicugna33freq'] = newChr1['vicugna33freq'].str.replace('[2-9]','0',0)
newChr1['vicugna33freq1st'] = newChr1['vicugna33freq'].str[0:1]
newChr1['vicugna33freq2nd'] = newChr1['vicugna33freq'].str[2:3]
firstAlleleVicugnaThirtyThree = sum(list(map(int,newChr1['vicugna33freq1st'].tolist())))
secondAlleleVicugnaThirtyThree = sum(list(map(int,newChr1['vicugna33freq2nd'].tolist())))
print(firstAlleleVicugnaThirtyThree,secondAlleleVicugnaThirtyThree)

# vicugna37
newChr1['vicugna37freq'] = newChr1['vicugna37'].str.replace('[:].*','',0)
newChr1['vicugna37freq'] = newChr1['vicugna37freq'].str.replace('[.]','0',0)
newChr1['vicugna37freq'] = newChr1['vicugna37freq'].str.replace('[2-9]','0',0)
newChr1['vicugna37freq1st'] = newChr1['vicugna37freq'].str[0:1]
newChr1['vicugna37freq2nd'] = newChr1['vicugna37freq'].str[2:3]
firstAlleleVicugnaThirtySeven= sum(list(map(int,newChr1['vicugna37freq1st'].tolist())))
secondAlleleVicugnaThirtySeven = sum(list(map(int,newChr1['vicugna37freq2nd'].tolist())))
print(firstAlleleVicugnaThirtySeven,secondAlleleVicugnaThirtySeven)

# vicugna38
newChr1['vicugna38freq'] = newChr1['vicugna38'].str.replace('[:].*','',0)
newChr1['vicugna38freq'] = newChr1['vicugna38freq'].str.replace('[.]','0',0)
newChr1['vicugna38freq'] = newChr1['vicugna38freq'].str.replace('[2-9]','0',0)
newChr1['vicugna38freq1st'] = newChr1['vicugna38freq'].str[0:1]
newChr1['vicugna38freq2nd'] = newChr1['vicugna38freq'].str[2:3]
firstAlleleVicugnaThirtyEight= sum(list(map(int,newChr1['vicugna38freq1st'].tolist())))
secondAlleleVicugnaThirtyEight = sum(list(map(int,newChr1['vicugna38freq2nd'].tolist())))
print(firstAlleleVicugnaThirtyEight,secondAlleleVicugnaThirtyEight)


# vicugna39
newChr1['vicugna39freq'] = newChr1['vicugna39'].str.replace('[:].*','',0)
newChr1['vicugna39freq'] = newChr1['vicugna39freq'].str.replace('[.]','0',0)
newChr1['vicugna39freq'] = newChr1['vicugna39freq'].str.replace('[2-9]','0',0)
newChr1['vicugna39freq1st'] = newChr1['vicugna39freq'].str[0:1]
newChr1['vicugna39freq2nd'] = newChr1['vicugna39freq'].str[2:3]
firstAlleleVicugnaThirtyNine= sum(list(map(int,newChr1['vicugna39freq1st'].tolist())))
secondAlleleVicugnaThirtyNine = sum(list(map(int,newChr1['vicugna39freq2nd'].tolist())))
print(firstAlleleVicugnaThirtyNine,secondAlleleVicugnaThirtyNine)


# vicugna40
newChr1['vicugna40freq'] = newChr1['vicugna40'].str.replace('[:].*','',0)
newChr1['vicugna40freq'] = newChr1['vicugna40freq'].str.replace('[.]','0',0)
newChr1['vicugna40freq'] = newChr1['vicugna40freq'].str.replace('[2-9]','0',0)
newChr1['vicugna40freq1st'] = newChr1['vicugna40freq'].str[0:1]
newChr1['vicugna40freq2nd'] = newChr1['vicugna40freq'].str[2:3]
firstAlleleVicugnaForty= sum(list(map(int,newChr1['vicugna40freq1st'].tolist())))
secondAlleleVicugnaForty = sum(list(map(int,newChr1['vicugna40freq2nd'].tolist())))
print(firstAlleleVicugnaForty,secondAlleleVicugnaForty)


BactrianCamelFreq = firstAlleleFreqCamel/secondAlleleFreqCamel
VicugnaEighteenFreq = firstAlleleVicugnaEighteen/secondAlleleVicugnaEighteen
VicugnaTwentyFreq = firstAlleleVicugnaTwenty/secondAlleleVicugnaTwenty
VicugnaThirtyThreeFreq = firstAlleleVicugnaThirtyThree/secondAlleleVicugnaThirtyThree
VicugnaThirtySevenFreq = firstAlleleVicugnaThirtySeven/secondAlleleVicugnaThirtySeven
VicugnaThirtyEightFreq = firstAlleleVicugnaThirtyEight/secondAlleleVicugnaThirtyEight
VicugnaThirtyNineFreq = firstAlleleVicugnaThirtyNine/secondAlleleVicugnaThirtyNine
VicugnaFortyFreq = firstAlleleVicugnaForty/secondAlleleVicugnaForty

frequencies = {'BactrianCamelFreq':BactrianCamelFreq,
'VicugnaEighteenFreq':VicugnaEighteenFreq,
'VicugnaTwentyFreq':VicugnaTwentyFreq,
'VicugnaThirtyThreeFreq':VicugnaThirtyThreeFreq,
'VicugnaThirtySevenFreq':VicugnaThirtySevenFreq,
'VicugnaThirtyEightFreq':VicugnaThirtyEightFreq,
'VicugnaThirtyNineFreq':VicugnaThirtyNineFreq,
'VicugnaFortyFreq':VicugnaFortyFreq}

pairwiseCompDict = dict()
for key,value in frequencies.items():
	for keyi,valuei in frequencies.items():
		if keyi != key:
			pairwiseCompDict[''+key+'-'+keyi]=abs(value-valuei)

# Set the c value for a threshold
c = 0.2
for key,value in pairwiseCompDict.items():
	if value > c:
		print(key,value)

