
import matplotlib.pyplot as plt
import csv
import numpy as np
import pandas as pd
import re

def computeAlleleFreq(dataframe,colname):
	dataframe[colname+'freq'] = dataframe[colname].str.replace('[:].*','',0)
	dataframe[colname+'freq'] = dataframe[colname+'freq'].str.replace('[.]','0',0)
	dataframe[colname+'freq'] = dataframe[colname+'freq'].str.replace('[2-9]','0',0)
	dataframe[colname+'freq1st'] = dataframe[colname+'freq'].str[0:1]
	dataframe[colname+'freq2nd'] = dataframe[colname+'freq'].str[2:3]
	firstAlleleFreq = sum(list(map(int,dataframe[colname+'freq1st'].tolist())))
	secondAlleleFreq = sum(list(map(int,dataframe[colname+'freq2nd'].tolist())))
	print(firstAlleleFreq,secondAlleleFreq)
	return firstAlleleFreq/secondAlleleFreq


def main():
	fname = input("File VCF input path:")
	print('File name: ',fname)
	dataframeVCF = pd.read_csv(fname, skiprows=range(0,43),sep='\t')
	computeColumns = input("List the columns you want to compare separated by commas: ")
	computeColumns=computeColumns.split(',')
	frequencies=dict()
	for column in computeColumns:
		frequencies[column] = computeAlleleFreq(dataframeVCF,column)
		# frequencies = {'Bactrian_Camel':BactrianCamelFreq,
		# 'Vicugna18':VicugnaEighteenFreq...

	pairwiseCompDict = dict()
	for key,value in frequencies.items():
		for keyi,valuei in frequencies.items():
			if keyi != key:
				pairwiseCompDict[''+key+'-'+keyi]=abs(value-valuei)

	# Set the c value for a threshold
	c = input("Set a threshold c value (float between 0 and 0.25):")
	c = float(c)
	for key,value in pairwiseCompDict.items():
		if value > c:
			print(key,value)

main()