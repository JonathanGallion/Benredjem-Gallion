import numpy as np
import scipy
from scipy import linalg
from numpy import dot

def creatematrix(filename,row,col):
	np.set_printoptions(threshold='nan')
	X=np.zeros((row,col))
	file=open(filename)
	rowcounter=0
	linecounter=0
	for line in file:
		if linecounter>1:
			colcounter=0
			parts=line.strip('\n').split('\t')
			'''print('File line',len(parts),parts)'''

			for i in parts[1:]:
				checkvalue=''
				try:
					float(i)
					checkvalue=True
				except ValueError:
					checkvalue=False
				print(i,checkvalue,rowcounter,colcounter)
				if checkvalue==True:
					X[rowcounter][colcounter]=i
				else:
					X[rowcounter][colcounter]='NaN'
				colcounter=colcounter+1
			rowcounter=rowcounter+1
		linecounter=linecounter+1
	#print('Xmatrix', X)
	return X
