import os,math,sys, numpy as np
np.set_printoptions(threshold=np.nan)
import sklearn
from sklearn.cluster import KMeans
from nmf_creatematrix_structure_random import creatematrix
from sklearn.decomposition import NMF
import nimfa
from NMF_python_sparse import nmf


def choose_patients():

	###Just like "Clustering_structure.py" This method is used to cluster structural data. However, this script creates random clustering using the same input data by first reading in the compound structure
	#input data, and then randomly mixing the values. Clustering is performed on these randomized matrices. The output is 1000 random frequency matrices.
	

	filename='../compoundStructure/Threecombined.txt'
	drugarray=[]
	countarray=[]
	assayarray=[]
	#assayarray is just used to count the number of columns in the dataset. In this case the columns aren't assays, but this doesn't matter. All we need is the length of assayarray
	linecounter=0
	for line in open(filename):
		if linecounter==1:
			print(line)

			parts=line.strip('\n').split('\t')
			print(parts)
			for i in parts[1:]:
				if i !='':
					try:
						assayarray.append(i)
					except:
						assayarray=[i,]
		linecounter=linecounter+1
	print('assayarray',assayarray)
	
	#Getting list of drugs and the order in which they occur in the Data# This array is used for the final output. Need compound list.
	for line in open(filename):
		parts=line.strip('\n').split('\t')
		if parts[0]!='':
			try:
				drugarray.append(parts[0])
			except:
				drugarray=[parts[0]]
	print('drugarray',drugarray)
	freqmatrix=np.zeros((len(drugarray),len(drugarray))) #Initialize final frequency matrix with 0's
	
	
	print('freqmatrix',freqmatrix.shape)
	print(len(drugarray),len(assayarray))
	row=len(drugarray)
	col=len(assayarray)
	print(row,col)
	numtrials=0
	while numtrials<1000: #Repeat clustering 1000 times
		row=len(drugarray)
		col=len(assayarray)
		print(row,col)
		X=creatematrix(filename,row,col)
		print(X)
		np.savetxt('../compoundStructure/Clusters/RandomMatrices/RandomInput_%s.txt'%(numtrials),X)
		freqmatrix=np.zeros((len(drugarray),len(drugarray)))

		for basis in range(2,8):
			for clustk in range(2,8):
				if basis==clustk:

					for i in range(0,numit):
						print(basis,clustk,i)
							
							
						W,H=nmf(X,basis)
						#print('Wmatrix',W)
						#print('Hmatrix',H)
						clusterdata=KMeans(n_clusters=clustk).fit_predict(W)
						for j in range (0,len(clusterdata)):
							for k in range (0,len(clusterdata)):
								'''print(j,k)'''
								if clusterdata[j]==clusterdata[k]:
									tempvalue=freqmatrix[j][k]
									freqmatrix[j][k]=tempvalue+1/(6*float(numit))
					
					#print(freqmatrix)
		outfile=open('../compoundStructure/Clusters/Random/Frequencymatrix_%s.txt'%(numtrials),'w')
		outfile.write('\t')
		for i in drugarray:
			outfile.write(i)
			outfile.write('\t')
		outfile.write('\n')
		row1,col1=freqmatrix.shape
		for i in range(0,row1):
			outfile.write(drugarray[i])
			outfile.write('\t')
			for j in range(0,col1):
				outfile.write(str(freqmatrix[i][j]))
				outfile.write('\t')
			outfile.write('\n')
		numtrials+=1
		
choose_patients()