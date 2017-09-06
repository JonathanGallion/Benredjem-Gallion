import os,math,sys, numpy as np
np.set_printoptions(threshold=np.nan)
import sklearn
from sklearn.cluster import KMeans
from nmf_creatematrix_structure import creatematrix
from sklearn.decomposition import NMF
import nimfa
from NMF_python_sparse import nmf


def choose_patients():
	##This method is a modification of the clustering used for the compound x assay data. It reads in a matrix which is composed of 3 compound x compound similarity matrices (ECFP6, FCFP6, and MDLkeys).
	# This code can be easily modified to cluster 1 or 2 similarity matrices, just modify the input matrix.
	numit=250
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
	freqmatrix=np.zeros((len(drugarray),len(drugarray)))
	
	
	print('freqmatrix',freqmatrix.shape)
	print(len(drugarray),len(assayarray))
	row=len(drugarray)
	col=len(assayarray)
	print(row,col)
	X=creatematrix(filename,row,col)#Making a matrix of the compound structure input data
	print(X)
	for basis in range(2,9):
		for clustk in range(2,9):
			if basis==clustk: #Clustering the input at different K's
				freqmatrix=np.zeros((len(drugarray),len(drugarray)))#Initialize temporary frequency matrix for this iteration

				for i in range(0,numit):
					print(basis,clustk,i)###Feature reduction, clustering, and frequency matrix creation performed in the same manner as the main clustering script. (IterativeClustering_SEM.py)
					W,H=nmf(X,basis)
					#print('Wmatrix',W)
					#print('Hmatrix',H)
					clusterdata=KMeans(n_clusters=clustk).fit_predict(W)
					for j in range (0,len(clusterdata)):
						for k in range (0,len(clusterdata)):
							'''print(j,k)'''
							if clusterdata[j]==clusterdata[k]:
								tempvalue=freqmatrix[j][k]
								freqmatrix[j][k]=tempvalue+1/float(numit)
								'''print('These two are the same cluster',j,k,freqmatrix[j][k])'''
					
				#print(freqmatrix)
				outfile=open('../compoundStructure/Clusters_3data/Frequencymatrix_Basis=%s_K=%s.txt'%(basis,clustk),'w')#Saving final frequency matrix
				outfile.write('\t')
				for i in drugarray:
					outfile.write(i)
					outfile.write('\t')
				outfile.write('\n')
				row,col=freqmatrix.shape
				for i in range(0,row):
					outfile.write(drugarray[i])
					outfile.write('\t')
					for j in range(0,col):
						outfile.write(str(freqmatrix[i][j]))
						outfile.write('\t')
					outfile.write('\n')
		
choose_patients()