import os,math,sys, numpy as np
import scipy
from scipy import linalg
from numpy import dot
from collections import defaultdict
from sklearn.cluster import KMeans
import glob
import matplotlib.pyplot as plt
from scipy.spatial import distance

###Once all the random clustering frequency matrices are created (see IterativeclusteringonSEM_Randomclustering.py) this method now loops over all those frequency matrix and identifies a single value which
#best separates compounds in the same cluster vs those in different clusters. This is accomplished by calculating the distance between every compound and separating these values into 2 values.
#The output is on graph for each data type. Threshold is determined from the graph manually.


def clustercompare():
	compoundlist_Mu=['Buprenorphine','DAMGO','Endomorphine1','Fentanyl','Loperamide','Meptazinol','Met-Enk','Morphine','Oxycodone','Tramadol','Cmp 1','Cmp 3','Cmp 4','Cmp 5','Cmp 6','Cmp 7','Cmp 8','Cmp 9','Cmp 10','Cmp 11','Cmp 12','Cmp 13','Cmp 14','Cmp 15','Cmp 16']
	compoundlist_Delta=['Buprenorphine','Fentanyl','Loperamide','Met-Enk','Morphine','Oxycodone','Cmp 1','Cmp 2','Cmp 3','Cmp 4','Cmp 5','Cmp 6','Cmp 7','Cmp 8','Cmp 9','Cmp 10','Cmp 11','Cmp 12','Cmp 13','Cmp 14','Cmp 15','Cmp 16']
	
	filetype='{0}'.format(sys.argv[1])#Specific the datatype to work on
	
	randomallarray=[]
	alldict=defaultdict(list)
	######Repeat above process for all randomly generated files
	countsamearray=[]
	if filetype='Structure':
		path='../compoundStructure/Clusters/Random/'
	else:
		path='../SEM/RandomClustering/%s/FrequencyMatrix/'%(filetype)
	filecounter=0
	for filename in glob.glob(os.path.join(path,'*Average*')):#Loop over all random frequency matrices
		if filecounter<200:
			file=open(filename)##Parsing over file to create a matrix in python
			#print(filename)
			rowcounter=0
			linecounter=0
			if filetype[1]=='D':#21 compounds in Delta Data
				randommatrix=np.zeros((21,21))
				posarray=range(0,21)

			else:
				randommatrix=np.zeros((25,25))#25 compounds in all other data types
			for line in file:
				if linecounter>0:
					colcounter=0
					parts=line.strip('\n').split('\t')
					#print('File line',len(parts),parts)

					for i in parts[1:-1]:
						randommatrix[rowcounter][colcounter]=i
						colcounter=colcounter+1
					rowcounter=rowcounter+1
				linecounter=linecounter+1
			randomcluster=KMeans(n_clusters=4).fit_predict(randommatrix)#Clustering the frequency matrix into 4 clusters. 4 Clusters were used in the actual data clustering as well, so the number of clusters should match between random and actual data. If # clusters changed in one, must be changed in both.
			threshdict=defaultdict(list)
			for i in posarray: #For every compound in matrix
				for j in posarray:#for Every compound in matrix
					if i!=j:
						index=str(randomcluster[i])+'-'+str(randomcluster[j])
						tempdistance=distance.euclidean(np.array(randommatrix[i]),np.array(randommatrix[j]))#Calculate the distance between compound 1 and compound 2 using the array (row) of similarity values from the frequency matrix
						#print('hmorarray',compoundlist_Mu[i],randommatrix[i])
						#print('hmorarray',compoundlist_Mu[j],randommatrix[i])
						#print('tempdistance',tempdistance)
						threshdict[index].append(tempdistance)#Save distance in a dictionary, indexed by comp1-comp2
			#print('threshdict',threshdict)
			for index in threshdict:
				for i in threshdict[index]:
					allarray.append(i)
				parts=index.split('-')
				if parts[1]==parts[0]:##Sort the distance dictionary (from above) into two groups based on whether the 2 compounds are in the same or different clusters. 
					for i in threshdict[index]:
						alldict['same'].append(i)
				else:
					for i in threshdict[index]:
						alldict['diff'].append(i)
				temparray=threshdict[index]
				meanvalue=np.mean(temparray)
				stdvalue=np.std(temparray)
				#print(index,meanvalue,stdvalue)
		filecounter+=1
	#plt.hist(allarray)
	#plt.show()
	print('alldict',alldict)
	plotdict=defaultdict(list)						
	for i in alldict:
		print(i,np.mean(alldict[i]),np.std(alldict[i]))
	#Plot the two distributions and Save
	bins=np.linspace(0,4.0,40)
	plt.hist(alldict['same'],bins,alpha=0.5,color='green', label='same')
	plt.hist(alldict['diff'],bins,alpha=0.5,color='red',label='diff')
	plt.title('Differences in distance betweeen compounds in the same and in different clusters')
	plt.legend()
	plt.xlabel('Euclidian Distance')
	plt.ylabel('Frequency')
	if filetype=='Structure':
		plt.savefig('../compoundStructure/Clusters/DistanceDifferences_Random.png')
	else:
		plt.savefig('../SEM/RandomClustering/DistanceDifferences_%s.png'%(filetype))
	#plt.show()
	
	
	
	
	
	
	

clustercompare()