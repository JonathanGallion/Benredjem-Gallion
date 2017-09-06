import os,math,sys, numpy as np
np.set_printoptions(threshold=np.nan)
from collections import defaultdict
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from scipy.stats import ks_2samp
from sklearn.cluster import KMeans



def Assayenrichment():
	######Input: Need Compound Frequency matrix and Phenotype frequency matrix already created for the dataset you want to use (e.g. hMor)
	######Input: Raw data in the specified location. This is the same data previously used for clustering
	#####Given these inputs this method uses KS-test to identify which parameters are separated into statistically significant groups.
	#First, the method creates groups of assays based on the frequency matrix and groups of compounds based on the frequency matrix. The K (number of clusters) can be adjusted in the method
	#Second, reads in raw assay data and sorts values based on the compound clusters
	#Finally, KS test using these groups for every assay in the input
	
	filetype='{0}'.format(sys.argv[1])# specify which data type to analyze.
	###Getting Assay clusters
	
	clusterassign=defaultdict(list)
	file=open('../SEM/%s/Assays/TotalFrequencyMatrix_Average.txt'%(filetype))#Read in Assay Frequency Matrix
	assaylist=[]
	assaymatrix=np.zeros((30,30))
	linecounter=0
	for line in file:
		parts=line.strip('\n').split('\t')
		if linecounter==0:
			for i in parts[1:-1]:
				assaylist.append(i)
		else:
			colcounter=0
			print(parts)
			for i in parts[1:-1]:
				assaymatrix[linecounter-1][colcounter]=i
				colcounter+=1
		linecounter+=1
	print(assaymatrix)
	print(assaylist)
	clusterdata=KMeans(n_clusters=3).fit_predict(assaymatrix)###Cluster Frequency matrix using Kmeans
	assaycounter=0
	clusterdict=defaultdict(list)
	print(len(assaylist))
	print(len(clusterdata))
	for assay in assaylist:
		print(assaycounter)
		clusterdict[str(clusterdata[assaycounter])].append(assay)
		assaycounter+=1
	print(clusterdata)
	print('clusterdict',clusterdict)##These are the assay clusters


	###getting compound clusters
	clusterdata=[]
	clusterassign=defaultdict(list)
	file=open('../SEM/%s/TotalFrequencyMatrix_Average.txt'%(filetype))#Read in Compound Frequency Matrix. Same data type (e.g. hMor) as the Assay matrix
	compoundlist=[]
	if filetype[1]=='D':#Initialize a matrix of 0's
		compoundmatrix=np.zeros((21,21))#Checking if the datatype is a type of delta receptor. Only 21 compounds were tested for Delta
	else:
		compoundmatrix=np.zeros((25,25))#25 Compounds tested for Mu
	linecounter=0
	for line in file:
		parts=line.strip('\n').split('\t')
		if linecounter==0:
			for i in parts[1:-1]:
				compoundlist.append(i)
		else:
			colcounter=0
			print(parts)
			for i in parts[1:-1]:
				compoundmatrix[linecounter-1][colcounter]=i
				colcounter+=1
		linecounter+=1
	clusterdata=KMeans(n_clusters=4).fit_predict(compoundmatrix)#Cluster the compound frequency matrix
	compcounter=0
	compoundclusterdict=defaultdict(list)
	print(len(compoundlist))
	print(len(clusterdata))
	print('clusterdata',clusterdata)
	for comp in compoundlist:
		print(compcounter,str(clusterdata[compcounter]))
		compoundclusterdict[str(clusterdata[compcounter])].append(comp)
		compcounter+=1
	print('compoundclusterdict',compoundclusterdict)##These are the compound clusters
	###Reading in the raw data
	datadict=defaultdict(list)
	file=open('../SEM/Data/%s/RawData_assays.txt'%(filetype))#Opening raw data and saving as a dictionary
	linecounter=0
	for line in file:
		parts=line.strip('\n').split('\t')
		if linecounter>1:
			comp=parts[0]
			for i in parts[1:]:
				#print(i)
				try:
					datadict[comp].append(float(i))
				except:
					datadict[comp].append(float(0))
		linecounter+=1
	print(datadict)
	
	finalassaylist=assaylist
	
	clustervaluedict=defaultdict(list)
	clustervaluedictall=defaultdict(list)
	try:
		os.mkdir('../SEM/%s/Assays/Assaycorrelations/'%(filetype))
	except:
		''
	for i in range(0,30):#There are 30 assay measures to loop over
		valuedict=defaultdict(list)
		endtuple=[]
		totallist=[]
		###Loop over compound clusters and then all compounds in each cluster. Sort raw data values into clusters based on the compound clusters.
		for clust in compoundclusterdict:	
			temptuple=[]
			for j in compoundclusterdict[clust]:
				#print(j,i)
				temptuple.append(datadict[j][i])
				totallist.append(datadict[j][i])
			print('temptuple',i,clust,temptuple)
			try:
				endtuple.append(temptuple)
			except:
				endtuple=[temptuple,]
		endtuple.append(totallist)#Adding the array of all raw values for that assay measure to the end of the tuple
		print('endtuple',endtuple)
		
		outfile=open('../SEM/%s/Assays/Assaycorrelations/Data_%s.txt'%(filetype,finalassaylist[i][:9]),'w')
		for part in endtuple:
			for smallerpart in part:
				outfile.write('%s\t'%(smallerpart))
			outfile.write('\n')
		parray=[] #Now compare each cluster of raw values to the array of all values using KS test
		t,p=ks_2samp(endtuple[0],endtuple[4])
		parray.append(round(p,4))
		t,p=ks_2samp(endtuple[1],endtuple[4])
		parray.append(round(p,4))
		t,p=ks_2samp(endtuple[2],endtuple[4])
		parray.append(round(p,4))
		t,p=ks_2samp(endtuple[3],endtuple[4])
		parray.append(round(p,4))
		
		paverage=np.mean(parray)
		
		clustervaluedict[finalassaylist[i]]=parray
		print('parray',parray,paverage)
		
		print('endtuple',endtuple)
		print(finalassaylist[i],np.mean(endtuple[0]),np.mean(endtuple[1]),np.mean(endtuple[2]),np.mean(endtuple[3]))
		plt.boxplot(endtuple)
		plt.title('Data distribution of %s -- %s\n%s'%(filetype,finalassaylist[i],parray))
		try:
			plt.savefig('../SEM/%s/Assays/Assaycorrelations/%s.png'%(filetype,finalassaylist[i][:9]))
		except:
			''
		plt.close()
	bins=np.linspace(min(totallist),max(totallist),11)
	plt.hist(endtuple[0],bins,alpha=.95,histtype='step',linewidth=4)
	plt.hist(endtuple[1],bins,alpha=.95,histtype='step',linewidth=4)
	plt.hist(endtuple[2],bins,alpha=.95,histtype='step',linewidth=4)
	plt.hist(endtuple[3],bins,alpha=.95,histtype='step',linewidth=4)
	#plt.show()
	plt.close()
	
	print('clustervaluedict',clustervaluedict)
	print('clusterdict',clusterdict)
	finaltuple=[]
	for clust in clusterdict:
		print('Start',clust,clusterdict[clust])
		temparray=[]
		for assay in clusterdict[clust]:
			for k in clustervaluedict[assay]:
				temparray.append(k)
				totallist.append(datadict[j][i])

		finaltuple.append(temparray)
	print('finaltuple',finaltuple)
	outfile=open('../SEM/%s/Assays/Assaycorrelations/PvalueArrays1.txt'%(filetype),'w')
	clustcounter=1
	for i in finaltuple:
	
		outfile.write('Cluster %s'%(str(clustcounter)))
		for j in i:
			outfile.write('\t')
			outfile.write(str(j))
		outfile.write('\n')
		clustcounter+=1
		
	
Assayenrichment()