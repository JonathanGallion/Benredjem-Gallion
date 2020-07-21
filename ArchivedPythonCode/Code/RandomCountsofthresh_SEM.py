import os,math,sys, numpy as np
import scipy
from scipy import linalg
from numpy import dot
from collections import defaultdict
from sklearn.cluster import KMeans
import glob
import matplotlib.pyplot as plt
from scipy.spatial import distance

def clustercompare():
	#This code is designed to calculate the similarity between the clustering result from 2 different data types (e.g. hMor vs hMor_Beta). The frequency matrices for both of the datatypes need to be calculated already
	#Additionally, This method uses the threshold values calculated in the script "RandomDefinethresh_SEM.py". These values are already hard coded into this script, so there is no need to recalculate.
	
	
	
	#Since the nomenclature and compounds tested vary between the Mu, Delta, and Structure data, these arrays are used to convert the naming
	compoundlist_Mu=['Buprenorphine','DAMGO','Endomorphine1','Fentanyl','Loperamide','Meptazinol','Met-Enk','Morphine','Oxycodone','Tramadol','Cmp 1','Cmp 3','Cmp 4','Cmp 5','Cmp 6','Cmp 7','Cmp 8','Cmp 9','Cmp 10','Cmp 11','Cmp 12','Cmp 13','Cmp 14','Cmp 15','Cmp 16']
	compoundlist_Delta=['Buprenorphine','Fentanyl','Loperamide','Met-Enk','Morphine','Oxycodone','Cmp 1','Cmp 3','Cmp 4','Cmp 5','Cmp 6','Cmp 7','Cmp 8','Cmp 9','Cmp 10','Cmp 11','Cmp 12','Cmp 13','Cmp 14','Cmp 15','Cmp 16']
	compoundlist_Mu_structure=['Buprenorphine','DAMGO','Endomorphin1','Fentanyl','Loperamide','Meptazinol','Met-Enkephalin','Morphine','Oxycodone','Tramadol','Compound-1','Compound-3','Compound-4','Compound-5','Compound-6','Compound-7','Compound-8','Compound-9','Compound-10','Compound-11','Compound-12','Compound-13','Compound-14','Compound-15','Compound-16']
	compoundlist_Mu_pfizer=['Buprenorphine','DAMGO','Endomorphin1','Fentanyl','Loperamide','Meptazinol','Met-Enkephalin','Morphine','Oxycodone','Compound 1','Compound 3','Compound 4','Compound 5','Compound 6','Compound 7','Compound 8','Compound 9','Compound 10','Compound 11','Compound 12','Compound 13','Compound 14','Compound 15','Compound 16']
	
	filetype='{0}'.format(sys.argv[1])#first data type you wish to compare
	filetype2='{0}'.format(sys.argv[2])#second datatype you wish to compare
	diffdict={}
	#These arrays are used to convert between naming in the various different datatype (e.g. convert from hMor to structure)
	morsamples=[0,3,4,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
	pfizersamples=[0,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
	dorsamples=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
	moronly=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
	

	###Below there are several paths the code may take depending on which filetypes are being compared. Since there are only 21 comps in Delta and 25 in everything else, need to ensure that the correct cmps are being compared
	#across the two data sets. For example, row 5 may be cmp1 in delta, but row 5 may be cmp7 in Mu, therefore need to correctly compare. Additionally, structural files used a different naming system, so correct for this as well
	simdict1=defaultdict(list)
	simdict2=defaultdict(list)
	valuedict1=[]
	valuedict2=[]
	#Parsing over the Frequency Matrix of the first datatype
	if filetype=='structure':#Different location for the Structure Clusters
		file=open('../compoundStructure/Clusters/FrequencyMatrix_Average.txt'%(filetype))
	else:
		file=open('../SEM/%s/TotalFrequencyMatrix_Average.txt'%(filetype))
	if filetype2[1]=='D' and filetype[1]!='D':
		print('deltareceptor')
		linecounter=0
		for line in file:
			if linecounter>0:
				parts=line.strip('\n').split('\t')
				partcounter=0
				for i in parts[1:]:
					if i!='' and partcounter in morsamples and linecounter-1 in morsamples:#Use the two arrays above to mask and correctly name the compounds
						#print(i)
						simdict1[parts[0]].append(float(i))
						valuedict1.append(float(i))
					partcounter+=1
			linecounter+=1
	
	else:
		linecounter=0
		for line in file:
			if linecounter>0:
				parts=line.strip('\n').split('\t')
				for i in parts[1:]:
					if i!='':
						#print(i)
						simdict1[parts[0]].append(float(i))
						valuedict1.append(float(i))
			linecounter+=1
	print('simdict1',simdict1)
	
	#Parsing over the frequency matrix of the second datatype
	if filetype2=='structure':
		file=open('../compoundStructure/Clusters/FrequencyMatrix_Average.txt')
	else:
		file=open('../SEM/%s/TotalFrequencyMatrix_Average.txt'%(filetype2))
	if filetype[1]=='D' and filetype2[1]!='D': #If the first file type is a Delta receptor and the second is not, need to handle differently since there are fewer compounds tested in delta, so have to ensure that the same compounds are being compared. 
		if filetype2=='Pfizer' or filetype2=='Pfizer_BRET':
			print('first option is not a delta receptor')
			linecounter=0
			for line in file:
				if linecounter>0:
					parts=line.strip('\n').split('\t')
					#print(parts)
					for i in parts[1:]:
						if i!='':
							#print(i)
							simdict2[parts[0]].append(float(i))
							valuedict2.append(float(i))

				linecounter+=1
		else:	
			linecounter=0
			for line in file:
				if linecounter>0:
					parts=line.strip('\n').split('\t')
					print(parts)
					partcounter=0
					for i in parts[1:]:
						if i!='' and partcounter in morsamples and linecounter-1 in morsamples:
							#print(i)
							simdict2[parts[0]].append(float(i))
							valuedict2.append(float(i))
						partcounter+=1

				linecounter+=1
	else:
		print('first option is not a delta receptor')
		linecounter=0
		for line in file:
			if linecounter>0:
				parts=line.strip('\n').split('\t')
				#print(parts)
				for i in parts[1:]:
					if i!='':
						#print(i)
						simdict2[parts[0]].append(float(i))
						valuedict2.append(float(i))

			linecounter+=1
	print('simdict2',simdict2)
	
	##End result of above are two dictionaries.
	
	##These are the thresholds calculated for each datatype in a different script. The threshold vary depending on the number of clusters used. Here we use 4 clusters, but 3 are also listed
	##A unique threshold exists for each datatype because the the specific values (range, max, min) along with the number of different assays will affect random clustering, predisposing clustering towards recurrent trends, even when random.
	
	####Threshold for K=3
	#threshdict={'4receptors':1.4,'hDor_Beta':0.95,'hDor_Gprot':1.55,'hDor_Gproteins':1.5,'hDor_Kir_CAMP':1.65,'hDor':1.3,'hMor_All4':1.6,'hMor_Beta':1.2,'hMor_CAMP_GAI2':1.7,'hMor_CAMP':2.1,'hMor_GAI2':1.9,'hMor_Gprot':1.7,'hMor_GRK2':1.7,'hMor_GRK6':1.5,'hMor_GRK26_GAI2':1.3,'hMor_GRK26':1.3,'hMor':1.5,'Pfizer_BRET':1.4,'Pfizer':1.6,'rMor':1.5}
	####Thresholds for K=4
	threshdict={'4receptors':1.2,'hDor_Beta':0.75,'hDor_Gprot':1.25,'hDor_Gproteins':1.2,'hDor_Kir_CAMP':1.25,'hDor':1.1,'hMor_All4':1.3,'hMor_Beta':0.9,'hMor_CAMP_GAI2':1.4,'hMor_CAMP':1.6,'hMor_GAI2':1.4,'hMor_Gprot':1.4,'hMor_GRK2':1.2,'hMor_GRK6':1.0,'hMor_GRK26_GAI2':1.1,'hMor_GRK26':1.0,'hMor':1.2,'Pfizer_BRET':1.25,'Pfizer':1.4,'rMor':1.2}
	
	thresh=threshdict[filetype]
	print('The threshold is',thresh)
	changearray=[0,0,0,0]
	###Calculating the difference between individual compounds between the hMor and other clusters	
	differencearray=[]
	
	###As above, there are multiple pathways within the code to ensure that the correct compounds are compared to one another when two different datatypes are used.
	
	if filetype[1]=='D' or filetype2[1]=='D':
		print('Taking pathway 1')
		if filetype2=='Pfizer' or filetype2=='Pfizer_BRET':
			thresh=threshdict[filetype2]
			for i in dorsamples:
				pfvalue=pfizersamples[i]
				for j in dorsamples:
					pfvalue2=pfizersamples[j]

					print(compoundlist_Delta[i],compoundlist_Delta[j])
					print(pfvalue,pfvalue2)
					print(compoundlist_Mu_pfizer[pfvalue],compoundlist_Mu_pfizer[pfvalue2])
					print('value1',simdict1[compoundlist_Delta[i]])
					print('value2',simdict1[compoundlist_Delta[j]])
					print(simdict2[compoundlist_Mu_pfizer[pfvalue]])
					print(simdict2[compoundlist_Mu_pfizer[pfvalue2]])
					
					tempdistance=distance.euclidean(np.array(simdict1[compoundlist_Delta[i]]),np.array(simdict1[compoundlist_Delta[j]]))
					tempdistance2=distance.euclidean(np.array(simdict2[compoundlist_Mu_pfizer[pfvalue]]),np.array(simdict2[compoundlist_Mu_pfizer[pfvalue2]]))
					#print(tempdistance,tempdistance2)
					if tempdistance<thresh:
						if tempdistance2<thresh:
							changearray[0]+=1
						else:
							changearray[1]+=1
					else:
						if tempdistance2<thresh:
							changearray[2]+=1
						else:
							changearray[3]+=1
		else:
			for i in compoundlist_Delta:
				for j in compoundlist_Delta:
					tempdistance=distance.euclidean(np.array(simdict1[i]),np.array(simdict1[j]))
					tempdistance2=distance.euclidean(np.array(simdict2[i]),np.array(simdict2[j]))
					#print(tempdistance,tempdistance2)
					if tempdistance<thresh:
						if tempdistance2<thresh:
							changearray[0]+=1
						else:
							changearray[1]+=1
					else:
						if tempdistance2<thresh:
							changearray[2]+=1
						else:
							changearray[3]+=1
		
						
						
	else:
		checkvalue=0
		if filetype[0]=='P' or filetype2[0]=='P':
			checkvalue=1
			if filetype[0]=='P' and filetype2[0]=='P':
				checkvalue=2
		if filetype[0]=='P' and filetype2[0]=='s':
			checkvalue=3
		else:
			''
		#checkvalue=1
		print('taking pathway 2, checkvalue=%s'%(checkvalue))
		if checkvalue==2:
			print('Both are Pfizer')
			for i in compoundlist_Mu_pfizer:
				for j in compoundlist_Mu_pfizer:
					print(i,j)
					print(simdict1[i])
					print(simdict1[j])
					print(simdict2[i])
					print(simdict1[j])
					tempdistance=distance.euclidean(np.array(simdict1[i]),np.array(simdict1[j]))
					tempdistance2=distance.euclidean(np.array(simdict2[i]),np.array(simdict2[j]))
					if tempdistance<thresh:
						if tempdistance2<thresh:
							changearray[0]+=1
						else:
							changearray[1]+=1
					else:
						if tempdistance2<thresh:
							changearray[2]+=1
						else:
							changearray[3]+=1
		elif checkvalue==1:

			print('Starting here')
			comp=0
			for i in compoundlist_Mu:
				comp1=0
				if i !='Tramadol':
					for j in compoundlist_Mu:
						if j!='Tramadol':
							print(i,j)
							print(comp,comp1)
							print(compoundlist_Mu_pfizer[comp],compoundlist_Mu_pfizer[comp1])
							print(simdict1[i])
							print(simdict1[j])
							print(simdict2[compoundlist_Mu_pfizer[comp]])
							print(simdict2[compoundlist_Mu_pfizer[comp1]])
							tempdistance=distance.euclidean(np.array(simdict1[i]),np.array(simdict1[j]))
							tempdistance2=distance.euclidean(np.array(simdict2[compoundlist_Mu_pfizer[comp]]),np.array(simdict2[compoundlist_Mu_pfizer[comp1]]))
							if tempdistance<thresh:
								if tempdistance2<thresh:
									changearray[0]+=1
								else:
									changearray[1]+=1
							else:
								if tempdistance2<thresh:
									changearray[2]+=1
								else:
									changearray[3]+=1
							comp1+=1
					comp+=1
		elif checkvalue==3:
			print('This is a pfizer structure combo')
			comp=0
			for i in compoundlist_Mu_structure:
				comp1=0
				if i !='Tramadol':
					for j in compoundlist_Mu_structure:
						if j!='Tramadol':
							print(i,j)
							print(comp,comp1)
							print(compoundlist_Mu_pfizer[comp],compoundlist_Mu_pfizer[comp1])
							print(simdict2[i])
							print(simdict2[j])
							print(simdict1[compoundlist_Mu_pfizer[comp]])
							print(simdict1[compoundlist_Mu_pfizer[comp1]])
							tempdistance2=distance.euclidean(np.array(simdict2[i]),np.array(simdict2[j]))
							#tempdistance=distance.euclidean(np.array(simdict2[compoundlist_Mu_pfizer[comp]]),np.array(simdict2[compoundlist_Mu_pfizer[comp1]]))
							tempdistance=distance.euclidean(np.array(simdict1[compoundlist_Mu_pfizer[comp]]),np.array(simdict1[compoundlist_Mu_pfizer[comp1]]))
							if tempdistance<thresh:
								if tempdistance2<thresh:
									changearray[0]+=1
								else:
									changearray[1]+=1
							else:
								if tempdistance2<thresh:
									changearray[2]+=1
								else:
									changearray[3]+=1
							comp1+=1
					comp+=1
		
		else:
			for i in compoundlist_Mu:
				for j in compoundlist_Mu:
					tempdistance=distance.euclidean(np.array(simdict1[i]),np.array(simdict1[j]))
					tempdistance2=distance.euclidean(np.array(simdict2[i]),np.array(simdict2[j]))
					if tempdistance<thresh:
						if tempdistance2<thresh:
							changearray[0]+=1
						else:
							changearray[1]+=1
					else:
						if tempdistance2<thresh:
							changearray[2]+=1
						else:
							changearray[3]+=1
	print('changearray',changearray)
	changearray=[changearray[0]+changearray[3],changearray[1]+changearray[2]]
	print('changearray',changearray)
	for i in changearray:
		print(i/float(sum(changearray)))
	
	
	
	

clustercompare()