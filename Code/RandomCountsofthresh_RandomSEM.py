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
	compoundlist_Mu=['Buprenorphine','DAMGO','Endomorphine1','Fentanyl','Loperamide','Meptazinol','Met-Enk','Morphine','Oxycodone','Tramadol','Cmp 1','Cmp 3','Cmp 4','Cmp 5','Cmp 6','Cmp 7','Cmp 8','Cmp 9','Cmp 10','Cmp 11','Cmp 12','Cmp 13','Cmp 14','Cmp 15','Cmp 16']
	compoundlist_Delta=['Buprenorphine','Fentanyl','Loperamide','Met-Enk','Morphine','Oxycodone','Cmp 1','Cmp 3','Cmp 4','Cmp 5','Cmp 6','Cmp 7','Cmp 8','Cmp 9','Cmp 10','Cmp 11','Cmp 12','Cmp 13','Cmp 14','Cmp 15','Cmp 16']
	compoundlist_Mu_structure=['Buprenorphine','DAMGO','Endomorphin1','Fentanyl','Loperamide','Meptazinol','Met-Enkephalin','Morphine','Oxycodone','Tramadol','Compound-1','Compound-3','Compound-4','Compound-5','Compound-6','Compound-7','Compound-8','Compound-9','Compound-10','Compound-11','Compound-12','Compound-13','Compound-14','Compound-15','Compound-16']
	compoundlist_Mu_pfizer=['Buprenorphine','DAMGO','Endomorphin1','Fentanyl','Loperamide','Meptazinol','Met-Enkephalin','Morphine','Oxycodone','Compound 1','Compound 3','Compound 4','Compound 5','Compound 6','Compound 7','Compound 8','Compound 9','Compound 10','Compound 11','Compound 12','Compound 13','Compound 14','Compound 15','Compound 16']
	filetype='{0}'.format(sys.argv[1])
	filetype2='{0}'.format(sys.argv[2])
	diffdict={}
	morsamples=[0,3,4,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
	pfizersamples=[0,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
	dorsamples=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
	moronly=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
	


	simdict1=defaultdict(list)
	simdict2=defaultdict(list)
	valuedict1=[]
	valuedict2=[]
	if filetype=='structure':
		file=open('../compoundStructure/Clusters_2data/FrequencyMatrix_Average.txt'%(filetype))
	else:
		file=open('../SEM/%s/FrequencyMatrix_SuperAverage.txt'%(filetype))
	if filetype2[1]=='D' and filetype[1]!='D':
		print('deltareceptor')
		linecounter=0
		for line in file:
			if linecounter>0:
				parts=line.strip('\n').split('\t')
				partcounter=0
				for i in parts[1:]:
					if i!='' and partcounter in morsamples and linecounter-1 in morsamples:
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
	#print('simdict1',simdict1)
	
		
	totalchangearray=[]
	path='../SEM/RandomClustering/%s/FrequencyMatrix/'%(filetype2)
	for filename in glob.glob(os.path.join(path,'FrequencyMatrix_*Average.txt')):
		simdict2=defaultdict(list)
		file=open(filename)
		print(filename)
		if filetype[1]=='D' and filetype2[1]!='D':
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
					if linecounter>=0:
						parts=line.strip('\n').split('\t')
						#print(parts)
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
		#print('simdict2',simdict2)
		#for i in simdict2:
			#print(i,len(simdict2[i]))
		####Threshold for K=3
		threshdict={'4receptors':1.4,'hDor_Beta':0.95,'hDor_Gprot':1.55,'hDor_Gproteins':1.5,'hDor_Kir_CAMP':1.65,'hDor':1.3,'hMor_All4':1.6,'hMor_Beta':1.2,'hMor_CAMP_GAI2':1.7,'hMor_CAMP':2.1,'hMor_GAI2':1.9,'hMor_Gprot':1.7,'hMor_GRK2':1.7,'hMor_GRK6':1.5,'hMor_GRK26_GAI2':1.3,'hMor_GRK26':1.3,'hMor':1.5,'Pfizer_BRET':1.4,'Pfizer':1.6,'rMor':1.5}
		####Thresholds for K=4
		threshdict={'4receptors':1.2,'hDor_Beta':0.75,'hDor_Gprot':1.25,'hDor_Gproteins':1.2,'hDor_Kir_CAMP':1.25,'hDor':1.1,'hMor_All4':1.3,'hMor_Beta':0.9,'hMor_CAMP_GAI2':1.4,'hMor_CAMP':1.6,'hMor_GAI2':1.4,'hMor_Gprot':1.4,'hMor_GRK2':1.2,'hMor_GRK6':1.0,'hMor_GRK26_GAI2':1.1,'hMor_GRK26':1.0,'hMor':1.2,'Pfizer_BRET':1.25,'Pfizer':1.4,'rMor':1.2}
		thresh=threshdict[filetype2]
		print('The threshold is',thresh)
		changearray=[0,0,0,0]
		###Calculating the difference between individual compounds between the hMor and other clusters	
		differencearray=[]
		if filetype[1]=='D' or filetype2[1]=='D':
			if filetype2=='Pfizer' or filetype2=='Pfizer_BRET':
		
				for i in dorsamples:
					pfvalue=pfizersamples[i]
					for j in dorsamples:
						print(i,j)
						pfvalue2=pfizersamples[j]
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
			else:
				''
			if checkvalue==2:
				print('Both are Pfizer')
				for i in compoundlist_Mu_pfizer:
					for j in compoundlist_Mu_pfizer:
						#print(i,j)
						#print(simdict1[i])
						#print(simdict1[j])
						#print(simdict2[i])
						#print(simdict1[j])
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
								'''print(i,j)
								print(comp,comp1)
								print(compoundlist_Mu_pfizer[comp],compoundlist_Mu_pfizer[comp1])
								print(simdict1[i])
								print(simdict1[j])
								print(simdict2[compoundlist_Mu_pfizer[comp]])
								print(simdict2[compoundlist_Mu_pfizer[comp1]])'''
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
		
			else:
				comp=0
				for i in compoundlist_Mu:
					comp1=0
					for j in compoundlist_Mu:
						'''print(comp,comp1)
						print(compoundlist_Mu[comp],compoundlist_Mu[comp1])
						print(simdict1[i])
						print(simdict1[j])
						print(simdict2[compoundlist_Mu[comp]])
						print(simdict2[compoundlist_Mu[comp1]])'''
						tempdistance=distance.euclidean(np.array(simdict1[i]),np.array(simdict1[j]))
						tempdistance2=distance.euclidean(np.array(simdict2[compoundlist_Mu[comp]]),np.array(simdict2[compoundlist_Mu[comp1]]))
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
		print('changearray',changearray)
		changearray=[changearray[0]+changearray[3],changearray[1]+changearray[2]]
		print('changearray',changearray)
		totalchangearray.append(changearray[1]/float(sum(changearray)))
		for i in changearray:
			print(i/float(sum(changearray)))
	
	
	plt.hist(totalchangearray)
	plt.title('Random Difference between hMor and %s'%(filetype2))
	plt.savefig('../SEM/RandomClustering/%s/DifferencesfromRandom_%s_%s.png'%(filetype2,filetype,filetype2))
	plt.close()
	plt.show()
	
	outfile=open('../SEM/RandomClustering/%s/DifferencesfromRandom_%s_%s.txt'%(filetype2,filetype,filetype2),'w')
	outfile.write('mean=%s'%(np.mean(totalchangearray)))
	outfile.write('\nstd=%s'%(np.std(totalchangearray)))
	for i in totalchangearray:
		outfile.write('\n%s'%(str(i)))

clustercompare()