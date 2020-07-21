import os,math,sys, numpy as np
import random
from collections import defaultdict


def makingrandomdistributionsofdifferences():
	##Assays needed for hMor figure
	##Assays for comparison between receptors
	assayarray=['4receptors','hMor','Pfizer','Pfizer_BRET']
	
	differencedict=defaultdict(list)
	poscounter=0
	for i in assayarray:
		for j in assayarray[poscounter+1:]:
			if i!=j:
				check='False'
				check1='False'
				#print('../SEM/RandomClustering/%s/DifferencesfromRandom_hDor_%s.txt'%(i,i))
				try:
					file=open('../compoundStructure/Clusters/DifferencesfromRandom_%s.txt'%(i))
					check='True'
				except:
					check='False'
					print('Could not find the file for this assay ',i)
				try:
					file1=open('../compoundStructure/Clusters/DifferencesfromRandom_%s.txt'%(j))
					check1='True'
				except:
					check1='False'
					print('Could not find the file for this assay2 ',j)
					
					
				if check=='True' and check1=='True':
					filearray=[]
					linecounter=0
					for line in file:
						if linecounter>1:
							filearray.append(float(line.strip('\n')))
						linecounter+=1
					print("The array for assay", i, len(filearray))
				
					filearray2=[]
					linecounter=0
					for line in file1:
						if linecounter>1:
							filearray2.append(float(line.strip('\n')))
						linecounter+=1
					print("The array for assay2", j, len(filearray2))
					
					differencedict[i+'-'+j]=[]
					
					numit=1000
					numcounter=0
					while numcounter<numit:
						index=i+'-'+j
						print(index)
						randiff=random.choice(filearray)-random.choice(filearray2)
						differencedict[index].append(randiff)
						numcounter+=1
		poscounter+=1
	for i in differencedict:
		print(i)
	outfile=open('../compoundStructure/Clusters/Differences_Significances.txt','w')			
	for i in differencedict:
		print('This is the index',i)
		#print(differencedict[i])
		outfile.write(i)
		for j in differencedict[i]:
			outfile.write('\t'+str(j))
		outfile.write('\n')
makingrandomdistributionsofdifferences()
