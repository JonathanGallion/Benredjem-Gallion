import os,math,sys, numpy as np
import random
from collections import defaultdict


def makingrandomdistributionsofdifferences():
	##Assays needed for hMor figure
	assayarray=['hDor_GAI2','hDor_CAMP','hDor_GRK2','hDor_GRK6','hDor_GRK26','hDor_GRK26_GAI2','hDor_CAMP_GAI2','hDor_Beta','hDor_Gprot','hDor_All4',]
	###Assays for hDor Figure
	assayarray=['hDor_Beta','hDor_Gprot','hDor_Kir_CAMP','hDor_Gproteins']
	##Assays for comparison between receptors
	assayarray=['hDor_Beta','hDor_Gprot','hDor_Kir_CAMP','hDor_Gproteins']
	
	differencedict=defaultdict(list)
	poscounter=0
	for i in assayarray:
		for j in assayarray[poscounter+1:]:
			if i!=j:
				check='False'
				check1='False'
				#print('../SEM/RandomClustering/%s/DifferencesfromRandom_hDor_%s.txt'%(i,i))
				try:
					file=open('../SEM/RandomClustering/%s/DifferencesfromRandom_hDor_%s.txt'%(i,i))
					check='True'
				except:
					check='False'
					print('Could not find the file for this assay ',i)
				try:
					file1=open('../SEM/RandomClustering/%s/DifferencesfromRandom_hDor_%s.txt'%(j,j))
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
					
					differencedict[i+'---'+j]=[]
					
					numit=1000
					numcounter=0
					while numcounter<numit:
						index=i+'-'+j
						print(index)
						randiff=random.choice(filearray)-random.choice(filearray2)
						differencedict[index].append(randiff)
						numcounter+=1
		poscounter+=1
	#print(differencedict['hDor_All4-hDor_Gprot'])
	#print(differencedict)
	outfile=open('../SEM/RandomClustering/hDorDifferences_Significances.txt','w')			
	for i in differencedict:
		print('This is the index',i)
		#print(differencedict[i])
		outfile.write(i)
		for j in differencedict[i]:
			outfile.write('\t'+str(j))
		outfile.write('\n')
makingrandomdistributionsofdifferences()
