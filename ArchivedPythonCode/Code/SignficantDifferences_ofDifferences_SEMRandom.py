import os,math,sys, numpy as np
import random
from collections import defaultdict


def makingrandomdistributionsofdifferences():
	###This method creates a random distribution representing the expected clustering difference of differences between to datasets (e.g. If hMor vs hMor_Beta= 0.05 and hMor vs hMor_Gprot =0.25 the difference is 0.20
	#This method provides a distribution to tell if 0.20 is significantly different. This is accomplished by comparing the differences from random of hMor vs hMor_Beta and hMor vs hMor_Gprot. These differences
	#should have already been calculated using "RandomCountsofthresh_RandomSEM.py". If by chance the datatypes for hMor_Beta and hMor_Gprot just happen to result in a difference of 0.20, this result is not significantly different.
	#The datatypes force this difference. However, if there is no innate bias the difference of differences should be close to 0. This method enables the calculation of a mean and std in order to compare random distribution to observed values.
	
	
	
	#This method loops over datatypes listed in the arrays below. A custom array can be made to suit the user
	
	#assayarray=['hMor_GAI2','hMor_CAMP','hMor_GRK2','hMor_GRK6','hMor_GRK26','hMor_GRK26_GAI2','hMor_CAMP_GAI2','hMor_Beta','hMor_Gprot','hMor_All4',]
	assayarray=['hDor_Beta','hDor_Gprot','hDor_Kir_CAMP','hDor_Gproteins']
	
	differencedict=defaultdict(list)
	poscounter=0
	for i in assayarray:
		for j in assayarray[poscounter+1:]:#For every combination of datatypes in the array
			if i!=j:
				check='False'
				check1='False'
				try:
					file=open('../SEM/RandomClustering/%s/DifferencesfromRandom_hDor_%s.txt'%(i,i))#Try opening the file which holds the array of differences for the first datatype
					check='True'
				except:
					check='False'
					print('Could not find the file for this assay ',i,' Have you already created the difference array?')#This file does not exist. Perhaps you need to run the code listed above first
				try:
					file1=open('../SEM/RandomClustering/%s/DifferencesfromRandom_hDor_%s.txt'%(j,j))#Try opening the file which holds the array of differences for the first datatype
					check1='True'
				except:
					check1='False'
					print('Could not find the file for this assay2 ',j)
					
				if check=='True' and check1=='True':#If both array files are found, proceed
					filearray=[]
					linecounter=0
					for line in file:#Parse and save as array
						if linecounter>1:
							filearray.append(float(line.strip('\n')))
						linecounter+=1
					print("The array for assay", i, len(filearray))
				
					filearray2=[]
					linecounter=0
					for line in file1:#Parse and save as array
						if linecounter>1:
							filearray2.append(float(line.strip('\n')))
						linecounter+=1
					print("The array for assay2", j, len(filearray2))
					index=i+'---'+j#Make an index for the type datatypes being compared
					differencedict[index]=[]
					
					numit=1000
					numcounter=0
					while numcounter<numit:#For the number of iterations listed above select a random difference value from array 1 and array 2, compare.
						#print(index)
						randiff=random.choice(filearray)-random.choice(filearray2)
						differencedict[index].append(randiff)#Save the random difference of differences in a dictionary
						numcounter+=1
		poscounter+=1
	#print(differencedict['hDor_All4-hDor_Gprot'])
	#print(differencedict)
	outfile=open('../SEM/RandomClustering/hDorDifferences_Significances.txt','w')#Save the random difference of differences			
	for i in differencedict:
		print('This is the index',i)
		#print(differencedict[i])
		outfile.write(i)
		for j in differencedict[i]:
			outfile.write('\t'+str(j))
		outfile.write('\n')
makingrandomdistributionsofdifferences()
