import os,math,sys, numpy as np
import csv
import networkx as nx
np.set_printoptions(threshold=np.nan)
import matplotlib
from collections import defaultdict
import glob
from sklearn.decomposition import NMF
from sklearn.cluster import KMeans
from NMF_python_sparse import nmf
import glob
from datetime import date
import time

'''This code is written and published by the lab of Olivier Lichtarge at Baylor College of Medicine for release with the manuscript:

In collaboration with Graciela Pineyro, Besma, and Michel Bouvier
'''
####This method reads in mutation dose response data with SEM values and WT response data.
#It samples the raw data within the bounds of the errors (convert SEM to stdev using #replicates) and then uses this sampled data to perform the clustering algorithm
#For each iteration, the raw sampled data results in frequency matrix. This repeats 1000 times and the end result is how often mutations cluster together when the SEM
#values are accounted for. 
#####must provide a filetype (hMor, hDor, GAI2, etc) 

def readinWT(filetype):#read in compound list, in order
	file=open('../SEM/Data/%s/compoundlist.txt'%(filetype))
	linecounter=0
	mutationlist=[]
	for line in file:
		mutationlist.append(line.strip('\n'))
	return(mutationlist)#returning a list of all compounds/mutations/entities based on the input data

def readreplicatecounts(filetype):#obtain the number of times each mutation was tested. This replicate data is used to convert SEM values to STD values
	replicatedict=defaultdict(list)
	file=open('../SEM/Data/%s/replicatecounts.txt'%(filetype))
	for line in file:
		parts=line.strip('\n').split('\t')
		mut=parts[0]
		for i in parts[1:]:
			replicatedict[mut].append(i)
	print('replicatedict',replicatedict)
	return replicatedict
	

def readinraw(mutationlist,filetype):#obtain the mean values and SEM values for each mutation. Returns a matrix of the raw data and a scalar indicating how many assays are used in the data
	file=open('../SEM/Data/%s/Rawdata.txt'%(filetype))
	for line in file:
		parts=line.strip('\n').split('\t')
		
	assaysize=len(parts)#The size of the data
	X=np.zeros((len(mutationlist),assaysize))#initializing matrix of 0's, these values are updated below
	rowcounter=0
	linecounter=0
	file=open('../SEM/Data/%s/Rawdata.txt'%(filetype))
	for line in file:#Parsing over the data final and appending the data matrix X
		colcounter=0
		parts=line.strip('\n').split('\t')
		'''print('File line',len(parts),parts)'''

		for i in parts:
			checkvalue=''
			try:#Checking if the value is a float. If it is a float, use the float value. If the value is not a float use 'NaN'
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
	print(np.shape(X))
	return(X,assaysize)
	
def sampleraw(rawmat,replicatedict,mutationlist,assaysize):#use SEM values to 'sample' random distribution around each mutation assay mean. Uses SEM values in sampling
	assayarray=range(0,assaysize,2)#Since every mean is paired with a SEM value, create an array that only samples the means (every other value)
	print('assayarray',assayarray)
	erroredmat=np.zeros((len(mutationlist),assaysize/2)) #28 mutations, assays (again, the assaysize is divided by 2 because there is a mean and SEM value
	replicatecountarray=range(2,assaysize/2,3)#create an array that is used to select the correct replicate from the replicate matrix. Note, the numbering is hard coded for 3 measures per assay
	print('replicatecountarray',replicatecountarray)
	for i in range(0,len(mutationlist)):
		colcounter=0
		for j in assayarray:
			stdev=0
			mu=0
			randomvalue=0
			#print(rawmat[i][j+1])
			#print(replicatecounts[i])
			checkvalue=''
			try:#Checking is SEM value is a float. If SEM is not a float, it means the mean value is NaN, therefore there is no sampling. the sampled matrix therefore has a NaN in this position
				float(rawmat[i][j+1])
				checkvalue=True
			except ValueError:
				checkvalue=False
			if checkvalue==False:
				randomvalue=rawmat[i][j]#Sampled matrix has NaN in this position
			print(i,j,rawmat[i][j+1],checkvalue)
			
			else:#If the SEM value is a float this necessitates that the mean value is also a number, we therefore have to sample from the mean and STD
				if rawmat[i][j+1]==0:#if the SEM value is 0 then mean is also 0 (for this data set) Therefore, skip entire process and use 0 as the value. No variation here.
					print('SEM value is 0')
					mu=rawmat[i][j]
					randomvalue=mu
					stdev=0
				else: #These are the measurements with SEM values need to sample
					replicatecount=0
					compound=mutationlist[i]
					for rep in replicatecountarray:#Getting the matched number of replicates from the raw data and replicate matrix using the assay created above
						if colcounter>rep:
							replicatecount+=1
					print(compound,replicatecount)
					replicates=replicatedict[compound][replicatecount]
					stdev=rawmat[i][j+1]*float(replicates)**0.5#converting SEM to stdev using number of replicates
					#stdev=0
					mu=rawmat[i][j]
					print(compound,colcounter,replicatecount,replicates,mutationlist[i],mu,rawmat[i][j+1],stdev)
					#randomvalue=mu
					randomvalue=np.random.normal(mu,stdev,1)#select a single value from the distribution around each mean value
					if stdev=='NaN' or stdev=='nan':
						print('This is not a number!')
			print(mutationlist[i],mu,rawmat[i][j+1],stdev,randomvalue)
			erroredmat[i,colcounter]=randomvalue
			colcounter+=1
	#print(erroredmat)
	return(erroredmat)
	
	
def standardizemat(erroredmat,mutationlist,assaysize):#In order to prevent any bias induced by scale differences between measures, we standardize each column between 0-1 based on the max and min of that column
	standmat=np.zeros((assaysize/2,len(mutationlist)))
	normtrans=np.transpose(erroredmat)
	print('transmatrixshape',np.shape(normtrans))
	rowcounter=0
	for i in normtrans:
		print(i)
		maxvalue=np.nanmax(i)
		minvalue=np.nanmin(i)
		print(maxvalue,minvalue)
		colcounter=0
		for j in i:
			if j!='NaN':
				tempvalue=(float(j)-minvalue)/float(maxvalue-minvalue)
				standmat[rowcounter][colcounter]=tempvalue
			else:
				standmat[rowcounter][colcounter]='NaN'
			
			colcounter+=1
		rowcounter+=1
	standmat=np.transpose(standmat)
	return(standmat)#output is a standardized matrix of the same dimensions as the sampled "erroredmat"


def Freqaverage(repl,mutationlist,filetype,today):#Method to average the frequency matrix over all K's and format the output with labels.
	rep=repl
	nmfmatrix=np.zeros(shape=(len(mutationlist),len(mutationlist)))#initialize matrix
	filecounter=0
	path='../SEM/%s/FrequencyMatrix/'%(filetype)
	for filename in glob.glob(os.path.join(path,'Frequencymatrix_%s_rep=%s*'%(today,rep))):#loop over all frequency matrices within this replicat
		file=open(filename)
		linecounter=0
		tempmatrix=np.loadtxt(filename)
		nmfmatrix=nmfmatrix+tempmatrix
		filecounter=filecounter+1
	linecounter=0
	outfile=open('../SEM/%s/FrequencyMatrix/FrequencyMatrix_rep=%s_Average.txt'%(filetype,rep),'w')#Creating output file
	outfile.write('\t')
	for gene in mutationlist:#appending the names of the compounds
		outfile.write(gene)
		outfile.write('\t')
	outfile.write('\n')
	rowcounter=0
	while rowcounter<len(mutationlist):#adding the averaged values to the file
		outfile.write(mutationlist[rowcounter])#Adding the row labels (compound names)
		outfile.write('\t')
		colcounter=0
		while colcounter<len(mutationlist):
			outfile.write(str((nmfmatrix[rowcounter][colcounter])/float(filecounter)))
			outfile.write('\t')
			colcounter=colcounter+1
		outfile.write('\n')
		rowcounter=rowcounter+1
	nmfmatrix=nmfmatrix/float(filecounter)
	return nmfmatrix
	
def main():
	##############
	##This is the main subroutine that calls all of the above scripts. Need to specify the type of data, e.g hMor or hDor or hMor_Beta etc.
	filetype='{0}'.format(sys.argv[1])
	try:#Creating directories if this is the first time this data type has been run. Files are relative to code location
		os.mkdir('../SEM/%s/'%(filetype))
	except:
		''
	try:
		os.mkdir('../SEM/%s/FrequencyMatrix/'%(filetype))#=
		os.mkdir('../SEM/%s/Standardizedmatrices/'%(filetype))
		os.mkdir('../SEM/%s/Randommatrices/'%(filetype))
	except:
		''
	today=date.fromtimestamp(time.time())
	print(today)
	mutationlist=readinWT(filetype)#getting the list of compounds/mutations/entities, etc
	replicatedict=readreplicatecounts(filetype)#Getting array of the number of replicates for each mutation
	totalfreqmatrix=np.zeros((len(mutationlist),len(mutationlist)))#initializing the final frequency matrix
	rawmat,assaysize=readinraw(mutationlist,filetype)#Obtaining mutant measurements with SEM
	replicates=4#This is where you define the number of total replicates. This will include new sampling from the SEM values within each replicate. 
	for repl in range(0,replicates):
		numit=4#Within each replicate you need to define the number of iterations to run the NMF/Kmeans method at each K
		erroredmat=sampleraw(rawmat,replicatedict,mutationlist,assaysize)#Created matrix of data sampled around the mean for each assay using the SEM values and the #replicates
		np.savetxt('../SEM/%s/Randommatrices/randommatrix_%s_%s.txt'%(filetype,today,str(repl)),erroredmat)

		#normmat=createnormal(WTmat,erroredmat)#use WT values to convert matrix values to normalized difference. (mut-wt)/(mut+wt)
		#np.savetxt('../Output/Normaldiffmatrices/Normdiffmatrix_%s.txt'%(str(repl)),normmat)
		standmat=standardizemat(erroredmat,mutationlist,assaysize)#standize the matrix values using the max and min within each assay column (value-min)/(max-min) where max and min are column specific
		np.savetxt('../SEM/%s/Standardizedmatrices/Standardmatrix_%s_%s.txt'%(filetype,today,str(repl)),standmat)
		freqmatrix=np.zeros((len(mutationlist),len(mutationlist)))
		##############################################
		###Begin clustering algorithm using the standardized sampled matrix you just created above
		nmfmatrix=np.zeros(shape=(len(mutationlist),len(mutationlist)))
		for basis in range(2,8):#Will perform nmf/kmeans method over each K between this range
			for clustk in range(2,8):
				if basis==clustk:#Only running if K=K for nnmf and Kmeans. This could be removed if you want all iterations 2,2; 2,3; 2,4;....8,7; 8,8
					freqmatrix=np.zeros((len(mutationlist),len(mutationlist)))
					for i in range(0,numit):
						print('This is replicate',repl,'and iteration=',i,' with basis and clustk=',basis)
							
							
						W,H=nmf(standmat,basis)#Creating basis vectors. ONly W is used in this script, but H is used in the Assay clustering script
						clusterdata=KMeans(n_clusters=clustk).fit_predict(W)
			
						for j in range (0,len(clusterdata)):#Quantifying if two compounds occur in the same cluster or not. Appending a an overall frequency matrix specific to this replicate.
															#So for each replicate, a frequency matrix created which calculates how often all compounds cluster together over the number of numit for Each K
							for l in range (0,len(clusterdata)):
								if clusterdata[j]==clusterdata[l]:
									tempvalue=freqmatrix[j][l]
									freqmatrix[j][l]=tempvalue+1/float(numit)#dividing by the number of iterations so the average will be between 0-1
				np.savetxt('../SEM/%s/FrequencyMatrix/Frequencymatrix_%s_rep=%s_basis=%s.txt'%(filetype,today,str(repl),str(basis)),freqmatrix)
		nmfmatrix=Freqaverage(repl,mutationlist,filetype,today)#averaging over all K
		np.savetxt('../SEM/%s/FrequencyMatrix/Frequencymatrix_%s_rep=%s_Average.txt'%(filetype,today,str(repl)),nmfmatrix)
		totalfreqmatrix=totalfreqmatrix+nmfmatrix/replicates#This line is a running average over all replicates. This is the final frequency matrix that is the culmination of this entire method
	outfile=open('../SEM/%s/TotalFrequencyMatrix_Average.txt'%(filetype),'w')#Saving the final frequency matrix in a refined format with labels
	linecounter=0
	outfile.write('\t')
	for gene in mutationlist:
		outfile.write(gene)
		outfile.write('\t')
	outfile.write('\n')
	rowcounter=0
	while rowcounter<len(mutationlist):
		outfile.write(mutationlist[rowcounter])
		outfile.write('\t')
		colcounter=0
		while colcounter<len(mutationlist):
			outfile.write(str((totalfreqmatrix[rowcounter][colcounter])))
			outfile.write('\t')
			colcounter=colcounter+1
		outfile.write('\n')
		rowcounter=rowcounter+1
		
		
		
		
main()