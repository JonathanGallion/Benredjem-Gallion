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

####This method reads in mutation dose response data with SEM values and WT response data to identify clusters of assays with similar performances across this dataset
#It samples the raw data within the bounds of the errors (convert SEM to stdev using #replicates) and then uses this sampled data to perform the clustering algorithm
#For each iteration, the raw sampled data results in cluster assignment. This repeats 1000 times and the end result is how often mutations cluster together when the SEM
#values are accounted for. 
#####must provide a filetype (hMor, hDor, GAI2, etc) 

def readinWT(filetype):#read in compound list, in order
	Assaylist=['Beta1--EC50','Beta1--Span','Beta1--LogR','Beta2--EC50','Beta2--Span','Beta2--LogR','GRK2--EC50','GRK2--Span','GRK2--LogR','GRK5--EC50','GRK5--Span','GRK5--LogR','GRK6--EC50','GRK6--Span','GRK6--LogR','EPAC--EC50','EPAC--Span','EPAC--LogR','Gai1--EC50','Gai1--Span','Gai1--LogR','Gai2--EC50','Gai2--Span','Gai2--LogR','Gao--EC50','Gao--Span','Gao--LogR','Kir--EC50','Kir--Span','Kir--LogR']
	print(len(Assaylist))
	file=open('../SEM/Data/%s/compoundlist.txt'%(filetype))
	linecounter=0
	mutationlist=[]
	for line in file:
		mutationlist.append(line.strip('\n'))
	return(mutationlist,Assaylist)

def readreplicatecounts(filetype):#obtain the number of times each mutation was tested
	replicatedict=defaultdict(list)
	file=open('../SEM/Data/%s/replicatecounts.txt'%(filetype))
	for line in file:
		parts=line.strip('\n').split('\t')
		mut=parts[0]
		for i in parts[1:]:
			replicatedict[mut].append(i)
	print(replicatedict)
	return replicatedict
	

def readinraw(mutationlist,filetype):#obtain the mean values and SEM values for each mutation
	#X=np.zeros((25,80))#Use when including hillC
	file=open('../SEM/Data/%s/Rawdata.txt'%(filetype))
	for line in file:
		parts=line.strip('\n').split('\t')
		
	assaysize=len(parts)
	X=np.zeros((len(mutationlist),assaysize))
	rowcounter=0
	linecounter=0
	file=open('../SEM/Data/%s/Rawdata.txt'%(filetype))
	for line in file:
		colcounter=0
		parts=line.strip('\n').split('\t')
		'''print('File line',len(parts),parts)'''

		for i in parts:
			checkvalue=''
			try:
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
	#assayarray=[0, 1, 3, 5, 7, 9, 11, 12, 14,16,18,20, 22, 23,25,  27, 29,  31,  33, 34, 36, 38, 40,  42, 44, 45, 47, 49, 51]#These are the columns of means, ignoring SEM columns
	assayarray=range(0,assaysize,2)
	print('assayarray',assayarray)
	erroredmat=np.zeros((len(mutationlist),assaysize/2)) #28 mutations, 40 assays
	#replicatecountarray=[4,8,12,16,20,24,28,32,36,40]
	replicatecountarray=range(3,assaysize/2,3)
	#print(range(0,52))
	for i in range(0,len(mutationlist)):
		colcounter=0
		for j in assayarray:
			stdev=0
			mu=0
			randomvalue=0
			#print(rawmat[i][j+1])
			#print(replicatecounts[i])
			checkvalue=''
			try:
				float(rawmat[i][j+1])
				checkvalue=True
			except ValueError:
				checkvalue=False
			if checkvalue==False:
				randomvalue=rawmat[i][j]
			print(i,j,rawmat[i][j+1],checkvalue)
			
			if checkvalue==False:
				randomvalue=rawmat[i][j]
			else:
				if rawmat[i][j+1]==0:#if the SEM value is 0 then mean is also 0 (for this data set) Therefore, skip entire process and use 0 as the value. No variation here.
					print('SEM value is 0')
					mu=rawmat[i][j]
					randomvalue=mu
					stdev=0
				else: #These are the measurements with SEM values need to sample
					replicatecount=0
					compound=mutationlist[i]
					for rep in replicatecountarray:
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
	
	
'''def createnormal(WTmat,erroredmat):
	normmat=np.zeros((25,30))
	for i in range(0,28):
		for j in range(0,29):
			wtvalue=WTmat[i][j]
			ranvalue=erroredmat[i][j]
			tempvalue=(ranvalue-float(wtvalue))/float(float(wtvalue)+ranvalue)#normalized difference value
			print (i,j,wtvalue,ranvalue,tempvalue)
			normmat[i][j]=tempvalue
	return(normmat)'''
	
def standardizemat(erroredmat,mutationlist,assaysize):
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
	return(standmat)



def Freqaverage(repl,Assaylist,filetype,today):
	rep=repl
	nmfmatrix=np.zeros(shape=(len(Assaylist),len(Assaylist)))
	filecounter=0
	path='../SEM/%s/Assays/FrequencyMatrix/'%(filetype)
	for filename in glob.glob(os.path.join(path,'Frequencymatrix_%s_rep=%s*'%(today,rep))):
		file=open(filename)
		linecounter=0
		tempmatrix=np.loadtxt(filename)
		nmfmatrix=nmfmatrix+tempmatrix
		filecounter=filecounter+1
	linecounter=0
	outfile=open('../SEM/%s/Assays/FrequencyMatrix/FrequencyMatrix_rep=%s_Average.txt'%(filetype,rep),'w')
	outfile.write('\t')
	for gene in Assaylist:
		outfile.write(gene)
		outfile.write('\t')
	outfile.write('\n')
	rowcounter=0
	while rowcounter<len(Assaylist):
		outfile.write(Assaylist[rowcounter])
		outfile.write('\t')
		colcounter=0
		while colcounter<len(Assaylist):
			outfile.write(str((nmfmatrix[rowcounter][colcounter])/float(filecounter)))
			outfile.write('\t')
			colcounter=colcounter+1
		outfile.write('\n')
		rowcounter=rowcounter+1
	nmfmatrix=nmfmatrix/float(filecounter)
	return nmfmatrix
	
def main():
	filetype='{0}'.format(sys.argv[1])
	try:
		os.mkdir('../SEM/%s/Assays/'%(filetype))
	except:
		''
	try:
		os.mkdir('../SEM/%s/Assays/FrequencyMatrix/'%(filetype))
		os.mkdir('../SEM/%s/Assays/Standardizedmatrices/'%(filetype))
		os.mkdir('../SEM/%s/Assays/Randommatrices/'%(filetype))
	except:
		''
	today=date.fromtimestamp(time.time())
	print(today)
	mutationlist,Assaylist=readinWT(filetype)
	replicatedict=readreplicatecounts(filetype)#Getting array of the number of replicates for each mutation
	totalfreqmatrix=np.zeros((len(Assaylist),len(Assaylist)))
	rawmat,assaysize=readinraw(mutationlist,filetype)#Obtaining mutant measurements with SEM
	replicates=4
	for repl in range(0,replicates):
		numit=4
		erroredmat=sampleraw(rawmat,replicatedict,mutationlist,assaysize)#Created matrix of data sampled around the mean for each assay using the SEM values and the #replicates
		np.savetxt('../SEM/%s/Assays/Randommatrices/randommatrix_%s_%s.txt'%(filetype,today,str(repl)),erroredmat)

		#normmat=createnormal(WTmat,erroredmat)#use WT values to convert matrix values to normalized difference. (mut-wt)/(mut+wt)
		#np.savetxt('../Output/Normaldiffmatrices/Normdiffmatrix_%s.txt'%(str(repl)),normmat)
		standmat=standardizemat(erroredmat,mutationlist,assaysize)#standize the matrix values using the max and min within each assay column (value-min)/(max-min) where max and min are column specific
		np.savetxt('../SEM/%s/Assays/Standardizedmatrices/Standardmatrix_%s_%s.txt'%(filetype,today,str(repl)),standmat)
		freqmatrix=np.zeros((len(Assaylist),len(Assaylist)))
		##############################################
		###Begin clustering algorithm using the standardized sampled matrix you just created above
		nmfmatrix=np.zeros(shape=(len(Assaylist),len(Assaylist)))
		for basis in range(2,8):
			for clustk in range(2,8):
				if basis==clustk:
					freqmatrix=np.zeros((len(Assaylist),len(Assaylist)))
					for i in range(0,numit):
						print('This is replicate',repl,'and iteration=',i,' with basis and clustk=',basis)
							
							
						W,H=nmf(standmat,basis)
						clusterdata=KMeans(n_clusters=clustk).fit_predict(np.transpose(H))
						for j in range (0,len(clusterdata)):
							for l in range (0,len(clusterdata)):
								if clusterdata[j]==clusterdata[l]:
									tempvalue=freqmatrix[j][l]
									freqmatrix[j][l]=tempvalue+1/float(numit)
				np.savetxt('../SEM/%s/Assays/FrequencyMatrix/Frequencymatrix_%s_rep=%s_basis=%s.txt'%(filetype,today,str(repl),str(basis)),freqmatrix)
		nmfmatrix=Freqaverage(repl,Assaylist,filetype,today)
		totalfreqmatrix=totalfreqmatrix+nmfmatrix/replicates
	
	
	#np.savetxt('../SEM/TotalFrequencymatrix_k=3_replicates=%s.txt'%(str(replicates)),totalfreqmatrix)
	outfile=open('../SEM/%s/Assays/TotalFrequencyMatrix_Average.txt'%(filetype),'w')
	linecounter=0
	outfile.write('\t')
	for gene in Assaylist:
		outfile.write(gene)
		outfile.write('\t')
	outfile.write('\n')
	rowcounter=0
	while rowcounter<len(Assaylist):
		outfile.write(Assaylist[rowcounter])
		outfile.write('\t')
		colcounter=0
		while colcounter<len(Assaylist):
			outfile.write(str((totalfreqmatrix[rowcounter][colcounter])))
			outfile.write('\t')
			colcounter=colcounter+1
		outfile.write('\n')
		rowcounter=rowcounter+1
		
		
		
		
main()