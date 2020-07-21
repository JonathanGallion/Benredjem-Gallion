import numpy as np
import scipy
from scipy import sparse
from scipy import linalg
from numpy import dot
from nmf_creatematrix import creatematrix
def nmf(X, numclusters, max_iter=10000, error_limit=1e-7, fit_error_limit=1e-7):
	
	'''print(X)'''
	X=sparse.csr_matrix(X)
    
	eps = 1e-5
	X = X.toarray()  # I am passing in a scipy sparse matrix
    # I am creating a mask matrix where mask[i][j]=1 where X is a number and mask=0 where X is Nan
	masktemp = np.isnan(X)
	rows, columns = masktemp.shape
	mask=np.zeros((rows,columns))
	for i in range(0, rows):
		for j in range(0, columns):
			if masktemp[i][j]==False:
				mask[i][j]=1
			else:
				mask[i][j]=0
    # initial matrices. Both A and Y are initially set as random values between [0,1].
	rows, columns = X.shape
	A = np.random.rand(rows, numclusters)
	#print('A step 1',A)
	A = np.maximum(A, eps)
	#print('A step 2',A)
	
	Y = np.random.rand(numclusters, columns)
	Y = np.maximum(Y, eps)
	
	###CREATING masked_X matrix so that any positition that is a not a number in X is ignored in the following calculations.
	rows, columns = X.shape
	masked_X=np.zeros((rows,columns))
	for i in range(0, rows):
		for j in range(0, columns):
			if mask[i][j]==0:
				masked_X[i][j]=0
			else:
				masked_X[i][j]=X[i][j]
		
	#print('X',X)
	#print('masked_X',masked_X)
	X_est_prev = dot(A, Y)
	for i in range(1, max_iter + 1):
		# ===== updates =====
		# Matlab: A=A.*(((W.*X)*Y')./((W.*(A*Y))*Y'));
		top = dot(masked_X, Y.T)
		bottom = (dot((mask * dot(A, Y)), Y.T)) + eps
		A *= top / bottom

		A = np.maximum(A, eps)
		#print 'A',  np.round(A, 2)

		# Matlab: Y=Y.*((A'*(W.*X))./(A'*(W.*(A*Y))));
		top = dot(A.T, masked_X)
		bottom = dot(A.T, mask * dot(A, Y)) + eps
		Y *= top / bottom
		Y = np.maximum(Y, eps)
		#print 'Y', np.round(Y, 2)


		# ==== evaluation ====
		if i % 5 == 0 or i == 1 or i == max_iter:
			'''print 'Iteration {}:'.format(i),'''
			X_est = dot(A, Y)
			err = mask * (X_est_prev - X_est)
			fit_residual = np.sqrt(np.sum(err ** 2))
			X_est_prev = X_est
			#print('X at the end',X)
			#print('X_est at the end',X_est)
			#print('X-X_est',(masked_X - X_est))
			#curRes = linalg.norm(mask * (X - X_est), ord='fro')
			curRes = linalg.norm(masked_X - X_est, ord='fro')
			'''print(i, 'fit residual', np.round(fit_residual, 6))
			print( 'total residual', np.round(curRes, 6))'''
			if curRes < error_limit or fit_residual < fit_error_limit:
				#if curRes < error_limit:
				#	print('I am stopping because error_limit is minimized',curRes)
				#else:
				#	print('I am stopping because fit_residual is minimized',fit_residual)
				break

	return(A,Y)