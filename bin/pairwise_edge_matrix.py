import argparse
import numpy as np
from numba import jit
from pathlib import Path
import math
import re
import scipy.sparse
import sys

# edge sim tolerance (start and end time)
startTol = .01

# edge score file extention
ext = ".txt" 

# indicies for MS1 feature file
intensCol = 1
pticCol = 3

# indicies for edge file
leftPeakIndex = 0
rightPeakIndex = 1
ticDiffIndex = 3

# look at fastmath and parallel
# and explicit parllel loops
@jit(nopython=True)
def fillInMatrix(edgeFile, leftFile, rightFile, nRow, \
				 lambda1, lambda2, lambda3, lambda4, \
				 alpha1, alpha2, alpha3):
	"""
	Assume that edge file is sorted by leftFileRT
	"""

	# diagonal values
	countTermArray = np.zeros(nRow,dtype=np.float32)
	intensTermArray = np.zeros(nRow,dtype=np.float32)
	pticTermArray = np.zeros(nRow,dtype=np.float32)
	diagScoreArray = np.zeros(nRow,dtype=np.float32)
	diagArray = np.zeros(nRow,dtype=np.uint32)

	countTermSum = np.float32(0.0); intensTermSum = np.float32(0.0);
	pticTermSum = np.float32(0.0)
	
	# non-diagonal values
	rowList = []; colList = []; 
	valList = [np.float64(x) for x in range(0)]
	edgeSimTermSum = 0.0

	for i in range(0,nRow):
		for j in range(i,nRow):
			edge1 = edgeFile[i]
			edge2 = edgeFile[j]

			edge1Left = int(edge1[leftPeakIndex])
			edge1Right = int(edge1[rightPeakIndex])

			edge2Left = int(edge2[leftPeakIndex])
			edge2Right = int(edge2[rightPeakIndex])

			# two left ms1 features must be within tolerance
			edge1LeftpTIC = leftFile[edge1Left,pticCol]
			edge2LeftpTIC = leftFile[edge2Left,pticCol]

			# two right MS1 features must be within tolerance
			edge1RightpTIC = rightFile[edge1Right,pticCol]
			edge2RightpTIC = rightFile[edge2Right,pticCol]
			
			if i == j: # edge1 == edge2
				# count term
				countTerm = 1.0
				countTermArray[i] = np.float32(countTerm * lambda1)
				countTermSum += np.float32(countTerm)

				# intensity term
				intensTerm = leftFile[edge1Left,intensCol] * \
							 rightFile[edge1Right,intensCol]
				intensTermArray[i] = np.float32(intensTerm * lambda2)
				intensTermSum += np.float32(intensTerm)

				# edge length term
				pticTerm = math.exp(-alpha1 * abs(leftFile[edge1Left,pticCol] - \
									rightFile[edge1Right,pticCol]))
				pticTermArray[i] = np.float32(pticTerm * lambda3)
				pticTermSum += np.float32(pticTerm)

				diagArray[i] = i
			elif (abs(edge1LeftpTIC - edge2LeftpTIC) <= startTol) and \
				 (abs(edge1RightpTIC - edge2RightpTIC) <= startTol):
				# edge shift term
				shiftTerm = math.exp(-alpha2 * abs(edge1[ticDiffIndex] - \
												   edge2[ticDiffIndex]))

				startTerm = math.exp(-alpha3 * abs((edge1LeftpTIC - \
													edge2LeftpTIC)))

				rowList.append(i)
				colList.append(j)
				valList.append(lambda4 * shiftTerm * startTerm)

				rowList.append(j)
				colList.append(i)
				valList.append(lambda4 * shiftTerm * startTerm)

				# sum twice for index i,j and j,i
				edgeSimTermSum += (shiftTerm * startTerm)
				edgeSimTermSum += (shiftTerm * startTerm)

			# since edges are ordered by left file RT 
			if (edge2LeftpTIC - edge1LeftpTIC) > startTol:
				break

	# normalize each array by cumulative sum
	if countTermSum != 0:
		countTermArray = countTermArray / countTermSum
	if intensTermSum != 0:
		intensTermArray = intensTermArray / intensTermSum
	if pticTermSum != 0:
		pticTermArray = pticTermArray / pticTermSum
	if edgeSimTermSum != 0:
		valList = [x/edgeSimTermSum for x in valList]
	
	postTermNorm = (lambda1 * countTermSum) + (lambda2 * intensTermSum) + \
				   (lambda3 * pticTermSum) + (lambda4 * edgeSimTermSum)
	diagScoreArray = countTermArray + intensTermArray + pticTermArray
	return(rowList,colList,valList,diagArray,diagScoreArray,postTermNorm)

leftFileIndex = 0
rightFileIndex = 1
ms1FeatureExt = "_ms1Peak.txt"
def getLeftRightFile(fullFileName,folderName):
	"""
	From the score file name get individual files names
	of the left and right input MS1 feature files
	Input: score file name 
	Input: folder that contains MS1 feature files
	Output: Left and right MS1 feature file name
	"""
	fullFileName_sp = fullFileName.split('___')
	leftFileName = folderName + "/" + fullFileName_sp[leftFileIndex] + ms1FeatureExt
	rightFileName = folderName + "/" + fullFileName_sp[rightFileIndex] + ms1FeatureExt
	return(leftFileName,rightFileName)


###############################################################################
def normalizeIntensity(inputFile):
	"""
	Input file is a MS1 feature file
	Normalize intensities by max value
	"""
	maxIntens = np.max(inputFile,axis=0)[intensCol]
	inputFile[:,intensCol] /= maxIntens


###############################################################################
def createEdgeSimMatrix(edgeFileName,peakFolderName,outputFolderName, \
						lambda1, lambda2, lambda3,lambda4, \
						alpha1, alpha2, alpha3):
	if Path(outputFolderName).is_dir() == False:
		raise Exception(outputFolderName + " does not exist")

	# determine new file name
	edgeFileName_basename =  Path(edgeFileName).stem
	newFileName = Path(outputFolderName) / Path(edgeFileName_basename + "___pairwise.npz")

	# check if output file already exists
	#if Path(newFileName).is_file():
	#	return

	edgeFile = np.loadtxt(edgeFileName,delimiter='\t',skiprows=1,ndmin=2)
	nRow = edgeFile.shape[0]

	if nRow != 0:
		# leftFile and rightFile are MS1 feature files
		leftFileName,rightFileName = getLeftRightFile(edgeFileName_basename,peakFolderName)
		leftFile = np.loadtxt(leftFileName,delimiter='\t',skiprows=1)
		rightFile = np.loadtxt(rightFileName,delimiter='\t',skiprows=1)

		normalizeIntensity(leftFile)
		normalizeIntensity(rightFile)

		# checks that the last edge is between MS1 features
		# that exist in the MS1 feature files
		assert(edgeFile[nRow-1][0] <= leftFile.shape[0])
		assert(edgeFile[nRow-1][1] <= rightFile.shape[0])

		rowList,colList,valList,diagRow,diagScore,postNormVal = \
			fillInMatrix(edgeFile, leftFile, rightFile, nRow,\
						 lambda1,lambda2,lambda3,lambda4,\
						 alpha1,alpha2,alpha3)

		rowList_tmp = np.array(rowList,dtype=np.uint32)
		del rowList
		colList_tmp = np.array(colList,dtype=np.uint32)
		del colList
		valList_tmp = np.array(valList,dtype=np.float32)
		del valList

		rowList_np = np.append(diagRow,rowList_tmp)
		del rowList_tmp
		colList_np = np.append(diagRow,colList_tmp)
		del colList_tmp
		valList_np = np.append(diagScore,valList_tmp)
		del valList_tmp
	else:
		rowList_np = np.array([])
		colList_np = np.array([])
		valList_np = np.array([])

	sparseMat = scipy.sparse.csr_matrix((valList_np, (rowList_np, colList_np)),shape=(nRow,nRow))
	sparseMat.eliminate_zeros()

	#print(edgeFileName_basename,nRow,rowList_np.size,postNormVal)
	scipy.sparse.save_npz(newFileName, sparseMat,compressed=False)
	return(edgeFileName_basename,nRow,rowList_np.size,postNormVal)
