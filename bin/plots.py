import argparse
import math
import numpy as np
import seaborn as sns
import pandas as pd
import re
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
from sklearn.manifold import MDS
from sklearn.metrics.pairwise import euclidean_distances
from sklearn import metrics
from pathlib import Path


###############################################################################
def getFileList(ms1_folder, metadataFileName):
	"""
	Reads metadata file to determine what file exists and what metadata we are
	interested in. Also finds all files in MS1 feature folder. Compares the
	two lists to see if anything is different. The metadata file determines the
	order of the two lists.

	Parameters
	----------
	ms1_folder: str, path
		Path of folder that contains MS1 feature files
	metadataFileName: str, path
		Path of tab delimited file that contains the filename and metadata label

	Returns
	------
	fileList: list
		List of MS1 files found in MS1 feature folder. This list is in the same
		order as metadataList.
		metadataList
	metadataList: list
		List of metadata labels. In same order as fileList.
	"""
	ext = "_ms1Peak"
	delim = "\t"

	# Find list of MS1 feature files
	ms1_file_dic = {}
	for f in Path(ms1_folder).glob("**/*_ms1Peak.txt"):
		f_remove_ext = re.sub(ext, "", str(f.stem))
		ms1_file_dic[f_remove_ext] = 1

	metadata_file_dic = {}
	metadata_file_list = []
	metadata_file = pd.read_csv(metadataFileName,sep='\t')
	for f in metadata_file['fileName']:
		cur_file = str(Path(f).stem)
		metadata_file_dic[cur_file] = 1
		metadata_file_list.append(cur_file)

	if ms1_file_dic.keys() != metadata_file_dic.keys():
		raise Exception("The number of MS1 files does not match number of \
lines in " + metadataFileName)
		
	return(metadata_file_list, list(metadata_file['metadataLabel']))


###############################################################################
def plotMDS(runMatrix, metadata_label, output_folder):
	"""
	Plots MDS on run similarity matrix. Plots a normal MDS.
	If species data also plots a second MDS to better visualize things
	Input1: matrix of pairwise run similarities
	Input2: list of labels (run order in run sim matrix)
	"""
    # MDS of run sim matrix
	dist = euclidean_distances(runMatrix)
	model = MDS(dissimilarity='precomputed',n_components=2,random_state=0,
			    normalized_stress="auto")
	out = model.fit_transform(dist)

    # normal MDS
	fig,ax = plt.subplots()
	ax = sns.scatterplot(x=out[:,0], y=out[:,1],hue=metadata_label)
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.tight_layout()
	fig.savefig(output_folder + "/mds.png")
	plt.close(fig)


###############################################################################
def plotHeatmap(inputRunMatrix, tick_label, output_folder):
	"""
	Plot two different heatmaps of the run similarity matrix
	One heatmap is normal coopraize score heatmap
	Second heatmap is where normalized by the diagonal
	Input1: run similarity matrix
	Input2: order of runs that run similarity matrix is
	"""
	# heatmap of run sim matrix
	inputRunMatrix = np.sqrt(inputRunMatrix)
	vmax = np.percentile(inputRunMatrix,95)
	vmin = np.amin(inputRunMatrix)
	
	fig,ax = plt.subplots()
	ax = sns.heatmap(inputRunMatrix,vmin=vmin,vmax=vmax, \
                     xticklabels=tick_label,yticklabels=tick_label)

	# square the color bar tick label to undo sqrt of sim matrix
	c_bar = ax.collections[0].colorbar
	ticLoc = c_bar.get_ticks()
	newTic = [int(x*x) for x in ticLoc]
	c_bar.set_ticks(ticLoc)
	c_bar.set_ticklabels(newTic)

	plt.tight_layout()
	fig.savefig(output_folder + "/heatmap.png")
	plt.close(fig)

#	fig = sns.clustermap(inputRunMatrix, xticklabels=tickLabel, 
#						 yticklabels=tickLabel, vmax=vmax)
#	fig.savefig("heatmap-dendogram.png")
#	plt.close()


###############################################################################
def fillInSimMatrixCooprize(ms1FileList,scoreFile,matrix):
	"""
	Fill in run similarity matrix based on cooprize output file
    Fill in run similarity matrix based on baseline output file
	Input1: List of ordered file names (same run order as matrix)
	Input2: output of baseline score file
	Input3: run similarity matrix to be filled in
	"""
	coopLeftFileIndex = 1
	coopRightFileIndex = 2
	scoreFileDelim = "___"
	numEdgeDic = {}
	with open(scoreFile,'r') as file1:
		for line1 in file1:
			if line1.startswith("filename___"):
				line1_sp = line1.split(scoreFileDelim)

				leftFile = line1_sp[coopLeftFileIndex]
				rightFile = line1_sp[coopRightFileIndex]

				if leftFile not in ms1FileList or \
				   rightFile not in ms1FileList:
					leftFileIndex = -1
					rightFileIndex = -1
				else:
					leftFileIndex = ms1FileList.index(leftFile)
					rightFileIndex = ms1FileList.index(rightFile)
			elif line1.startswith("Loaded raw SPSSD "):
				line1_sp = line1.split(' ')
				token = line1_sp[3]
				token_sp = token.split("x")
				numEdges = float(token_sp[0])
			elif line1.startswith("Summary valuation:"):
				# to get # selected edges
				line1_sp = line1.split(',')
				token = line1_sp[0]
				token_sp = token.split('=')
				numSelectedEdge = float(token_sp[1])

				# to get score
				line1_sp = line1.split("=")
				lastIndex = len(line1_sp) - 1
				score = float(line1_sp[lastIndex].strip())

				if leftFileIndex == -1 and \
				   rightFileIndex == -1:
					continue
				matrix[leftFileIndex,rightFileIndex] = score
				matrix[rightFileIndex,leftFileIndex] = score 

				# fill in numEdgeDic
				dicKey = str(leftFileIndex) + "-" + str(rightFileIndex)
				numEdgeDic[dicKey] = numEdges
				dicKey = str(rightFileIndex) + "-" + str(leftFileIndex)
				numEdgeDic[dicKey] = numEdges
	return(numEdgeDic)


###############################################################################
def assertDiagonal(curRunMatrix):
	"""
	Checks that all diagonals are non-zero. This is for case we didn't score
	a run but is present in metadata file
	"""
	nRow = len(curRunMatrix)

	# Ignore if only two runs as a 2x2 matrix does not have a diagonal
	if nRow == 2:
		return
	for i in range(0,nRow):
		assert(curRunMatrix[i][i] != 0), "Metadata file contains file that was not present in score file. Remove from metadata file. %s line" %(i)
	return


##############################################################################
def postNormBySetE(inputRunMatrix, setEScores, ms1FileList):
	"""
	Allows for post normalization of coopraize scores by some value
	found in setEScores file
	"""
	ext = "_ms1Peak.txt"
	# currently assumes that post norm value is found in 4th column of file
	with open(setEScores,'r') as normFile:
		for line1 in normFile:
			line1_sp = line1.strip().split('\t')

			curFile = line1_sp[0]
			curFile_sp = curFile.split("___")

			leftFile = curFile_sp[0]
			rightFile = curFile_sp[1]

			leftFileIndex = ms1FileList.index(leftFile)
			rightFileIndex = ms1FileList.index(rightFile)
			
			normVal = float(line1_sp[3])

			if leftFileIndex == rightFileIndex:
				inputRunMatrix[leftFileIndex][rightFileIndex] = \
				inputRunMatrix[leftFileIndex][rightFileIndex] * normVal
			else:
				inputRunMatrix[leftFileIndex][rightFileIndex] = \
				inputRunMatrix[leftFileIndex][rightFileIndex] * normVal

				inputRunMatrix[rightFileIndex][leftFileIndex] = \
				inputRunMatrix[rightFileIndex][leftFileIndex] * normVal


###############################################################################			
def createRunSimMatrix(ms1PeakFolderName, scoreFileName, metadataFileName, \
					   edgeCountFileName, output_folder):
	"""
	Main driver script
	"""
	fileList, metadataList = getFileList(ms1PeakFolderName,metadataFileName)
	runMatrix = np.zeros((len(fileList),len(fileList)))

	# create run similarity matrix
	edgeDic = fillInSimMatrixCooprize(fileList, scoreFileName, runMatrix)
	
	# assert diagonal is non-zero
	assertDiagonal(runMatrix)

	# post-normalize each value in runMatrix
	# by value of last col in file
	postNormBySetE(runMatrix,edgeCountFileName,fileList)

	plotHeatmap(runMatrix, metadataList, output_folder)
	plotMDS(runMatrix, metadataList, output_folder)

	np.savetxt(output_folder + "/output_score_matrix.txt", runMatrix, \
			   delimiter='\t',fmt='%f')


###############################################################################
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='From baseline metric output '+
	"Create a pairwise run similarity matrix. Plot heatmap and " +
	"MDS of this matrix. Metadata file determines order of runs " +
	"in output files.")
	parser.add_argument("ms1PeakFolder", help='Folder containing MS1 peaks')
	parser.add_argument("baselineOutput", help='coopraize score file')
	parser.add_argument("metadataFile", help='Metadata file for LB data')
	parser.add_argument('edgeCountFile', help='pairwise edge count file. ' +
	"This file contains the score used for post-normalization")
	parser.add_argument("output", help="output folder")
	args = parser.parse_args()
	createRunSimMatrix(args.ms1PeakFolder, args.baselineOutput, args.metadataFile,\
					   args.edgeCountFile, args.output)
