import argparse
import glob
import os
from datetime import datetime
import subprocess
import re

# default values used for 20200416-LB-samples
# mzTolerance = .01 
# ticTolerance = 0.5

# clean file name needs to be changed
def cleanFileName(fileName):
	fileName = os.path.basename(fileName)
	fileName = re.sub('_ms1Peak.txt','',fileName)
	return(fileName)

def createEdges(inputFolderName,outputFolderName,binaryPath,mzTolerance,ticTolerance):
	if os.path.isdir(outputFolderName) == False:
		os.mkdir(outputFolderName)

	fileList = [] # list of file names
	for f in glob.glob(inputFolderName+'/*_ms1Peak.txt'):
		fileList.append(f)

	fileList.sort()
	newFileCnt = 0
	for i in range(0,len(fileList)):
		for j in range(i,len(fileList)):
			leftFile = fileList[i]
			rightFile = fileList[j]
			outFile = outputFolderName + '/' + \
					  cleanFileName(leftFile) + "___" + \
					  cleanFileName(rightFile) + "___score.txt"
			mzTol = mzTolerance
			ticTol = ticTolerance

			if os.path.isfile(outFile) == False:
				subprocess.call([binaryPath,leftFile,rightFile,outFile,str(mzTol),str(ticTol)])

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Determine set of edges between each pairwise input of MS1 feature files. MS1 feature files come from OpenMS')
	parser.add_argument("inputFolder",help="Folder containing sorted MS1 peak files")
	parser.add_argument("outputFolder",help='Folder to put output')
	parser.add_argument("createEdgeBinary",help='Path to createEdge binary')
	parser.add_argument("mzTol",help='m/z tolerance for valid edge',type=float)
	parser.add_argument("ticTol",help='TIC tolerance for valid edge',type=float)
	args = parser.parse_args()
	createEdges(args.inputFolder,args.outputFolder,args.createEdgeBinary,args.mzTol,args.ticTol)
