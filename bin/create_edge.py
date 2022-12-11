import argparse
import glob
import os
import subprocess
import re
from pathlib import Path

# clean file name needs to be changed
def cleanFileName(fileName):
	fileName = str(Path(fileName).stem)
	fileName = re.sub('_ms1Peak','',fileName)
	return(fileName)

def create_edges(inputFolderName,outputFolderName,binaryPath, mz_tol,\
				 tic_tol):
	if os.path.isdir(outputFolderName) == False:
		os.mkdir(outputFolderName)

	fileList = [] # list of file names
	for f in glob.glob(inputFolderName+'/*_ms1Peak.txt'):
		fileList.append(f)

	# TODO need better way to check if two files have been compared
	# ie x___y___score.txt is the same as y___x___.score.txt
	fileList.sort()
	for i in range(0,len(fileList)):
		for j in range(i,len(fileList)):
			leftFile = fileList[i]
			rightFile = fileList[j]
			outFile = outputFolderName + '/' + \
					  cleanFileName(leftFile) + "___" + \
					  cleanFileName(rightFile) + "___score.txt"

			if Path(outFile).is_file():
				continue
			subprocess.call([binaryPath, leftFile, rightFile, outFile,
							str(mz_tol), str(tic_tol)])
