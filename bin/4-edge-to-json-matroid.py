import argparse
import json
import os.path
import re

fileExt=".txt"
outputFolder="matroidFile"

leftFeatureCol = 0
rightFeatureCol = 1

def jsonHelperPartitionMatroid(dic,limit,name):
	blockArray = []
	for key in dic:
		curBlock = {"block":dic[key],
					"limit":limit}
		blockArray.append(curBlock)
	singleMatroid = {
						"partition-matroid": {
							"name":name,
							"blocks":blockArray
						}
					}
	return(singleMatroid)

def makeJson(leftDic,rightDic,limit=1,name="bipartite graph"):
	leftMatroid = jsonHelperPartitionMatroid(leftDic,limit,"edges incident to left nodes,V")
	rightMatroid = jsonHelperPartitionMatroid(rightDic,limit,"edges incident to right nodes,U")

	json_object = {
					"intersection-of-matroids": {
						"name":name,
						"comment":"",
						"partition-matroids":[leftMatroid,rightMatroid]
					}
				  }	
	return(json_object)

def createJsonMatroid(inputFileName):
	if os.path.isdir(outputFolder) == False:
		os.mkdir(outputFolder)
	
	newFileName = os.path.basename(inputFileName)
	newFileName = re.sub(fileExt,'',newFileName)
	newFileName = outputFolder+ "/" + newFileName + "___matroid.json"
	if os.path.exists(newFileName):
		return

	# dic shows which edges are assciated with each MS1 feature
	leftFeatureEdgeDic = {} # key is feature index. value is list of edge indicies
	rightFeatureEdgeDic = {}
	with open(inputFileName,'r') as file1:
		header = file1.readline().strip()
		edgeIndex = 0
		for line1 in file1:
			line1_sp = line1.split('\t')
			
			leftFeatureIndex = line1_sp[leftFeatureCol]
			rightFeatureIndex = line1_sp[rightFeatureCol]

			if leftFeatureIndex not in leftFeatureEdgeDic:
				leftFeatureEdgeDic[leftFeatureIndex] = [edgeIndex]
			else:
				tmpList = leftFeatureEdgeDic[leftFeatureIndex]
				tmpList.append(edgeIndex)
				leftFeatureEdgeDic[leftFeatureIndex] = tmpList

			if rightFeatureIndex not in rightFeatureEdgeDic:
				rightFeatureEdgeDic[rightFeatureIndex] = [edgeIndex]
			else:
				tmpList = rightFeatureEdgeDic[rightFeatureIndex]
				tmpList.append(edgeIndex)
				rightFeatureEdgeDic[rightFeatureIndex] = tmpList
			edgeIndex += 1

	json_object = makeJson(leftFeatureEdgeDic,rightFeatureEdgeDic)
	with open(newFileName,'w') as newFile:
		newFile.write(json.dumps(json_object,indent=2))
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='From edge file create json matroid constraint')
	parser.add_argument("inputFile",help='Input edge file')
	args = parser.parse_args()
	createJsonMatroid(args.inputFile)
