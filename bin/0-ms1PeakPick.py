from pyopenms import *
import argparse
import bisect
from os.path import basename
import re

ext1 = '\.mzML\.gz'
ext2 = '\.mzML'
min_pTIC = .05
max_pTIC = .95

# Prepare data loading (save memory by only
# loading MS1 spectra into memory)
options = PeakFileOptions()
options.setMSLevels([1])

# calc TIC from MS1 scans
def calcTIC(exp):
	tic = 0.0
	for scan in exp:
		if scan.getMSLevel() == 1:
			mz, i = scan.get_peaks()
			tic += sum(i)
	return(tic)

def peakPick(fileName):
	fh = MzMLFile()
	fh.setOptions(options)

	# Load data
	input_map = MSExperiment()
	fh.load(fileName, input_map)

	# calc TIC
	totalTic = calcTIC(input_map)

	# convert TIC to pTIC
	rtList = []
	pTicList = []
	sumTIC = 0.0
	for scan in input_map:
		if scan.getMSLevel() == 1:
			mz, i = scan.get_peaks()

			rtList.append(scan.getRT())
			pTicList.append(sumTIC/totalTic)
			sumTIC += sum(i)
	print("Finished converting TIC to pTIC")
	print()

	input_map.updateRanges()
	ff = FeatureFinder()
	ff.setLogType(LogType.CMD) # progress log

	# Run the feature finder
	name = "centroided"
	features = FeatureMap()
	seeds = FeatureMap()
	params = FeatureFinder().getParameters(name)
	ff.run(name, input_map, features, params, seeds)
	features.setUniqueIds()
	fh = FeatureXMLFile()
	#fh.store("output.featureXML", features)
	print("Found", features.size(), "features")

	# get info for each feature 
	featureList = []
	for f in features:
		curIntens = f.getIntensity()
		curMz = round(f.getMZ(),4)
		curRt = round(f.getRT(),4)
		curCharge = f.getCharge()

		curIndex = bisect.bisect_left(rtList,curRt)
		leftSideRT = rtList[curIndex-1]
		rightSideRT = rtList[curIndex]

		leftSideDiff = curRt - leftSideRT
		rightSideDiff = rightSideRT - curRt

		pTIC = 0 
		if leftSideDiff > rightSideDiff:
			pTIC = round(pTicList[curIndex],4)
		else:
			pTIC = round(pTicList[curIndex-1],4)

		if pTIC >= min_pTIC and pTIC <= max_pTIC:
			feature = (curMz,curIntens,curRt,pTIC,curCharge)
			featureList.append(feature)

	# sort feature list by m/z
	featureList.sort(key=lambda x:x[0])

	# print MS1 peak file
	newFileName = re.sub(ext1,'',basename(fileName))
	newFileName = re.sub(ext2,'',newFileName)+"_ms1Peak.txt"
	with open(newFileName,'w') as newFile:
		newFile.write("mz\tintensity\tRT\tpTIC\tcharge\n")
		for item in featureList:
			printLine = '\t'.join(str(x) for x in item)
			newFile.write(printLine+'\n')


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Take as input a mzML file. Run OpenMS for MS1 peak detection.')
	parser.add_argument('mzMLFile',help='mzML file')
	args = parser.parse_args()
	peakPick(args.mzMLFile)
