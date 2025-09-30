from pyopenms import *
from pathlib import Path
import bisect
import logging

LOGGER = logging.getLogger(__name__)

min_pTIC = .05
max_pTIC = .95

# Save memory by only loading MS1 spectra into memory
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

def peakPick(file_name, folder_loc, top_n):
	"""Performs MS1 feature detection on input file and saves output. TODO fill
	in more.

	Parameters
	----------
	file_name : str
		Name of mzML file to convert to MS1 feature file.
	folder_loc : str
		Location of folder to save output file.
	top_n : int
		Top N most intense MS1 features to save.

	Returns
	-------
	"""
	fh = MzMLFile()
	fh.setOptions(options)

	# Load data
	input_map = MSExperiment()
	fh.load(file_name, input_map)

	# calc TIC
	totalTic = calcTIC(input_map)

	# convert TIC to pTIC
	rtList = []
	pTicList = []
	sumTIC = 0.0
	LOGGER.info("Converting TIC to pTIC")
	for scan in input_map:
		if scan.getMSLevel() == 1:
			mz, i = scan.get_peaks()

			rtList.append(scan.getRT())
			pTicList.append(sumTIC/totalTic)
			sumTIC += sum(i)

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
	LOGGER.info("Found %s features", features.size())

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
			feature = (curMz, curIntens, curRt, pTIC, curCharge)
			featureList.append(feature)

	# TODO the version that we ran for the paper on calculated pTIC
	# on features that were kept (ie denom only contained top N intensity)
	# This is kind of odd and doesn't feel right Need to test if this is better.

	# sort feature list by intensity
	featureList.sort(key=lambda x:x[1], reverse=True)

	# keep N most intense
	intens_features = featureList[0:top_n]

	# sort feature list by m/z
	intens_features.sort(key=lambda x:x[0])

	# print MS1 peak file
	newFileName = folder_loc + "/" + str(Path(file_name).stem) +\
"_ms1Peak.txt"
	with open(newFileName,'w') as newFile:
		newFile.write("mz\tintensity\tRT\tpTIC\tcharge\n")
		for item in intens_features:
			printLine = '\t'.join(str(x) for x in item)
			newFile.write(printLine + '\n')
