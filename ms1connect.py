import argparse
from pathlib import Path
from bin import ms1_feature_detection

def ms1Connect(mzml_folder, ms1_folder, top_n):
	'''
	Main script for MS1Connect
	'''
	# Perform MS1 feature detection on mzML files
	# Keeps top N most intense MS1 features per file
	# Writes each output file to disk
	print(top_n)
	for f in Path(mzml_folder).glob("**/*mzML"):
		if Path(ms1_folder + "/" + f.stem + "_ms1Peak.txt").is_file():
			continue

		ms1_feature_detection.peakPick(str(f), ms1_folder, top_n)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Runs MS1Connect on a set of \
mzML files.")
	parser.add_argument("mzml", help="Folder containing mzML files")
	parser.add_argument("ms1", help="Folder containing MS1 feature files")
	parser.add_argument("--topN", help="Keep top N most intense MS1 features.\
	Default=4000", type=int, default=4000)
	args = parser.parse_args()
	ms1Connect(args.mzml, args.ms1, args.topN)
