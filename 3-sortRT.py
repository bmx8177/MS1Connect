import pandas as pd
import argparse
from pathlib import Path

def sortByRT(edgeFileName):
	file1 = pd.read_csv(edgeFileName,sep='\t')
	file1 = file1.sort_values(['leftFileRT','leftFileIndex'])

	newFileName = Path(edgeFileName).name
	file1.to_csv(newFileName,sep='\t',index=False)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='sort edge file by leftFileRT')
	parser.add_argument("edgeFile",help='edge file')
	args = parser.parse_args()
	sortByRT(args.edgeFile)
