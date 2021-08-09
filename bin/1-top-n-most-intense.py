import pandas as pd
import argparse
from pathlib import Path

def keep_topN_most_intense_peaks(curFileName,numKeep):
	curFile = pd.read_csv(curFileName,sep='\t')
	curFile = curFile.sort_values(by=['intensity'],ascending=False)
	curFile_topN = curFile[:numKeep]

	curFile_topN = curFile_topN.sort_values(by=['RT'])
	totalIntens = curFile_topN['intensity'].sum()
	curFile_topN['pTIC'] = curFile_topN['intensity'].cumsum() / totalIntens

	curFile_topN = curFile_topN.sort_values(by=['mz'])

	newFileName = Path(curFileName).name
	curFile_topN.to_csv(newFileName,index=False,sep='\t')

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Determines how similar two MS1 runs are based on proportion of MS1 peaks in common')
	parser.add_argument("ms1File",help='MS1 feature file')
	parser.add_argument('topNIntens',help='Keep top N most intense peaks',type=int)
	args = parser.parse_args()
	keep_topN_most_intense_peaks(args.ms1File,args.topNIntens)
