import argparse
import subprocess
from pathlib import Path
from bin import ms1_feature_detection
from bin import create_edge
from bin import edge_to_json_matroid
from bin import pairwise_edge_matrix
from bin import plots


def ms1Connect(mzml_folder, ms1_folder, edge_folder, matroid_folder,
			   edge_sim_folder, output_folder, top_n, mz_tol, tic_tol,
			   metadata_file, lambda1, lambda2, lambda3, lambda4, alpha, beta,
			   gamma):
	'''Main script for MS1Connect.

	Parameters
	----------
	mzml_folder : str, path
		Path of folder that contains mzML files
	ms1_folder : str, path
		Path of folder to put MS1 files. These files will be generated by
		MS1Connect.
	edge_folder : str, path
		Path of folder to put edges files. These files will be generated by
		MS1Connect.
	matroid_folder : str, path
		Path of folder to put matroid files. These files will be generated by
		MS1Connect. Folder will be created if it does not exist.
	edge_sim_folder : str, path
		Path of folder to put sparse edge similarity matricies. These files will
		be generated by MS1Connect. Folder will be created if it does not exist.
	output_folder : str, path
		Path of folder to put output files and plots. This folder wil be created
		if it does not exist.
	top_n : int
		Number of high intensity MS1 features to keep for use in edge generation.
	mz_tol : float
		The m/z tolerance (in ppm) that two MS1 features need to be within in
		order to generate an edge.
	tic_tol : float
		The normalized retention time tolerance that two MS1 features need to be
		within in order to generate an edge.

	Returns
	-------
	'''
	if Path(ms1_folder).exists() == False:
		Path(ms1_folder).mkdir()

	# Perform MS1 feature detection on mzML files
	# Keeps top N most intense MS1 features per file
	# Writes each output file to disk
	for f in Path(mzml_folder).glob("**/*mzML"):
		if Path(ms1_folder + "/" + f.stem + "_ms1Peak.txt").is_file():
			continue
		ms1_feature_detection.peakPick(str(f), ms1_folder, top_n)

	# Generate set of edges from each pair of runs
	create_edge.create_edges(ms1_folder, edge_folder, "bin/createEdge", mz_tol,
		tic_tol)

	# Generate matroid file for each edge file
	if Path(matroid_folder).is_dir() == False:
		Path(matroid_folder).mkdir()
	for f in Path(edge_folder).glob("**/*___score.txt"):
		edge_to_json_matroid.createJsonMatroid(f, matroid_folder)

	# Generate sparse edge similarity matrix for each edge file
	if Path(edge_sim_folder).is_dir() == False:
		Path(edge_sim_folder).mkdir()
	with open("pairwise-edge.log.txt", 'w') as file1, \
		 open("coopraize.log.txt", 'w') as file2:
		#for f in Path(edge_folder).glob("**/*___score.txt"):
		a = Path(edge_folder).glob("**/*___score.txt")
		b = list(a)
		b.sort()
		for f in b:
			fileName, nRow, rowListSize, postNormVal = \
				pairwise_edge_matrix.createEdgeSimMatrix(str(f), 
				ms1_folder, edge_sim_folder, lambda1, lambda2, lambda3, lambda4,
				alpha, beta, gamma)
			file1.write(fileName + '\t' + str(nRow) + '\t' + str(rowListSize) +
						'\t' + str(postNormVal) + '\n')

			npz_file = fileName + "___pairwise.npz"
			matroid_file = fileName + "___matroid.json"

			# log file can be directly generated from coopraize using the below
			# -flogfilename /output/coopraiz_log.txt
			cmd = "singularity exec --bind ./:/input/ --bind " +\
			"./:/output/ bin/coopraiz-singularity " +\
			"/submarine/build/opic-coopraiz -spssdfilename /input/" +\
			edge_sim_folder + "/" + npz_file + " -imjson /input/" +\
			matroid_folder + "/" + matroid_file + " -cloglevel info " +\
			"-ctrl-logsolution -flogtruncate false"
		
			file2.write("filename___" + f.stem + "\n")	
			result = subprocess.run(cmd, shell=True, capture_output=True)
			file2.write(result.stdout.decode())

	if Path(output_folder).is_dir() == False:
		Path(output_folder).mkdir()
	plots.createRunSimMatrix(ms1_folder, "coopraize.log.txt", metadata_file,
							 "pairwise-edge.log.txt", output_folder)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Runs MS1Connect on a set of \
mzML files.")
	parser.add_argument("mzml", help="Folder containing mzML files")
	parser.add_argument("ms1", help="Folder containing MS1 feature files")
	parser.add_argument("edge", help="Folder containing edge files")
	parser.add_argument("matroid", help="Folder containing matroid files.\
						Folder will be created if it does not exist.")
	parser.add_argument("edgeSimMatrix", help='Folder containing sparse edge\
						similarity matricies. Folder will be created if\
						it does not exist.')
	parser.add_argument("metadata",help='Metadata file')
	parser.add_argument("output",help='Folder to place outputs')
	parser.add_argument("--topN", help="Keep top N most intense MS1 features.\
	Default=4000", type=int, default=4000)
	parser.add_argument("--mzTol", help='m/z tolerance in ppm to create an edge.\
 	Default=4', type=float, default=4)
	parser.add_argument("--ticTol",help='TIC tolerance to create an edge.\
	Default=1', type=float, default=1.0)
	parser.add_argument("--lambda1",help='lambda 1 hyperparameter. Default=0.0',
						default=0, type=float)
	parser.add_argument("--lambda2",help='lambda 2 hyperparameter. Default=0.1',
						default=0.1, type=float)
	parser.add_argument("--lambda3",help='lambda 3 hyperparameter. Default=0.0',
						default=0.0, type=float)
	parser.add_argument("--lambda4",help='lambda 4 hyperparameter. Default=0.9',
						default=0.9, type=float)
	parser.add_argument("--alpha",help='alpha hyperparameter. Default=0.0',
						default=0.0, type=float)
	parser.add_argument("--beta",help='beta hypereparameter. Default=1e-5',
						default=0.00001, type=float)
	parser.add_argument("--gamma",help='gamma hyperparameter. Default=1.0',
						default=1.0, type=float)
	args = parser.parse_args()
	ms1Connect(args.mzml, args.ms1, args.edge, args.matroid, args.edgeSimMatrix,
			   args.output, args.topN, args.mzTol, args.ticTol, args.metadata,
			   args.lambda1, args.lambda2, args.lambda3, args.lambda4,
			   args.alpha, args.beta, args.gamma)
