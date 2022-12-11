#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <chrono>
#include <assert.h>
#include <algorithm>
using namespace std;

// Inputs
// 1) File containing MS1 features. Columns are mz,
//    intensity, RT, pTIC, and charge. Sorted by mz
// 2) Same type of input as 1)
// 3) Output file name
// 4) m/z tolerance for edge creation (in ppm)
// 5) pTIC tolerance for edge creation

// # of columns in input 1 and 2 (note 0 index)
int numCols = 5;
int pTicCol = 3;
int chargeCol= 4;

struct edge_struct {
	int leftIndex;
	int rightIndex;
	double mzDiff;
	double ticDiff;
	double leftFileRT;
};

bool cmp_edge(const edge_struct& edge1, const edge_struct& edge2)
{
	if (edge1.leftFileRT < edge2.leftFileRT) {
		return true;
	}
	else if (edge1.leftFileRT == edge2.leftFileRT) {
		if (edge1.leftIndex < edge2.leftIndex)
			return true;
		else
			return false;
	}
	else {
		return false;
	}
}

std::string readFileToMatrix(const std::string& fileName,
                      std::vector<double>& lines) {

	std::ifstream inFile(fileName);
	std::string header;	
	std::string line;

	std::getline(inFile,header);
	while (std::getline(inFile,line)) {
		std::stringstream ss(line);
		while (getline(ss,line,'\t')) {
			lines.push_back(std::stof(line));
		}
	}
	return header;
}

double calcPpmDiff(double mass1, double mass2) {
	// Calculate the difference in parts-per-million between two peptide masses.
	// mass 1 is mass of precursor 1
	// mass 2 is mass of precursor 2
	// window is ppm window
	double ppmDiff = (1000000 * (mass1 - mass2)) / (0.5 * (mass1 + mass2));
	return ppmDiff;
}

int main(int argc, char * argv[])
{
	// check that argument count is correct
	if (argc != 6) {
		cout << "Not enough arguments" << std::endl;
		exit(0);
    }

	// get arguments
	std::string fileName1 = argv[1]; // filename of input 1
	std::string fileName2 = argv[2]; // filename of input 2
	std::string outFileName = argv[3]; // filename of output file
	double mzTol = stof(argv[4]); // mz tolerance (in ppm)
	double ticTol = stof(argv[5]); // pTIC tolerance

	std::vector<double> file1Vector;
	std::vector<double> file2Vector;
	std::vector<edge_struct> edgeVector;

	std::string header;
	header = readFileToMatrix(fileName1,file1Vector);
	readFileToMatrix(fileName2,file2Vector);

	std::cout << fileName1 << "\t" << fileName2 << std::endl;

	// calc number of rows
	int numRowsFile1 = file1Vector.size() / numCols;
	int numRowsFile2 = file2Vector.size() / numCols;

	int startIndex = 0; int endIndex = 0; int edgeCnt = 0;
	double leftMz; double rightMz;
	int breakLoop = 0; // becomes 1 when right file index goes past EOF
	int leftCharge; int rightCharge;
	double mzDiff; double ticDiff; double ppmDiff;
	for (int i=0;i<numRowsFile1;i++) {
		// mz is first column
		leftMz = file1Vector[numCols*i];
		rightMz = file2Vector[numCols*startIndex];
		ppmDiff = calcPpmDiff(rightMz,leftMz);

		if (ppmDiff > mzTol) {
			continue;
		}

		// get updated start index of file2
		for (int j=startIndex;j<numRowsFile2;j++) {
			rightMz = file2Vector[numCols*j];
			ppmDiff = calcPpmDiff(leftMz,rightMz);

			if (ppmDiff < mzTol) {
				startIndex = j;
				break;
			}

			if (j == numRowsFile2 - 1 ) { // reached EOF
				breakLoop = 1;
			}
		}
		
		if (breakLoop == 1) { // right file start index past EOF
			break;
		}

		// get update end index of file 2
		for (int j=endIndex; j<numRowsFile2; j++) {
			rightMz = file2Vector[numCols*j];
			ppmDiff = calcPpmDiff(rightMz, leftMz);
			if (ppmDiff >= mzTol) {
				endIndex = j;
				break;
			}

			if (j == numRowsFile2-1) {
				endIndex=j;
			}
		}

		for (int k=startIndex; k<=endIndex; k++) {
			// header line is skipped when reading in file
			mzDiff = calcPpmDiff(leftMz,file2Vector[numCols*k]);

			// need b/c of edges case
			if (abs(mzDiff) > mzTol) {
				continue;
			}

			ticDiff = file1Vector[numCols*i + pTicCol] - file2Vector[numCols*k + pTicCol];
			if (abs(ticDiff) > ticTol) {
				continue;
			}

			leftCharge = file1Vector[numCols*i + chargeCol];
			rightCharge = file2Vector[numCols*k + chargeCol];
			if (leftCharge != rightCharge) {
				continue;
			}

			edge_struct newEdge;
			newEdge.leftIndex = i;
			newEdge.rightIndex = k;
			newEdge.mzDiff = mzDiff;
			newEdge.ticDiff = ticDiff;
			newEdge.leftFileRT = file1Vector[numCols*i + pTicCol];
			edgeVector.push_back(newEdge);
		}
	}

	// sort edges by retention time of the feaure in the left file
	// and write edges to output file
	std::ofstream outFile(outFileName);
	outFile << "leftFileIndex\trightFileIndex\tmzDiff\tticDiff\tleftFileRT" <<
			   std::endl;
	std::sort(edgeVector.begin(), edgeVector.end(), cmp_edge);
	for (auto edge=edgeVector.begin(); edge!=edgeVector.end(); ++edge) {
		outFile << (*edge).leftIndex << '\t'; // left file line index
		outFile << (*edge).rightIndex << '\t'; // right file line index
		outFile << (*edge).mzDiff << '\t';
		outFile << (*edge).ticDiff <<  '\t';
		outFile << (*edge).leftFileRT;
		outFile << std::endl;
		edgeCnt += 1;
	}
	outFile.close();

	//std::cout << "nrow: " << numRowsFile1 << '\t' << numRowsFile2 << std::endl;
	//std::cout << "Total # edges: " << edgeCnt << std::endl;
	//std::cout << "Done" << std::endl << std::endl;
	return 0;
}
