//AUTHOR:       Oliver Marx - ojm40@cam.ac.u/k

#include <fstream>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <string>
#include <vector>
#include <array>
#include <sstream>

// Libs
#include <libconfig.h++>
#include <armadillo>

#include "BoundaryPatch.H"
#include "Face.H"
#include "Cell.H"
#include "fvMatrix.H"
#include "fvMesh.H"
#include "fvMeshParser.H"

void outputState(std::ofstream& file, arma::mat x, arma::vec T)
{
	for(unsigned int i=0;i<T.size();i++)
	{
		file << x(i,0) << " "
		     << x(i,1) << " "
		     << T(i) << " "
		     << std::endl;
		if((i+1)%300==0)
			file << "\n";
	}
	file << "\n\n";
}

void initializeState(arma::vec& T, arma::vec& gamma, arma::vec& S,
		arma::mat& u, int& tRes, fvMesh& thisMesh, int testNum)
{
	const char* configFileName = "./stg/config.cfg";
	libconfig::Config cfg;
	
	try{
		cfg.readFile(configFileName);
	}
	catch(const libconfig::FileIOException &fioex){
		std::cerr << "I/O error while reading file." << std::endl;
	}
	catch (libconfig::ParseException &e){
        /*inform user about the parse exception*/
        std::cerr << "Parse error at " << e.getFile() << ":" << e.getLine()
                  << " - " << e.getError() << std::endl;
	}
	
	const libconfig::Setting& root = cfg.getRoot();

	try{
		const libconfig::Setting& tests = root["simulation"]["tests"];
		const libconfig::Setting& test = tests[testNum];

		double T0, gamma0, ux, uy, uz;

		if(!(test.lookupValue("T0", T0)
			&& test.lookupValue("gamma0", gamma0)
			&& test.lookupValue("ux", ux)
			&& test.lookupValue("uy", uy)
			&& test.lookupValue("uz", uz)
			&& test.lookupValue("tRes", tRes)))
			std::cout << "Settings for test " << testNum+1 << " read in."
				<< std::endl;

	unsigned int nCells = gamma.size();

	for(unsigned int i=0;i<nCells;i++)
	{
		if(thisMesh.allCells()[i].getCellCentroid()[0]<=150.0)
			T(i)     = T0;
		else
			T(i)     = T0+1;
		gamma(i) = gamma0;
		u(i,0)   = ux;
		u(i,1)   = uy;
		u(i,2)   = uz;
	}

	// Initialize Boundary Conditions from settings
	const libconfig::Setting& boundaries = tests[testNum]["inits"];

	int count = boundaries.getLength();

	for(int n=0;n<count;n++)
	{
		const libconfig::Setting& boundary = boundaries[n];
		std::string bpname;
		int type;
		patch p;
		double value;

		if(!(
			boundary.lookupValue("patch", bpname)
			&& boundary.lookupValue("type", type)
			&& boundary.lookupValue("value", value)
		    ))
			continue;

		switch (type)
		{
			case 1 :
			{
				p = fixedValue;
				break;
			}
			case 2 : 
			{
				p = fixedGrad;
				break;
			}
			case 3 :
			{
				p = empty;
				break;
			}
			case 4 :
			{
				p = mixed;
				break;
			}
		}

		std::vector<BoundaryPatch>& bpArr = thisMesh.allBPs();

		for(unsigned int n=0;n<bpArr.size();n++)
		{
			if(bpArr[n].getBoundaryPatchName()==bpname)
			{
				bpArr[n].setType(p);
				bpArr[n].setValue(value);
				break;
			}
		}
	}

	const libconfig::Setting& sources = tests[testNum]["sources"];

	count = sources.getLength();

	for(int n=0;n<count;n++)
	{
		const libconfig::Setting& source = sources[n];
		std::array<double,3> pos;
		double value;
		double px, py, pz;

		if(!(
			source.lookupValue("px", px)
			&& source.lookupValue("py", py)
			&& source.lookupValue("pz", pz)
			&& source.lookupValue("value", value)
		    ))
			continue;
		pos = {px, py, pz};

		std::vector<Face> faceArr = thisMesh.allFaces();
		for(unsigned int i=0;i<faceArr.size();i++)
		{
			if(pos==faceArr[i].getFaceCentroid())
				S[faceArr[i].getOwner()] = value;
		}
	}



	} catch(const libconfig::SettingNotFoundException &nfex){
		//Ignore
	}

}

// File reading will be done within the main method here
// Command line reading comes from https://stackoverflow.com/a/868894
int main(int argc, char* argv[]){

	fvMesh thisMesh = fvMesh();

	fvMeshParser parser = fvMeshParser("/home/ojm40/Documents/MPhil_SIMPLE/cavity/constant/polyMesh/");
	parser.readPointsFromFile(thisMesh);
	parser.readFacesFromFile(thisMesh);
	parser.readCellsFromFile(thisMesh);
	parser.readBoundariesFromFile(thisMesh);

	thisMesh.calculateFaceDeltaCoeffs();
	thisMesh.calculateFaceCellDistanceRatios();

	std::cout << thisMesh.displayMeshDetails() << std::endl;
	//std::cout << thisMesh.displayVolumesAndAreas() << std::endl;
	//std::cout << thisMesh.displayCentroids() << std::endl;
	//std::cout << thisMesh.displayBoundaryFaces() << std::endl;
	//std::cout << thisMesh.displayCellNeighbors() << std::endl;
	//std::cout << thisMesh.displayFaceOwnerNeighbor() << std::endl;

	std::ofstream outStream("/home/ojm40/Documents/MPhil_SIMPLE/dat/2D_test.dat");

	const unsigned int nCells = thisMesh.allCells().size();

	// Create storage arrays
	arma::vec T(nCells);
	arma::vec gamma(nCells);
	arma::vec S(nCells);

	// Array of vectors is a n x 3 matrix
	arma::mat x(nCells, 3);
	arma::mat u(nCells, 3);

	int testNum=0;
	bool _test=false;

	if(argc==2)
		testNum = strtol(argv[1], nullptr, 0);
	else if(argc==1)
		std::cout << "No test number parsed. EXIT" << std::endl;
	else if(argc==3)
	{
		std::string arg = argv[2];
		testNum = strtol(argv[3], nullptr, 0);
		if(arg=="-t")
			_test=true;
	}

	int tRes=0;
	
	initializeState(T, gamma, S, u, tRes, thisMesh, testNum);

	fvMatrix thisMatrix = fvMatrix(nCells);

	thisMesh.mergeT(T);
	thisMesh.copyX(x);

	outputState(outStream, x, T);

	double Tb = 0.0;

	double t1 = 1.0;
	double dt = t1 / tRes;

	int numLoops = 1;

	double t = 0.0;
	do{
		t=t+dt;

		fvMatrix::discretizeRateofChange(thisMatrix, T, thisMesh.allCells(), dt);
		fvMatrix::discretizeDiffusion(thisMatrix, gamma, thisMesh.allFaces());
		fvMatrix::discretizeBoundaryConditions(thisMatrix, T, thisMesh.allBPs(), gamma, u, thisMesh.allFaces());
		fvMatrix::discretizeConvectionUpwind(thisMatrix, u, thisMesh.allFaces());

		//if(_test) thisMatrix.printDiscretization();
		//thisMatrix.printDiscretization();

		thisMatrix.solveLinearSystem(T);

		outputState(outStream, x, T);

		fvMatrix::resetMatrix(thisMatrix, thisMesh.allCells());

		std::cout << numLoops++ << std::endl;

	} while (t<t1);

	outputState(outStream, x, T);

        return 0;
}
