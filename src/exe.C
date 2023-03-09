//AUTHOR:       Oliver Marx - ojm40@cam.ac.u/k

#include <cstddef>
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
#include "fvSimulation.H"
#include "VectorUtils.H"


// File reading will be done within the main method here
// Command line reading comes from https://stackoverflow.com/a/868894
int main(int argc, char* argv[]){

	fvMesh thisMesh = fvMesh();

	fvMeshParser parser = fvMeshParser("./cavity/constant/polyMesh/");

	parser.readPointsFromFile(thisMesh);
	parser.readFacesFromFile(thisMesh);
	parser.readCellsFromFile(thisMesh);
	parser.readBoundariesFromFile(thisMesh);

	thisMesh.calculateFaceDeltaCoeffs();
	thisMesh.calculateFaceCellDistanceRatios();

	//outputMeshDetails();

	std::ofstream outStream("./dat/2D_test.dat");

	const int nCells = thisMesh.allCells().size();

	int testNum=0;

	if(argc==2)
		testNum = strtol(argv[1], nullptr, 0);
	else if(argc==1)
		std::cout << "No test number parsed. EXIT" << std::endl;


	// Create storage arrays
	arma::mat x(nCells, 3);

	std::array<arma::vec,2> u;
	u[0] = arma::vec(nCells);
	u[1] = arma::vec(nCells);

	arma::vec P(nCells, 1.0);

	arma::vec gamma(nCells);

	int nFaces = thisMesh.allFaces().size();
	arma::vec F(nFaces);

	int tRes=0;
	double t1 = 0.0, Re = 0.0, L = 0.0;
	
	initializeState(u, gamma, t1, tRes, Re, L, thisMesh, testNum);

	fvMatrix uMat = fvMatrix(thisMesh);
	fvMatrix pMat = fvMatrix(thisMesh);

	// Copy cell centroid positions into n x 3 matrix (x)
	thisMesh.copyX(x);

	outputState(outStream, nCells, x, u, P);

	double dt = t1 / tRes;

	// Init residual storage
	double magR, norm=-1.0, alphaP=0.2, alphaU=0.8;
	arma::vec resP(nCells);

	int numLoops = 1;

	double t = 0.0;
	do{
		t=t+dt;

		int nIterations = 0;

		calculateFaceFluxes(F, thisMesh, u);

		// Start of SIMPLE LOOP
		do{
			nIterations++;

			// Calculate kinematic viscosity
			calculateKinematicViscosity(gamma, u, L, Re);

			// Discretize the momentum equation
			discretizeMomentumEqn(uMat, u, gamma, F, thisMesh, dt, alphaU);
			//std::cout << "Momentum Discretization:" <<std::endl;
			//uMat.printDiscretization();

			// Solve the momentum equation
			uMat.solveLinearSystem(u[0], true);
			uMat.solveLinearSystem(u[1], false);
			//std::cout << "X Velocity:" <<std::endl;
			//u[0].print();
			//std::cout << "Y Velocity:" <<std::endl;
			//u[1].print();

			// Calculate uncorrected face fluxes from velocity field
			calculateFaceFluxes(F, thisMesh, u);

			// Discretize Pressure equation
			const arma::vec POld = P;
			//std::cout << "POld:" << std::endl;
			//POld.print();

			discretizePressureEqn(pMat, P, uMat, F, u, thisMesh);
			std::cout << "Pressure Discretization:" <<std::endl;
			pMat.printDiscretization();

			pMat.solveLinearSystem(P, true);

			//std::cout << "Pnew:" << std::endl;
			//P.print();

			//explicitUnderRelax(P, POld, alphaP);

			std::cout << "FOld:" << std::endl;
			F.print();
			correctF(F, uMat, P, thisMesh);
			//std::cout << "FNew:" << std::endl;
			//F.print();
			//std::cout <<"||F|| " << arma::norm(F) << std::endl;

			//std::cout << "PUR:" << std::endl;
			//P.print();

			correctU(uMat, P, u, thisMesh);
			//std::cout << "X Velocity:" <<std::endl;
			//u[0].print();
			//std::cout << "Y Velocity:" <<std::endl;
			//u[1].print();

			fvMatrix::calculateUResidual(resP, pMat, P, thisMesh, norm);
			//std::cout << "P Residual:" <<std::endl;
			//resP.print();

			magR = arma::norm(resP);
			std::cout << "||r|| " << magR << std::endl;
		}
		while (magR > 1e-8 && nIterations<10);
		// End of SIMPLE LOOP

		outputState(outStream, nCells, x, u, P);

		fvMatrix::resetMatrix(uMat, nCells);

		std::cout << (numLoops++) 
			  //<< ": t=" << t 
			  << "\n SIMPLE Iters: " << nIterations
			  << std::endl;
	}
	while (t<t1);

	outputState(outStream, nCells, x, u, P);

        return 0;
}

