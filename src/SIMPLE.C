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

	const int nCells = thisMesh.allCells().size();

	int testNum=0, totalIterations=0;

	if(argc==3)
	{
		testNum = strtol(argv[1], nullptr, 0);
		totalIterations = strtol(argv[2], nullptr, 0);
	}
	else if(argc<3)
		std::cout << "Incorrect number of arguments parsed. EXIT" << std::endl;

	std::ofstream outStream("./dat/"+std::to_string(totalIterations)+".dat");


	// Create storage arrays
	arma::mat x(nCells, 3);

	arma::mat report(7, 4);

	std::array<arma::vec,2> u;
	u[0] = arma::vec(nCells);
	u[1] = arma::vec(nCells);

	arma::vec P(nCells, 1.0);

	arma::vec gamma(nCells);

	int nFaces = thisMesh.allFaces().size();
	arma::vec F(nFaces);

	int tRes=0;
	double t1 = 0.0, Re = 0.0, L = 0.0, alphaU=0.0, alphaP=0.0;
	
	initializeState(u, gamma, t1, tRes, Re, L, alphaU, alphaP, thisMesh, testNum);

	fvMatrix uMat = fvMatrix(thisMesh);
	fvMatrix pMat = fvMatrix(thisMesh);

	// Copy cell centroid positions into n x 3 matrix (x)
	thisMesh.copyX(x);

	outputState(outStream, nCells, x, u, P, F, thisMesh);

	double dt = t1 / tRes;

	// Init residual storage
	double magR, norm=-1.0;
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

			std::cout << "\n-------------------------------SIMPLE Iteration #"
				  << nIterations << "-------------------------------\n"
				  << std::endl;

			// Calculate kinematic viscosity
			calculateKinematicViscosity(gamma, u, L, Re);

			// Discretize the momentum equation
			fvMatrix::resetMatrix(uMat, nCells);
			discretizeMomentumEqn(uMat, u, P, gamma, F, thisMesh, dt, alphaU);

			// Solve the momentum equation
			uMat.solveLinearSystem(u[0], true);
			uMat.solveLinearSystem(u[1], false);

			// Calculate uncorrected face fluxes from velocity field
			calculateFaceFluxes(F, thisMesh, u);
			std::cout << "\nMomentum predictor step completed..." << std::endl;
			divU(F, thisMesh);

			// Store old pressure and print it
			const arma::vec POld = P;

			// Discretize Pressure equation
			fvMatrix::resetMatrix(pMat, nCells);
			discretizePressureEqn(pMat, P, uMat, F, u, thisMesh);

			if(nIterations==1)
				buildReport(report, uMat, pMat, thisMesh);

			// Solve Pressure Equation
			pMat.solveLinearSystem(P, true);

			// Correct face fluzes with new pressure field
			correctF(F, uMat, P, thisMesh);
			std::cout << "\nPressure corrector step completed..." << std::endl;
			divU(F, thisMesh);

			// Explicitly underrelax pressure
			explicitUnderRelax(P, POld, alphaP);

			// Correct velocities with new pressure field
			correctU(uMat, P, u, thisMesh);

			//Calculate and print pressure matrix residual
			fvMatrix::calculateUResidual(resP, pMat, P, thisMesh, norm);

			//outputState(outStream, nCells, x, u, P, F, thisMesh);

		}
		while (nIterations<totalIterations);
		// End of SIMPLE LOOP

		std::cout 
			  //<< (numLoops++) 
			  //<< ": t=" << t 
			  << "\n SIMPLE Iters: " << nIterations
			  << std::endl;
	}
	while (t<t1);

	printReport(report);

	outputState(outStream, nCells, x, u, P, F, thisMesh);

        return 0;
}

