//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#ifndef SIM_HEADER
#define SIM_HEADER

#include <iostream>
#include <array>
#include <string>
#include <vector>
#include <armadillo>

#include "Cell.H"
#include "BoundaryPatch.H"
#include "fvMesh.H"
#include "fvMatrix.H"

void outputState(std::ofstream& file, unsigned int nCells, const arma::mat& x,
		const std::array<arma::vec,2>& u, const arma::vec& gradP);

void initializeState(std::array<arma::vec, 2>& u, arma::vec& gamma, double& t1,
		int& tRes, double& Re, double& L, fvMesh& thisMesh, int testNum);

void outputMeshDetails(fvMesh&);

void discretizeMomentumEqn(fvMatrix& thisMatrix,
		std::array<arma::vec,2>& u, arma::vec& P, arma::vec& gamma,
		arma::vec& F, fvMesh& thisMesh, double dt, double alphaU);

void discretizePressureEqn(fvMatrix& pMat, arma::vec& P, fvMatrix& uMat,
		arma::vec& F, std::array<arma::vec,2>& u, fvMesh& thisMesh);

void calculateFaceFluxes(arma::vec& F, fvMesh& thisMesh, std::array<arma::vec,2>& u);

void correctF(arma::vec& F, fvMatrix& uMat, arma::vec& P, fvMesh& thisMesh);

void correctU(fvMatrix& uMat, arma::vec& P,
		std::array<arma::vec,2>& u, fvMesh& thisMesh);

void calculateGradP(std::array<arma::vec,2>& gradP, arma::vec& P, fvMesh& thisMesh);

void divUCell(arma::vec& F, fvMesh& thisMesh);

void explicitUnderRelax(arma::vec& x, const arma::vec& xOld, const double& alphax);

void calculateKinematicViscosity(arma::vec& gamma, const std::array<arma::vec,2> u,
		double L, double Re);

#endif
