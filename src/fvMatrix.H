//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#ifndef MATRIX_HEADER
#define MATRIX_HEADER

#include <iostream>
#include <array>
#include <string>
#include <vector>
#include <armadillo>

#include "Cell.H"
#include "BoundaryPatch.H"
#include "fvMesh.H"

class fvMatrix{
        public:
		// Constructor
                fvMatrix(fvMesh&);

		// Public member functions
		static void discretizeRateofChange(fvMatrix&, const std::array<arma::vec,2>& T,
				fvMesh&, const double&);

		static void discretizeLaplacian(fvMatrix&, const arma::vec& gamma,
				fvMesh&, bool plus);

		static void discretizeConvectionUpwind(fvMatrix& m, fvMesh& thisMesh,
				const arma::vec& F, std::array<arma::vec,2>& T);

		static void discretizeGradP(fvMatrix& m, arma::vec& P, fvMesh& thisMesh);

		static void implicitUnderRelax(fvMatrix& m, std::array<arma::vec,2>& u, fvMesh& thisMesh,
				double alphaU);

		static void discretizeContinuity(fvMatrix& m, fvMesh& thisMesh, 
				const arma::vec& F);

		static void calculateUResidual(arma::vec& res, fvMatrix& m,
				arma::vec& x, fvMesh& thisMesh, double& norm);

		static void resetMatrix(fvMatrix&, int nCells);

		// Member function solver
		void solveLinearSystem(arma::Col<double>&, bool);
		arma::vec recipD() const;

		void printDiscretization();

		double& A(int, int);
		double A(int, int) const;

		double& bx(int);
		double  bx(int) const;

		double& by(int);
		double  by(int) const;
        private:
		// Private member variables
		int nCells;
		arma::mat m_A;
		arma::vec m_bx;
		arma::vec m_by;

		//arma::sp_mat m_D;
		//arma::sp_mat m_U;
		//arma::sp_mat m_L;
};

#endif

