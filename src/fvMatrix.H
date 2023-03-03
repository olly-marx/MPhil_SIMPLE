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

class fvMatrix{
        public:
		// Constructor
                fvMatrix(int);

		// Public member functions
		static void discretizeRateofChange(fvMatrix&, const arma::vec& T,
				const std::vector<Cell>&, const double&);

		static void discretizeDiffusion(fvMatrix&, const arma::vec& gamma,
				const std::vector<Face>&);

		static void discretizeBoundaryConditions(fvMatrix& m, const arma::vec& T,
				const std::vector<BoundaryPatch>& bpArr, const arma::vec& gamma,
				const arma::mat& u, const std::vector<Face>& faceArr);

		static void discretizeConvectionUpwind(fvMatrix&, const arma::mat& u,
				const std::vector<Face>&);

		static void discretizeConvectionCentral(fvMatrix&, const arma::mat& u,
				const std::vector<Face>&);

		static void resetMatrix(fvMatrix&, const std::vector<Cell>&);

		// Member function solver
		void solveLinearSystem(arma::Col<double>&);

		void printDiscretization();

		double& A(int, int);
		double A(int, int) const;

		double& b(int);
		double  b(int) const;
        private:
		// Private member variables
		arma::Mat<double> m_A;
		arma::Col<double> m_b;
};

#endif

