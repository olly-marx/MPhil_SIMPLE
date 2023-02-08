//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#ifndef MATRIX_HEADER
#define MATRIX_HEADER

#include <iostream>
#include <array>
#include <string>
#include <vector>

#include "Cell.H"
#include "BoundaryPatch.H"

class fvMatrix{
        public:
		// Constructor
                fvMatrix(int);

		// Public member functions
		static void discretizeRateofChange(fvMatrix&, const std::vector<double>&,
				const std::vector<Cell>&, const double&);

		static void discretizeDiffusion(fvMatrix&, const std::vector<double>&,
				const std::vector<Face>&);

		static void discretizeDirichlet(fvMatrix& m, const std::vector<double>& T,
				const std::vector<BoundaryPatch>& bpArr, const std::vector<Face>& faceArr,
				const double& Tb);

		std::vector<double>& operator[](int);
		std::vector<double>  operator[](int) const;

		double& b(int);
		double  b(int) const;
        private:
		// Private member variables
		std::vector<std::vector<double>> m_a;
		std::vector<double>              m_b;
};

#endif

