//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#ifndef BOUNDARY_HEADER
#define BOUNDARY_HEADER

#include <iostream>
#include <array>
#include <string>
#include <vector>

#include "Face.H"

class BoundaryPatch{
        public:
		// Constructor
                BoundaryPatch(int, int, std::string);

		// Public member functions
		std::string          getBoundaryPatchName() const;
		int                  getBoundaryPatchStartFace() const;
		int                  getBoundaryPatchLength() const;
        private:
		// Private member variables
		std::string          m_boundaryPatchName;
                int                  m_startFace;
                int                  m_length;

		// Private members funtions
};

#endif

