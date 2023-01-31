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
                BoundaryPatch(const std::vector<Face>&, std::vector<int>, std::string);

		// Public member functions
		std::string          getBoundaryPatchName();
                std::vector<int>     getBoundaryFaceIndices();
        private:
		// Private member variables
		std::string          m_boundaryPatchName;
                std::vector<int>     m_boundaryFaces;

		// Private members funtions
};

#endif

