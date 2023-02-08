//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#ifndef MESHPARSER_HEADER
#define MESHPARSER_HEADER

#include <iostream>
#include <array>
#include <string>
#include <vector>

#include "fvMesh.H"

class fvMeshParser{
        public:
		// Constructor
                fvMeshParser(std::string);

		void readPointsFromFile(fvMesh&);
		void readFacesFromFile(fvMesh&);
		void readCellsFromFile(fvMesh&);
		void readBoundariesFromFile(fvMesh&);


        private:
		// Private member variables
		std::string m_fileDirectory;
};

#endif
