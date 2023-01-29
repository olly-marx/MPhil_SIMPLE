//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#ifndef CELL_HEADER
#define CELL_HEADER

#include <iostream>
#include <array>
#include <string>
#include <vector>

#include "Face.H"

class Cell{
        public:
		// Constructor
                Cell(const std::vector<Face>&, std::vector<int>);

		// Public member functions
                std::vector<int>     getCellFaceIndices();
		std::array<double,3> getCellCentroid() const;
		double               getCellVolume() const;

        private:
		// Private member variables
                std::vector<int>     m_faces;
                std::array<double,3> m_cellCentroid;
                double               m_cellVolume;

		// Private member functions
                std::array<double,3> calcCellCentroid(const std::vector<Face>&);
                double               calcCellVolume(const std::vector<Face>&);
};

#endif
