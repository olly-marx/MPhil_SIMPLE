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
                Cell(std::vector<int>, int);

		// Public member functions
                std::vector<int>     getCellFaceIndices() const;
		std::array<double,3> getCellCentroid() const;
		double               getCellVolume() const;
		int                  getCellId();

		double               getT() const;
		void                 setT(double);
		double               getP() const;
		void                 setP(double);
		
		void addNeighborFace(int);
		void completeCell(const std::vector<Face>&);

        private:
		// Private member variables
                std::vector<int>     m_faces;
                std::array<double,3> m_cellCentroid;
                double               m_cellVolume;
		int                  m_cellId;
		
		double               T;
		double               p;

		// Private member functions
                std::array<double,3> calcCellCentroid(const std::vector<Face>&);
                double               calcCellVolume(const std::vector<Face>&);
};

#endif
