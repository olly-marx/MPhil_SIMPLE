//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#ifndef MESH_HEADER
#define MESH_HEADER

#include <iostream>
#include <array>
#include <string>
#include <vector>

#include "Face.H"
#include "Cell.H"

class fvMesh{
        public:
                fvMesh();
		void addPoint(std::array<double,3>);
		void addFace(Face);
		void addCell(Cell);
		const Face& getFace(int);
		const std::array<double,3>& getPoint(int);
		const std::vector<std::array<double,3>>& allPoints();
		const std::vector<Face>& allFaces();
		const std::vector<Cell>& allCells();
		std::string getMeshDetails();
        private:
                std::vector<std::array<double,3>> m_points;
                std::vector<Face> m_faces;
		std::vector<Cell> m_cells;
};

#endif

