//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#ifndef MESH_HEADER
#define MESH_HEADER

#include <iostream>
#include <array>
#include <string>
#include <vector>

#include <armadillo>

#include "BoundaryPatch.H"
#include "Face.H"
#include "Cell.H"

class fvMesh{
        public:
		//Constructors
                fvMesh();

		//Mesh construction
		void addPoint(std::array<double,3>);
		void addFace(Face);
		void addCell(Cell);
		void addBoundaryPatch(BoundaryPatch);

		//Mesh fetching
		const Face& getFace(int);
		const std::array<double,3>&        getPoint(int);
		std::vector<std::array<double,3>>& allPoints();
		std::vector<Face>&                 allFaces();
		std::vector<Cell>&                 allCells();
		std::vector<BoundaryPatch>&        allBPs();

		//Mesh utilities
		void                               copyT(arma::vec&) const;
		void                               mergeT(arma::vec&);
		void                               copyX(arma::mat&) const;

		//Random mesh detail methods
		std::array<int, 4>                 getMeshDetails();
		double                             cellVolume(int);
		std::array<double,3>               faceAreaVector(int);
		std::array<double,3>               cellCentroid(int);
		std::array<double,3>               faceCentroid(int);
		std::vector<std::array<int,2>>     boundaryFaceOwner(BoundaryPatch);
		std::vector<int>                   cellNeighbors(int);
		std::array<int,2>                  faceOwnerNeighbor(int);

		std::string displayMeshDetails();
		std::string displayVolumesAndAreas();
		std::string displayCentroids();
		std::string displayBoundaryFaces();
		std::string displayCellNeighbors();
		std::string displayFaceOwnerNeighbor();

		//Mesh finalizing methods
		void                               calculateFaceCellDistanceRatios();
		void                               calculateFaceDeltaCoeffs();

        private:
                std::vector<std::array<double,3>> m_points;
                std::vector<Face> m_faces;
		std::vector<Cell> m_cells;
		std::vector<BoundaryPatch> m_boundaryPatches;
};

#endif

