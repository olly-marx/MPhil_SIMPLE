//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk
#include <iostream>
#include <array>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>

#include "Face.H"
#include "VectorUtils.H"
#include "Cell.H"


// Face constructor, takes a vector of points to define a face of n Points
Cell::Cell(std::vector<Face> &faceArr, std::vector<int> indices, int cellId){
	m_cellId = cellId;

        m_faces = indices;

        //Implement cell centroid calculation
        m_cellCentroid = calcCellCentroid(faceArr);
	//std::cout << m_cellCentroid[0] << " " << m_cellCentroid[1] << " " << m_cellCentroid[2] << std::endl;

        //Implement tesselation calc
        m_cellVolume = calcCellVolume(faceArr);
	//std::cout << "V: " << std::to_string(m_cellVolume) << std::endl;

	//Set Owner and Neighbor cells based on the sign of the fP and faceVec
	for(int i : m_faces){
		std::array<double,3> fArea = faceArr[i].getFaceAreaVector(),
			             fP = diff( faceArr[i].getFaceCentroid() , m_cellCentroid );
		double direction = dot( fArea , fP );
		if(direction > 0) faceArr[i].setOwner(m_cellId);
		else if (direction < 0) faceArr[i].setNeighbor(m_cellId);
	}
}

// Return vector of points forming face
std::vector<int> Cell::getCellFaceIndices(){
        return m_faces;
}

// Return const cell centroid vector
std::array<double,3> Cell::getCellCentroid() const{
	return m_cellCentroid;
}

// Return const cell volume double
double Cell::getCellVolume() const{
	return m_cellVolume;
}

int Cell::getCellId(){
	return m_cellId;
}

// Algorithm to take points defining face and calculate position of face
// centroid in 3D
// Defined as the vector average of the points defining a faces
std::array<double,3> Cell::calcCellCentroid(const std::vector<Face> &faceArr){

        // Init empty result array
        std::array<double,3> result = {0.0, 0.0, 0.0};
        // Loop over points of face
        for(int i : m_faces){

                // Add current point to result
		result = sum(result, faceArr[i].getFaceCentroid());
        }

        // divide by number of points in face to take average
        result = scalarMult(1.0/m_faces.size(), result);

        return result;
}

// Algorithm to tesselate face points and calculate face area
double Cell::calcCellVolume(const std::vector<Face> &faceArr){
        //Init zero result
	double result = 0.0;
	
	for(unsigned int i=0;i<m_faces.size();i++){
		Face f = faceArr[m_faces[i]];
		result += mod ( scalarMult ( mod( f.getFaceAreaVector() ) , diff( f.getFaceCentroid() , m_cellCentroid ) ) );
	}
        return result / 3;
}
