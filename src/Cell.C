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
Cell::Cell(std::vector<int> indices, int cellId){
	m_cellId = cellId;

        m_faces = indices;
}

// Return vector of points forming face
std::vector<int> Cell::getCellFaceIndices() const
{
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

void Cell::addNeighborFace(int faceId)
{
	m_faces.push_back(faceId);
}

void Cell::completeCell(const std::vector<Face>& faceArr)
{
	m_cellCentroid = calcCellCentroid(faceArr);
	m_cellVolume = calcCellVolume(faceArr);
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
		result = result + faceArr[i].getFaceCentroid();
        }

        // divide by number of points in face to take average
        result = (1.0/m_faces.size()) * result;

        return result;
}

// Algorithm to tesselate face points and calculate face area
double Cell::calcCellVolume(const std::vector<Face> &faceArr){
        //Init zero result
	double result = 0.0;
	
	for(unsigned int i=0;i<m_faces.size();i++){
		Face f = faceArr[m_faces[i]];
		result += mod ( mod( f.getFaceAreaVector() ) * ( f.getFaceCentroid() - m_cellCentroid ) );
	}
        return result / 3;
}
