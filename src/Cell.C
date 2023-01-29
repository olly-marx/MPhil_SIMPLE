//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk
#include <iostream>
#include <array>
#include <vector>
#include <string>
#include <cstring>

#include "Face.H"
#include "VectorUtils.H"
#include "Cell.H"

// Face constructor, takes a vector of points to define a face of n Points
Cell::Cell(const std::vector<Face> &faceArr, std::vector<int> indices){
        m_faces = indices;

        //Implement cell centroid calculation
        m_cellCentroid = calcCellCentroid(faceArr);
	std::cout << m_cellCentroid[0] << " " << m_cellCentroid[1] << " " << m_cellCentroid[2] << std::endl;

        //Implement tesselation calc
        m_cellVolume = calcCellVolume(faceArr);
	std::cout << "V: " << std::to_string(m_cellVolume) << std::endl;
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

// Algorithm to take points defining face and calculate position of face
// centroid in 3D
// Defined as the vector average of the points defining a faces
std::array<double,3> Cell::calcCellCentroid(const std::vector<Face> &faceArr){

        // Init empty result array
        std::array<double,3> result = {0.0, 0.0, 0.0};
        // Loop over points of face
        for(int i=0;i<m_faces.size();i++){

                // Add current point to result
		result = sum(result, faceArr[m_faces[i]].getFaceCentroid());
        }

        // divide by number of points in face to take average
        result = scalarMult(1.0/m_faces.size(), result);

        return result;
}

// Algorithm to tesselate face points and calculate face area
double Cell::calcCellVolume(const std::vector<Face> &faceArr){
        //Init zero result
        double result = 0.0;
	
	for(int i=0;i<m_faces.size();i++){
		Face f = faceArr[m_faces[i]];
		result += mod( diff( f.getFaceCentroid() , m_cellCentroid ) ) 
				* mod( f.getFaceAreaVector() );
	}
        return result / 3;
}
