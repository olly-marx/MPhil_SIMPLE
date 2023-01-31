//AUTHOR:       Oliver Marx - ojm40@cam.ac.u/k

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>

#include "fvMesh.H"
#include "Face.H"
#include "Cell.H"
#include "VectorUtils.H"

fvMesh::fvMesh(){
	std::vector<std::array<double,3>> m_points;	
	std::vector<Face> m_faces;	
	std::vector<Cell> m_cells;
}

void fvMesh::addPoint(std::array<double,3> p){
	m_points.push_back(p);
}

void fvMesh::addFace(Face f){
	m_faces.push_back(f);
}

void fvMesh::addCell(Cell c){
	m_cells.push_back(c);
}

void fvMesh::addBoundaryPatch(BoundaryPatch bp){
	m_boundaryPatches.push_back(bp);
}

const Face& fvMesh::getFace(int index){
	return m_faces[index];
}

const std::array<double,3>& fvMesh::getPoint(int index){
	return m_points[index];
}

std::vector<std::array<double,3>>& fvMesh::allPoints(){
	return m_points;
}

std::vector<Face>& fvMesh::allFaces(){
	return m_faces;
}

std::vector<Cell>& fvMesh::allCells(){
	return m_cells;
}

std::string fvMesh::getMeshDetails(){
	std::string out = "fvMesh Details:\n";
	
	out += "# Points: " + std::to_string(m_points.size()) +" \n";
	out += "# Faces: " + std::to_string(m_faces.size()) +" \n";
	out += "# Cells: " + std::to_string(m_cells.size()) +" \n";
	out += "# Boundaries: " + std::to_string(m_boundaryPatches.size()) +" \n";
	return out;
}

std::string fvMesh::displayCentroids(){
	std::string out = "Cell Volumes: \n";
	for(Cell c : m_cells) out += std::to_string(c.getCellId()) + " " 
		+ std::to_string(c.getCellCentroid()[0]) + " " 
		+ std::to_string(c.getCellCentroid()[1]) + " " 
		+ std::to_string(c.getCellCentroid()[2]) + "\n";
	out += "\nFace Centroids:\n";
	for(Face f : m_faces) out += std::to_string(f.getFaceId()) + " " 
		+ printV(f.getFaceAreaVector()) + std::to_string(f.getOwner()) + " " 
		+ std::to_string(f.getNeighbor()) + "\n";
	return out;
}

std::string fvMesh::displayVolumesAndAreas(){
	std::string out = "Cell Volumes: \n";
	for(Cell c : m_cells) out += std::to_string(c.getCellId()) + " " 
		+ std::to_string(c.getCellVolume()) + "\n";
	out += "\nFace Areas:\n";
	for(Face f : m_faces) out += std::to_string(f.getFaceId()) + " " 
		+ printV(f.getFaceAreaVector()) + " -- Mag:"
		+ std::to_string( mod(f.getFaceAreaVector() ) ) + "\n";
	return out;
}

std::string fvMesh::displayBoundaryFaces(){
	std::string out = "Boundary Faces (by patch)\n";
	for(BoundaryPatch bp : m_boundaryPatches){
		out+=bp.getBoundaryPatchName() + "\n";
		std::vector<int> patchFaces = bp.getBoundaryFaceIndices();
		for(unsigned int i=0;i<bp.getBoundaryFaceIndices().size();i++){
			Face f = m_faces[patchFaces[i]];
			out += std::to_string(f.getFaceId()) + " : Owner " + std::to_string(f.getOwner()) + "\n";
		}
	}
	return out;
}

std::string fvMesh::displayCellNeighbors(){
	std::string out = "Cell Neighbors:\n";
	for(Cell c : m_cells){
		std::vector<int> cellFaces = c.getCellFaceIndices();
		out += std::to_string(c.getCellId()) + " : ";
		for(unsigned int i=0;i<cellFaces.size();i++){
			int owner = m_faces[cellFaces[i]].getOwner();
			int neighbor = m_faces[cellFaces[i]].getNeighbor();
			if(c.getCellId()==owner && neighbor>=0) out += std::to_string(neighbor) + " ";
			else if(c.getCellId()==neighbor) out += std::to_string(owner) + " ";
		}
		out += "\n";
	}
	return out;
}

std::string fvMesh::displayFaceOwnerNeighbor(){
	std::string out = "Face Owner-Neighbor Pairs:\n";
	for(Face f : m_faces) out += std::to_string(f.getFaceId()) + " " 
		+ "O -> " + std::to_string(f.getOwner()) + " " 
		+ "N -> " + (f.getNeighbor()==-1 ? "b" : std::to_string(f.getNeighbor())) + "\n";
	return out;
}
