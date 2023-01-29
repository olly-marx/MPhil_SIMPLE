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

const Face& fvMesh::getFace(int index){
	return m_faces[index];
}

const std::array<double,3>& fvMesh::getPoint(int index){
	return m_points[index];
}

const std::vector<std::array<double,3>>& fvMesh::allPoints(){
	return m_points;
}

const std::vector<Face>& fvMesh::allFaces(){
	return m_faces;
}

const std::vector<Cell>& fvMesh::allCells(){
	return m_cells;
}

std::string fvMesh::getMeshDetails(){
	std::string out = "fvMesh Details:\n";
	
	out += "# Points: " + std::to_string(m_points.size()) +" \n";
	out += "# Faces: " + std::to_string(m_faces.size()) +" \n";
	out += "# Cells: " + std::to_string(m_cells.size()) +" \n";
	out += "# Boundaries: " ;//+ std::to_string(m_points.size()) +" \n";
	return out;
}
