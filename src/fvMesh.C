//AUTHOR:       Oliver Marx - ojm40@cam.ac.u/k

#include <cstddef>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>

#include "BoundaryPatch.H"
#include "fvMesh.H"
#include "Face.H"
#include "Cell.H"
#include "VectorUtils.H"

/* -------------------------------------------------------------------*/         
/*                                                                    */         
/*                    Mesh Initialization Methods                     */         
/*                                                                    */         
/* -------------------------------------------------------------------*/

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

std::vector<BoundaryPatch>& fvMesh::allBPs(){
	return m_boundaryPatches;
}

void fvMesh::copyX(arma::mat& X) const
{
	for(unsigned int i=0;i<m_cells.size();i++)
	{
		for(int n=0;n<3;n++)
		{
			X(i,n) = m_cells[i].getCellCentroid()[n];
		}
	}
}

int fvMesh::getNInternalFaces() const
{
	return m_internalFaces;
}

void fvMesh::setNInternalFaces(int n)
{
	m_internalFaces = n;
}

int fvMesh::getNCells() const
{
	return m_nCells;
}

void fvMesh::setNCells(int n)
{
	m_nCells = n;
}

std::array<int, 4> fvMesh::getMeshDetails(){
	std::array<int, 4> details;
	details[0] = m_points.size();
	details[1] = m_faces.size();
	details[2] = m_cells.size();
	details[3] = m_boundaryPatches.size();
	return details;
}

double fvMesh::cellVolume(int c) const
{
	return m_cells[c].getCellVolume();
}

std::array<double,3> fvMesh::faceAreaVector(int f) const
{
	return m_faces[f].getFaceAreaVector();
}

std::array<double,3> fvMesh::cellCentroid(int c) const
{
	return m_cells[c].getCellCentroid();
}

std::array<double,3> fvMesh::faceCentroid(int f) const
{
	return m_faces[f].getFaceCentroid();
}

std::vector<std::array<int,2>> fvMesh::boundaryFaceOwner(BoundaryPatch bp){
	std::vector<std::array<int,2>> faceCellPairs;
	for(int f=0;f<bp.getBoundaryPatchLength();f++)
	{
		faceCellPairs.push_back({m_faces[bp.getBoundaryPatchStartFace()+f].getFaceId()
				, m_faces[bp.getBoundaryPatchStartFace()+f].getOwner()});
	}
	return faceCellPairs;
}

std::vector<int> fvMesh::cellNeighbors(int c) const
{
	std::vector<int> allCellNeighbors;
	std::vector<int> cellFaces = m_cells[c].getCellFaceIndices();
	for(unsigned int i=0;i<cellFaces.size();i++){
		int owner = m_faces[cellFaces[i]].getOwner();
		int neighbor = m_faces[cellFaces[i]].getNeighbor();
		if(c==owner && neighbor>=0) allCellNeighbors.push_back(neighbor);
		else if(c==neighbor) allCellNeighbors.push_back(owner);
	}
	return allCellNeighbors;
}

std::array<int,2> fvMesh::faceOwnerNeighbor(int f) const
{
	return {m_faces[f].getOwner(), m_faces[f].getNeighbor()};
}

void fvMesh::calculateFaceCellDistanceRatios()
{
	for(std::size_t i=0;i<m_faces.size();i++)
	{
		Face& f = m_faces[i];
		double fx;

		const Cell& c = m_cells[f.getOwner()];

		if(!f.isBoundary())
		{
			const Cell& n = m_cells[f.getNeighbor()];
			
			fx = mod( m_faces[i].getFaceCentroid() - n.getCellCentroid() ) /
				mod( n.getCellCentroid() - c.getCellCentroid() );

			f.setfx(fx);
		}
	}
} 

void fvMesh::calculateFaceDeltaCoeffs()
{
	for(std::size_t i=0;i<m_faces.size();i++)
	{
		Face& f = m_faces[i];
		double delta;

		const Cell& c = m_cells[f.getOwner()];

		if(f.isBoundary())
		{
			delta = 1.0 / mod( f.getFaceCentroid() - c.getCellCentroid() );

			f.setDelta(delta);
		}
		else
		{
			const Cell& n = m_cells[f.getNeighbor()];
			
			delta = 1.0 / mod( n.getCellCentroid() - c.getCellCentroid() );

			f.setDelta(delta);
		}
	}
} 

std::string fvMesh::displayMeshDetails(){
	std::string out = "fvMesh Details:\n";
	
	out += "# Points: " + std::to_string(m_points.size()) +" \n";
	out += "# Faces: " + std::to_string(m_faces.size()) +" \n";
	out += "# Cells: " + std::to_string(m_cells.size()) +" \n";
	out += "# Boundaries: " + std::to_string(m_boundaryPatches.size()) +" \n";
	return out;
}

std::string fvMesh::displayCentroids(){
	std::string out = "Cell Centroids: \n";
	for(Cell c : m_cells) out += std::to_string(c.getCellId()) + " " 
		+ std::to_string(c.getCellCentroid()[0]) + " " 
		+ std::to_string(c.getCellCentroid()[1]) + " " 
		+ std::to_string(c.getCellCentroid()[2]) + "\n";
	out += "\nFace Centroids:\n";
	for(Face f : m_faces) out += std::to_string(f.getFaceId()) + " " 
		+ printV(f.getFaceCentroid()) + std::to_string(f.getOwner()) + " " 
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
		for(int i=0;i<bp.getBoundaryPatchLength();i++){
			Face f = m_faces[bp.getBoundaryPatchStartFace()+i];
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
