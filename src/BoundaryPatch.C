//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#include <iostream>
#include <array>
#include <string>
#include <vector>

#include "Face.H"
#include "BoundaryPatch.H"

// Constructor
BoundaryPatch::BoundaryPatch(const std::vector<Face>& faceArr, std::vector<int> patchFaceIndices, std::string patchName){
	m_boundaryPatchName = patchName;
	
	m_boundaryFaces = patchFaceIndices;

	for(unsigned int i=0;i<patchFaceIndices.size();i++){
		Face b = faceArr[patchFaceIndices[i]];
		int faceOwner = b.getOwner();
		int faceNeighbor = b.getNeighbor();

		if(faceOwner==-1 && faceNeighbor>=0){
			b.setOwner(faceNeighbor);
			// For now have a special case where the cell index for
			// a boundary is -1
			b.setNeighbor(-1);
		}
	}
}

// Public member functions
std::string BoundaryPatch::getBoundaryPatchName(){
	return m_boundaryPatchName;
}

std::vector<int> BoundaryPatch::getBoundaryFaceIndices(){
	return m_boundaryFaces;
}
