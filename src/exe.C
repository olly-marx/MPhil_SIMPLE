//AUTHOR:       Oliver Marx - ojm40@cam.ac.u/k

#include <fstream>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <string>
#include <vector>
#include <sstream>

#include "BoundaryPatch.H"
#include "Face.H"
#include "Cell.H"
#include "fvMesh.H"
#include "fvMeshParser.H"

// File reading will be done within the main method here
// Command line reading comes from https://stackoverflow.com/a/868894
int main(int argc, char* argv[]){

	fvMesh thisMesh = fvMesh();

	fvMeshParser parser = fvMeshParser("/home/ojm40/Documents/MPhil_SIMPLE/cavity/constant/polyMesh/");
	parser.readPointsFromFile(thisMesh);
	parser.readFacesFromFile(thisMesh);
	parser.readCellsFromFile(thisMesh);
	parser.readBoundariesFromFile(thisMesh);


	std::cout << thisMesh.displayMeshDetails() << std::endl;
	std::cout << thisMesh.displayVolumesAndAreas() << std::endl;
	std::cout << thisMesh.displayCentroids() << std::endl;
	std::cout << thisMesh.displayBoundaryFaces() << std::endl;
	std::cout << thisMesh.displayCellNeighbors() << std::endl;
	std::cout << thisMesh.displayFaceOwnerNeighbor() << std::endl;

        return 0;
}
