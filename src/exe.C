//AUTHOR:       Oliver Marx - ojm40@cam.ac.u/k

#include <fstream>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <string>
#include <vector>
#include <sstream>

#include "Face.H"
#include "fvMesh.H"

// File reading will be done within the main method here
// Command line reading comes from https://stackoverflow.com/a/868894
int main(int argc, char* argv[]){
        std::ifstream meshFile("./stg/mesh.txt");

        if(!meshFile){
                std::cerr << "File could not be opened" << std::endl;
                exit(EXIT_FAILURE);
	}

	fvMesh thisMesh = fvMesh();

	// Read in the number of expected points
        int numPoints;
	meshFile >> numPoints;

	std::cout << numPoints << std::endl;
	// Loop over the number of points skipping ( and )
	for(int i=0;i<numPoints+3;i++){
		// Read next line
		std::string line;
		std::getline(meshFile, line);

		// Skip the first and last lines because of the parenthesis
		if(i>0 && i<numPoints+1){
			// Find the position of the parenthesis for this point
			size_t start = line.find_first_of("("),
				end = line.find_first_of(")"); 
			// Read point data into triplet
			std::array<double,3> point;
			std::istringstream iss(line.substr(start+1, end-start-1));
			iss >> point[0] >> point[1] >> point[2];
			thisMesh.addPoint(point);
		}
	}
	
	int numFaces;
	std::string temp;
	std::getline(meshFile, temp);
	numFaces = std::stoi(temp);

	// Loop over the number of faces skipping ( and )
	for(int i=0;i<numFaces+2;i++){
		// Read next line
		std::string line;
		std::getline(meshFile, line);

		std::cout << line << std::endl;

		// Skip the first and last lines because of the parenthesis
		if(i>0 && i<numFaces+1){
			// Find the position of the parenthesis for this point
			size_t start = line.find_first_of("("),
				end = line.find_first_of(")"); 
			// Read point data into triplet
			std::vector<int> vertices;
			int len = std::stoi(line.substr(start-1,start));
			vertices.resize(len);
			std::istringstream iss(line.substr(start+1, end-start-1));
			for(int j=0;j<len;j++) iss >> vertices[j];
			Face f = Face(thisMesh.allPoints(), vertices);
			thisMesh.addFace(f);
		}
	}
	std::cout << thisMesh.getMeshDetails() << std::endl;

        return 0;
}
