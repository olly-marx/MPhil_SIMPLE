//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#include <cstddef>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <string>
#include <vector>
#include <sstream>

#include "VectorUtils.H"
#include "fvMesh.H"
#include "fvMeshParser.H"

fvMeshParser::fvMeshParser(std::string dir)
{
	m_fileDirectory = dir;
}

void fvMeshParser::readPointsFromFile(fvMesh& m)
{
	std::cout << "Reading points..." << std::endl;

	std::string file {"points"};
        std::ifstream meshFile(m_fileDirectory+file);

        if(!meshFile){
                std::cerr << file << " : file could not be opened" << std::endl;
                exit(EXIT_FAILURE);
	}

	std::string trash;
	while(std::getline(meshFile, trash) && !trash.empty());
	std::getline(meshFile, trash);

	// Read in the number of expected points
	int numPoints;
	std::string temp;
	std::getline(meshFile, temp);
	numPoints = std::stoi(temp);

	// Loop over the number of points skipping ( and )
	for(int i=0;i<numPoints+2;i++){
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
			m.addPoint(point);
		}
	}

	std::cout << "Done!" << std::endl;

	meshFile.close();
}

void fvMeshParser::readFacesFromFile(fvMesh& m)
{
	std::cout << "Reading faces..." << std::endl;

	std::string file {"faces"};
        std::ifstream meshFile(m_fileDirectory+file);

        if(!meshFile){
                std::cerr << file << " : file could not be opened" << std::endl;
                exit(EXIT_FAILURE);
	}

	std::string trash;
	while(std::getline(meshFile, trash) && !trash.empty());
	std::getline(meshFile, trash);
	
	int numFaces;
	std::string temp;
	std::getline(meshFile, temp);
	numFaces = std::stoi(temp);

	// Loop over the number of faces skipping ( and )
	for(int i=0;i<numFaces+2;i++){
		// Read next line
		std::string line;
		std::getline(meshFile, line);

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
			Face f = Face(m.allPoints(), vertices, i-1);
			m.addFace(f);
		}
	}

	std::cout << "Done!" << std::endl;
	meshFile.close();
}

void fvMeshParser::readCellsFromFile(fvMesh& m)
{
	std::cout << "Reading owners..." << std::endl;

	std::string file {"owner"};
        std::ifstream meshFile(m_fileDirectory+file);

        if(!meshFile){
                std::cerr << file << " file could not be opened" << std::endl;
                exit(EXIT_FAILURE);
	}

	std::string trash;
	while(std::getline(meshFile, trash) && !trash.empty());
	std::getline(meshFile, trash);

	// Get the number of faces in owner file
	int numFaces;
	std::string temp;
	std::getline(meshFile, temp);
	numFaces = std::stoi(temp);

	// The highest cell owner index is also the number of cell in domain -1
	std::vector<int> owners;
	int numCells = 0;
	// For each cell store the faces attached to the cell
	std::vector<std::vector<int>> cellsfaces;

	// Loop over the number of cells skipping ( and )
	for(int i=0;i<numFaces+2;i++){
		// Read next line
		std::string line;
		std::getline(meshFile, line);

		// Skip the first and last lines because of the parenthesis
		if(i>0 && i<numFaces+1)
		{
			// Read next face owner index
			int cellId = std::stoi(line);

			// push owner owners array
			owners.push_back(cellId);

			// check if the index is largest and set num Cells
			// accordingly
			if(cellId > numCells)
			{
				numCells = cellId;
			}

			// Set current face owner
			m.allFaces()[i-1].setOwner(cellId);
		}
	}

	std::cout << "Done!" << std::endl;

	meshFile.close();

	cellsfaces.resize(numCells+1);

	// push faces owned by cells onto cell face array
	// here j is the face number
	for(std::size_t j=0;j<owners.size();j++)
	{
		// depending on cell index stored in owners push j onto the
		// cell's face array
		cellsfaces[owners[j]].push_back(j);
	}

	for(std::size_t k=0;k<cellsfaces.size();k++)
	{
		Cell c = Cell(cellsfaces[k], k);
		m.addCell(c);
	}

	std::cout << "Reading neighbours..." << std::endl;

	file = "neighbour";
        meshFile.open(m_fileDirectory+file);

        if(!meshFile){
                std::cerr << file << " file could not be opened" << std::endl;
                exit(EXIT_FAILURE);
	}

	while(std::getline(meshFile, trash) && !trash.empty());
	std::getline(meshFile, trash);

	std::getline(meshFile, temp);
	numFaces = std::stoi(temp);

	// Loop over the number of internal faces skipping ( and )
	for(int i=0;i<numFaces+2;i++){
		// Read next line
		std::string line;
		std::getline(meshFile, line);

		// Skip the first and last lines because of the parenthesis
		if(i>0 && i<numFaces+1)
		{
			int cellId = std::stoi(line);

			m.allFaces()[i-1].setNeighbor(cellId);
			m.allFaces()[i-1].makeInternalFace();

			m.allCells()[cellId].addNeighborFace(i-1);
		}
	}

	for(std::size_t i=0;i<m.allCells().size();i++)
	{
		m.allCells()[i].completeCell(m.allFaces());
	}

	std::cout << "Done" << std::endl;
	
	meshFile.close();
}

void fvMeshParser::readBoundariesFromFile(fvMesh& m)
{
	std::cout << "Reading boundaries..." << std::endl;

	std::string file {"boundary"};
        std::ifstream meshFile(m_fileDirectory+file);

        if(!meshFile){
                std::cerr << file << " file could not be opened" << std::endl;
                exit(EXIT_FAILURE);
	}

	std::string trash;
	while(std::getline(meshFile, trash) && !trash.empty());

	int numBoundaryPatches;
	std::string temp;
	std::getline(meshFile, temp);
	numBoundaryPatches = std::stoi(temp);

	//skip the initial backet
	std::getline(meshFile, temp);
	// Loop over the number of patches skipping ( and )
	for(int i=0;i<numBoundaryPatches;i++){
		std::vector<int> faces;
		int len;
		std::string patchName;
		for(int j=0;j<7;j++){
			// Read next line
			std::string line;
			std::getline(meshFile, line);
			std::istringstream iss(line);

			switch (j) {
				case 0: {
					iss >> patchName;
					break;
				}
				case 4: {
					std::string item;
					iss >> item;
					if(item=="nFaces")
					{
						iss >> len;
						faces.resize(len);
					}
					break;
				}
				case 5: {
					int start;
					std::string item;
					iss >> item;
					if(item=="startFace")
					{
						iss >> start;
					}
					for(int k=0;k<len;k++) iss >> faces[k];
					BoundaryPatch bp = BoundaryPatch(start, len, patchName);
					m.addBoundaryPatch(bp);
					break;
				}
			}	
		}
	}

	std::cout << "Done!" << std::endl;

	meshFile.close();
}

