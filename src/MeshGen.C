//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#include <iostream>
#include <array>
#include <string>
#include <vector>
#include <cmath>

#include "MeshGen.H"

std::string generateMesh(std::string name){
	double xlen, ylen, zlen;
	int xnum, ynum, znum;
	if(name == "1D")
	{
		xlen = 100;
		ylen = 0.01;
		zlen = 0.01;
		xnum = 100;
		ynum = 1;
		znum = 1;
	}

	double dx = xlen / xnum,
	       dy = ylen / ynum,
	       dz = zlen / znum;

	const int nPoints = (xnum+1) * (ynum+1) * (znum+1);
	const int nFaces = (xnum+1)*(ynum)*(znum) + (xnum)*(ynum+1)*(znum) + (xnum)*(ynum)*(znum+1);
	const int nCells = xnum * ynum * znum;
	
	std::vector<std::array<double,3>> pointArray;
	std::vector<std::vector<int>> faceArray;
	std::vector<std::vector<int>> cellArray;

	for(int k=0;k<znum+1;k++)
	{
		for(int j=0;j<ynum+1;j++)
		{
			for(int i=0;i<xnum+1;i++)
			{
				pointArray.push_back({i*dx, j*dy, k*dz});
			}
		}
	}

}
