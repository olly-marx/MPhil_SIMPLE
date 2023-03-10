//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#include <iostream>
#include <array>
#include <string>
#include <vector>
#include <armadillo>
#include <libconfig.h++>

#include "Cell.H"
#include "BoundaryPatch.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "VectorUtils.H"
#include "fvSimulation.H"

void discretizeMomentumEqn(fvMatrix& uMat,
		std::array<arma::vec,2>& u, arma::vec& P, arma::vec& gamma,
		arma::vec& F, fvMesh& thisMesh, double dt, double alphaU)
{
	//fvMatrix::discretizeRateofChange(uMat, u, thisMesh, dt);
	fvMatrix::discretizeLaplacian(uMat, gamma, thisMesh, false);
	fvMatrix::discretizeConvectionUpwind(uMat, thisMesh, F, u);
	//fvMatrix::discretizeGradP(uMat, P, thisMesh);
	//fvMatrix::implicitUnderRelax(uMat, u, thisMesh, alphaU);
}

void discretizePressureEqn(fvMatrix& pMat,
		arma::vec &P, fvMatrix &uMat, arma::vec &F,
		std::array<arma::vec, 2> &u, fvMesh& thisMesh)
{
	const arma::vec rD = uMat.recipD();
	fvMatrix::discretizeLaplacian(pMat, rD, thisMesh, true);
	fvMatrix::discretizeContinuity(pMat, thisMesh, F);
}

void calculateFaceFluxes(arma::vec& F, fvMesh& thisMesh, std::array<arma::vec,2>& u)
{
	std::vector<Face>& faceArr = thisMesh.allFaces();
	std::vector<BoundaryPatch>& bpArr = thisMesh.allBPs();

	for(int i=0;i<thisMesh.getNInternalFaces();i++)
	{
		const Face& f = faceArr[i];
			
		int owner = f.getOwner();
		int neighbor = f.getNeighbor();

		std::array<double,3> Sf = f.getFaceAreaVector();
		double fx =   f.getfx();

		std::array<double,3> uO = {u[0](owner), u[1](owner), 0.0};
		std::array<double,3> uN = {u[0](neighbor), u[1](neighbor), 0.0};
		std::array<double,3> uF = fx*uO + (1.0-fx)*uN;
		F(i) = dot(Sf, uF);
	}

	for(std::size_t i=0;i<bpArr.size();i++)
	{
		const int length = bpArr[i].getBoundaryPatchLength();
		const int start = bpArr[i].getBoundaryPatchStartFace();

		for(int j=0;j<length;j++)
		{
			const Face& f = faceArr[start+j];

			int owner = f.getOwner();

			std::array<double,3> Sf = f.getFaceAreaVector();
			std::array<double,3> gb = bpArr[i].getBoundaryPatchValue();

			if(bpArr[i].getBoundaryPatchType()==fixedValue)
			{
				F(start+j) = dot(Sf, gb);
			}
		}
	}
}

void correctF(arma::vec& F, fvMatrix& uMat, arma::vec& P, fvMesh& thisMesh)
{
	// Store gradP cell centroid values
	std::array<arma::vec,2> gradP;
	gradP[0] = arma::vec(P.size());
	gradP[1] = arma::vec(P.size());

	// Calculate gradP
	calculateGradP(gradP, P, thisMesh);
	
	std::cout << "P" << std::endl;
	P.print();
	std::cout << "gradP X" << std::endl;
	gradP[0].print();
	std::cout << "gradP Y" << std::endl;
	gradP[1].print();

	// Calculate 1/aP cell values
	const arma::vec rD = uMat.recipD();

	// Loop over internal faces
	std::vector<Face>& faceArr = thisMesh.allFaces();
	std::vector<BoundaryPatch>& bpArr = thisMesh.allBPs();

	for(int i=0;i<thisMesh.getNInternalFaces();i++)
	{
		std::array<double,3> Sf = faceArr[i].getFaceAreaVector();
		double fx = faceArr[i].getfx();
		
		int o = faceArr[i].getOwner();
		int n = faceArr[i].getNeighbor();

		// gradP face interpolation
		std::array<double,3> gradPf;
		gradPf[0] = fx*gradP[0](o) + (1.0-fx)*gradP[0](n);
		gradPf[1] = fx*gradP[1](o) + (1.0-fx)*gradP[1](n);
		gradPf[2] = 0.0;

		// 1/aP face reconstruction
		double rDf = fx*rD(o) + (1.0-fx)*rD(n);

		std::cout << "\no " << o
			  << "\nn " << n
			  << "\ngradPf " << gradPf[0]
			  << " "<< gradPf[1]
			  << "\nSf " << Sf[0] << " " << Sf[1]
			  << "\nrD " << rD[o] << " " << rD(n) << " " << rDf
			  << "\nF " << F(i);
		//Finally, correct F
		F(i) = F(i) - rDf*dot(Sf, gradPf);
		std::cout << "\nFcorr " << F(i) << std::endl;
	}

	// Now loop over boundary patches
	for(std::size_t i=0;i<bpArr.size();i++)
	{
		const int length = bpArr[i].getBoundaryPatchLength();
		const int start = bpArr[i].getBoundaryPatchStartFace();

		for(int j=0;j<length;j++)
		{
			if(bpArr[i].getBoundaryPatchType()==fixedValue)
			{
				const Face& f = faceArr[start+j];

				int o = f.getOwner();

				std::array<double,3> Sf = f.getFaceAreaVector();
				std::array<double,3> gradPf;
				gradPf[0] = gradP[0](o);
				gradPf[1] = gradP[1](o);
				gradPf[2] = 0.0;

				// 1/aP face reconstruction
				double rDf = rD(o);

				// Correct F
				//F(i) = F(i) - rDf*dot(Sf, gradPf);
			}
		}
	}
}

void correctU(fvMatrix& uMat, arma::vec& P,
		std::array<arma::vec,2>& u, fvMesh& thisMesh)
{
	// Store gradP cell centroid values
	std::array<arma::vec,2> gradP;
	gradP[0] = arma::vec(P.size());
	gradP[1] = arma::vec(P.size());

	// Calculate gradP
	calculateGradP(gradP, P, thisMesh);

	// Calculate 1/aP cell values
	const arma::vec rD = uMat.recipD();

	//Finally, loop through cells and correct cell centre velocities
	for(std::size_t k=0;k<thisMesh.allCells().size();k++)
	{
		u[0](k) = u[0](k) - rD(k)*gradP[0](k);
		u[1](k) = u[1](k) - rD(k)*gradP[1](k);
	}

}

void explicitUnderRelax(arma::vec &x, const arma::vec &xOld, const double &alphax)
{
	for(int i=0;i<x.size();i++)
	{
		x(i) = xOld(i) + alphax*(x(i)-xOld(i));
	}
}

void calculateGradP(std::array<arma::vec,2>& gradP, arma::vec& P, fvMesh& thisMesh)
{
	// Loop over internal faces
	std::vector<Face>& faceArr = thisMesh.allFaces();
	std::vector<BoundaryPatch>& bpArr = thisMesh.allBPs();

	for(int i=0;i<thisMesh.getNInternalFaces();i++)
	{
		std::array<double,3> Sf = faceArr[i].getFaceAreaVector();
		double fx = faceArr[i].getfx();
		
		int o = faceArr[i].getOwner();
		int n = faceArr[i].getNeighbor();

		// gradP calculations
		double Pf = fx*P(o) + (1.0-fx)*P(n);

		gradP[0](o) += Sf[0]*Pf;
		gradP[1](o) += Sf[1]*Pf;

		gradP[0](n) -= Sf[0]*Pf;
		gradP[1](n) -= Sf[1]*Pf;
	}
}

void divUCell(arma::vec& F, fvMesh& thisMesh)
{
	std::vector<double> Fluxes;
	Fluxes.resize(9);
	std::vector<Face>& faceArr = thisMesh.allFaces();
	for(int i=0;i<faceArr.size();i++)
	{
		int o = faceArr[i].getOwner();
		int n = faceArr[i].getNeighbor();
		
		Fluxes[o] += F(i);
		if(n>=0)
			Fluxes[n] -= F(i);
		else
			break;
	}

		std::cout << "Fluxes:" 
			  << "\n0 " << (fabs(Fluxes[0])<1e-18 ? 0.0 : Fluxes[0])
			  << "\n1 " << (fabs(Fluxes[1])<1e-18 ? 0.0 : Fluxes[1])
			  << "\n2 " << (fabs(Fluxes[2])<1e-18 ? 0.0 : Fluxes[2])
			  << "\n3 " << (fabs(Fluxes[3])<1e-18 ? 0.0 : Fluxes[3])
			  << "\n4 " << (fabs(Fluxes[4])<1e-18 ? 0.0 : Fluxes[4])
			  << "\n5 " << (fabs(Fluxes[5])<1e-18 ? 0.0 : Fluxes[5])
			  << "\n6 " << (fabs(Fluxes[6])<1e-18 ? 0.0 : Fluxes[6])
			  << "\n7 " << (fabs(Fluxes[7])<1e-18 ? 0.0 : Fluxes[7])
			  << "\n8 " << (fabs(Fluxes[8])<1e-18 ? 0.0 : Fluxes[8])
			  << std::endl;
}

void outputState(std::ofstream& file, unsigned int nCells, const arma::mat& x,
		const std::array<arma::vec,2>& u, const arma::vec& P)
{
	for(int i=0;i<nCells;i++)
	{
		file << x(i,0) << " "
		     << x(i,1) << " "
		     << u[0](i) << " "
		     << u[1](i) << " "
		     << P(i) << " "
		     << std::endl;
		if((i+1)%40==0)
			file << "\n";
	}
	file << "\n\n";
}

void initializeState(std::array<arma::vec, 2>& u, arma::vec& gamma, double& t1,
		int& tRes, double& Re, double& L, fvMesh& thisMesh, int testNum)
{
	const char* configFileName = "./stg/config.cfg";
	libconfig::Config cfg;
	
	try{
		cfg.readFile(configFileName);
	}
	catch(const libconfig::FileIOException &fioex){
		std::cerr << "I/O error while reading file." << std::endl;
	}
	catch (libconfig::ParseException &e){
        /*inform user about the parse exception*/
        std::cerr << "Parse error at " << e.getFile() << ":" << e.getLine()
                  << " - " << e.getError() << std::endl;
	}
	
	const libconfig::Setting& root = cfg.getRoot();

	try{
		const libconfig::Setting& tests = root["simulation"]["tests"];
		const libconfig::Setting& test = tests[testNum];

		double gamma0, ux, uy, uz;

		if(!(test.lookupValue("gamma0", gamma0)
			&& test.lookupValue("ux", ux)
			&& test.lookupValue("uy", uy)
			&& test.lookupValue("uz", uz)
			&& test.lookupValue("Re", Re)
			&& test.lookupValue("L", L)
			&& test.lookupValue("t1", t1)
			&& test.lookupValue("tRes", tRes)))
			std::cout << "Settings for test " << testNum+1 << " read in."
				<< std::endl;

		unsigned int nCells = gamma.size();

		for(unsigned int i=0;i<nCells;i++)
		{
			u[0](i)   = ux;
			u[1](i)   = uy;
		}
		calculateKinematicViscosity(gamma, u, L, Re);

		// Initialize Boundary Conditions from settings
		const libconfig::Setting& boundaries = tests[testNum]["inits"];

		int count = boundaries.getLength();

		for(int n=0;n<count;n++)
		{
			const libconfig::Setting& boundary = boundaries[n];
			std::string bpname;
			int type;
			patch p;
			std::array<double,3> value;

			if(!(
				boundary.lookupValue("patch", bpname)
				&& boundary.lookupValue("type", type)
				&& boundary.lookupValue("vx", value[0])
				&& boundary.lookupValue("vy", value[1])
				&& boundary.lookupValue("vz", value[2])
			    ))
				continue;

			switch (type)
			{
				case 1 :
				{
					p = fixedValue;
					break;
				}
				case 2 : 
				{
					p = fixedGrad;
					break;
				}
				case 3 :
				{
					p = empty;
					break;
				}
				case 4 :
				{
					p = mixed;
					break;
				}
			}

			std::vector<BoundaryPatch>& bpArr = thisMesh.allBPs();

			for(unsigned int n=0;n<bpArr.size();n++)
			{
				if(bpArr[n].getBoundaryPatchName()==bpname)
				{
					bpArr[n].setType(p);
					bpArr[n].setValue(value);
					break;
				}
			}
		}
	} catch(const libconfig::SettingNotFoundException &nfex){
		//Ignore
	}

}

void calculateKinematicViscosity(arma::vec& gamma, const std::array<arma::vec,2> u,
		double L, double Re)
{
	for(std::size_t i=0;i<gamma.size();i++)
	{
		gamma(i) = L/Re;
	}
}

void outputMeshDetails(fvMesh& thisMesh)
{

	std::cout << thisMesh.displayMeshDetails() << std::endl;
	std::cout << thisMesh.displayVolumesAndAreas() << std::endl;
	std::cout << thisMesh.displayCentroids() << std::endl;
	std::cout << thisMesh.displayBoundaryFaces() << std::endl;
	std::cout << thisMesh.displayCellNeighbors() << std::endl;
	std::cout << thisMesh.displayFaceOwnerNeighbor() << std::endl;
}
