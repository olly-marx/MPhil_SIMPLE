//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#include <iomanip>
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
	fvMatrix::implicitUnderRelax(uMat, u, thisMesh, alphaU);
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
	// Calculate 1/aP cell values
	const arma::vec rD = uMat.recipD();

	// Loop over internal faces
	std::vector<Face>& faceArr = thisMesh.allFaces();
	std::vector<BoundaryPatch>& bpArr = thisMesh.allBPs();

	for(int i=0;i<thisMesh.getNInternalFaces();i++)
	{
		std::array<double,3> Sf = faceArr[i].getFaceAreaVector();
		double fx = faceArr[i].getfx();
		double dx = faceArr[i].getDelta();
		
		int o = faceArr[i].getOwner();
		int n = faceArr[i].getNeighbor();

		// 1/aP face reconstruction
		double rDf = fx*rD(o) + (1.0-fx)*rD(n);
		//Finally, correct F
		F(i) = F(i) - rDf*mod(Sf)*dx*(P(n)-P(o));
	}
}

void correctU(fvMatrix& uMat, arma::vec& P,
		std::array<arma::vec,2>& u, fvMesh& thisMesh)
{
	// Store gradP cell centroid values
	std::array<arma::vec,2> gradP;
	gradP[0] = arma::vec(thisMesh.getNCells());
	gradP[1] = arma::vec(thisMesh.getNCells());

	// Calculate gradP
	calculateGradP(gradP, P, thisMesh);

	// Calculate 1/aP cell values
	const arma::vec rD = uMat.recipD();

	//Finally, loop through cells and correct cell centre velocities
	for(int k=0;k<thisMesh.getNCells();k++)
	{
		//if(k==1561)
		//{
		//	std::cout << std::setw(10) 
		//		  << "\n1/aP = " << rD(k)
		//		  << "\ngradP = ( " << gradP[0](k) << " , " << gradP[1](k) << " )"
		//		  << "\nu = ( " << u[0](k) << " , " << u[1](k) << " )"
		//		  << "\nu_corr = ( " << u[0](k) - rD(k)*gradP[0](k) << " , " << u[1](k) - rD(k)*gradP[1](k) << " )"
		//		  << std::endl;
		//}
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
		double df = faceArr[i].getDelta();
		
		int o = faceArr[i].getOwner();
		int n = faceArr[i].getNeighbor();

		// gradP calculations
		double Pf = fx*P(o) + (1.0-fx)*P(n);

		gradP[0](o) += Sf[0]*Pf*df;
		gradP[1](o) += Sf[1]*Pf*df;

		gradP[0](n) -= Sf[0]*Pf*df;
		gradP[1](n) -= Sf[1]*Pf*df;
	}

	for(std::size_t i=0;i<bpArr.size();i++)
	{
		const int length = bpArr[i].getBoundaryPatchLength();
		const int start = bpArr[i].getBoundaryPatchStartFace();

		for(int j=0;j<length;j++)
		{
			const Face& f = faceArr[start+j];

			int o = f.getOwner();

			std::array<double,3> Sf = f.getFaceAreaVector();
			double df = faceArr[i].getDelta();

			if(bpArr[i].getBoundaryPatchType()==fixedValue)
			{
				gradP[0](o) += Sf[0]*P(o)*df;
				gradP[1](o) += Sf[1]*P(o)*df;
			}
		}
	}
}

void divU(arma::vec& F, fvMesh& thisMesh)
{
	arma::vec localFluxes(thisMesh.getNCells());

	std::vector<Face>& faceArr = thisMesh.allFaces();
	for(int i=0;i<faceArr.size();i++)
	{
		int o = faceArr[i].getOwner();
		int n = faceArr[i].getNeighbor();
		
		localFluxes(o) += F(i);
		if(n>=0)
		{
			localFluxes(n) -= F(i);
		}
	}

	double localAve = 0.0;
	for(int i=0;i<thisMesh.getNCells();i++)
	{
		localAve += localFluxes(i);
	}
	localAve = localAve / thisMesh.getNCells();

	std::cout << "\n----------- divU -----------\n"
		  << "\n local flux in internal cell : " << localFluxes(823) << "\n"
		  << "\n local flux in top (lid) cell : " << localFluxes(1570) << "\n"
		  << "\n local ave flux : " << localAve << "\n"
		  << "\n----------------------------\n" << std::endl;
}

void outputState(std::ofstream& file, unsigned int nCells, const arma::mat& x,
		const std::array<arma::vec,2>& u, const arma::vec& P,
		const arma::vec& F, fvMesh& thisMesh)
{
	arma::vec localFluxes(thisMesh.getNCells());
	std::vector<Face>& faceArr = thisMesh.allFaces();
	for(int i=0;i<faceArr.size();i++)
	{
		int o = faceArr[i].getOwner();
		int n = faceArr[i].getNeighbor();
		
		localFluxes(o) += F(i);
		if(n>=0)
			localFluxes(n) -= F(i);
	}
	for(int i=0;i<nCells;i++)
	{
		file << x(i,0) << " "
		     << x(i,1) << " "
		     << u[0](i) << " "
		     << u[1](i) << " "
		     << localFluxes(i) << " "
		     << P(i) << " "
		     << std::endl;
		if((i+1)%40==0)
			file << "\n";
	}
	file << "\n\n";
}

void buildReport(arma::mat& report, const fvMatrix& uMat, const fvMatrix& pMat,
		const fvMesh& thisMesh)
{
	int internalCell = 823;
	std::vector<int> intNeighbors = thisMesh.cellNeighbors(internalCell);
	int topCell = 1570;
	std::vector<int> topNeighbors = thisMesh.cellNeighbors(topCell);

	report(0,0) = uMat.A(internalCell,internalCell);
	report(1,0) = uMat.A(internalCell, intNeighbors[0]);
	report(2,0) = uMat.A(internalCell, intNeighbors[1]);
	report(3,0) = uMat.A(internalCell, intNeighbors[2]);
	report(4,0) = uMat.A(internalCell, intNeighbors[3]);
	report(5,0) = uMat.bx(internalCell);
	report(6,0) = uMat.by(internalCell);

	report(0,1) = uMat.A(topCell,topCell);
	report(1,1) = uMat.A(topCell, topNeighbors[0]);
	report(2,1) = uMat.A(topCell, topNeighbors[1]);
	report(3,1) = uMat.A(topCell, topNeighbors[2]);
	report(4,1) = uMat.bx(topCell);
	report(5,1) = uMat.by(topCell);

	report(0,2) = pMat.A(internalCell,internalCell);
	report(1,2) = pMat.A(internalCell, intNeighbors[0]);
	report(2,2) = pMat.A(internalCell, intNeighbors[1]);
	report(3,2) = pMat.A(internalCell, intNeighbors[2]);
	report(4,2) = pMat.A(internalCell, intNeighbors[3]);
	report(5,2) = pMat.bx(internalCell);

	report(0,3) = pMat.A(topCell,topCell);
	report(1,3) = pMat.A(topCell, topNeighbors[0]);
	report(2,3) = pMat.A(topCell, topNeighbors[1]);
	report(3,3) = pMat.A(topCell, topNeighbors[2]);
	report(4,3) = pMat.bx(topCell);
}

void printReport(arma::mat &report)
{
	std::cout << "#################### REPORT ####################\n\n";
	std::cout << "Initial momentum equation coefficients for an internal cell\n\n";

	std::cout << std::setw(5) << "Ap = "   << report(0,0) << "\n" 
		                  << "An_1 = " << report(1,0) << "\n"
		                  << "An_2 = " << report(2,0) << "\n"
		                  << "An_3 = " << report(3,0) << "\n"
		                  << "An_4 = " << report(4,0) << "\n"
		                  << "bx = "   << report(5,0) << "\n"
		                  << "by = "   << report(6,0) << "\n"
				  << std::endl;

	std::cout << "Initial momentum equation coefficients for a top (lid) cell\n\n";
	std::cout << std::setw(5) << "Ap = "   << report(0,1) << "\n" 
		                  << "An_1 = " << report(1,1) << "\n"
		                  << "An_2 = " << report(2,1) << "\n"
		                  << "An_3 = " << report(3,1) << "\n"
		                  << "bx = "   << report(4,1) << "\n"
		                  << "by = "   << report(5,1) << "\n"
				  << std::endl;

	std::cout << "Initial pressure equation coefficients for an internal cell\n\n";
	std::cout << std::setw(5) << "Ap = "   << report(0,2) << "\n" 
		                  << "An_1 = " << report(1,2) << "\n"
		                  << "An_2 = " << report(2,2) << "\n"
		                  << "An_3 = " << report(3,2) << "\n"
		                  << "An_4 = " << report(4,2) << "\n"
		                  << "b = "    << report(5,2) << "\n"
				  << std::endl;

	std::cout << "Initial pressure equation coefficients for an top (lid) cell\n\n";
	std::cout << std::setw(5) << "Ap = "   << report(0,3) << "\n" 
		                  << "An_1 = " << report(1,3) << "\n"
		                  << "An_2 = " << report(2,3) << "\n"
		                  << "An_3 = " << report(3,3) << "\n"
		                  << "b = "    << report(4,3) << "\n"
				  << std::endl;
}

void initializeState(std::array<arma::vec, 2>& u, arma::vec& gamma, double& t1,
		int& tRes, double& Re, double& L, double& alphaU, double& alphaP,
		fvMesh& thisMesh, int testNum)
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
			&& test.lookupValue("alphaU", alphaU)
			&& test.lookupValue("alphaP", alphaP)
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
