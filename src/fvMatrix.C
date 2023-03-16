//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <array>
#include <string>
#include <vector>
#include <cmath>

#include "Cell.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "VectorUtils.H"
#include "BoundaryPatch.H"

fvMatrix::fvMatrix(fvMesh& mesh)
{
	nCells = mesh.getNCells();
	m_bx = arma::vec(nCells);
	m_by = arma::vec(nCells);
	m_A  = arma::mat(nCells, nCells);
}

// Public member functions
void fvMatrix::discretizeRateofChange(fvMatrix& m, const std::array<arma::vec,2>& T,
		fvMesh& thisMesh, const double& dt)
{
	const std::vector<Cell>& cellArr = thisMesh.allCells();

	for(std::size_t i=0;i<cellArr.size();i++)
	{
		double ap = cellArr[i].getCellVolume() / dt;

		m.A(i,i) = ap;

		double TPx = T[0](i);
		double TPy = T[1](i);

		double bx = cellArr[i].getCellVolume() * TPx / dt;
		double by = cellArr[i].getCellVolume() * TPy / dt;

		m.bx(i) = bx;
		m.by(i) = by;
	}
}

void fvMatrix::discretizeLaplacian(fvMatrix& m, const arma::vec& gamma,
		fvMesh& thisMesh, bool pressure)
{
	const std::vector<BoundaryPatch>& bpArr = thisMesh.allBPs();
	const std::vector<Face>& faceArr = thisMesh.allFaces();

	for(int i=0;i<thisMesh.getNInternalFaces();i++)
	{
		const Face& f = faceArr[i];
			
		int owner = f.getOwner();
		int neighbor = f.getNeighbor();

		double Sf =   mod (f.getFaceAreaVector());
		double df =   f.getDelta();
		double fx =   f.getfx();

		double aN = (fx*gamma(owner) + (1.0-fx)*gamma(neighbor)) * df * Sf;

		m.A(owner,owner) += pressure ? -aN : aN;
		m.A(neighbor,neighbor) += pressure ? -aN : aN;
		m.A(owner,neighbor) += pressure ? aN : -aN;
		m.A(neighbor,owner) += pressure ? aN : -aN;
	}
	
	for(std::size_t i=0;i<bpArr.size();i++)
	{
		const int length = bpArr[i].getBoundaryPatchLength();
		const int start = bpArr[i].getBoundaryPatchStartFace();

		for(int j=0;j<length;j++)
		{
			const Face& f = faceArr[start+j];

			int owner = f.getOwner();

			std::array<double,3> gb = bpArr[i].getBoundaryPatchValue();

			if(bpArr[i].getBoundaryPatchType()==fixedValue && !pressure)
			{
				// Dirichlet Diffusion
				std::array<double,3> Sf = f.getFaceAreaVector();
				double               df = f.getDelta();

				double bx =   gamma(owner) * mod(Sf) * gb[0] * df;
				double by =   gamma(owner) * mod(Sf) * gb[1] * df;
				double aP =   gamma(owner) * mod(Sf) * df;

				m.bx(owner) += bx;
				m.by(owner) += by;
				m.A(owner, owner) += aP;
			}
			else if(bpArr[i].getBoundaryPatchType()==fixedValue && pressure)
			{
				// Neumann Diffusion for pressure
				//std::array<double,3> Sf = f.getFaceAreaVector();
				//double               df = f.getDelta();

				//double bx =   gamma(owner) * mod(Sf) * P(owner) * df;
				//double by =   gamma(owner) * mod(Sf) * P(owner) * df;
				//double aP =   gamma(owner) * mod(Sf) * df;

				//m.bx(owner) += pressure ? -bx : bx;
				//m.by(owner) += pressure ? -by : by;
				//m.A(owner, owner) += pressure ? -aP : aP;
			}
			else if(bpArr[i].getBoundaryPatchType()==fixedGrad)
			{
				// Neumann Diffusion
				m.bx(owner) += gb[0];
				m.by(owner) += gb[1];
			}
		}
	}
}

void fvMatrix::discretizeConvectionUpwind(fvMatrix& m, fvMesh& thisMesh,
		const arma::vec& F, std::array<arma::vec,2>& T)
{
	const std::vector<BoundaryPatch>& bpArr = thisMesh.allBPs();
	const std::vector<Face>& faceArr = thisMesh.allFaces();

	for(int i=0;i<thisMesh.getNInternalFaces();i++)
	{
		const Face& f = faceArr[i];
			
		int owner = f.getOwner();
		int neighbor = f.getNeighbor();

		double fx =   f.getfx();

		double aP = std::max(F(i), 0.0);
		double aN = std::min(F(i), 0.0);
		
		m.A(owner,owner) += aP;
		m.A(owner,neighbor) += aN;

		aP = std::max(-F(i), 0.0);
		aN = std::min(-F(i), 0.0);

		m.A(neighbor,neighbor) += aP;
		m.A(neighbor,owner) += aN;
	}
	
	for(std::size_t i=0;i<bpArr.size();i++)
	{
		const int length = bpArr[i].getBoundaryPatchLength();
		const int start = bpArr[i].getBoundaryPatchStartFace();

		for(int j=0;j<length;j++)
		{
			const Face& f = faceArr[start+j];

			int owner = f.getOwner();

			std::array<double,3> gb = bpArr[i].getBoundaryPatchValue();

			if(bpArr[i].getBoundaryPatchType()==fixedValue)
			{
				// Dirichlet Convection
				double bx = std::min(F(start+j),0.0) * gb[0];
				double by = std::min(F(start+j),0.0) * gb[1];

				// NOT SURE IF THIS SHOULD BE NEGATIVE
				m.bx(owner) -= bx;
				m.by(owner) -= by;
			}
			else if(bpArr[i].getBoundaryPatchType()==fixedGrad)
			{
				// Neumann Convection
				double df = 1.0/f.getDelta();
				double Tbx = T[0](owner) + df*gb[0];
				double Tby = T[1](owner) + df*gb[1];
				
				double bx  = std::max(F(start+j), 0.0) * Tbx;
				double by  = std::max(F(start+j), 0.0) * Tby;

				// FOR CONSISTENCY THIS IS NEGATIVE
				m.bx(owner) -= bx;
				m.by(owner) -= by;
			}
		}
	}
}

void fvMatrix::discretizeContinuity(fvMatrix& m, fvMesh& thisMesh, 
		const arma::vec& F)
{
	std::vector<Face>& faceArr = thisMesh.allFaces();
	std::vector<BoundaryPatch>& bpArr = thisMesh.allBPs();

	for(int i=0;i<faceArr.size();i++)
	{
		int o = faceArr[i].getOwner();
		int n = faceArr[i].getNeighbor();
		
		m.bx(o) += F(i);
		if(n>=0)
			m.bx(n) -= F(i);
	}
}

void fvMatrix::calculateUResidual(arma::vec& res, fvMatrix& m,
		arma::vec& x, fvMesh& thisMesh, double& norm)
{
	for(std::size_t i=0;i<thisMesh.allCells().size();i++)
	{
		res(i) += m.bx(i);
		res(i) -= m.A(i,i)*x(i);
	}

	const std::vector<Face>& faceArr = thisMesh.allFaces();

	for(int i=0;i<thisMesh.getNInternalFaces();i++)
	{
		int o = faceArr[i].getOwner();
		int n = faceArr[i].getNeighbor();

		res(o) -= m.A(o,n)*x(n);
		res(n) -= m.A(n,o)*x(o);
	}

	if(norm==-1.0)
		norm = arma::norm(res);
}

void fvMatrix::implicitUnderRelax(fvMatrix &m, std::array<arma::vec,2> &u, fvMesh &thisMesh, double alphaU)
{
	int nCells = thisMesh.allCells().size();
	for(int i=0;i<nCells;i++)
	{
		m.A(i,i) = m.A(i,i)/alphaU;
		m.bx(i) += (1.0-alphaU)*m.A(i,i)*u[0](i)/alphaU;
		m.by(i) += (1.0-alphaU)*m.A(i,i)*u[1](i)/alphaU;
	}
}

void fvMatrix::resetMatrix(fvMatrix& m, int nCells)
{
	for(int i=0;i<nCells;i++)
	{
		for(int j=0;j<nCells;j++)
		{
			m.A(i,j) = 0.0;
		}
		m.bx(i) = 0.0;
		m.by(i) = 0.0;
	}
}

void fvMatrix::solveLinearSystem(arma::vec& T, bool x)
{
	//std::cout << "Solving " << (x ? "bx" : "by") << std::endl;
	T = arma::solve(m_A, (x ? m_bx : m_by), arma::solve_opts::force_approx);
}

arma::vec fvMatrix::recipD() const
{
	arma::vec rD(nCells);
	for(int i=0;i<nCells;i++)
		rD(i) = 1.0 / m_A(i,i);
	return rD;
}

double& fvMatrix::A(int i, int j)
{
	return m_A(i,j);
}

double fvMatrix::A(int i, int j) const
{
	return m_A(i,j);
}

double& fvMatrix::bx(int i)
{
	return m_bx(i);
}

double fvMatrix::bx(int i) const
{
	return m_bx(i);
}

double& fvMatrix::by(int i)
{
	return m_by(i);
}

double fvMatrix::by(int i) const
{
	return m_by(i);
}

void fvMatrix::printDiscretization()
{
	for(int j=0;j<m_bx.size();j++)
	{
		std::cout<< std::setprecision(2)  << "|";
		for(int i=0;i<m_bx.size();i++)
		{
			std::cout << std::setw(10) << m_A(i,j);
		}
		std::cout << "| | x" << j << " | = | " 
				<< std::setw(10) << m_bx(j) << "\t| , | " 
				<< m_by(j) << " |" 
				<< std::endl;
	}
}
