//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <array>
#include <string>
#include <vector>
#include <cmath>

#include "Cell.H"
#include "fvMatrix.H"
#include "VectorUtils.H"
#include "BoundaryPatch.H"

fvMatrix::fvMatrix(int nCells)
{
	m_A = arma::mat(nCells, nCells);
	m_b = arma::vec(nCells);

	std::cout << "Mat len " << m_A.size() << std::endl;
}

// Public member functions
void fvMatrix::discretizeRateofChange(fvMatrix& m, const arma::vec& T,
		const std::vector<Cell>& cellArr, const double& dt)
{
	for(std::size_t i=0;i<cellArr.size();i++)
	{
		double ap = cellArr[i].getCellVolume() / dt;

		//std::cout << "cell " << i << " : " 
		//	<< "\nV  " << cellArr[i].getCellVolume() 
		//	<< "\ndt " << dt 
		//	<< "\naP " << ap 
		//	<< std::endl;

		m.A(i,i) = ap;

		double TP = T(i);

		double b = cellArr[i].getCellVolume() * TP / dt;

		//std::cout << "cell " << i << " : " 
		//	<< "\nT " << T(i) 
		//	<< "\nV " << cellArr[i].getCellVolume() 
		//	<< "\ndt " << dt 
		//	<< "\nb " << b 
		//	<< std::endl;

		m.b(i) = b;
	}
}

void fvMatrix::discretizeDiffusion(fvMatrix& m, const arma::vec& gamma,
		const std::vector<Face>& faceArr)
{
	for(std::size_t i=0;i<faceArr.size();i++)
	{
		const Face& f = faceArr[i];
			
		int owner = f.getOwner();
		int neighbor = f.getNeighbor();

		double Sf =   mod (f.getFaceAreaVector());
		double df =   f.getDelta();

		double aN = gamma(owner) * df * Sf;

		//std::cout << "face " << i 
		//	<< " : gamma " << gamma(owner) 
		//	<< "\ndelta " << faceArr[i].getDelta() 
		//	<< "\narea " << mod(faceArr[i].getFaceAreaVector()) 
		//	<< "\naN " << aN 
		//	<< "\no " << owner 
		//	<< "\nn " << neighbor 
		//	<< std::endl;
 

		if(neighbor!=-1)
		{
			//std::cout << "B m["<<owner<<"]["<<owner<<"] = " << m.A(owner, owner) << std::endl;
			m.A(owner,owner) += aN;
			m.A(neighbor,neighbor) += aN;
			m.A(owner,neighbor) = -aN;
			m.A(neighbor,owner) = -aN;
			//std::cout << "A m["<<owner<<"]["<<owner<<"] = " << m.A(owner, owner) << std::endl;
		}
	}
}

void fvMatrix::discretizeConvectionCentral(fvMatrix& m, 
		const arma::mat& u, const std::vector<Face>& faceArr)
{
	for(std::size_t i=0;i<faceArr.size();i++)
	{
		const Face& f = faceArr[i];
			
		int owner = f.getOwner();
		int neighbor = f.getNeighbor();

		if(neighbor!=-1)
		{
			std::array<double,3> Sf = f.getFaceAreaVector();
			std::array<double,3> ui = {u(owner,0), u(owner,1), u(owner,2)};
			double F = dot(Sf, ui);

			double fx =   f.getfx();

			double aP = fx*F;
			double aN = (1.0 - fx)*F;

			m.A(owner,owner) += aP;
			m.A(owner,neighbor) += aN;

			m.A(neighbor,neighbor) -= aP;
			m.A(neighbor,owner) -= aN;
		}
	}
}

void fvMatrix::discretizeConvectionUpwind(fvMatrix& m, 
		const arma::mat& u, const std::vector<Face>& faceArr)
{
	for(std::size_t i=0;i<faceArr.size();i++)
	{
		const Face& f = faceArr[i];
			
		int owner = f.getOwner();
		int neighbor = f.getNeighbor();

		if(neighbor!=-1)
		{
			std::array<double,3> Sf = f.getFaceAreaVector();
			std::array<double,3> ui = {u(owner,0), u(owner,1), u(owner,2)};
			double F = dot(Sf, ui);

			double fx =   f.getfx();

			double aP = std::max(F, 0.0);
			double aN = std::min(F, 0.0);
			
			m.A(owner,owner) += aP;
			m.A(owner,neighbor) += aN;

			aP = std::max(-F, 0.0);
			aN = std::min(-F, 0.0);

			m.A(neighbor,neighbor) += aP;
			m.A(neighbor,owner) += aN;
		}
	}
}

void fvMatrix::discretizeBoundaryConditions(fvMatrix& m, const arma::vec& T,
		const std::vector<BoundaryPatch>& bpArr, const arma::vec& gamma,
		const arma::mat& u, const std::vector<Face>& faceArr)
{
	for(std::size_t i=0;i<bpArr.size();i++)
	{
		const int length = bpArr[i].getBoundaryPatchLength();
		const int start = bpArr[i].getBoundaryPatchStartFace();

		for(int j=0;j<length;j++)
		{
			const Face& f = faceArr[start+j];

			int owner = f.getOwner();

			std::array<double,3> Sf = f.getFaceAreaVector();
			std::array<double,3> ui = {u(owner,0), u(owner,1), u(owner,2)};
			double F = dot(Sf, ui);

			double gb = bpArr[i].getBoundaryPatchValue();

			if(bpArr[i].getBoundaryPatchType()==fixedValue)
			{
				// Dirichlet Diffusion
				double df = f.getDelta();

				double b  =   gamma(owner) * mod(Sf) * gb * df;
				double aP =   gamma(owner) * mod(Sf) * df;

				m.b(owner) += b;
				m.A(owner, owner) += aP;

				// Dirichlet Convection
				b = std::max(F,0.0) * gb;

				// NOT SURE IF THIS SHOULD BE NEGATIVE
				m.b(owner) -= b;
			}
			else if(bpArr[i].getBoundaryPatchType()==fixedGrad)
			{

				// Neumann Diffusion
				m.b(owner) += gb;

				// Neumann Convection
				double df = 1.0/f.getDelta();
				double Tb = T(owner) + df*gb;
				
				double b  = std::max(F, 0.0) * Tb;

				// FOR CONSISTENCY THIS IS NEGATIVE
				m.b(owner) -= b;
			}
		}
	}
}

void fvMatrix::resetMatrix(fvMatrix& m, const std::vector<Cell>& cellArr)
{
	for(std::size_t i=0;i<cellArr.size();i++)
	{
		for(std::size_t j=0;j<cellArr.size();j++)
		{
			m.A(i,j) = 0.0;
		}
		m.b(i) = 0.0;
	}
}

void fvMatrix::solveLinearSystem(arma::vec& T)
{
	T = arma::solve(m_A, m_b);
}

double& fvMatrix::A(int i, int j)
{
	return m_A(i,j);
}

double fvMatrix::A(int i, int j) const
{
	return m_A(i,j);
}

double& fvMatrix::b(int i)
{
	return m_b(i);
}

double fvMatrix::b(int i) const
{
	return m_b(i);
}

void fvMatrix::printDiscretization()
{
	for(int j=0;j<m_b.size();j++)
	{
		std::cout << " | ";
		for(int i=0;i<m_b.size();i++)
		{
			std::cout << m_A(i,j) << " ";
		}
		std::cout << " | | x" << j << " | = | " << m_b[j] << " |" << std::endl;
	}
}
