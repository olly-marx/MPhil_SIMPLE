//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#include <cstddef>
#include <iostream>
#include <array>
#include <string>
#include <vector>

#include "Cell.H"
#include "fvMatrix.H"
#include "VectorUtils.H"
#include "BoundaryPatch.H"

fvMatrix::fvMatrix(int nCells)
{
	m_a.resize(nCells);

	for(int i=0;i<nCells;i++)
	{
		m_a[i].resize(nCells);
	}

	m_b.resize(nCells);
}

// Public member functions
void fvMatrix::discretizeRateofChange(fvMatrix& m, const std::vector<double>& T,
		const std::vector<Cell>& cellArr, const double& dt)
{
	for(std::size_t i=0;i<cellArr.size();i++)
	{
		m[i][i] += cellArr[i].getCellVolume() / dt;

		m.b(i) += cellArr[i].getCellVolume() * T[i] / dt;
	}
}

void fvMatrix::discretizeDiffusion(fvMatrix& m, const std::vector<double>& gamma,
		const std::vector<Face>& faceArr)
{
	for(std::size_t i=0;i<faceArr.size();i++)
	{
		int owner = faceArr[i].getOwner();
		int neighbor = faceArr[i].getNeighbor();

		double aN = gamma[i] * faceArr[i].getDelta() * mod( faceArr[i].getFaceAreaVector() );

		m[owner][owner] -= aN;
		m[owner][neighbor] = aN;
		m[neighbor][owner] = aN;
	}
}

void fvMatrix::discretizeDirichlet(fvMatrix& m, const std::vector<double>& T,
		const std::vector<BoundaryPatch>& bpArr, const std::vector<Face>& faceArr,
		const double& Tb)
{
	for(std::size_t i=0;i<bpArr.size();i++)
	{
		const int length = bpArr[i].getBoundaryPatchLength();
		const int start = bpArr[i].getBoundaryPatchStartFace();

		for(int j=0;j<length;j++)
		{
			const Face& f = faceArr[start+j];

			int owner = f.getOwner();

			double b = Tb * (1 / f.getDelta());
			double aP = - T[owner] * (1 / f.getDelta());

			m.b(owner) += b;
			m[owner][owner] += aP;
		}
	}
}

std::vector<double>& fvMatrix::operator[](int i)
{
	return m_a[i];
}

std::vector<double>  fvMatrix::operator[](int i) const
{
	return m_a[i];
}

double& fvMatrix::b(int i)
{
	return m_b[i];
}

double fvMatrix::b(int i) const
{
	return m_b[i];
}
