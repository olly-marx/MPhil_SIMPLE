//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#include <iostream>
#include <array>
#include <string>
#include <vector>

#include "Face.H"
#include "BoundaryPatch.H"

// Constructor
BoundaryPatch::BoundaryPatch(int startFace, int length, std::string patchName)
{
	m_boundaryPatchName = patchName;
	m_startFace = startFace;
	m_length = length;
}

// Public member functions
std::string BoundaryPatch::getBoundaryPatchName() const
{
	return m_boundaryPatchName;
}

int BoundaryPatch::getBoundaryPatchStartFace() const
{
	return m_startFace;
}

int BoundaryPatch::getBoundaryPatchLength() const
{
	return m_length;
}

patch BoundaryPatch::getBoundaryPatchType() const
{
	return m_type;
}

std::array<double,3> BoundaryPatch::getBoundaryPatchValue() const
{
	return m_value;
}

void BoundaryPatch::setType(patch p)
{
	m_type = p;
}

void BoundaryPatch::setValue(std::array<double,3> v)
{
	m_value = v;
}
