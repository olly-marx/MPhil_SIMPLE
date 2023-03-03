//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#ifndef BOUNDARY_HEADER
#define BOUNDARY_HEADER

#include <iostream>
#include <array>
#include <string>
#include <vector>

#include "Face.H"

enum patch { fixedValue=1, fixedGrad=2, empty=3, mixed=4 };

class BoundaryPatch{
        public:
		// Constructor
                BoundaryPatch(int, int, std::string);

		// Public member functions
		std::string getBoundaryPatchName() const;
		int         getBoundaryPatchStartFace() const;
		int         getBoundaryPatchLength() const;
		patch       getBoundaryPatchType() const;
		double      getBoundaryPatchValue() const;
		void        setType(patch);
		void        setValue(double);

        private:
		// Private member variables
		std::string m_boundaryPatchName;
                int         m_startFace;
                int         m_length;
		patch       m_type;
		double      m_value;

		// Private members funtions
};

#endif

