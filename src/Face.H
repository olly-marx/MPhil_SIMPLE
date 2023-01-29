//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#ifndef FACE_HEADER
#define FACE_HEADER

#include <iostream>
#include <array>
#include <string>
#include <vector>

class Face{
        public:
		// Constructor
                Face(const std::vector<std::array<double,3>>&, std::vector<int>);

		// Public member functions
                std::vector<int> getFaceVertexIndices();
		std::array<double,3> getFaceCentroid() const;
		std::array<double,3> getFaceAreaVector() const;
        private:
		// Private member variables
                std::vector<int> m_vertices;
                std::array<double,3> m_faceCentroid;
                std::array<double,3> m_faceAreaVector;

		// Private members funtions
                std::array<double,3> calcFaceCentroid(const std::vector<std::array<double,3>>&);
                std::array<double,3> calcFaceAreaVector(const std::vector<std::array<double,3>>&);
};

#endif

