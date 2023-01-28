//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#ifndef FACE_HEADER
#define FACE_HEADER

#include <iostream>
#include <array>
#include <string>
#include <vector>

class Face{
        public:
                Face(const std::vector<std::array<double,3>>&, std::vector<int>);
                std::vector<int> getFaceVertexIndices();
                std::array<double,3> calcFaceCentroid(const std::vector<std::array<double,3>>&);
                std::array<double,3> calFaceAreaVector(const std::vector<std::array<double,3>>&);
        private:
                std::vector<int> m_vertices;
                std::array<double,3> m_faceCentroid;
                std::array<double,3> m_faceAreaVector;
};

#endif

