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
                Face(const std::vector<std::array<double,3>>&, std::vector<int>, int);

		// Public member functions
		int                  getFaceId();
                std::vector<int>     getFaceVertexIndices();
		std::array<double,3> getFaceCentroid() const;
		std::array<double,3> getFaceAreaVector() const;
		int                  getOwner() const;
		int                  getNeighbor() const;
		void                 setOwner(int);
		void                 setNeighbor(int);
		void                 makeInternalFace();
		bool                 isBoundary() const;

		double               getfx() const;
		void                 setfx(double fx);
		double               getDelta() const;
		void                 setDelta(double fx);
        private:
		// Private member variables
		int                  m_faceId;
                std::vector<int>     m_vertices;
                std::array<double,3> m_faceCentroid;
                std::array<double,3> m_faceAreaVector;
		int                  m_owner;
		int                  m_neighbor;
		bool                 m_isBoundary;
		int                  m_boundaryPatchId;
		
		double               m_fx;
		double               m_delta;

		// Private members funtions
                void                 calcFaceCentroid(const std::vector<std::array<double,3>>&);
                void                 calcFaceAreaVector(const std::vector<std::array<double,3>>&);
};

#endif

