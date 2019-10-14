#pragma once

#include <iostream>

#include "pxr/pxr.h"

#include "pxr/usd/sdf/layer.h"
#include "pxr/usd/sdf/path.h"
#include "pxr/usd/usd/stage.h"
#include "pxr/usd/usdGeom/mesh.h"
#include "pxr/base/vt/array.h"

#include "pxr/base/gf/range3f.h"

#include <Eigen/Dense>

#define THREE_DIMENSION
//#define TWO_DIMENSION

PXR_NAMESPACE_USING_DIRECTIVE

template<class T, int d> // d is the dimension of the mesh elements, e.g. 3 for triangles, 4 for quads
struct AnimatedMesh
{
    using Vectord = Eigen::Matrix< T , d , 1>;
    
    SdfLayerRefPtr m_layer;
    UsdStageRefPtr m_stage;
#ifdef THREE_DIMENSION
    UsdGeomMesh m_boundaryMesh;
#endif
#ifdef TWO_DIMENSION
    UsdGeomMesh m_mesh;
#endif
    UsdAttribute m_pointsAttribute;

    GfRange3f m_extent;
    int m_lastFrame;
    
    std::vector<std::array<int, d+1>> m_meshElements;
#ifdef THREE_DIMENSION
    std::vector<std::array<int, d>> m_meshBoundaryElements;
#endif
    std::vector<Vectord> m_particleX;

    AnimatedMesh()
        :m_lastFrame(-1)
    {}
    
    void initializeUSD(const std::string filename)
    {
        // Create the layer to populate.
        m_layer = SdfLayer::CreateNew(filename);

        // Create a UsdStage with that root layer.
        m_stage = UsdStage::Open(m_layer);
    }

#ifdef THREE_DIMENSION
    void initializeTopology()
    {
        // Create a mesh for this surface
        m_boundaryMesh = UsdGeomMesh::Define(m_stage, SdfPath("/MeshBoundarySurface"));

        // Detect boundary mesh elements

        // Boundary Mesh Creation
        //   Stage 1: First, compute all candidate triangle faces

        std::vector<std::array<int, 3>> m_candidateBoundaryElements;
        for (const auto& tetElement : m_meshElements){
            m_candidateBoundaryElements.push_back(std::array<int, 3>{ tetElement[0], tetElement[2], tetElement[1] });
            m_candidateBoundaryElements.push_back(std::array<int, 3>{ tetElement[0], tetElement[1], tetElement[3] });
            m_candidateBoundaryElements.push_back(std::array<int, 3>{ tetElement[0], tetElement[3], tetElement[2] });
            m_candidateBoundaryElements.push_back(std::array<int, 3>{ tetElement[1], tetElement[2], tetElement[3] });
        }

        //   Stage 2: Mark occurences of triangles

        std::map<std::array<int, 3>, int> triangleOccurences;
        for (const auto& triElement : m_candidateBoundaryElements){
            std::array<int, 3> sortedTriangle = triElement;
            std::sort(sortedTriangle.begin(), sortedTriangle.end());
            auto search = triangleOccurences.find(sortedTriangle);
            if(search != triangleOccurences.end())
                search->second++;
            else
                triangleOccurences.insert({sortedTriangle, 1});
        }

        //   Stage 3: Iterate over candidate triangles again, use only those that appeared exactly once

        for (const auto& triElement : m_candidateBoundaryElements){
            std::array<int, 3> sortedTriangle = triElement;
            std::sort(sortedTriangle.begin(), sortedTriangle.end());
            auto search = triangleOccurences.find(sortedTriangle);
            if(search == triangleOccurences.end())
                throw std::logic_error("triangle not found");
            else if(search->second == 1)
                m_meshBoundaryElements.push_back(triElement);
            else if(search->second != 2)
                throw std::logic_error("triangle belonging to more than 2 tetrahedra found");
        }

        std::cout << "Boundary mesh contains " << m_meshBoundaryElements.size() << " triangles" << std::endl;
        
        // Create appropriate buffers for vertex counts and indices, and populate them
        VtIntArray faceVertexCounts, faceVertexIndices;
        for (const auto& element : m_meshBoundaryElements) {
            faceVertexCounts.push_back(element.size());
            for (const auto& vertex : element)
                faceVertexIndices.push_back(vertex);
        }
        
        // Now set the attributes
        m_boundaryMesh.GetFaceVertexCountsAttr().Set(faceVertexCounts);
        m_boundaryMesh.GetFaceVertexIndicesAttr().Set(faceVertexIndices);
    }
#endif

#ifdef TWO_DIMENSION
    void initializeTopology()
    {
        // Create a mesh for this surface
        m_mesh = UsdGeomMesh::Define(m_stage, SdfPath("/MeshSurface"));
        
        // Create appropriate buffers for vertex counts and indices, and populate them
        VtIntArray faceVertexCounts, faceVertexIndices;
        for (const auto& element : m_meshElements) {
            faceVertexCounts.push_back(element.size());
            for (const auto& vertex : element)
                faceVertexIndices.push_back(vertex);
        }
        
        // Now set the attributes
        m_mesh.GetFaceVertexCountsAttr().Set(faceVertexCounts);
        m_mesh.GetFaceVertexIndicesAttr().Set(faceVertexIndices);
    }
#endif
    
    void initializeParticles()
    {
        // Grab the points (Positions) attribute, and indicate it is time-varying
#ifdef THREE_DIMENSION
        m_pointsAttribute = m_boundaryMesh.GetPointsAttr();
#endif
#ifdef TWO_DIMENSION
        m_pointsAttribute = m_mesh.GetPointsAttr();
#endif
        m_pointsAttribute.SetVariability(SdfVariabilityVarying);
    }
    
    void writeFrame(const int frame)
    {
        std::cout << "Writing frame " << frame << " ..." << std::endl;
        
        // Check that there are any particles to write at all
        if (m_particleX.empty())
            throw std::logic_error("Empty array of input vertices");

        // Check that frames have been written in sequence
        if(frame != m_lastFrame+1)
            throw std::logic_error("Non-consequtive frame sequence requested in writeFrame()");
        m_lastFrame = frame;

        // Update extent
        for (const auto& pt : m_particleX)
#ifdef THREE_DIMENSION
            m_extent.UnionWith(GfVec3f(pt[0], pt[1], pt[2]));
#endif
#ifdef TWO_DIMENSION
            m_extent.UnionWith(GfVec3f(pt[0],pt[1],T()));
#endif

        // Copy particleX into VtVec3fArray for Usd
        VtVec3fArray usdPoints(m_particleX.size());
        for(int p = 0; p < m_particleX.size(); p++)
#ifdef THREE_DIMENSION
            usdPoints[p] = GfVec3f(m_particleX[p][0], m_particleX[p][1], m_particleX[p][2]);
#endif
#ifdef TWO_DIMENSION
            usdPoints[p] = GfVec3f(m_particleX[p][0],m_particleX[p][1],T());
#endif
        
        // Write the points attribute for the given frame
        m_pointsAttribute.Set(usdPoints, (double) frame);
    }
    
    void writeUSD()
    {
        // Set up the timecode
        m_stage->SetStartTimeCode(0.);
        m_stage->SetEndTimeCode((double) m_lastFrame);

        // Set the effective extent
        VtVec3fArray extentArray(2);
        extentArray[0] = m_extent.GetMin();
        extentArray[1] = m_extent.GetMax();
#ifdef THREE_DIMENSION
        m_boundaryMesh.GetExtentAttr().Set(extentArray);
#endif
#ifdef TWO_DIMENSION
        m_mesh.GetExtentAttr().Set(extentArray);
#endif
        // Save USD file
        m_stage->GetRootLayer()->Save();
        std::cout << "USD file saved!" << std::endl;
    }
};

