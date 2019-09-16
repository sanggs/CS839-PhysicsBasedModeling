///////////////////////////////////////////////////////////////////////////////
//                   ALL STUDENTS COMPLETE THESE SECTIONS
// Title:            Assignment 1 : Build an animated USD
// Files:            main.cpp
// Semester:         (CS839) Fall 2019
//
// Author:           Sangeetha Grama Srinivasan
// Email:            sgsrinivasa2@wisc.edu
// CS Login:         sgsrinivasa2@cs.wisc.edu
// Lecturer's Name:  Prof. Eftychios Sifakis
//
//////////////////// STUDENTS WHO GET HELP FROM OTHER THAN THEIR PARTNER //////
//                   must fully acknowledge and credit those sources of help.
//                   Instructors and TAs do not have to be credited here,
//                   but tutors, roommates, relatives, strangers, etc do.
//
// Online sources:   * This code uses the code given at
//                   https://github.com/uwgraphics/PhysicsBasedModeling-Demos
//                   as the template and develops a simple animation using it.
//                   * The explanation provided in this link has been referred
//                   to for coming up with the math behind this simulation:
//                   http://mathforum.org/library/drmath/view/55327.html
//////////////////////////// 80 columns wide //////////////////////////////////


#include "pxr/pxr.h"


#include "pxr/usd/sdf/layer.h"
#include "pxr/usd/sdf/path.h"
#include "pxr/usd/usd/stage.h"
#include "pxr/usd/usdGeom/mesh.h"
#include "pxr/base/vt/array.h"

#include "pxr/base/gf/range3f.h"

#include <iostream>
#include <cmath>

PXR_NAMESPACE_USING_DIRECTIVE

template<class T, int d> // d is the dimension of the mesh elements, e.g. 3 for triangles, 4 for quads
struct AnimatedMesh
{
    SdfLayerRefPtr m_layer;
    UsdStageRefPtr m_stage;
    UsdGeomMesh m_mesh;
    UsdAttribute m_pointsAttribute;

    GfRange3f m_extent;
    int m_lastFrame;
    
    std::vector<std::array<int, d>> m_meshElements;
    std::vector<GfVec3f> m_particleX;

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

    void initializeParticles()
    {
        // Grab the points (Positions) attribute, and indicate it is time-varying
        m_pointsAttribute = m_mesh.GetPointsAttr();
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
            m_extent.UnionWith(pt);

        // Copy particleX into VtVec3fArray for Usd
        VtVec3fArray usdPoints;
        usdPoints.assign(m_particleX.begin(), m_particleX.end());
        
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
        m_mesh.GetExtentAttr().Set(extentArray);

        // Save USD file
        m_stage->GetRootLayer()->Save();
        std::cout << "USD file saved!" << std::endl;
    }
};

template<class T>
struct LatticeMesh : public AnimatedMesh<T, 4>
{
    using Base = AnimatedMesh<T, 4>;
    using Base::m_meshElements;
    using Base::m_particleX;
    using Base::initializeUSD;
    using Base::initializeTopology;
    using Base::initializeParticles;

    std::array<int, 2> m_cellSize; // dimensions in grid cells
    T m_gridDX;
    int m_nFrames;

    void initialize()
    {
        initializeUSD("assignment1.usda");

        // Create a Cartesian lattice topology
        for(int cell_i = 0; cell_i < m_cellSize[0]; cell_i++)
        for(int cell_j = 0; cell_j < m_cellSize[1]; cell_j++)
            m_meshElements.emplace_back(
                std::array<int, 4>{
                    gridToParticleID(cell_i  , cell_j  ), 
                    gridToParticleID(cell_i+1, cell_j  ),
                    gridToParticleID(cell_i+1, cell_j+1),
                    gridToParticleID(cell_i  , cell_j+1)
                }
            );
        initializeTopology();

        // Also initialize the associated particles
        for(int node_i = 0; node_i <= m_cellSize[0]; node_i++) {
            double u = (double) node_i / (double) m_cellSize[0];
            for(int node_j = 0; node_j <= m_cellSize[1]; node_j++) {
                double v = (double) node_j / (double) m_cellSize[0];
                m_particleX.emplace_back(v*v, m_gridDX * (T)node_j, u);
            }
        }
        initializeParticles();

    }


    void prepareFrame(const int frame)
    {
        constexpr double pi = 3.1415926535897932;
        double angle = ((double) frame / (double) m_nFrames) * (pi*5);
        
        // Animate frame to show the rolling of a surface into a conical form
        for(int node_i = 0, p = 0; node_i <= m_cellSize[0]; node_i++){
            double u = (double) node_i / (double) m_cellSize[0];
            for(int node_j = 0; node_j <= m_cellSize[1]; node_j++) {
                double v = (double) node_j / (double) m_cellSize[0];
                m_particleX[p][0] = (v * v) * std::cos(u * angle) - (u) * std::sin(u * angle);
                m_particleX[p][1] = m_gridDX * (T)node_j;
                m_particleX[p][2] = (u) * std::cos(u * angle) + (v * v) * std::sin(u * angle);
                p++;
            }
        }
    }


private:
    inline int gridToParticleID(const int i, const int j) { return i * (m_cellSize[1]+1) + j; }
};

int main(int argc, char *argv[])
{
    LatticeMesh<float> simulationMesh;
    simulationMesh.m_cellSize = { 40, 40 };
    simulationMesh.m_gridDX = 0.020;
    simulationMesh.m_nFrames = 100;

    // Initialize the simulation example
    simulationMesh.initialize();
    
    // Output the initial shape of the surface
    simulationMesh.writeFrame(0);

    // Construct the cone, output results at each frame
    for(int frame = 1; frame <= simulationMesh.m_nFrames; frame++){
        simulationMesh.prepareFrame(frame);
        simulationMesh.writeFrame(frame);
    }
    
    // Write the entire timeline to USD
    simulationMesh.writeUSD();

    return 0;
}

