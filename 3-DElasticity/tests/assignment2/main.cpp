/*
 ************** Assignment 2 ******************
 * Title:            Linear Elasticity
 * Files:            main.cpp
 * Semester:         (CS839) Fall 2019
 
 * Author:           Sangeetha Grama Srinivasan
 * Email:            sgsrinivasa2@wisc.edu
 * CS Login:         sgsrinivasa2
 * Lecturer's Name:  Prof. Eftychios Sifakis
 * Description:      This assignment simulates Linear Elasticity in 2D and 3D scenarios.
                     The macro THREE_DIMENSION simulates linear elasticity for a 3D "X".
                     The macro TWO_DIMENSION simulated linear elasticity for a 2D membrane.
                     The assignment is based on the Demos provided in class: This can
                     be found at https://github.com/uwgraphics/PhysicsBasedModeling-Demos
                     The 2D example is the same as the MembraneStretch1 example
                     given in the above link.
 NOTE: For 2D/3D elasticity, the macro has to be changed in main.cpp and AnimatedMesh.h files.
 */

#include "FiniteElementMesh.h"

#include <map>

#define THREE_DIMENSION
//#define TWO_DIMENSION

template<class T, int d>
struct LatticeMesh : public FiniteElementMesh<T, d>
{
    using Base = FiniteElementMesh<T, d>;
    
    // from AnimatedTetrahedonMesh
    using Base::m_meshElements;
    using Base::m_particleX;
    using Base::initializeUSD;
    using Base::initializeTopology;
    using Base::initializeParticles;
    using Vectord = typename Base::Vectord;
    
    // from FiniteElementMesh
    using Base::m_particleV;
    using Base::initializeUndeformedConfiguration;
    using Base::m_stepEndTime;
    
    std::array<int, d> m_cellSize; // dimensions in grid cells
    T m_gridDX;
    
    std::vector<std::array<int, d>> m_activeCells; // Marks the "active" cells in the lattice
    std::map<std::array<int, d>, int> m_activeNodes; // Maps the "active" nodes to their particle index
    
    std::vector<Vectord> m_particleUndeformedX;
    std::vector<int> m_leftHandleIndices;
    std::vector<int> m_rightHandleIndices;
    Vectord m_leftHandleVelocity;
    Vectord m_rightHandleVelocity;
    
    LatticeMesh():Base(1.e2, 1., 4., .05)
    {
#ifdef THREE_DIMENSION
        m_leftHandleVelocity  = Vectord(-.2, 0.0, 0.0);
        m_rightHandleVelocity = Vectord( .2, 0.0, 0.0);
#endif
#ifdef TWO_DIMENSION
        m_leftHandleVelocity  = Vectord(-.2, 0.0);
        m_rightHandleVelocity = Vectord( .2, 0.0);
#endif
    }
#ifdef THREE_DIMENSION
    void initialize3DX() {
        // Activates cells in the shape of a 3-D X.
        for(int cell_i = 0; cell_i < m_cellSize[0]; cell_i++) {
            
            std::set<int> yValues;
            int cell_j;
            cell_j = cell_i;
            yValues.insert(cell_j);
            cell_j = m_cellSize[1] - 1 - cell_i;
            yValues.insert(cell_j);
            
            if(cell_i > m_cellSize[1]/2 - 1 ) {
                
                std::cout << cell_i << " is greater than " << m_cellSize[1]/2 << std::endl;
                
                cell_j = m_cellSize[1] - 1 - cell_i - 1;
                yValues.insert(cell_j);
                cell_j = cell_i + 1;
                yValues.insert(cell_j);
            }
            else {
                std::cout << cell_i << " is less than or equal to " << m_cellSize[1]/2 << std::endl;
                cell_j = m_cellSize[1] - 1 - cell_i + 1;
                yValues.insert(cell_j);
                cell_j = cell_i - 1;
                yValues.insert(cell_j);
            }
            
            for(auto it = yValues.begin(); it != yValues.end(); it++) {
                
                int cell_j = *it;
                if(cell_j >= 0 && cell_j < m_cellSize[1] ) {
                    std::set<int> zValues;
                    int cell_k;
                    cell_k = cell_j;
                    zValues.insert(cell_k);
                    cell_k = m_cellSize[2] - 1 - cell_j;
                    zValues.insert(cell_k);
                    
                    if(cell_j > m_cellSize[1]/2) {
                        cell_k = m_cellSize[2] - 1 - cell_j + 1;
                        zValues.insert(cell_k);
                        cell_k = cell_j - 1;
                        zValues.insert(cell_k);
                    } else {
                        cell_k = m_cellSize[2] - 1 - cell_j - 1;
                        zValues.insert(cell_k);
                        cell_k = cell_j + 1;
                        zValues.insert(cell_k);
                    }
                    
                    for(auto it = zValues.begin(); it != zValues.end(); it++) {
                        m_activeCells.push_back(std::array<int, 3>{cell_i, cell_j, *it});
                    }
                    
                }
                
            }
        }
    }

    void initialize()
    {
        initializeUSD("assignment2.usda");
        
        initialize3DX();
        if(m_activeCells.size() == 0)
            throw std::logic_error("3D-X not initialised correctly");
        
        std::cout << "Created a model including " << m_activeCells.size() << " lattice cells" <<std::endl;
        
        // Create (uniquely numbered) particles at the node corners of active cells
        // Particle number is the size of the particle array
        for(const auto& cell: m_activeCells){
            std::cout << "Cell: " << cell[0] << " " << cell[1] << " "<< cell[2] << std::endl;
            std::array<int, 3> node;
            for(node[0] = cell[0]; node[0] <= cell[0]+1; node[0]++)
                for(node[1] = cell[1]; node[1] <= cell[1]+1; node[1]++)
                    for(node[2] = cell[2]; node[2] <= cell[2]+1; node[2]++){
                        auto search = m_activeNodes.find(node);
                        if(search == m_activeNodes.end()){ // Particle not yet created at this lattice node location -> make one
                            m_activeNodes.insert({node, m_particleX.size()});
                            m_particleX.emplace_back(m_gridDX * T(node[0]), m_gridDX * T(node[1]), m_gridDX * T(node[2]));
                        }
                    }
        }
        std::cout << "Model contains " << m_particleX.size() << " particles" << std::endl;
        
        // Make tetrahedra out of all active cells (6 tetrahedra per cell)
        for(const auto& cell: m_activeCells){
            std::cout << cell[0] << " " << cell[1] << " "<< cell[2] << std::endl;
            int vertexIndices[2][2][2];
            for(int i = 0; i <= 1; i++)
                for(int j = 0; j <= 1; j++)
                    for(int k = 0; k <= 1; k++){
                        std::array<int, 3> node{cell[0] + i, cell[1] + j, cell[2] + k};
                        
                        
                        auto search = m_activeNodes.find(node);
                        if(search != m_activeNodes.end())
                            vertexIndices[i][j][k] = search->second;
                        else
                            throw std::logic_error("particle at cell vertex not found");
                    }
            std::cout << vertexIndices[0][0][0] << std::endl;
            std::cout << vertexIndices[1][1][1] << std::endl;
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][0][0], vertexIndices[1][1][0], vertexIndices[1][1][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][0][0], vertexIndices[1][1][1], vertexIndices[1][0][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][0][1], vertexIndices[1][1][1], vertexIndices[0][0][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][1][1], vertexIndices[0][1][1], vertexIndices[0][0][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][1][1], vertexIndices[0][1][0], vertexIndices[0][1][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][1][0], vertexIndices[0][1][0], vertexIndices[1][1][1]});
        }
        
        // Perform the USD-specific initialization of topology & particles
        // (this will also create a boundary *surface* to visualuze
        
        initializeTopology();
        initializeParticles();
        
        // Check particle indexing in mesh
        
        for(const auto& element: m_meshElements)
            for(const auto vertex: element)
                if(vertex < 0 || vertex >= m_particleX.size())
                    throw std::logic_error("mismatch between mesh vertex and particle array");
        
        // Also resize the velocities to match
        m_particleV.resize(m_particleX.size(), Vectord::Zero());
        
        // Initialize rest shape matrices and particle mass
        initializeUndeformedConfiguration();
        
        // Also record rest shape
        m_particleUndeformedX = m_particleX;
        
        // Identify particles on left and right handles
        std::vector<std::array<int, 3>> leftCells;
        std::vector<std::array<int, 3>> rightCells;
        for(auto element: m_activeNodes) {
            std::array<int, 3> particle = element.first;
            if(particle[0] == 0) {
                leftCells.push_back(element.first);
                std::array<int, 3> adjParticle = {particle[0]+1, particle[1], particle[2]};
                leftCells.push_back(adjParticle);
            }
            else if (particle[0] == m_cellSize[0]) {
                rightCells.push_back(element.first);
                std::array<int, 3> adjParticle = {particle[0]-1, particle[1], particle[2]};
                rightCells.push_back(adjParticle);
            }
        }
        
        for(auto element : leftCells) {
            auto search = m_activeNodes.find(element);
            if(search != m_activeNodes.end())
                m_leftHandleIndices.push_back(search->second);
            else
                throw std::logic_error("Left particle at cell vertex not found");
        }
        for(auto element : rightCells) {
            auto search = m_activeNodes.find(element);
            if(search != m_activeNodes.end())
                m_rightHandleIndices.push_back(search->second);
            else
                throw std::logic_error("Right particle at cell vertex not found");
        }
        
    }
#endif
#ifdef TWO_DIMENSION
    inline int gridToParticleID(const int i, const int j) const { return i * (m_cellSize[1]+1) + j; }
    
    void initialize()
    {
        initializeUSD("assignment2.usda");
        
        // Create a Cartesian lattice topology
        for(int cell_i = 0; cell_i < m_cellSize[0]; cell_i++)
            for(int cell_j = 0; cell_j < m_cellSize[1]; cell_j++){
                m_meshElements.emplace_back(
                                            std::array<int, 3>{
                                                gridToParticleID(cell_i  , cell_j  ),
                                                gridToParticleID(cell_i+1, cell_j  ),
                                                gridToParticleID(cell_i+1, cell_j+1)
                                            }
                                            );
                m_meshElements.emplace_back(
                                            std::array<int, 3>{
                                                gridToParticleID(cell_i  , cell_j  ),
                                                gridToParticleID(cell_i+1, cell_j+1),
                                                gridToParticleID(cell_i  , cell_j+1)
                                            }
                                            );
            }
        initializeTopology();
        
        // Also initialize the associated particles
        for(int node_i = 0; node_i <= m_cellSize[0]; node_i++)
            for(int node_j = 0; node_j <= m_cellSize[1]; node_j++)
                m_particleX.emplace_back(m_gridDX * (T)node_i, m_gridDX * (T)node_j);
        initializeParticles();
        
        // Check particle indexing in mesh
        for(const auto& element: m_meshElements)
            for(const auto vertex: element)
                if(vertex < 0 || vertex >= m_particleX.size())
                    throw std::logic_error("mismatch between mesh vertex and particle array");
        
        // Also resize the velocities to match
        m_particleV.resize(m_particleX.size(), Vectord::Zero());
        
        // Initialize rest shape matrices and particle mass
        initializeUndeformedConfiguration();
        
        // Also record rest shape
        m_particleUndeformedX = m_particleX;
        
        // Identify particles on left and right handles
        for(int node_j = 0; node_j <= m_cellSize[1]; node_j++){
            m_leftHandleIndices.push_back(gridToParticleID(0, node_j));
            m_rightHandleIndices.push_back(gridToParticleID(m_cellSize[0], node_j));
        }
    }
#endif
    //Clears the computed elastic forces for handle particles
    //making a symmetric positive definite matrix for CG
    void clearConstrainedParticles(std::vector<Vectord>& x) override
    {
        for(const auto v: m_leftHandleIndices)
            x[v] = Vectord::Zero();
        for(const auto v: m_rightHandleIndices)
            x[v] = Vectord::Zero();
    }
    
    //The handles are moved for a duration of effectiveTime
    void setBoundaryConditions() override
    {
        T effectiveTime = std::min<T>(m_stepEndTime, 1.0);
        
        for(const auto v: m_leftHandleIndices){
            m_particleX[v] = m_particleUndeformedX[v] + effectiveTime * m_leftHandleVelocity;
            m_particleV[v] = m_leftHandleVelocity;
        }
        for(const auto v: m_rightHandleIndices){
            m_particleX[v] = m_particleUndeformedX[v] + effectiveTime * m_rightHandleVelocity;
            m_particleV[v] = m_rightHandleVelocity;
        }
    }
    
};

int main(int argc, char *argv[])
{
#ifdef THREE_DIMENSION
    LatticeMesh<float, 3> simulationMesh;
    simulationMesh.m_cellSize = { 22, 22, 22};
    simulationMesh.m_gridDX = 0.05;
    simulationMesh.m_nFrames = 1000;
    simulationMesh.m_subSteps = 1;
    simulationMesh.m_frameDt = 0.1;
    
    // Initialize the simulation example
    simulationMesh.initialize();
    
    // Output the initial shape of the mesh
    simulationMesh.writeFrame(0);
    
    // Perform the animation, output results at each frame
    for(int frame = 1; frame <= simulationMesh.m_nFrames; frame++){
        simulationMesh.simulateFrame(frame);
        simulationMesh.writeFrame(frame);
    }
    
    // Write the entire timeline to USD
    simulationMesh.writeUSD();
#endif
    
#ifdef TWO_DIMENSION
    LatticeMesh<float, 2> simulationMesh;
    simulationMesh.m_cellSize = { 40, 40 };
    simulationMesh.m_gridDX = 0.025;
    simulationMesh.m_nFrames = 50;
    simulationMesh.m_subSteps = 1;
    simulationMesh.m_frameDt = 0.1;
    
    // Initialize the simulation example
    simulationMesh.initialize();

    // Output the initial shape of the mesh
    simulationMesh.writeFrame(0);
    
    // Perform the animation, output results at each frame
    for(int frame = 1; frame <= simulationMesh.m_nFrames; frame++){
        simulationMesh.simulateFrame(frame);
        simulationMesh.writeFrame(frame);
    }
    
    // Write the entire timeline to USD
    simulationMesh.writeUSD();
#endif
    
    return 0;

}

