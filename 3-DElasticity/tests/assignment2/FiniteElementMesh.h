#pragma once

#include "AnimatedMesh.h"

#include "CGVectorWrapper.h"
#include "CGSystemWrapper.h"
#include "ConjugateGradient.h"

#include <Eigen/Dense>

template<class T, int d>
struct FiniteElementMesh : public AnimatedMesh<T, d>
{
    using Base = AnimatedMesh<T, d>;
    using Base::m_meshElements;
    using Base::m_particleX;
    using Vectord = typename Base::Vectord;
    using Matrixdd = Eigen::Matrix< T , d , d>;

    int m_nFrames;
    int m_subSteps;
    T m_frameDt;
    T m_stepDt;
    T m_stepEndTime;

    const T m_density;
    const T m_mu;
    const T m_lambda;
    const T m_rayleighCoefficient;

    std::vector<T> m_particleMass;
    std::vector<Vectord> m_particleV;
    std::vector<Matrixdd> m_DmInverse;
    std::vector<T> m_restVolume;
    
    FiniteElementMesh(const T density, const T mu, const T lambda, const T rayleighCoefficient)
        :m_density(density), m_mu(mu), m_lambda(lambda), m_rayleighCoefficient(rayleighCoefficient)
    {}

    void initializeUndeformedConfiguration()
    {
        // Initialize rest shape and particle mass (based on constant density)
        m_particleMass.resize(m_particleX.size(), T()); // Initialize all particle masses to zero
        for(const auto& element: m_meshElements)
        {
            Matrixdd Dm;
            for(int j = 0; j < d; j++)
                Dm.col(j) = m_particleX[element[j+1]]-m_particleX[element[0]];
            T restVolume = .5 * Dm.determinant();
            if(restVolume < 0)
                throw std::logic_error("Inverted element");
            m_DmInverse.emplace_back(Dm.inverse());
            m_restVolume.push_back(restVolume);
            T elementMass = m_density * restVolume;
            for(const int v: element)
                m_particleMass[v] += (1./(float)(d+1)) * elementMass;
        }
    }
    
    void addElasticForce(std::vector<Vectord>& f) const
    {
        for(int e = 0; e < m_meshElements.size(); e++)
        {
            const auto& element = m_meshElements[e];

            // Linear Elasticity
            Matrixdd Ds;
            for(int j = 0; j < d; j++)
                Ds.col(j) = m_particleX[element[j+1]]-m_particleX[element[0]];
            Matrixdd F = Ds * m_DmInverse[e];

            Matrixdd strain = .5 * (F + F.transpose()) - Matrixdd::Identity();
            Matrixdd P = 2. * m_mu * strain + m_lambda * strain.trace() * Matrixdd::Identity();

            Matrixdd H = -m_restVolume[e] * P * m_DmInverse[e].transpose();
            
            for(int j = 0; j < d; j++){
                f[element[j+1]] += H.col(j);
                f[element[0]] -= H.col(j);
            }
        }
    }
    //Multiply with K without actually constructing K (K=Stiffness matrix)
    void addProductWithStiffnessMatrix(std::vector<Vectord>& w, std::vector<Vectord>& f, const T scale) const
    {
        for(int e = 0; e < m_meshElements.size(); e++)
        {
            const auto& element = m_meshElements[e];

            // Linear Damping
            Matrixdd Ds_dot;
            for(int j = 0; j < d; j++)
                Ds_dot.col(j) = w[element[j+1]]-w[element[0]];
            Matrixdd F_dot = Ds_dot * m_DmInverse[e];

            Matrixdd strain_rate = .5 * (F_dot + F_dot.transpose());
            Matrixdd P_damping = scale * (2. * m_mu * strain_rate + m_lambda * strain_rate.trace() * Matrixdd::Identity());

            Matrixdd H_damping = m_restVolume[e] * P_damping * m_DmInverse[e].transpose();
            
            for(int j = 0; j < d; j++){
                f[element[j+1]] += H_damping.col(j);
                f[element[0]] -= H_damping.col(j);
            }
        }
    }
    
    // Computes LHS, RHS in the Linear Elasticity formula
    // and uses CG to solve for dx
    void simulateSubstep()
    {
        using FEMType = FiniteElementMesh<T, d>;

        const int nParticles = m_particleX.size();

        // Save last time step velocities

        std::vector<Vectord> lastV = m_particleV;
        
        // Construct initial guess for next-timestep
        //   Velocities -> Same as last timestep
        //   Positions -> Using Forward Euler
        
        for(int p = 0; p < nParticles; p++)
            m_particleX[p] += m_stepDt * m_particleV[p];

        // Overwrite boundary conditions with desired values

        setBoundaryConditions();
        
        // Solve for everything else using Conjugate Gradients

        std::vector<Vectord> dx(nParticles, Vectord::Zero());
        std::vector<Vectord> rhs(nParticles, Vectord::Zero());
        std::vector<Vectord> q(nParticles, Vectord::Zero());
        std::vector<Vectord> s(nParticles, Vectord::Zero());
        std::vector<Vectord> r(nParticles, Vectord::Zero());
        CGVectorWrapper<Vectord> dxWrapper(dx);
        CGVectorWrapper<Vectord> rhsWrapper(rhs);
        CGVectorWrapper<Vectord> qWrapper(q);
        CGVectorWrapper<Vectord> sWrapper(s);
        CGVectorWrapper<Vectord> rWrapper(r);
        CGSystemWrapper<Vectord, FEMType> systemWrapper(*this);
        
        addElasticForce(rhs);
        for(int p = 0; p < nParticles; p++)
            rhs[p] += (m_particleMass[p] / m_stepDt) * (lastV[p] - m_particleV[p]);
        addProductWithStiffnessMatrix(m_particleV, rhs, -m_rayleighCoefficient);
        clearConstrainedParticles(rhs);
        
        ConjugateGradient<T>::Solve(systemWrapper,
            dxWrapper, rhsWrapper, qWrapper, sWrapper, rWrapper,
            1e-4, 50);

        // Apply corrections to positions and velocities

        const T oneOverDt = T(1.) / m_stepDt;
        for(int p = 0; p < nParticles; p++){
            m_particleX[p] += dx[p];
            m_particleV[p] += oneOverDt * dx[p];
        }
    }

    void simulateFrame(const int frame)
    {
        m_stepDt = m_frameDt / (T) m_subSteps;

        for(int step = 1; step <= m_subSteps; step++){
            m_stepEndTime = m_frameDt * (T) (frame-1) + m_stepDt * (T) step;
            simulateSubstep();
        }
    }

    virtual void clearConstrainedParticles(std::vector<Vectord>& x) {}
    virtual void setBoundaryConditions() {}
};

