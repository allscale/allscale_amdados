//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {
namespace app {

//=================================================================================================
// Interface to any solver, e.g. for advection-diffusion equation.
//=================================================================================================
class ISolver
{
public:
    // Destructor.
    virtual ~ISolver() {}

    // Function initializes the solver given configuration parameters.
    virtual void InitSolver(const utils::Configuration & conf) = 0;

    // Function runs the solver given configuration parameters.
    virtual void RunSolver(const utils::Configuration & conf) = 0;
};

// Function creates an instance of a solver using Euler time integration (EI),
// finite-difference discretization (FD) and Kalman Filter (KF) for data assimilation.
std::unique_ptr<ISolver> CreateSolver_EI_FD_KF(const utils::Configuration & conf);

} // namespace app
} // namespace amdados

