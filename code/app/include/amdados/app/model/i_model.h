//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {
namespace app {

//=================================================================================================
// Interface to any model that propagates state according to a 2D physical model in hand
// (e.g. advection-diffusion equation).
//=================================================================================================
template<size_t SizeX, size_t SizeY>
class IModel
{
public:
    static const size_t PROBLEM_SIZE = SizeX * SizeY;
    using vector_t = utils::Vector<PROBLEM_SIZE>;
    using matrix_t = utils::Matrix<PROBLEM_SIZE, PROBLEM_SIZE>;

    // Destructor.
    virtual ~IModel() {}

    // Function initializes the model matrix, then updates the state and covariance
    // one time step ahead using the current ones and the model matrix. If 'M' is the model
    // matrix, then:  state <- M * state;  covar <- M * covar * M^T.
    // \param  flow_x      abscissa of the flow vector.
    // \param  flow_y      ordinate of the flow vector.
    // \param  time_delta  time integration step.
    // \param  step_size   space discretization step.
    // \param  state       in: current state vector; out: estimated new state vector.
    // \param  covar       in: current covariance; out: estimated new covariance.
    virtual void Update(double flow_x, double flow_y, double time_delta, double step_size,
                        vector_t & state, matrix_t & covar) = 0;
};

} // namespace app
} // namespace amdados

