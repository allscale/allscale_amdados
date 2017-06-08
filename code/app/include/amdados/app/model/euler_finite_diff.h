//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {
namespace app {

//=================================================================================================
// Euler method for advection-diffusion equation discretized by finite difference approximation.
//=================================================================================================
template<size_t SizeX, size_t SizeY>
class EulerFiniteDifferenceModel : public IModel<SizeX, SizeY>
{
public:
    using base_model_t = IModel<SizeX, SizeY>;
    using base_model_t::PROBLEM_SIZE;
    using typename base_model_t::vector_t;
    using typename base_model_t::matrix_t;
    using sp_matrix_t = utils::SpMatrix<PROBLEM_SIZE, PROBLEM_SIZE>;
    using triplet_t = utils::Triplet;

public:
//-------------------------------------------------------------------------------------------------
// Constructor.
//-------------------------------------------------------------------------------------------------
IBM_NOINLINE
EulerFiniteDifferenceModel(const utils::Configuration & conf) : mConf(conf)
{
    mM.reset(new sp_matrix_t());
    mTriplets.reset(new triplet_arr_t());
}

//-------------------------------------------------------------------------------------------------
// Destructor.
//-------------------------------------------------------------------------------------------------
IBM_NOINLINE
virtual ~EulerFiniteDifferenceModel()
{
    mM.reset();             // this is not necessary but for debugging
    mTriplets.reset();
}

//-------------------------------------------------------------------------------------------------
//! @copydoc IModel::Update()
//-------------------------------------------------------------------------------------------------
IBM_NOINLINE
virtual void Update(double flow_x, double flow_y, double time_delta, double step_size,
                    vector_t & state, matrix_t & covar)
{
    MakeModelMatrix(flow_x, flow_y, time_delta, step_size);
    UpdateState(state);
    UpdateCovariance(covar);
}

private:
    using triplet_arr_t = std::vector<triplet_t>;

    const utils::Configuration &   mConf;       ///< reference to external parameter handler
    std::unique_ptr<sp_matrix_t>   mM;          ///< model matrix
    std::unique_ptr<triplet_arr_t> mTriplets;   ///< triplets for constructing the model matrix

private:
//-------------------------------------------------------------------------------------------------
// Function computes the model matrix, so that the transition from the field 'u_t' at time 't'
// to the field 'u_{t+1}' at time 't+1' (both fields are unrolled into vectors) can be written
// as follows:  u_{t+1} = mM * u_t.
//-------------------------------------------------------------------------------------------------
IBM_NOINLINE
void MakeModelMatrix(double flow_x, double flow_y, double time_delta, double step_size)
{
    using namespace amdados::app::utils;

    Sub2Ind<SizeX,SizeY> sub2ind;   // converts 2D index (x,y) to plain index

    mTriplets->clear();
    mTriplets->reserve(5 * PROBLEM_SIZE);

    const double c = mConf.asDouble("diffusion_coef");
    const double diffmult = (time_delta * c) / std::pow(step_size,2);
    const double advmult_x = (flow_x * time_delta) / (2*step_size);
    const double advmult_y = (flow_y * time_delta) / (2*step_size);
    const int NX = static_cast<int>(SizeX);
    const int NY = static_cast<int>(SizeY);

    // Finite difference approximation.
    //
    // u_{t+1}(x,y) = u_t(x,y)
    //      + (dt*c/ds^2) * (u_t(x+1,y) + u_t(x-1,y) + u_t(x,y+1) + u_t(x,y-1) - 4*u_t(x,y)
    //      - (dt*flow_x/(2*ds)) * (u_t(x+1,y) - u_t(x-1,y))
    //      - (dt*flow_y/(2*ds)) * (u_t(x,y+1) - u_t(x,y-1)).
    //
    for (int x = 1; x < NX-1; ++x) {
        for (int y = 1; y < NY-1; ++y) {
            int r = sub2ind(x,y);
            mTriplets->push_back(triplet_t(r, sub2ind(x+1,y), diffmult - advmult_x));
            mTriplets->push_back(triplet_t(r, sub2ind(x-1,y), diffmult + advmult_x));
            mTriplets->push_back(triplet_t(r, sub2ind(x,y+1), diffmult - advmult_y));
            mTriplets->push_back(triplet_t(r, sub2ind(x,y-1), diffmult + advmult_y));
            mTriplets->push_back(triplet_t(r, r, 1 - 4*diffmult));
        }
    }

    // Boundary: x == 0.
    for (int x = 0, y = 1; y < NY-1; ++y) {
        int r = sub2ind(x,y);
        mTriplets->push_back(triplet_t(r, sub2ind(x+1,y), diffmult - 2*advmult_x));
        mTriplets->push_back(triplet_t(r, sub2ind(x  ,y), diffmult + 2*advmult_x));
        mTriplets->push_back(triplet_t(r, sub2ind(x,y+1), diffmult - advmult_y));
        mTriplets->push_back(triplet_t(r, sub2ind(x,y-1), diffmult + advmult_y));
        mTriplets->push_back(triplet_t(r, r, 1 - 4*diffmult));
    }

    // Boundary: x == NX-1.
    for (int x = NX-1, y = 1; y < NY-1; ++y) {
        int r = sub2ind(x,y);
        mTriplets->push_back(triplet_t(r, sub2ind(x  ,y), diffmult - 2*advmult_x));
        mTriplets->push_back(triplet_t(r, sub2ind(x-1,y), diffmult + 2*advmult_x));
        mTriplets->push_back(triplet_t(r, sub2ind(x,y+1), diffmult - advmult_y));
        mTriplets->push_back(triplet_t(r, sub2ind(x,y-1), diffmult + advmult_y));
        mTriplets->push_back(triplet_t(r, r, 1 - 4*diffmult));
    }

    // Boundary: y == 0.
    for (int y = 0, x = 1; x < NX-1; ++x) {
        int r = sub2ind(x,y);
        mTriplets->push_back(triplet_t(r, sub2ind(x+1,y), diffmult - advmult_x));
        mTriplets->push_back(triplet_t(r, sub2ind(x-1,y), diffmult + advmult_x));
        mTriplets->push_back(triplet_t(r, sub2ind(x,y+1), diffmult - 2*advmult_y));
        mTriplets->push_back(triplet_t(r, sub2ind(x,y  ), diffmult + 2*advmult_y));
        mTriplets->push_back(triplet_t(r, r, 1 - 4*diffmult));
    }

    // Boundary: y == NY-1.
    for (int y = NY-1, x = 1; x < NX-1; ++x) {
        int r = sub2ind(x,y);
        mTriplets->push_back(triplet_t(r, sub2ind(x+1,y), diffmult - advmult_x));
        mTriplets->push_back(triplet_t(r, sub2ind(x-1,y), diffmult + advmult_x));
        mTriplets->push_back(triplet_t(r, sub2ind(x,y  ), diffmult - 2*advmult_y));
        mTriplets->push_back(triplet_t(r, sub2ind(x,y-1), diffmult + 2*advmult_y));
        mTriplets->push_back(triplet_t(r, r, 1 - 4*diffmult));
    }

    // At the corner points.
    const int corner_x[4] = {0, 0, NX-1, NX-1};
    const int corner_y[4] = {0, NY-1, NY-1, 0};
    for (int k = 0; k < 4; ++k) {
        int x = corner_x[k];
        int y = corner_y[k];
        int r = sub2ind(x,y);
        mTriplets->push_back(triplet_t(r, sub2ind(std::min(x+1,NX-1),y), diffmult - 2*advmult_x));
        mTriplets->push_back(triplet_t(r, sub2ind(std::max(x-1,   0),y), diffmult + 2*advmult_x));
        mTriplets->push_back(triplet_t(r, sub2ind(x,std::min(y+1,NY-1)), diffmult - 2*advmult_y));
        mTriplets->push_back(triplet_t(r, sub2ind(x,std::max(y-1,   0)), diffmult + 2*advmult_y));
        mTriplets->push_back(triplet_t(r, r, 1 - 4*diffmult));
    }

    // Assemble the model matrix.
    assert(mTriplets->capacity() == 5 * PROBLEM_SIZE);      // no reallocation had happened
    mM->SetFromTriplets(*mTriplets, MemoryPolicy::RETAIN_MEMORY);
}

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
void UpdateState(utils::Vector<PROBLEM_SIZE> & state)
{
    std::unique_ptr<vector_t> old_state(new vector_t());    // TODO: pool of free vectors
    *old_state = state;
    SparseMulVector(state, *mM, *old_state);                // state = M * old_state
}

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
void UpdateCovariance(utils::Matrix<PROBLEM_SIZE, PROBLEM_SIZE> & covar)
{
    std::unique_ptr<matrix_t> old_covar(new matrix_t());    // TODO: pool of free matrices
    DenseMulSparseTr(*old_covar, covar, *mM);
    SparseMulDense(covar, *mM, *old_covar);                 // covar = M * old_covar * M^T
}

}; // class EulerFiniteDifferenceModel

} // namespace app
} // namespace amdados

