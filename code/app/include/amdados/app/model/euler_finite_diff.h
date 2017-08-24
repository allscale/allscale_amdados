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
template<int SizeX, int SizeY>
class EulerFiniteDifferenceModel : public IModel<SizeX, SizeY>
{
public:
    using base_model_t = IModel<SizeX, SizeY>;
    using base_model_t::PROBLEM_SIZE;
    using typename base_model_t::sp_matrix_t;
    using triplet_t = utils::Triplet;

public:
//-------------------------------------------------------------------------------------------------
// Constructor.
//-------------------------------------------------------------------------------------------------
IBM_NOINLINE
EulerFiniteDifferenceModel(const utils::Configuration & conf) : mConf(conf), mM(), mTriplets()
{
}

//-------------------------------------------------------------------------------------------------
// Destructor.
//-------------------------------------------------------------------------------------------------
IBM_NOINLINE
virtual ~EulerFiniteDifferenceModel()
{
}

//-------------------------------------------------------------------------------------------------
//! @copydoc IModel::Update()
//-------------------------------------------------------------------------------------------------
IBM_NOINLINE
virtual const sp_matrix_t & ModelMatrix(double flow_x, double flow_y,
                                        double time_delta, double space_delta, double t)
{
    (void)t;
    MakeModelMatrix(flow_x, flow_y, time_delta, space_delta);
    return mM;
}

private:
    using triplet_arr_t = std::vector<triplet_t>;

    const utils::Configuration & mConf;     ///< reference to external parameter handler
    sp_matrix_t                  mM;        ///< model matrix
    triplet_arr_t                mTriplets; ///< triplets for constructing the model matrix

private:
//-------------------------------------------------------------------------------------------------
// Function computes the model matrix, so that the transition from the field 'u_t' at time 't'
// to the field 'u_{t+1}' at time 't+1' (both fields are unrolled into vectors) can be written
// as follows:  u_{t+1} = mM * u_t.
//-------------------------------------------------------------------------------------------------
IBM_NOINLINE
void MakeModelMatrix(double flow_x, double flow_y, double time_delta, double space_delta)
{
    using namespace amdados::app::utils;

    Sub2Ind<SizeX,SizeY> sub2ind;   // converts 2D index (x,y) to plain index

    mTriplets.resize(5 * PROBLEM_SIZE);
    triplet_t * pTriplets = mTriplets.data();
    int         count = 0;

    const double c = mConf.asDouble("diffusion_coef");
    const double diffmult = (time_delta * c) / std::pow(space_delta,2);
    const double advmult_x = (flow_x * time_delta) / (2*space_delta);
    const double advmult_y = (flow_y * time_delta) / (2*space_delta);
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
            assert_true(count + 5 <= 5 * PROBLEM_SIZE);
            int r = sub2ind(x,y);
            pTriplets[count++] = triplet_t(r, sub2ind(x+1,y), diffmult - advmult_x);
            pTriplets[count++] = triplet_t(r, sub2ind(x-1,y), diffmult + advmult_x);
            pTriplets[count++] = triplet_t(r, sub2ind(x,y+1), diffmult - advmult_y);
            pTriplets[count++] = triplet_t(r, sub2ind(x,y-1), diffmult + advmult_y);
            pTriplets[count++] = triplet_t(r, r, 1 - 4*diffmult);
        }
    }

    // Boundary: x == 0.
    for (int x = 0, y = 1; y < NY-1; ++y) {
        assert_true(count + 5 <= 5 * PROBLEM_SIZE);
        int r = sub2ind(x,y);
        pTriplets[count++] = triplet_t(r, sub2ind(x+1,y), diffmult - 2*advmult_x);
        pTriplets[count++] = triplet_t(r, sub2ind(x  ,y), diffmult + 2*advmult_x);
        pTriplets[count++] = triplet_t(r, sub2ind(x,y+1), diffmult - advmult_y);
        pTriplets[count++] = triplet_t(r, sub2ind(x,y-1), diffmult + advmult_y);
        pTriplets[count++] = triplet_t(r, r, 1 - 4*diffmult);
    }

    // Boundary: x == NX-1.
    for (int x = NX-1, y = 1; y < NY-1; ++y) {
        assert_true(count + 5 <= 5 * PROBLEM_SIZE);
        int r = sub2ind(x,y);
        pTriplets[count++] = triplet_t(r, sub2ind(x  ,y), diffmult - 2*advmult_x);
        pTriplets[count++] = triplet_t(r, sub2ind(x-1,y), diffmult + 2*advmult_x);
        pTriplets[count++] = triplet_t(r, sub2ind(x,y+1), diffmult - advmult_y);
        pTriplets[count++] = triplet_t(r, sub2ind(x,y-1), diffmult + advmult_y);
        pTriplets[count++] = triplet_t(r, r, 1 - 4*diffmult);
    }

    // Boundary: y == 0.
    for (int y = 0, x = 1; x < NX-1; ++x) {
        assert_true(count + 5 <= 5 * PROBLEM_SIZE);
        int r = sub2ind(x,y);
        pTriplets[count++] = triplet_t(r, sub2ind(x+1,y), diffmult - advmult_x);
        pTriplets[count++] = triplet_t(r, sub2ind(x-1,y), diffmult + advmult_x);
        pTriplets[count++] = triplet_t(r, sub2ind(x,y+1), diffmult - 2*advmult_y);
        pTriplets[count++] = triplet_t(r, sub2ind(x,y  ), diffmult + 2*advmult_y);
        pTriplets[count++] = triplet_t(r, r, 1 - 4*diffmult);
    }

    // Boundary: y == NY-1.
    for (int y = NY-1, x = 1; x < NX-1; ++x) {
        assert_true(count + 5 <= 5 * PROBLEM_SIZE);
        int r = sub2ind(x,y);
        pTriplets[count++] = triplet_t(r, sub2ind(x+1,y), diffmult - advmult_x);
        pTriplets[count++] = triplet_t(r, sub2ind(x-1,y), diffmult + advmult_x);
        pTriplets[count++] = triplet_t(r, sub2ind(x,y  ), diffmult - 2*advmult_y);
        pTriplets[count++] = triplet_t(r, sub2ind(x,y-1), diffmult + 2*advmult_y);
        pTriplets[count++] = triplet_t(r, r, 1 - 4*diffmult);
    }

    // At the corner points.
    const int corner_x[4] = {0, 0, NX-1, NX-1};
    const int corner_y[4] = {0, NY-1, NY-1, 0};
    for (int k = 0; k < 4; ++k) {
        assert_true(count + 5 <= 5 * PROBLEM_SIZE);
        int x = corner_x[k];
        int y = corner_y[k];
        int r = sub2ind(x,y);
        pTriplets[count++] = triplet_t(r, sub2ind(std::min(x+1,NX-1),y), diffmult - 2*advmult_x);
        pTriplets[count++] = triplet_t(r, sub2ind(std::max(x-1,   0),y), diffmult + 2*advmult_x);
        pTriplets[count++] = triplet_t(r, sub2ind(x,std::min(y+1,NY-1)), diffmult - 2*advmult_y);
        pTriplets[count++] = triplet_t(r, sub2ind(x,std::max(y-1,   0)), diffmult + 2*advmult_y);
        pTriplets[count++] = triplet_t(r, r, 1 - 4*diffmult);
    }

    // Assemble the model matrix.
    assert_true(count == 5 * PROBLEM_SIZE);
    mM.SetFromTriplets(mTriplets, MemoryPolicy::RETAIN_MEMORY);
}

}; // class EulerFiniteDifferenceModel

} // namespace app
} // namespace amdados

