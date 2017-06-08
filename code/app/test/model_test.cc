//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#include <gtest/gtest.h>
#include "allscale/utils/assert.h"
#include "amdados/app/amdados_grid.h"
#include "amdados/app/utils/common.h"
#include "amdados/app/utils/amdados_utils.h"
#include "amdados/app/utils/matrix.h"
#include "amdados/app/utils/sparse_matrix.h"
#include "amdados/app/utils/configuration.h"
#include "amdados/app/model/i_model.h"
#include "amdados/app/model/euler_finite_diff.h"

using namespace amdados::app::utils;

static Configuration gConf;         // a single global configuration

// Tolerance on relative error.
const double TOL = std::sqrt(std::numeric_limits<double>::epsilon());

//-------------------------------------------------------------------------------------------------
// Straightforward state propagation using finite difference discretization of advection-diffusion
// equation and Euler time integration method.
//-------------------------------------------------------------------------------------------------
template<size_t SizeX, size_t SizeY>
void EulerFiniteDifferenceDirect(Matrix<SizeX,SizeY> & state,
                                 double flow_x, double flow_y, double dt, double ds)
{
    const double c = gConf.asDouble("diffusion_coef");
    const double diffmult = (dt * c) / std::pow(ds,2);
    const double advmult_x = (flow_x * dt) / (2*ds);
    const double advmult_y = (flow_y * dt) / (2*ds);
    const int NX = static_cast<int>(SizeX);
    const int NY = static_cast<int>(SizeY);

    std::unique_ptr< Matrix<SizeX,SizeY> > new_state(new Matrix<SizeX,SizeY>());
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            double sum = state[{i,j}];

            int ip = std::min(i+1,NX-1);    int jp = std::min(j+1,NY-1);
            int im = std::max(i-1,0);       int jm = std::max(j-1,0);

            // Diffusion term.
            sum += (state[{ip,j}] +
                    state[{im,j}] +
                    state[{i,jp}] +
                    state[{i,jm}] - 4 * state[{i,j}]) * diffmult;

            double coef_x = (((0 < i) && (i+1 < NX)) ? 1 : 2) * advmult_x;
            double coef_y = (((0 < j) && (j+1 < NY)) ? 1 : 2) * advmult_y;

            // Advection term.
            sum += -coef_x * (state[{ip,j}] - state[{im,j}]);
            sum += -coef_y * (state[{i,jp}] - state[{i,jm}]);

            (*new_state)[{i,j}] = sum;
        }
    }
    state = *new_state;
}

//-------------------------------------------------------------------------------------------------
// Test straightforward state propagation against matrix based approach implemented in
// EulerFiniteDifferenceModel<..> class.
//-------------------------------------------------------------------------------------------------
template<size_t SizeX, size_t SizeY>
void TestEulerFiniteDifference(double & max_rel_diff)
{
    std::mt19937                           gen(std::time(nullptr));
    std::uniform_real_distribution<double> distrib(0.0, 2.0 * M_PI);
    const double                           angle = distrib(gen);

    const double flow_x = std::cos(angle);
    const double flow_y = std::sin(angle);
    const double time_delta = 0.1;
    const double step_size = 0.1;

    using model_t = amdados::app::EulerFiniteDifferenceModel<SizeX, SizeY>;
    using vector_t = typename model_t::vector_t;
    using matrix_t = typename model_t::matrix_t;
    using field_t = Matrix<SizeX, SizeY>;

    model_t                   model(gConf);
    std::unique_ptr<matrix_t> covar(new matrix_t());
    std::unique_ptr<vector_t> state(new vector_t()), state2(new vector_t());

    MakeIdentityMatrix(*covar);
    MakeRandomVector(*state);
    *state2 = *state;

    // Make several "iterations".
    for (int iterNo = 0; iterNo < 3; ++iterNo) {
        // Time integration by means of model matrix.
        model.Update(flow_x, flow_y, time_delta, step_size, *state, *covar);

        // Direct time integration on the grid.
        {
            std::unique_ptr<field_t> state_field(new field_t());
            Reshape1Dto2D<SizeX,SizeY>(*state_field, *state2);
            EulerFiniteDifferenceDirect(*state_field, flow_x, flow_y, time_delta, step_size);
            Reshape2Dto1D<SizeX,SizeY>(*state2, *state_field);
        }
    }

    max_rel_diff = std::max(max_rel_diff, NormVecDiff(*state, *state2) / NormVec(*state));
    EXPECT_LT(max_rel_diff, TOL) << "Relative error of exceeded tolerance";
}

//-------------------------------------------------------------------------------------------------
// Function tests the matrix library written by means of Allscale API.
//-------------------------------------------------------------------------------------------------
TEST(ModelTests, Basic)
{
    // Read configuration settings.
    gConf.ReadConfigFile("../../amdados_unittest.conf");
    gConf.PrintParameters();
    MakeDirectory(gConf.asCString("output_dir"));

    std::string filename = gConf.asString("output_dir") + "/model_test.log";
    std::fstream fout(filename, std::ios::out | std::ios::trunc);
    assert_true(fout.good()) << "failed to oped the summary file: " << filename << std::endl;

    {
        double max_rel_err = 0.0;
        for (int testNo = 0; testNo < 3; ++testNo) {
            TestEulerFiniteDifference<37,67>(max_rel_err);
            TestEulerFiniteDifference<43,53>(max_rel_err);
            TestEulerFiniteDifference<73,11>(max_rel_err);
            TestEulerFiniteDifference<57,67>(max_rel_err);
            TestEulerFiniteDifference<89,97>(max_rel_err);
            TestEulerFiniteDifference<97,17>(max_rel_err);
            TestEulerFiniteDifference<37,67>(max_rel_err);
            TestEulerFiniteDifference<33,41>(max_rel_err);
            TestEulerFiniteDifference<89,47>(max_rel_err);
            TestEulerFiniteDifference<97,41>(max_rel_err);
        }
        fout << "TestEulerFiniteDifference(): max. relative error: " << max_rel_err << std::endl;
    }
    fout << std::endl << std::flush;
}

