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
    using sp_matrix_t = utils::SpMatrix<PROBLEM_SIZE, PROBLEM_SIZE>;

    // Destructor.
    virtual ~IModel() {}

    // Function initializes and returns the model matrix, then updates the state and covariance
    // one time step ahead using the current ones and the model matrix. If 'M' is the model
    // matrix, then:  state <- M * state;  covar <- M * covar * M^T.
    // \param  flow_x       abscissa of the flow vector.
    // \param  flow_y       ordinate of the flow vector.
    // \param  time_delta   time integration step.
    // \param  space_delta  space discretization step.
    // \param  t            current time.
    virtual const sp_matrix_t & ModelMatrix(double flow_x, double flow_y,
                                            double time_delta, double space_delta, double t) = 0;
};

#if 0
//-------------------------------------------------------------------------------------------------
// Function updates the state and covariance: state <- M * state; covar <- M * covar * M^T.
// The function is slow (but safe) because temporary objects are allocated on the heap.
//-------------------------------------------------------------------------------------------------
template<size_t MSIZE>
IBM_NOINLINE
void UpdateStateSlow(const utils::Matrix<MSIZE,MSIZE> & M,
                           utils::Vector<MSIZE>       & state,
                           utils::Matrix<MSIZE,MSIZE> & covar)
{
    using vector_t = utils::Vector<MSIZE>;
    using matrix_t = utils::Matrix<MSIZE, MSIZE>;

    std::unique_ptr<vector_t> tmp_state(new vector_t());
    *tmp_state = state;
    utils::MatVecMult(state, M, *tmp_state);                // state <- M * state

    std::unique_ptr<matrix_t> tmp_covar(new matrix_t());
    *tmp_covar = covar;
    utils::MatMultTr(*tmp_covar, covar, M);
    utils::MatMult(covar, M, *tmp_covar);                   // covar <- M * old_covar * M^T
}

//-------------------------------------------------------------------------------------------------
// Function updates the state and covariance: state <- M * state; covar <- M * covar * M^T.
// The function is slow (but safe) because temporary objects are allocated on the heap.
//-------------------------------------------------------------------------------------------------
template<size_t MSIZE>
IBM_NOINLINE
void UpdateStateSlow(const utils::SpMatrix<MSIZE,MSIZE> & M,
                           utils::Vector<MSIZE>         & state,
                           utils::Matrix<MSIZE,MSIZE>   & covar)
{
    using vector_t = utils::Vector<MSIZE>;
    using matrix_t = utils::Matrix<MSIZE, MSIZE>;

    std::unique_ptr<vector_t> tmp_state(new vector_t());
    *tmp_state = state;
    utils::MatVecMult(state, M, *tmp_state);                // state <- M * state

    std::unique_ptr<matrix_t> tmp_covar(new matrix_t());
    *tmp_covar = covar;
    utils::MatMultTr(*tmp_covar, covar, M);
    utils::MatMult(covar, M, *tmp_covar);                   // covar <- M * old_covar * M^T
}
#endif

//-------------------------------------------------------------------------------------------------
// Function updates the state and covariance:
//      state <- M * state
//      covar <- M * covar * M^T
// The function is slow (but safe) because temporary objects are allocated on the heap
// avoiding potential stack overflow.
//-------------------------------------------------------------------------------------------------
template<typename MODEL, typename STATE, typename COVAR>
IBM_NOINLINE
void UpdateState(const MODEL & M, STATE & state, COVAR & covar)
{
#if 1
    std::unique_ptr<STATE> tmp_state(new STATE(state));
    /**tmp_state = state;*/
    utils::MatVecMult(state, M, *tmp_state);                // state <- M * state

    std::unique_ptr<COVAR> tmp_covar(new COVAR(covar));
    /**tmp_covar = covar;*/
    utils::MatMultTr(*tmp_covar, covar, M);
    utils::MatMult(covar, M, *tmp_covar);                   // covar <- M * old_covar * M^T
#else   // this is fast:
    STATE tmp_state{state};
    utils::MatVecMult(state, M, tmp_state);                 // state <- M * state
    COVAR tmp_covar{covar};
    utils::MatMultTr(tmp_covar, covar, M);
    utils::MatMult(covar, M, tmp_covar);                    // covar <- M * old_covar * M^T
#endif
}

} // namespace app
} // namespace amdados

