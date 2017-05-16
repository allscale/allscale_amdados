#include <gtest/gtest.h>

#include "amdados/app/utils/cholesky.h"

using namespace amdados::app::utils;

//-------------------------------------------------------------------------------------------------
// Function tests the Cholesky decomposition and the linear systems solvers given problem size.
//-------------------------------------------------------------------------------------------------
template<int PROBLEM_SIZE>
void TestCholeskyGivenProblemSize(int numTests)
{
	std::random_device               rd;
	std::mt19937                     gen(rd());
	std::uniform_real_distribution<> distrib(0.0, 1.0);

	const int NCOLS = PROBLEM_SIZE % 37 + 7;
	using matrix_t = allscale::utils::grid<double, PROBLEM_SIZE, PROBLEM_SIZE>;
	using matrix_other_t = allscale::utils::grid<double, PROBLEM_SIZE, NCOLS>;
	using vector_t = allscale::utils::grid<double, PROBLEM_SIZE>;

	matrix_t       A, tmpA;
	matrix_other_t B, X;
	vector_t       b, x;

	//std::cout << "Testing Cholesky decomposition, problem size: " << PROBLEM_SIZE
	//	<< std::endl << std::endl;
	for(int testNo = 0; testNo < numTests; ++testNo) {
		double norm_b = 0.0, norm_B = 0.0, diff_norm = 0.0;

		// Create a random, square, symmetric, positive-definite matrix A: A = tmpA * tmpA.
		for(int i = 0; i < PROBLEM_SIZE; ++i) {
			for(int j = 0; j < PROBLEM_SIZE; ++j) {
				tmpA[{i, j}] = distrib(gen);
			}
		}
		Symmetrize(tmpA);
		for(int r = 0; r < PROBLEM_SIZE; ++r) {
			for(int c = 0; c < PROBLEM_SIZE; ++c) {
				double v = 0.0;
				for(int i = 0; i < PROBLEM_SIZE; ++i) {
					v += tmpA[{r, i}] * tmpA[{i, c}];
				}
				A[{r, c}] = v;
			}
		}
		Symmetrize(A);

		// Create the right-hand side matrix and vector.
		for(int r = 0; r < PROBLEM_SIZE; ++r) {
			for(int c = 0; c < NCOLS; ++c) {
				double v = distrib(gen);
				B[{r, c}] = v;
				norm_B += v * v;
			}
			double v = distrib(gen);
			b[{r}] = v;
			norm_b += v * v;
		}
		norm_b = std::sqrt(norm_b);
		norm_B = std::sqrt(norm_B);

		// Solve linear systems.
		Cholesky<PROBLEM_SIZE> chol;
		chol.ComputeDecomposition(A);
		chol.Solve(x, b);
		chol.BatchSolve(X, B);

		// Compute |A*x - b| and print the relative error.
		diff_norm = 0.0;
		for(int r = 0; r < PROBLEM_SIZE; ++r) {
			double v = 0.0;
			for(int i = 0; i < PROBLEM_SIZE; ++i) {
				v += A[{r, i}] * x[{i}];
			}
			v -= b[{r}];
			diff_norm += v * v;
		}
		diff_norm = std::sqrt(diff_norm);
		EXPECT_LT(diff_norm / norm_b, 1.0e-5) << "Relative error |A*x - b|/|b| exceeded tolerance";
		//std::cout << "Relative error |A*x - b|/|b|: " << (diff_norm / norm_b) << std::endl;

		// Compute |A*X - B| and print the relative error.
		diff_norm = 0.0;
		for(int r = 0; r < PROBLEM_SIZE; ++r) {
			for(int c = 0; c < NCOLS; ++c) {
				double v = 0.0;
				for(int i = 0; i < PROBLEM_SIZE; ++i) {
					v += A[{r, i}] * X[{i, c}];
				}
				v -= B[{r, c}];
				diff_norm += v * v;
			}
		}
		diff_norm = std::sqrt(diff_norm);
		EXPECT_LT(diff_norm / norm_b, 1.0e-5) << "Relative error |A*x - b|/|b| exceeded tolerance";
		//std::cout << "Relative error |A*X - B|/|B|: " << (diff_norm / norm_B) << std::endl;
	}
	//std::cout << std::endl;
}

//-------------------------------------------------------------------------------------------------
// Function tests the Cholesky decomposition and the linear systems solvers
// for different matrix sizes.
//-------------------------------------------------------------------------------------------------
TEST(Cholesky, Basic) {
	
	const int numTests = 10;
	TestCholeskyGivenProblemSize<11>(numTests);
	TestCholeskyGivenProblemSize<17>(numTests);
	TestCholeskyGivenProblemSize<23>(numTests);
	TestCholeskyGivenProblemSize<37>(numTests);
	TestCholeskyGivenProblemSize<41>(numTests);
	TestCholeskyGivenProblemSize<53>(numTests);
	TestCholeskyGivenProblemSize<67>(numTests);
	TestCholeskyGivenProblemSize<73>(numTests);
	TestCholeskyGivenProblemSize<89>(numTests);
	TestCholeskyGivenProblemSize<97>(numTests);

}
