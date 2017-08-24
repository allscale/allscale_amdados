#include <gtest/gtest.h>
#include "amdados/app/amdados_grid.h"
#include "allscale/utils/assert.h"
#include "amdados/app/utils/common.h"
#include "amdados/app/utils/amdados_utils.h"
#include "amdados/app/utils/matrix.h"

using namespace amdados::app;
using namespace amdados::app::utils;

template<int SizeX, int SizeY>
void Test1()
{
    Matrix<SizeX, SizeY> covar;
    covar(0,0) = 0;
    std::cout << covar(0,0) << std::endl;
    std::cout << "Test1 succeeded" << std::endl;
}

template<int SizeX, int SizeY>
void Test2()
{
    const int PROBLEM_SIZE {SizeX * SizeY};
    auto covar = new Matrix<PROBLEM_SIZE, PROBLEM_SIZE>();
    //(*covar)(-1,0) = 0;
    std::cout << (*covar)(0,0) << std::endl;
    delete covar;
    std::cout << "Test2 succeeded" << std::endl;
}

template<int SizeX, int SizeY>
void Test3()
{
    const int PROBLEM_SIZE {SizeX * SizeY};
    Matrix<PROBLEM_SIZE, PROBLEM_SIZE> covar;
    std::cout << covar(0,0) << std::endl;
    std::cout << "Test3 succeeded" << std::endl;
}

TEST(FailTests, Basic)
{
    Test1<37,67>();     // fine
    Test2<37,67>();     // fine
    //Test3<37,67>();     // segfault: likely stack overflow!
}

