#include <gtest/gtest.h>
#include "amdados/app/amdados_grid.h"

template<size_t SizeX, size_t SizeY>
void Test1()
{
    allscale::utils::grid<double, SizeX, SizeY> covar;
    covar[{0,0}] = 0;
    std::cout << covar[{0,0}] << std::endl;
    std::cout << "Test1 succeeded" << std::endl;
}

template<size_t SizeX, size_t SizeY>
void Test2()
{
    const size_t PROBLEM_SIZE {SizeX * SizeY};
    auto covar = new allscale::utils::grid<double, PROBLEM_SIZE, PROBLEM_SIZE>();
    std::cout << (*covar)[{0,0}] << std::endl;
    delete covar;
    std::cout << "Test2 succeeded" << std::endl;
}

template<size_t SizeX, size_t SizeY>
void Test3()
{
    const size_t PROBLEM_SIZE {SizeX * SizeY};
    allscale::utils::grid<double, PROBLEM_SIZE, PROBLEM_SIZE> covar;
    std::cout << covar[{0,0}] << std::endl;
    std::cout << "Test3 succeeded" << std::endl;
}

TEST(FailTests, Basic)
{
    Test1<37,67>();     // fine
    Test2<37,67>();     // fine
    //Test3<37,67>();     // segfault: likely stack overflow!
}

