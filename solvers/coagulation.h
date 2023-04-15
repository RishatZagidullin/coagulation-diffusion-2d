#pragma once
#include <stdio.h>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
namespace solvers
{

    template <typename T>
    class Coagulation
    {
    public:
        T L2(const int &N, const int &i, const T *n);
        T L1(const int &N, const int &i, const T *n);

        T dt;
        int size;
        T kernel;

        T K(const int & u, const int &v) {
            return kernel;
        }

        Coagulation(int size, T dt, T kernel);
        void iteration(T * data, T * data_new, int coef);
    };
}
