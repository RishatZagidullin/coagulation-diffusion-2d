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
        T h;

        T K(const int & u, const int &v, const T h) {
            //double u1=(u+1.0) * h;
            //double v1=(v+1.0) * h;
            //double c = 1.;
            double result = 2.0;//(c/u1+c/v1)*(u1+v1);
            return result;
        }

        Coagulation(int size, T h, T dt);
        void iteration(T * data, T * data_new, int coef);
    private:
    };
}
