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
#include "coagulation/tensor_train.h"

namespace solvers
{
    class TKernel: public TMatrix
    {
    public:
        TKernel (const int &m, const int &n) : TMatrix(m , n) {}
        double value (const int &i, const int &j) override {return 5.;}
    };

    TCross_Parallel_v1 default_crossed_kernel(const double & tolerance, 
                                              const int & size);

    template <typename T>
    class Coagulation
    {
    public:
        //T *L1_res_vec, *L2_res_vec;
        fftw_complex *ub, *vb;
        fftw_plan *plan_v, *plan_u, *plan_i;
        int R_value;
        int M;
        int N;
        
    
        T L2(const int &N, const int &i, const T *n);
        T L1(const int &N, const int &i, const T *n);

        T dt;
        int size;
        TCross_Parallel_v1& kernel;

        Coagulation(int size, T dt, TCross_Parallel_v1 & kernel) :
                    size(size), dt(dt), kernel(kernel)
        {
            R_value = kernel.get_rank();
            M = kernel.get_rows_number();
            N = kernel.get_columns_number();
            //L2_res_vec = new T [M];
            double V_value = kernel.get_columns_number();
            ub = (fftw_complex *) fftw_malloc(R_value * V_value * 
                                               sizeof(fftw_complex));
            vb = (fftw_complex *) fftw_malloc(R_value * V_value * 
                                               sizeof(fftw_complex));
            plan_v = (fftw_plan *) fftw_malloc(R_value * 
                                               sizeof(fftw_plan));
            plan_u = (fftw_plan *) fftw_malloc(R_value * 
                                               sizeof(fftw_plan));
            plan_i = (fftw_plan *) fftw_malloc(R_value * 
                                               sizeof(fftw_plan));
            for (int i = 0; i < R_value; i++)
            {
                plan_v[i] = fftw_plan_dft_1d(size, vb+i*size, vb+i*size, 
                                           FFTW_FORWARD, FFTW_ESTIMATE);
                plan_u[i] = fftw_plan_dft_1d(size, ub+i*size, ub+i*size, 
                                           FFTW_FORWARD, FFTW_ESTIMATE);
                plan_i[i] = fftw_plan_dft_1d(size, ub+i*size, ub+i*size,
                                          FFTW_BACKWARD, FFTW_ESTIMATE);	
            }

        };
        ~Coagulation();
        void iteration_direct(T * data, T * data_new, int coef);
        void iteration_low_rank(T * data, T * data_new, int coef);
    };
}
