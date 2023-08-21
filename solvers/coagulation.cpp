#include "coagulation.h"
#include <cmath>

namespace solvers
{
    
    TCross_Parallel_v1 default_crossed_kernel(const double & tolerance, const int & size)
    {
        TCross_Parallel_v1_Parameters parameters;
        parameters.tolerance = tolerance;
        parameters.maximal_iterations_number = 0;
        TKernel kernel(size, size);
        TCross_Parallel_v1 crossed_kernel;
        crossed_kernel.Approximate(&kernel, parameters);
        return crossed_kernel;
    }

    template<typename T>
    void Coagulation<T>::iteration_low_rank(T * data, T * data_new,
                                            int coef)
    {
        T * n_k = new T [size];
        T * L1_res_vec;// = new T [M];
        T * L2_res_vec = new T [M];
        for (int m = 0; m < size; m++)
            n_k[m] = data[coef*m];

        kernel.matvec(n_k, L2_res_vec);
        L1_res_vec = kernel.smol_conv_discrete(n_k, ub, vb, plan_v, 
                                               plan_u, plan_i);
        //if (L2_res_vec[1] > 0)
        //    std::cout << L1_res_vec[1] << " " << L2_res_vec[1] << "\n"; 
        #pragma omp parallel for
        for (int m = 0; m < size; m++)
        {
            int ind = coef*m;
            data_new[ind] = ( L1_res_vec[m] * 0.5 
                             - n_k[m] * L2_res_vec[m] )
                             * dt + n_k[m];
            if (data_new[ind] < 0.0) data_new[ind] = 0.0;
        }
        delete [] L2_res_vec;
        free(L1_res_vec);
        delete [] n_k;
    }

    template<typename T>
    void Coagulation<T>::iteration_direct(T * data, T * data_new,
                                          int coef)
    {
        T * n_k = new T [size];
        for (int m = 0; m < size; m++)
            n_k[m] = data[coef*m];

        //#pragma omp parallel for
        for (int m = 0; m < size; m++)
        {
            int ind = coef*m;
            data_new[ind] = ( L1(size, m, n_k) * 0.5 
                             - n_k[m] * L2(size, m, n_k) )
                             * dt + n_k[m];
             if (data_new[ind] < 0.0) data_new[ind] = 0.0;
        }
        //if (L2(size, 1, n_k)> 0)
        //    std::cout << L1(size, 1, n_k) << " " << L2(size, 1, n_k) << "\n"; 
        delete [] n_k;
    }

    template<typename T>
    T Coagulation<T>::L1(const int &N, const int &i, const T *n)
    {
        T l1 = 0;
        for (int i1 = 0; i1 < i; i1++)
            l1 += n[i1] * n[i-i1-1] * kernel.value(i - i1 - 1, i1);
        return l1;
    }

    template<typename T>
    T Coagulation<T>::L2(const int &N, const int &i, const T *n)
    {
        T l2 = 0;
        for (int i1 = 0; i1 < N; i1++)
            l2 += n[i1] * kernel.value(i, i1);
        return l2;
    }
    
    template<typename T>
    Coagulation<T>::~Coagulation()
    {
        //delete [] L2_res_vec;
        for (int i = 0; i < R_value; i++)
        {
            fftw_destroy_plan(plan_v[i]);
            fftw_destroy_plan(plan_u[i]);
            fftw_destroy_plan(plan_i[i]);	
        }
        fftw_free(vb);
        fftw_free(ub);
        fftw_free(plan_u);
        fftw_free(plan_v);
        fftw_free(plan_i);
    }


    //template class Coagulation<double>;
    template class Coagulation<float>;
}
