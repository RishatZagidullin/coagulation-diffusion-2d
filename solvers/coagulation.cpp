#include "coagulation.h"
#include <cmath>

namespace solvers
{
    template<typename T>
    Coagulation<T>::Coagulation(int size, T dt, T kernel)
    {
        this->dt = dt;
        this->size = size;
        this->kernel = kernel;
    }

    template<typename T>
    void Coagulation<T>::iteration(T * data, T * data_new, int coef)
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
        delete [] n_k;
    }

    template<typename T>
    T Coagulation<T>::L1(const int &N, const int &i, const T *n)
    {
        T l1 = 0;
        for (int i1 = 0; i1 < i; i1++)
            l1 += n[i1] * n[i-i1-1] * K(i - i1 - 1, i1);
        return l1;
    }

    template<typename T>
    T Coagulation<T>::L2(const int &N, const int &i, const T *n)
    {
        T l2 = 0;
        for (int i1 = 0; i1 < N; i1++)
            l2 += n[i1] * K(i, i1);
        return l2;
    }


    template class Coagulation<double>;
    //template class Coagulation<float>;
}
