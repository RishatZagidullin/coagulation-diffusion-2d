#include "diffusion_reg.h"

namespace solvers
{
    template<typename T>
    void Diffusion3d_reg<T>::init_data()
    {
        this->vels = new Vector3d<T> [N*M*K];
        int c = 0;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                for (int k = 0; k < K; k++)
                {
                    vels[c].x = 0.0;
                    vels[c].y = 0.0;
                    vels[c].z = 0.0;
                    c++;
                }
            }
        }
    }

    //keep in mind that a,b,c are altered after the execution
    template<typename T>
    void Diffusion3d_reg<T>::solve_matrix (int n, T *a, T *c, 
                                       T *b, T *f, T *x)
    {
        T m;
        for (int i = 1; i < n; i++)
        {
            m = a[i] / c[i-1];
            c[i] = c[i] - m * b[i-1];
            f[i] = f[i] - m * f[i-1];
        }
        x[n-1] = f[n-1]/c[n-1];

        for (int i = n - 2; i >= 0; i--)
        {
            x[i] = ( f[i] - b[i] * x[i+1] ) / c[i];
        }
    }

    template<typename T>
    void Diffusion3d_reg<T>::iteration(T* data, T* data_new, int coef, bool boundary)
    {
        if (coef == 1)
            data[source.x+source.y*K+source.z*K*M] = 1.0;
        //=======================N=================================
        T * a = new T [N];
        T * b = new T [N];
        T * c = new T [N];
        T * rhs = new T [N];

        T * output = new T [N];

        for (int j = 1; j < M-1; j++)
        {
            for (int k = 0; k < K; k++)
            {
                for (int i = 1; i < N-1; i++)
                {
                    a[i] = -dt * D * 0.5 / dd / dd *
                            exp(vels[k+j*K+(i-1)*K*M].z/D*dd/2.0);
                    b[i] = -dt * D * 0.5 / dd / dd *
                            exp(-vels[k+j*K+(i+1)*K*M].z/D*dd/2.0);
                    c[i] = 1. + dt * 
                          (D * 0.5 * exp(vels[k+j*K+i*K*M].z/D*dd/2.0)
                          + D * 0.5 * exp(-vels[k+j*K+i*K*M].z/D*dd/2.0) 
                          ) / (dd*dd);
                    rhs[i] = 1. * data[k+j*K+i*K*M] + 
                             0.5 * dt * D / dd / dd *
                             (data[k+j*K+(i+1)*K*M]*
                              exp(-vels[k+j*K+(i+1)*K*M].z/D*dd/2.0)
                              -data[k+j*K+(i)*K*M]*
                              exp(-vels[k+j*K+(i)*K*M].z/D*dd/2.0)
                              -data[k+j*K+(i)*K*M]*
                              exp(vels[k+j*K+(i)*K*M].z/D*dd/2.0)
                              +data[k+j*K+(i-1)*K*M]*
                              exp(vels[k+j*K+(i-1)*K*M].z/D*dd/2.0)
                             )+
                             dt * D / dd / dd *
                             (data[k+(j+1)*K+i*K*M]*
                              exp(-vels[k+(j+1)*K+i*K*M].y/D*dd/2.0)
                              -data[k+(j)*K+i*K*M]*
                              exp(-vels[k+(j)*K+i*K*M].y/D*dd/2.0)
                              -data[k+(j)*K+i*K*M]*
                              exp(vels[k+(j)*K+i*K*M].y/D*dd/2.0)
                              +data[k+(j-1)*K+i*K*M]*
                              exp(vels[k+(j-1)*K+i*K*M].y/D*dd/2.0)
                             )/*+
                             dt * D / dd / dd *
                             (data[(k+1)+j*K+i*K*M]*
                              exp(-vels[(k+1)+j*K+i*K*M].x/D*dd/2.0)
                              -data[(k)+j*K+i*K*M]*
                              exp(-vels[(k)+j*K+i*K*M].x/D*dd/2.0)
                              -data[(k)+j*K+i*K*M]*
                              exp(vels[(k)+j*K+i*K*M].x/D*dd/2.0)
                              +data[(k-1)+j*K+i*K*M]*
                              exp(vels[(k-1)+j*K+i*K*M].x/D*dd/2.0)
                             )*/;
                }
                a[0] = 0.0;
                b[0] = 1.0;
                c[0] = -1.0;
                a[N-1] = -1.0;
                b[N-1] = 0.0;
                c[N-1] = 1.0;
                rhs[0] = 0.0;
                rhs[N-1] = 0.0;
                solve_matrix(N, a, c, b, rhs, output);
                for (int i = 0; i < N; i++)
                {
                    data_new[k+j*K+i*K*M] = output[i];
                }
            }
        }
        delete [] a;
        delete [] b;
        delete [] c;
        delete [] rhs;
        delete [] output;
        //=======================M=================================
        a = new T [M];
        b = new T [M];
        c = new T [M];
        rhs = new T [M];

        output = new T [M];

        for (int i = 1; i < N-1; i++)
        {
            for (int k = 0; k < K; k++)
            {
                for (int j = 1; j < M-1; j++)
                {
                    a[j] = -dt * D * 0.5 / dd / dd *
                            exp(vels[k+(j-1)*K+i*K*M].y/D*dd/2.0);
                    b[j] = -dt * D * 0.5 / dd / dd *
                            exp(-vels[k+(j+1)*K+i*K*M].y/D*dd/2.0);
                    c[j] = 1. + dt * 
                          (D * 0.5 * exp(vels[k+j*K+i*K*M].y/D*dd/2.0)
                          + D * 0.5 * exp(-vels[k+j*K+i*K*M].y/D*dd/2.0) 
                          ) / (dd*dd);
                    rhs[j] = 1. * data_new[k+j*K+i*K*M] + 
                             0.5 * dt * D / dd / dd *
                             (-data[k+(j+1)*K+i*K*M]*
                              exp(-vels[k+(j+1)*K+i*K*M].y/D*dd/2.0)
                              +data[k+(j)*K+i*K*M]*
                              exp(-vels[k+(j)*K+i*K*M].y/D*dd/2.0)
                              +data[k+(j)*K+i*K*M]*
                              exp(vels[k+(j)*K+i*K*M].y/D*dd/2.0)
                              -data[k+(j-1)*K+i*K*M]*
                              exp(vels[k+(j-1)*K+i*K*M].y/D*dd/2.0)
                             );
                }
                a[0] = 0.0;
                b[0] = 1.0;
                c[0] = -1.0;
                a[M-1] = -1.0;
                b[M-1] = 0.0;
                c[M-1] = 1.0;
                rhs[0] = 0.0;
                rhs[M-1] = 0.0;
                solve_matrix(M, a, c, b, rhs, output);
                for (int j = 0; j < M; j++)
                {
                    data_new[k+j*K+i*K*M] = output[j];
                }
            }
        }
        delete [] a;
        delete [] b;
        delete [] c;
        delete [] rhs;
        delete [] output;
        //=======================K=================================
        /*a = new T [K];
        b = new T [K];
        c = new T [K];
        rhs = new T [K];

        output = new T [K];

        for (int i = 1; i < N-1; i++)
        {
            for (int j = 1; j < M-1; j++)
            {
                for (int k = 1; k < K-1; k++)
                {
                    a[k] = -dt * D * 0.5 / dd / dd *
                            exp(vels[(k-1)+j*K+i*K*M].x/D*dd/2.0);
                    b[k] = -dt * D * 0.5 / dd / dd *
                            exp(-vels[(k+1)+j*K+i*K*M].x/D*dd/2.0);
                    c[k] = 1. + dt * 
                          (D * 0.5 * exp(vels[k+j*K+i*K*M].x/D*dd/2.0)
                          + D * 0.5 * exp(-vels[k+j*K+i*K*M].x/D*dd/2.0) 
                          ) / (dd*dd);
                    rhs[k] = 1. * data_new[k+j*K+i*K*M] + 
                             0.5 * dt * D / dd / dd *
                             (-data[(k+1)+j*K+i*K*M]*
                              exp(-vels[(k+1)+j*K+i*K*M].x/D*dd/2.0)
                              +data[(k)+j*K+i*K*M]*
                              exp(-vels[(k)+j*K+i*K*M].x/D*dd/2.0)
                              +data[(k)+j*K+i*K*M]*
                              exp(vels[(k)+j*K+i*K*M].x/D*dd/2.0)
                              -data[(k-1)+j*K+i*K*M]*
                              exp(vels[(k-1)+j*K+i*K*M].x/D*dd/2.0)
                             );
                }
                a[0] = 0.0;
                b[0] = 1.0;
                c[0] = -1.0;
                a[K-1] = -1.0;
                b[K-1] = 0.0;
                c[K-1] = 1.0;
                rhs[0] = 0.0;
                rhs[K-1] = 0.0;
                solve_matrix(K, a, c, b, rhs, output);
                for (int k = 0; k < K; k++)
                {
                    data_new[k+j*K+i*K*M] = output[k];
                }
            }
        }
        delete [] a;
        delete [] b;
        delete [] c;
        delete [] rhs;
        delete [] output;*/
        //=========================================================
        return;
    }

    template<typename T>
    Diffusion3d_reg<T>::~Diffusion3d_reg()
    {
        delete [] vels;
    }

    template class Diffusion3d_reg<double>;
    template class Diffusion3d_reg<float>;
}
