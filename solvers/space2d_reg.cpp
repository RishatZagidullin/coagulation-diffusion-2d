#include "space2d_reg.h"

namespace solvers
{
    template<typename T>
    void Space2d_reg<T>::init_data(T x, T y)
    {
        this->vels = new Vector3d<T> [N*M];
        int c = 0;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                vels[c].x = x;
                vels[c].y = y;
                vels[c].z = 0;
                c++;
            }
        }
    }

    //keep in mind that a,b,c are altered after the execution
    template<typename T>
    void Space2d_reg<T>::solve_matrix (int n, T *a, T *c, 
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
    void Space2d_reg<T>::diffusion_step(T* data, T* data_new, bool is_monomer)
    {
        //=======================N=================================
        T * a = new T [N];
        T * b = new T [N];
        T * c = new T [N];
        T * rhs = new T [N];

        T * output = new T [N];

        for (int j = 1; j < M-1; j++)
        {
            for (int i = 1; i < N-1; i++)
            {
                a[i] = -dt * D_y * 0.5 / dy / dy;
                b[i] = -dt * D_y * 0.5 / dy / dy;
                c[i] = 1. + dt * D_y / (dy*dy);
                rhs[i] = 1. * data[j+i*M] + 
                             0.5 * dt * D_y / dy / dy *
                             (data[j+(i+1)*M]-data[j+(i)*M]
                              -data[j+(i)*M]+data[j+(i-1)*M]
                             )+
                             dt * D_x / dx / dx *
                             (data[(j+1)+i*M]-data[(j)+i*M]
                              -data[(j)+i*M]+data[(j-1)+i*M]
                             );
            }
            a[0] = 0.0;
            b[0] = 1.0;
            c[0] = -1.0;
            a[N-1] = -1.0;
            b[N-1] = 0.0;
            c[N-1] = 1.0;
            rhs[0] = 0.0;
            rhs[N-1] = 0.0;
            if (j >= source.x-radius_x && j <= source.x+radius_x)
            {
                if (source.y-radius_y >=0)
                {
                    a[source.y-radius_y] = 0.0;
                    b[source.y-radius_y] = 1.0;
                    c[source.y-radius_y] = -1.0;
                    rhs[source.y-radius_y] = (is_monomer ? -J*dy*0.5 : -rhs[source.y-radius_y]*dy*0.5);
                }
                if (source.y+radius_y <=N-1)
                {
                    a[source.y+radius_y] = -1.0;
                    b[source.y+radius_y] = 0.0;
                    c[source.y+radius_y] = 1.0;
                    rhs[source.y+radius_y] = (is_monomer ? J*dy*0.5 : rhs[source.y+radius_y]*dy*0.5);
                }
                for (int i = source.y - radius_y+1; i < source.y + radius_y; i++)
                {
                    a[i] = 0.0;
                    b[i] = 0.0;
                    c[i] = 1.0;
                    rhs[i] = (is_monomer ? J : rhs[i]);
                }
            }
            solve_matrix(N, a, c, b, rhs, output);
            for (int i = 0; i < N; i++)
            {
                data_new[j+i*M] = output[i];
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
            for (int j = 1; j < M-1; j++)
            {
                a[j] = -dt * D_x * 0.5 / dx / dx;
                b[j] = -dt * D_x * 0.5 / dx / dx;
                c[j] = 1. + dt * D_x / (dx*dx);
                rhs[j] = 1. * data_new[j+i*M] + 
                         0.5 * dt * D_x / dx / dx *
                         (-data[(j+1)+i*M]+data[(j)+i*M]
                          +data[(j)+i*M]-data[(j-1)+i*M]
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
            if (i >= source.y-radius_y && i <= source.y+radius_y)
            {
                if (source.x-radius_x >=0)
                {
                    a[source.x-radius_x] = 0.0;
                    b[source.x-radius_x] = 1.0;
                    c[source.x-radius_x] = -1.0;
                    rhs[source.x-radius_x] = (is_monomer ? -J*dx*0.5 : -rhs[source.x-radius_x]*dx*0.5);
                }
                if (source.x+radius_x <=M-1)
                {
                    a[source.x+radius_x] = -1.0;
                    b[source.x+radius_x] = 0.0;
                    c[source.x+radius_x] = 1.0;
                    rhs[source.x+radius_x] = (is_monomer ? J*dx*0.5 : rhs[source.x+radius_x]*dx*0.5);
                }
                for (int j = source.x - radius_x+1; j < source.x + radius_x; j++)
                {
                    a[j] = 0.0;
                    b[j] = 0.0;
                    c[j] = 1.0;
                    rhs[j] = (is_monomer ? J : rhs[j]);
                }
            }
            solve_matrix(M, a, c, b, rhs, output);
            for (int j = 0; j < M; j++)
            {
                data_new[j+i*M] = output[j];
            }
        }
        delete [] a;
        delete [] b;
        delete [] c;
        delete [] rhs;
        delete [] output;
        //=========================================================
        return;
    }


    template<typename T>
    int Space2d_reg<T>::find_inds(T src_pos, int size, T dd)
    {
        int ind = -1;
        for (int i = 0; i < size; i++)
        {
            if (src_pos <= 0.0 || src_pos >= (size-1)*dd)
            {
                ind = -1;
                break;
            }
            else if (src_pos < i*dd)
            {
                //between i and i-1
                //std::cout << fabs(src_pos - i*dd) << " " << fabs(src_pos-(i-1)*dd) << "\n";
                ind = fabs(src_pos-i*dd) > fabs(src_pos-(i-1)*dd) ? i-1 : i;
                break;
            }
            else if (src_pos == i*dd)
            {
                ind = i;
                break;
            }
        }
        return ind;
    }

    template<typename T>
    void Space2d_reg<T>::advection_step(T* data, T* data_new)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                int idx = j+i*M;
                
                Vector3d<T> cur_pos = {j*dx, i*dy, 0};
                Vector3d<T> src_pos = {cur_pos.x-vels[idx].x*dt, 
                                       cur_pos.y-vels[idx].y*dt, 
                                       cur_pos.z-vels[idx].z*dt};
                int h = find_inds(src_pos.y, N, dy);
                int w = find_inds(src_pos.x, M, dx);
                if (h!=-1 && w!=-1) data_new[idx] = data[w+h*M];
            }
        }
    }

    template<typename T>
    Space2d_reg<T>::~Space2d_reg()
    {
        delete [] vels;
    }

    template class Space2d_reg<double>;
    //template class Space2d_reg<float>;
}
