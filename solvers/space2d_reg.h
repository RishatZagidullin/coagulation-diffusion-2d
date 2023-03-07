#pragma once
#include <iostream>
#include <vector>
#include "../geometry/vector3d.h"

namespace solvers
{
    template <typename T>
    class Space2d_reg
    {
    private:
        Vector3d<int> source;
        //N->z, M->y, K->x
        int N, M;
        T dx, dy, dt, D_x, D_y;
        Vector3d<T> * vels;
        void init_data(T x, T y);
        void solve_matrix (int n, T *a, T *c, T *b, T *f, T *x);

        int find_inds(T src_pos, int size, T dd);

    public:
        void add_source(T * data, int radius_x, int radius_y);
        void diffusion_step(T * data, T * data_new);
        void advection_step(T * data, T * data_new);
        Space2d_reg(T D_x, T D_y, int N, int M, T dx, T dy, T dt, Vector3d<int> source, T x = 0., T y = 0.):
                    D_x(D_x), D_y(D_y), N(N), M(M), dx(dx), dy(dy), dt(dt), source(source)
                    {
                        init_data(x, y);
                    };
        ~Space2d_reg();
    };
}
