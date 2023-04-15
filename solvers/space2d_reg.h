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
        int radius_x;
        int radius_y;
        int N, M;
        T dx, dy, dt, D_x, D_y, J;
        Vector3d<T> * vels;
        void init_data(T x, T y);
        void solve_matrix (int n, T *a, T *c, T *b, T *f, T *x);

        int find_inds(T src_pos, int size, T dd);

    public:
        void diffusion_step(T * data, T * data_new, bool is_monomer);
        void advection_step(T * data, T * data_new);
        Space2d_reg(T D_x, T D_y, T J, int N, int M, T dx, T dy, T dt, Vector3d<int> source, T x, T y, int radius_x, int radius_y):
                    D_x(D_x), D_y(D_y), J(J), N(N), M(M), dx(dx), dy(dy), dt(dt), source(source), radius_x(radius_x), radius_y(radius_y)
                    {
                        std::cout << "radius x: " << radius_x << " raidus y: " << radius_y << "\n";
                        std::cout << "source x: " << source.x << " source y: " << source.y << "\n";
                        init_data(x, y);
                    };
        ~Space2d_reg();
    };
}
