#include "utils/utils.h"
#include "solvers/diffusion_reg.h"
#include "solvers/coagulation.h"
#include "geometry/projection_reg.h"
#include "ppm/ppm.h"
#include <time.h>
#include <sys/time.h>

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int main(int argc, char ** argv)
{
    double start = get_wall_time();

    int N = 100;
    int M = 100;
    int K = 1;
    double dd = 0.05;
    grid_params proj_grid = grid_params(K, M, N, -0.1, -2.0, -2.0, dd);
    projection_reg proj(proj_grid);

    double dt = 0.01;
    Vector3d<int> source{0, std::atoi(argv[1]), std::atoi(argv[2])};
    int S = 32;
    solvers::Diffusion3d_reg<double> eqn(0.1, N, M, K, dd, dt, source);
    solvers::Coagulation<double> coag(S, 1.0, dt);

    int TIME_MAX = 100;
    std::cout <<"Preprocessing time: "<<get_wall_time()-start<<"\n";
    start = get_wall_time();

    double * data = new double [N*M*K*S];
    double * data_new = new double [N*M*K*S];

    for (int time = 0; time < TIME_MAX; time++)
    {
        for (int i = 0; i < S; i++)
            eqn.iteration(data+i*N*M*K, data_new+i*N*M*K, i==0 ? 1 : 0, false);

        if(std::atoi(argv[3]) == 1)
            for (int i = 0; i < N*M*K; i++)
                coag.iteration(data_new+i, data+i, N*M*K);
        else
            std::swap(data, data_new);

        for (int i = 0; i < S; i++)
        {
            std::string name = std::string("./imgs/") + 
                               std::to_string(i) + 
                               std::string("_res_") + 
                               std::to_string(time) + 
                               std::string(".ppm");
            create_ppm(N, M, data+i*N*M*K, name);
        }

        printProgress((double)time/TIME_MAX);
    }
    std::cout <<"\nComputation time: "<< get_wall_time()-start<<"\n";
    delete [] data;
    delete [] data_new;
    return 0;
}



