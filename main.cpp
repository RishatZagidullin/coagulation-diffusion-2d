#include "utils/utils.h"
#include "solvers/space2d_reg.h"
#include "solvers/coagulation.h"
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

    double hours = 48.0;
    double kms = 1.0;
    double size_dim = kms * 1000.0;
    double time_dim = hours *60*60;
    double diffusion_dim = 1e-2; // in m^2 s^-1
    double vel_dim = std::stod(argv[4]); // in m^1 s^-1

    double angle = std::stod(argv[5])/180.*3.14159265358979311600;
    double x_vel = cos(-angle);
    double y_vel = sin(-angle);
    std::cout << "x_vel: " << x_vel << " y_vel: " << y_vel << "\n";
    int N = (int) (1./fabs(y_vel)*200);
    int M = (int) (1./fabs(x_vel)*200);
    double dy = 1./N;
    double dx = 1./M;

    int TIME_MAX = (int) (time_dim/((size_dim*dy)/(fabs(y_vel)*vel_dim)));
    int TIME_MAX_OTHER = (int) (time_dim/((size_dim*dx)/(fabs(x_vel)*vel_dim)));
    
    if (TIME_MAX == 0 && TIME_MAX_OTHER == 0)
    {
        std::cout << "NO VELOCITY, THIS WILL PROBABLY BREAK\n";
    }
    if (TIME_MAX == 0)
    {
        TIME_MAX = TIME_MAX_OTHER;
        dy = dx;
        N = M;
    }
    if (TIME_MAX_OTHER == 0)
    {
        TIME_MAX_OTHER = TIME_MAX;
        dx = dy;
        M = N;
    }
    
    //TIME_MAX = TIME_MAX < 100 ? 100 : TIME_MAX;
    std::cout << "X: " << M << " Y: " << N << "\n";
    std::cout << "TOTAL ITERATIONS: " << TIME_MAX << "\n";
    double dt = 1./TIME_MAX;
    Vector3d<int> source{(int) (std::stod(argv[1])*M), (int) (std::stod(argv[2])*N), 0};
    int S = 1;

    double diffusion_undim_x = diffusion_dim * pow(M*dx,2)/pow(size_dim,2) * time_dim/(TIME_MAX*dt);
    double diffusion_undim_y = diffusion_dim * pow(N*dy,2)/pow(size_dim,2) * time_dim/(TIME_MAX*dt);

    std::cout << "size in meters: " << size_dim << "\n";
    std::cout << "time in hours: " << time_dim/60/60 << "\n";
    std::cout << "diffusion coef (in m^2 s^-1): " << diffusion_dim << "\n";
    //std::cout << "undimensional diffusion coefs: " << diffusion_undim << "\n";
    double vel_undim_x = vel_dim * M*dx/size_dim * time_dim/(TIME_MAX*dt);
    double vel_undim_y = vel_dim * N*dy/size_dim * time_dim/(TIME_MAX*dt);
    std::cout << "velocity magnitude (in m^1 s^-1): " << vel_dim << "\n";
    //std::cout << "undimensional velocity magnitude: " << vel_undim << "\n";

    std::cout << "angle (in radians): " << angle/3.14159265358979311600 << " pi\n";

    //return 0;
    solvers::Space2d_reg<double> eqn(diffusion_undim_x, diffusion_undim_y, N, M, dx, dy, dt, source, vel_undim_x*x_vel, vel_undim_y*y_vel);
    solvers::Coagulation<double> coag(S, 1.0, dt);

    
    std::cout <<"Preprocessing time: "<<get_wall_time()-start<<"\n";
    start = get_wall_time();

    double * data = new double [N*M*S];
    double * data_new = new double [N*M*S];

    for (int i = 0; i < N*M*S; i++)
    {
        data[i] = 0.0;
        data_new[i] = 0.0;
    }

    int img_num = 0;
    for (int time = 0; time < TIME_MAX; time++)
    {
        for (int i = 0; i < S; i++)
        {
            if (i==0) eqn.add_source(data, M*0.005, N*0.005);
            eqn.diffusion_step(data+i*N*M, data_new+i*N*M);
            eqn.advection_step(data_new+i*N*M, data+i*N*M);
            //std::swap(data, data_new);
        }

        if(std::atoi(argv[3]) == 1)
        {
            for (int i = 0; i < N*M; i++)
                coag.iteration(data+i, data_new+i, N*M);
            std::swap(data, data_new);
        }

        if (time % (int) (TIME_MAX/100.) == 0)
        {
            for (int i = 0; i < S; i++)
            {
                std::string name = std::string("./imgs/") + 
                               std::to_string(i) + 
                               std::string("_res_") + 
                               std::to_string(img_num) + 
                               std::string(".ppm");
                create_ppm(N, M, data+i*N*M, name);
            }
            img_num++;
        }

        printProgress((double)time/TIME_MAX);
    }
    std::cout <<"\nComputation time: "<< get_wall_time()-start<<"\n";
    delete [] data;
    delete [] data_new;
    return 0;
}



