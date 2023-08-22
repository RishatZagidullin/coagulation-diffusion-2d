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

    double hours = 6.0;
    double kms = 10.0;
    double source_len_kms = 0.1;//0.01;
    double size_dim = kms * 1000.0;
    double time_dim = hours *60*60;
    double diffusion_dim = 1e0; // in m^2 s^-1
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
        std::cout << "NO VELOCITY, ONLY DIFFUSION\n";
        N = 200;
        M = 200;
        dy = 1./N;
        dx = 1./M;
        TIME_MAX = 400;
    }
    else if (TIME_MAX == 0)
    {
        TIME_MAX = TIME_MAX_OTHER;
        dy = dx;
        N = M;
    }
    else if (TIME_MAX_OTHER == 0)
    {
        TIME_MAX_OTHER = TIME_MAX;
        dx = dy;
        M = N;
    }

    double dt = 1./TIME_MAX;

    std::cout << "X: " << M << " Y: " << N << "\n";
    std::cout << "TOTAL ITERATIONS: " << TIME_MAX << "\n";

    Vector3d<int> source{(int) (std::stod(argv[1])*M), (int) (std::stod(argv[2])*N), 0};
    int S = 500;

    double diffusion_undim_x = diffusion_dim * pow(M*dx,2)/pow(size_dim,2) * time_dim/(TIME_MAX*dt);
    double diffusion_undim_y = diffusion_dim * pow(N*dy,2)/pow(size_dim,2) * time_dim/(TIME_MAX*dt);

    std::cout << "size in meters: " << size_dim << "\n";
    std::cout << "time in hours: " << time_dim/60/60 << "\n";
    std::cout << "diffusion coef (in m^2 s^-1): " << diffusion_dim << "\n";
    std::cout << "undimensional diffusion coefs: " << diffusion_undim_x << " " << diffusion_undim_y << "\n";
    double vel_undim_x = vel_dim * M*dx/size_dim * time_dim/(TIME_MAX*dt);
    double vel_undim_y = vel_dim * N*dy/size_dim * time_dim/(TIME_MAX*dt);
    std::cout << "velocity magnitude (in m^1 s^-1): " << vel_dim << "\n";
    std::cout << "undimensional velocity magnitudes: " << vel_undim_x << " " << vel_undim_y << "\n";

    std::cout << "angle (in radians): " << angle/3.14159265358979311600 << " pi\n";

    double J_dim = std::stod(argv[6]); // in сm^-3 s^-1
    
    double J = J_dim > 100 ? 100 : J_dim;
    double J_scale = J_dim > 100 ? J_dim/100 : 1;
    std::cout << "source term magnitude (in cm^-3): " << J_dim << "\n";

    int radius_x = (int)(source_len_kms / (kms/M) );
    int radius_y = (int)(source_len_kms / (kms/N) );
    
    if (radius_x == 0 || radius_y == 0)
    {
        std::cout << "source radius is too small, try increasing it\n";
        std::cout << "radius x: " << radius_x << " raidus y: " << radius_y << "\n";
        return -1;
    }
    solvers::Space2d_reg<double> eqn(diffusion_undim_x, diffusion_undim_y, J, N, M, dx, dy, dt, source,
                                     vel_undim_x*x_vel, vel_undim_y*y_vel, radius_x, radius_y);
    auto kernel = solvers::default_crossed_kernel(1e-4, S);
    solvers::Coagulation<double> coag(S, dt, kernel);

    
    std::cout << "Preprocessing time: " << get_wall_time() - start << "\n";
    start = get_wall_time();

    double * data = new double [N*M*S];
    double * data_new = new double [N*M*S];

    for (int i = 0; i < N*M*S; i++)
    {
        data[i] = 0.0;
        data_new[i] = 0.0;
    }

    int img_num = 0;
    double maxi = 0.0;
    for (int time = 0; time < TIME_MAX; time++)
    {
        if(std::atoi(argv[3]) == 1)
        {
            //#pragma omp parallel for
            for (int i = 0; i < N*M; i++)
                //coag.iteration_direct(data+i, data_new+i, N*M);
                coag.iteration_low_rank(data+i, data_new+i, N*M);
            std::swap(data, data_new);
        }

        #pragma omp parallel for
        for (int i = 0; i < S; i++)
        {
            bool is_monomer = (i==0);
            eqn.diffusion_step(data+i*N*M, data_new+i*N*M, is_monomer);
            if (vel_dim != 0.0)
                eqn.advection_step(data_new+i*N*M, data+i*N*M);
        }
        if (vel_dim == 0.0)
            std::swap(data, data_new);

        if (time == 0)
        {
            for (int i = 0; i < S; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    for (int k = 0; k < N; k++)
                    {
                        if (data[j+k*M+i*N*M] > maxi)
                            maxi = data[j+k*M+i*N*M];
                    }
                }
            }
            std::cout << "max concentration (in cm^-3): " << maxi * J_scale << "\n";
            std::cout << "max undimensional concentration: " << maxi << "\n";
            std::ofstream concentration;
            concentration.open("concentration.txt");
            concentration << maxi * J_scale;
            concentration.close();
        }

        if (time == TIME_MAX-1)
        {
            #pragma omp parallel for
            for (int i = 0; i < S; i++)
            {
                std::string name = std::string("./imgs/") + 
                               std::to_string(i) + 
                               std::string("_res.ppm");
                create_ppm(N, M, data+i*N*M, name, maxi);
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



