#pragma once
#include <iostream>
#include <fstream>

//d->x, w->y, h->z
struct grid_params{
    int d, w, h;
    double left_d, left_w, left_h, dd;
    grid_params(int d, int w, int h,
                double left_d, double left_w,
                double left_h, double dd) : d(d), w(w), h(h),
                    left_d(left_d), left_w(left_w), left_h(left_h),
                    dd(dd){}
};

class projection_reg {
public:
    grid_params grid;
    int max_v, max_i;
    float * vertices;
    unsigned int * indices;
    projection_reg(grid_params & grid) : grid(grid)
    {
        init_data();
    }
    void out_to_file(std::string vfile, std::string ffile);
    ~projection_reg();
private:
    void init_data();
};

void projection_reg::out_to_file(std::string vfile, std::string ffile)
{
    std::ofstream v_reader;
    v_reader.open(vfile);
    v_reader << grid.d << " " << grid.w << " " << grid.h << " " << grid.dd << "\n";
    std::ofstream f_reader;
    f_reader.open(ffile);
    for (int i = 0; i < grid.h; i++)
        for (int j = 0; j < grid.w; j++)
            v_reader << grid.left_h+i*grid.dd << " " 
                     << grid.left_w+j*grid.dd << " 0.0\n";

    for (int i = 0; i < grid.h-1; i++){
        int offset = i*grid.w;
        for (int j = 0; j < grid.w-1; j++){
            f_reader << j+offset << " " << j+1+offset << " "
                     << j+grid.w+offset << "\n";
            f_reader << j+1+offset << " " << j+grid.w+offset
                     << " " << j+grid.w+1+offset << "\n";
        }
    }
    for (int i = 0; i < grid.w; i++)
        for (int j = 0; j < grid.d; j++)
            v_reader << "0.0 " << grid.left_w+i*grid.dd << " "
                     << grid.left_d+j*grid.dd << std::endl;

    for (int i = 0; i < grid.w-1; i++){
        int offset = i*grid.d + grid.h*grid.w;
        for (int j = 0; j < grid.d-1; j++){
            f_reader << j+offset << " " << j+1+offset << " "
                     << j+grid.d+offset << "\n";
            f_reader << j+1+offset << " " << j+grid.d+offset
                     << " " << j+grid.d+1+offset << "\n";
        }
    }

    for (int i = 0; i < grid.h; i++)
        for (int j = 0; j < grid.d; j++)
            v_reader << grid.left_h+i*grid.dd << " 0.0 "
                     << grid.left_d+j*grid.dd << std::endl;

    for (int i = 0; i < grid.h-1; i++){
        int offset = i*grid.d + grid.h*grid.w + grid.w*grid.d;
        for (int j = 0; j < grid.d-1; j++){
            f_reader << j+offset << " " << j+1+offset << " "
                     << j+grid.d+offset << "\n";
            f_reader << j+1+offset << " " << j+grid.d+offset
                     << " " << j+grid.d+1+offset << "\n";
        }
    }

    v_reader.close();
    f_reader.close();
}

void projection_reg::init_data()
{
    max_v = (grid.h*grid.w+grid.h*grid.d+grid.w*grid.d);
    max_i = ((grid.h-1)*(grid.w-1)+(grid.h-1)*(grid.d-1)
                                  +(grid.w-1)*(grid.d-1));
    vertices = new float [max_v*4];
    indices = new unsigned int [max_i*6];

    for (int i = 0; i < grid.h; i++)
        for (int j = 0; j < grid.w; j++)
        {
            vertices[4*(i*grid.w+j)] = grid.left_h+i*grid.dd; 
            vertices[4*(i*grid.w+j)+1] = grid.left_w+j*grid.dd;
            vertices[4*(i*grid.w+j)+2] = 0.0;
            vertices[4*(i*grid.w+j)+3] = 0.0;
        }

    for (int i = 0; i < grid.h-1; i++)
    {
        int offset = i*grid.w;
        for (int j = 0; j < grid.w-1; j++)
        {
            indices[6*(i*(grid.w-1)+j)] = j+offset;
            indices[6*(i*(grid.w-1)+j)+1] = j+1+offset;
            indices[6*(i*(grid.w-1)+j)+2] = j+grid.w+offset;
            indices[6*(i*(grid.w-1)+j)+4] = j+1+offset;
            indices[6*(i*(grid.w-1)+j)+3] = j+grid.w+offset;
            indices[6*(i*(grid.w-1)+j)+5] = j+grid.w+1+offset;
        }
    }
    int global_offset = grid.h*grid.w;
    for (int i = 0; i < grid.w; i++)
        for (int j = 0; j < grid.d; j++)
        {
            vertices[4*(i*grid.d+j+global_offset)] = 0.0; 
            vertices[4*(i*grid.d+j+global_offset)+1] = 
                                        grid.left_w+i*grid.dd;
            vertices[4*(i*grid.d+j+global_offset)+2] =
                                        grid.left_d+j*grid.dd;
            vertices[4*(i*grid.d+j+global_offset)+3] = 0.0;
        }
    global_offset = (grid.h-1)*(grid.w-1);
    for (int i = 0; i < grid.w-1; i++)
    {
        int offset = i*grid.d + grid.h*grid.w;
        for (int j = 0; j < grid.d-1; j++)
        {
            indices[6*(i*(grid.d-1)+j+global_offset)] = j+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+1] = j+1+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+2] = 
                                                     j+grid.d+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+4] = j+1+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+3] = 
                                                     j+grid.d+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+5] = 
                                                   j+grid.d+1+offset;
        }
    }
    global_offset = grid.h*grid.w + grid.w*grid.d;
    for (int i = 0; i < grid.h; i++)
        for (int j = 0; j < grid.d; j++)
        {
            vertices[4*(i*grid.d+j+global_offset)] = 
                                        grid.left_h+i*grid.dd; 
            vertices[4*(i*grid.d+j+global_offset)+1] = 0.0;
            vertices[4*(i*grid.d+j+global_offset)+2] =
                                        grid.left_d+j*grid.dd;
            vertices[4*(i*grid.d+j+global_offset)+3] = 0.0;
        }
    global_offset = (grid.h-1)*(grid.w-1) + (grid.w-1)*(grid.d-1);
    for (int i = 0; i < grid.h-1; i++)
    {
        int offset = i*grid.d + grid.h*grid.w + grid.w*grid.d;
        for (int j = 0; j < grid.d-1; j++)
        {
            indices[6*(i*(grid.d-1)+j+global_offset)] = j+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+1] = j+1+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+2] = 
                                                    j+grid.d+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+4] = j+1+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+3] = 
                                                    j+grid.d+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+5] = 
                                                  j+grid.d+1+offset;
        }
    }

    return;
}

projection_reg::~projection_reg()
{
    delete [] vertices;
    delete [] indices;
}
