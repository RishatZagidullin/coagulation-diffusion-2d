#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <omp.h>

int create_ppm(int N, int M, float const * const data, std::string filename, double maxi);
