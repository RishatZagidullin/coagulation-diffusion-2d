#pragma once
#include <sstream>
#include <vector>
#include <string>
#include <iterator>

#include <fstream>

template<typename T>
std::vector<T> split(const std::string& line)
{
    std::istringstream is(line);
    return std::vector<T>(std::istream_iterator<T>(is),
                          std::istream_iterator<T>());
}

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}
