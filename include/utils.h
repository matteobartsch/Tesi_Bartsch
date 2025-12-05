#ifndef __UTILS_H__
#define __UTILS_H__

#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <chrono> 
#include <filesystem>
#include <limits>
#include <algorithm>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    for (size_t i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i + 1 != v.size()) os << " ";
    }
    return os;
}

auto noise(std::size_t nrow, std::size_t ncol, double sigma, std::mt19937& gen) {
    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(nrow, ncol);
    std::normal_distribution<> __noise(0.0, sigma);
    for (std::size_t i = 0; i < nrow; ++i){
      for( std::size_t j = 0; j < ncol; ++j ){
       res(i, j) = __noise(gen); 
      }
    }
    return res;
}

std::string executable_name(const char* argv0) {
    std::string s(argv0);
    size_t pos = s.find_last_of("/\\");
    if (pos == std::string::npos)
        return s; 
    return s.substr(pos + 1);
}

#endif // __UTILS_H__