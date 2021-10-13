#pragma once

#include<vector>
#include <stdexcept>
#include <algorithm>

std::vector<double> operator+(const std::vector<double>& lhs, const std::vector<double>& rhs);
std::vector<double> operator-(const std::vector<double>& lhs, const std::vector<double>& rhs);
std::vector<double> operator*(const double& n, const std::vector<double>& rhs);
std::vector<double> operator*(const std::vector<double>& lhs, const double& n);
std::vector<double> operator/(const std::vector<double>& lhs, const double& n);
std::vector<double> abs(const std::vector<double>& v);
double max(const std::vector<double>& v);
double sign(double& a);
