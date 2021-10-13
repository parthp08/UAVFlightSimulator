#pragma once

#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>
#include <iostream>
#include "utils.h"

std::vector<double> SIMPLEX(std::function<double (std::vector<double>)> cost_func, const std::vector<double>& x0,
	const double xtol = 1e-3, const double ftol = 1e-3, const int maxiter = 1000);
