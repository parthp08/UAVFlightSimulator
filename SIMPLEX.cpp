#include "SIMPLEX.h"

// private functions
namespace {
	std::vector<size_t> sort_indexes(const std::vector<double>& v) {

		// initialize original index locations
		std::vector<size_t> idx(v.size());
		std::iota(idx.begin(), idx.end(), 0);

		std::stable_sort(idx.begin(), idx.end(),
			[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

		return idx;
	}

	template <typename T>
	std::vector<T> take(const std::vector<T>& v, const std::vector<size_t>& idx) {
		// init sorted array
		std::vector<T> sorted_v;
		sorted_v.resize(v.size());
		for (int i = 0; i < (int)idx.size(); i++) {
			sorted_v[i] = v[idx[i]];
		}
		return sorted_v;
	}

	bool check_convergence(const std::vector<std::vector<double>>& points, const std::vector<double>& sol_arr,
		const double& nstates, const double& xtol, const double& ftol) {

		std::vector<double> max_point_vec;
		std::vector<double> max_sol_vec;
		for (int i = 1; i < nstates + 1; i++) {
			max_point_vec.push_back(max(abs(points[i] - points[0])));
			max_sol_vec.push_back(std::abs(sol_arr[0] - sol_arr[i]));
		}

		return ((max(max_point_vec) <= xtol) && (max(max_sol_vec) <= ftol));
	}
} // end namespace

std::vector<double> SIMPLEX(std::function<double (std::vector<double>)> cost_func, const std::vector<double>& x0, 
	const double xtol, const double ftol, const int maxiter) {

	const double n = (double)x0.size(); // no of states

	const double delta = 0.05; // adding 5%
	const double delta0 = 0.00025; // replace 0s with this

	const double rho = 1;  // reflect
	const double chi = 2;  // extend
	const double psi = 0.5;  // contract
	const double sigma = 0.5;  // shrink

	// simplex points and solution array
	std::vector<std::vector<double>> points;
	std::vector<double> sol_arr;
	points.push_back(x0);
	sol_arr.push_back(cost_func(x0));
	for (int i = 0; i < n; i++) {
		points.push_back(x0);
		if (x0[i] != 0) {
			points[i + 1][i] *= (1 + delta);
		}
		else {
			points[i + 1][i] = delta0;
		}
		sol_arr.push_back(cost_func(points[i + 1]));
	}

	// sort points and sol array
	std::vector<size_t> sorted_indx = sort_indexes(sol_arr);
	sol_arr = take(sol_arr, sorted_indx);
	points = take(points, sorted_indx);

	int niter = 1; // no of iteration
	bool shrink = false;

	std::vector<double> xm, xr, xs, xc, xcc;
	double fxr, fxs, fxc, fxcc;

	while (niter < maxiter) {

		// check convergence
		if (check_convergence(points, sol_arr, n, xtol, ftol)) {
			break;
		}

		// increase the no of iteration
		niter += 1;

		// reflection point
		xm = points[0];
		for (int i = 1; i < n; i++) {
			xm = xm + points[i];
		}
		xm = xm / n;
		xr = (1 + rho) * xm - rho * points[n];
		fxr = cost_func(xr);

		if (fxr < sol_arr[0]) {
			// expansion point
			xs = (1 + rho * chi) * xm - rho * chi * points[n];
			fxs = cost_func(xs);

			if (fxs < fxr) {
				points[n] = xs;
				sol_arr[n] = fxs;
			}
			else {
				points[n] = xr;
				sol_arr[n] = fxr;
			}
		}
		else {
			if (fxr >= sol_arr[n - 1]) {
				// contract point
				if (fxr <= sol_arr[n]) {
					// contract outside
					xc = (1 + psi * rho) * xm - psi * rho * points[n];
					fxc = cost_func(xc);
					if (fxc <= fxr) {
						points[n] = xc;
						sol_arr[n] = fxc;
					}
					else {
						shrink = true;
					}
				}
				else {
					// contract inside
					xcc = (1 - psi) * xm + psi * points[n];
					fxcc = cost_func(xcc);
					if (fxcc <= sol_arr[n - 1]) {
						points[n] = xcc;
						sol_arr[n] = fxcc;
					}
					else {
						shrink = true;
					}
				}
			}
			else {
				points[n] = xr;
				sol_arr[n] = fxr;
			}
		}

		// shrink points
		if (shrink) {
			for (int i = 1; i < n + 1; i++) {
				points[i] = points[0] + 0.5 * (points[i] - points[0]);
				sol_arr[i] = cost_func(points[i]);
			}
			shrink = false;
		}

		// sorting
		sorted_indx = sort_indexes(sol_arr);
		sol_arr = take(sol_arr, sorted_indx);
		points = take(points, sorted_indx);

	} // end while

	std::cout << "Exiting SIMPLEX function: \n" <<
		"\tfval = " << sol_arr[0] <<
		"\n\titerations = " << niter << std::endl;

	return points[0];
}

