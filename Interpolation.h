#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <vector>

// Lagrange Interpolation Method
double lagrangeInterpolation(const std::vector<double>& x, const std::vector<double>& y, double xi);

// Helper functions for Newton interpolation
std::vector<double> delta_finder(const std::vector<double>& x, const std::vector<double>& y, int n);
double helper_of_newton(std::vector<double> x, double xi);
double newton_forward(std::vector<double> first_deltas, std::vector<double> x, double y0, double xi);
double newton_backward(std::vector<double> last_deltas, std::vector<double> x, double yn, double xi);

// Newton Interpolation Method
double NewtonInterpolation(const std::vector<double>& x, const std::vector<double>& y, double xi);

#endif // INTERPOLATION_H
