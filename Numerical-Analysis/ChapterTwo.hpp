#include <iostream>
#include <string>
#include <cmath>
#include <vector>







// Lagrange Interpolation Method
// Interpolates a function based on a set of known points
double lagrangeInterpolation(const std::vector<double>& x, const std::vector<double>& y, double xi) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("Input vectors must be non-empty and of the same size");
    }

    double result = 0.0;
    int n = x.size();

    for (int i = 0; i < n; i++) {
        double term = y[i];
        for (int j = 0; j < n; j++) {
            if (j != i) {
                if (x[i] == x[j]) {
                    throw std::runtime_error("Input x values must be distinct");
                }
                term *= (xi - x[j]) / (x[i] - x[j]);
            }
        }
        result += term;
    }

    return result;
}

// Calculate factorial
int factorial(int n) {
    if (n <= 1) return 1;
    return n * factorial(n - 1);
}

// Newton's Forward Difference Interpolation
// Uses forward differences and works best when interpolating near the beginning of the data set
double newtonForwardInterpolation(const std::vector<double>& x, const std::vector<double>& y, double xi) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("Input vectors must be non-empty and of the same size");
    }

    int n = x.size();

    // Check if x values are equally spaced
    double h = x[1] - x[0];
    for (int i = 2; i < n; i++) {
        if (std::abs((x[i] - x[i - 1]) - h) > 1e-10) {
            throw std::runtime_error("Newton's forward difference requires equally spaced x values");
        }
    }

    // Calculate the parameter u = (xi - x[0]) / h
    double u = (xi - x[0]) / h;

    // Create forward difference table
    std::vector<std::vector<double>> forwardDiff(n, std::vector<double>(n));

    // Fill first column with y values
    for (int i = 0; i < n; i++) {
        forwardDiff[i][0] = y[i];
    }

    // Calculate forward differences
    for (int j = 1; j < n; j++) {
        for (int i = 0; i < n - j; i++) {
            forwardDiff[i][j] = forwardDiff[i + 1][j - 1] - forwardDiff[i][j - 1];
        }
    }

    // Apply Newton's forward formula
    double result = forwardDiff[0][0];
    double uTerm = 1.0;

    for (int i = 1; i < n; i++) {
        uTerm *= (u - i + 1) / i;
        result += uTerm * forwardDiff[0][i];
    }

    return result;
}

// Newton's Backward Difference Interpolation
// Uses backward differences and works best when interpolating near the end of the data set
double newtonBackwardInterpolation(const std::vector<double>& x, const std::vector<double>& y, double xi) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("Input vectors must be non-empty and of the same size");
    }

    int n = x.size();

    // Check if x values are equally spaced
    double h = x[1] - x[0];
    for (int i = 2; i < n; i++) {
        if (std::abs((x[i] - x[i - 1]) - h) > 1e-10) {
            throw std::runtime_error("Newton's backward difference requires equally spaced x values");
        }
    }

    // Calculate the parameter u = (xi - x[n-1]) / h
    double u = (xi - x[n - 1]) / h;

    // Create backward difference table
    std::vector<std::vector<double>> backwardDiff(n, std::vector<double>(n));

    // Fill first column with y values
    for (int i = 0; i < n; i++) {
        backwardDiff[i][0] = y[i];
    }

    // Calculate backward differences
    for (int j = 1; j < n; j++) {
        for (int i = n - 1; i >= j; i--) {
            backwardDiff[i][j] = backwardDiff[i][j - 1] - backwardDiff[i - 1][j - 1];
        }
    }

    // Apply Newton's backward formula
    double result = backwardDiff[n - 1][0];
    double uTerm = 1.0;

    for (int i = 1; i < n; i++) {
        uTerm *= (u + i - 1) / i;
        result += uTerm * backwardDiff[n - 1][i];
    }

    return result;
}