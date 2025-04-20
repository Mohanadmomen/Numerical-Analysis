#include <iostream>
#include <string>
#include <cmath>
#include "FunctionHandler.hpp"
#include <stdexcept>

#pragma once 





// Bisection Method for finding roots of a function
// Takes a parser, expression, interval bounds, and tolerance


// ... Assume ExpressionParser is properly defined elsewhere ...

double bisectionMethod(ExpressionParser& parser, const std::string& expression,
    double a, double b, double tol = 1e-6, int maxIter = 100) {
    // Evaluate initial endpoints
    parser.setVariable("x", a);
    double fa = parser.evaluate(expression);
    parser.setVariable("x", b);
    double fb = parser.evaluate(expression);

    // Check if root is at endpoints
    if (std::abs(fa) < tol) return a;
    if (std::abs(fb) < tol) return b;

    // Validate interval signs
    if (fa * fb > 0) {
        throw std::runtime_error("Function must have opposite signs at interval endpoints");
    }

    int iter = 0;
    double c, fc;

    while ((b - a) > tol && iter < maxIter) {
        c = (a + b) / 2;  // Midpoint calculation
        parser.setVariable("x", c);
        fc = parser.evaluate(expression);

        // Early exit if root found
        if (std::abs(fc) < tol) break;

        // Update interval boundaries
        if (fa * fc < 0) {
            b = c;
            fb = fc;
        }
        else {
            a = c;
            fa = fc;
        }

        iter++;
    }

    // Return final midpoint (best approximation)
    return (a + b) / 2.0;
}

// Secant Method for finding roots of a function
// Takes a parser, expression, two initial guesses, and tolerance
double secantMethod(ExpressionParser& parser, const std::string& expression, double x0, double x1, double tol = 1e-6, int maxIter = 100) {
    parser.setVariable("x", x0);
    double f0 = parser.evaluate(expression);

    parser.setVariable("x", x1);
    double f1 = parser.evaluate(expression);

    double x2;
    int iter = 0;

    while (std::abs(f1) > tol && iter < maxIter && std::abs(x1 - x0) > tol) {
        // Check if division by zero might occur
        if (std::abs(f1 - f0) < 1e-10) {
            throw std::runtime_error("Division by near-zero value in secant method");
        }

        // Calculate the next approximation
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0);

        // Update values for next iteration
        x0 = x1;
        f0 = f1;
        x1 = x2;

        // Evaluate function at new point
        parser.setVariable("x", x1);
        f1 = parser.evaluate(expression);

        iter++;
    }

    return x1;
}

// Secant Method with detailed output
// Same as above but returns iteration details
std::vector<std::tuple<int, double, double, double>> secantMethodDetailed(
    ExpressionParser& parser, const std::string& expression, double x0, double x1, double tol = 1e-6, int maxIter = 100) {

    std::vector<std::tuple<int, double, double, double>> iterations;

    parser.setVariable("x", x0);
    double f0 = parser.evaluate(expression);

    parser.setVariable("x", x1);
    double f1 = parser.evaluate(expression);

    double x2;

    for (int iter = 0; iter < maxIter; iter++) {
        // Store iteration details
        iterations.push_back(std::make_tuple(iter, x0, x1, f1));

        // Check if we found the root or reached tolerance
        if (std::abs(f1) < tol || std::abs(x1 - x0) < tol) {
            break;
        }

        // Check if division by zero might occur
        if (std::abs(f1 - f0) < 1e-10) {
            throw std::runtime_error("Division by near-zero value in secant method");
        }

        // Calculate the next approximation
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0);

        // Update values for next iteration
        x0 = x1;
        f0 = f1;
        x1 = x2;

        // Evaluate function at new point
        parser.setVariable("x", x1);
        f1 = parser.evaluate(expression);
    }

    return iterations;
}