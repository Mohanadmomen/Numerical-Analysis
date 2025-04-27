#include <iostream>
#include <string>
#include <cmath>
#include "FunctionHandler.hpp"
#include <stdexcept>

#pragma once 




std::pair<double, double> findRootInterval(ExpressionParser& parser, const std::string& expression) {
   
    auto evaluate = [&](double x) {
        parser.setVariable("x", x);
        return parser.evaluate(expression);
        };

    double a = 0.00000001;
    double b = 1;

    double prevvalue;
    double nextvalue;
    while (true) {
        prevvalue = evaluate(a);
        nextvalue = evaluate(b);
        if (nextvalue >= 0 && prevvalue < 0)
            break;
        a++;
        b++;
    }


    return { a, b };
}
// Main bisection method that takes external function to find interval
double bisectionMethod(ExpressionParser& parser, const std::string& expression,
    double tol = 1e-6, int maxIter = 100) {

    // Find suitable interval containing a root
    auto [a, b] = findRootInterval(parser, expression);

    // Evaluate at endpoints
    parser.setVariable("x", a);
    double fa = parser.evaluate(expression);
    parser.setVariable("x", b);
    double fb = parser.evaluate(expression);

    // Check if root is at endpoints
    if (std::abs(fa) < tol) return a;
    if (std::abs(fb) < tol) return b;

    // Standard bisection method
    int iter = 0;
    double c, fc;
    while ((b - a) > tol && iter < maxIter) {
        c = (a + b) / 2.0;  // Midpoint calculation
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

// Secant method implementation for finding roots of equations
double secantMethod(ExpressionParser& parser, const std::string& expression,
    double tol = 1e-6, int maxIter = 100) {
    // Find suitable interval containing a root
    auto [a, b] = findRootInterval(parser, expression);

    // Initialize with the interval boundaries
    double x_prev = 2;  // x₀
    double x_curr = 3;  // x₁

    // Evaluate function at initial points
    parser.setVariable("x", x_prev);
    double f_prev = parser.evaluate(expression);
    
    parser.setVariable("x", x_curr);
    double f_curr = parser.evaluate(expression);

    // Check if either initial point is already a root
    if (std::abs(f_prev) < tol) return x_prev;
    if (std::abs(f_curr) < tol) return x_curr;

    // Iterate using secant method
    int iter = 0;
    double x_next, f_next;

    while (std::abs(f_curr) > tol && iter < maxIter) {
        // Calculate next approximation using secant formula
        x_next = x_curr - (f_curr * (x_curr - x_prev)) / (f_curr - f_prev);

        // Evaluate function at new point
        parser.setVariable("x", x_next);
        f_next = parser.evaluate(expression);

        // Early exit if root found
        if (std::abs(f_next) < tol) break;

        // Update values for next iteration
        x_prev = x_curr;
        f_prev = f_curr;
        x_curr = x_next;
        f_curr = f_next;

        iter++;
    }

    // Return final approximation
    return x_curr;
}



