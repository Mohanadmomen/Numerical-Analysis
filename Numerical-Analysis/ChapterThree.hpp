#include <iostream>
#include <string>
#include <cmath>
#include "FunctionHandler.hpp"

#pragma once 


// Trapezoidal method for numerical integration
double trapezoidalIntegration(ExpressionParser& parser, const std::string& expression, double a, double b, int n) {
    if (n <= 0) {
        throw std::runtime_error("Number of intervals must be positive");
    }

    if (a >= b) {
        throw std::runtime_error("Lower bound must be less than upper bound");
    }

    double h = (b - a) / n;
    double sum = 0.0;

    // f(a) and f(b) are multiplied by 1/2
    parser.setVariable("x", a);
    double fa = parser.evaluate(expression);

    parser.setVariable("x", b);
    double fb = parser.evaluate(expression);

    sum = 0.5 * (fa + fb);

    // Sum of f(a+ih) for i=1 to n-1
    for (int i = 1; i < n; i++) {
        parser.setVariable("x", a + i * h);
        sum += parser.evaluate(expression);
    }

    return h * sum;
}
// Simpson's 1/3 rule for numerical integration
// Requires n to be even
double simpsonOneThird(ExpressionParser& parser, const std::string& expression, double a, double b, int n) {
    if (n <= 0 || n % 2 != 0) {
        throw std::runtime_error("Number of intervals must be positive and even");
    }

    if (a >= b) {
        throw std::runtime_error("Lower bound must be less than upper bound");
    }

    double h = (b - a) / n;
    double sum = 0.0;

    // f(a) and f(b)
    parser.setVariable("x", a);
    double fa = parser.evaluate(expression);

    parser.setVariable("x", b);
    double fb = parser.evaluate(expression);

    sum = fa + fb;

    // Sum of 4*f(a+ih) for i=1,3,5,...,n-1 (odd indices)
    for (int i = 1; i < n; i += 2) {
        parser.setVariable("x", a + i * h);
        sum += 4 * parser.evaluate(expression);
    }

    // Sum of 2*f(a+ih) for i=2,4,6,...,n-2 (even indices)
    for (int i = 2; i < n; i += 2) {
        parser.setVariable("x", a + i * h);
        sum += 2 * parser.evaluate(expression);
    }

    return (h / 3) * sum;
}

// Simpson's 3/8 rule for numerical integration
// Requires n to be divisible by 3
double simpsonThreeEighth(ExpressionParser& parser, const std::string& expression, double a, double b, int n) {
    if (n <= 0 || n % 3 != 0) {
        throw std::runtime_error("Number of intervals must be positive and divisible by 3");
    }

    if (a >= b) {
        throw std::runtime_error("Lower bound must be less than upper bound");
    }

    double h = (b - a) / n;
    double sum = 0.0;

    // f(a) and f(b)
    parser.setVariable("x", a);
    double fa = parser.evaluate(expression);

    parser.setVariable("x", b);
    double fb = parser.evaluate(expression);

    sum = fa + fb;

    // Sum of 3*f(a+ih) for i=1,2,4,5,7,8,...,n-2,n-1 (indices not divisible by 3)
    for (int i = 1; i < n; i++) {
        if (i % 3 != 0) {
            parser.setVariable("x", a + i * h);
            sum += 3 * parser.evaluate(expression);
        }
    }

    // Sum of 2*f(a+ih) for i=3,6,9,...,n-3 (indices divisible by 3)
    for (int i = 3; i < n; i += 3) {
        parser.setVariable("x", a + i * h);
        sum += 2 * parser.evaluate(expression);
    }

    return (3 * h / 8) * sum;
}