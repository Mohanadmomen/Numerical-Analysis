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
    std::string DetectVariable = "x";

    if (expression.find('y') != std::string::npos) {
        DetectVariable = "y";
    }

    // f(a) and f(b) are multiplied by 1/2
    parser.setVariable(DetectVariable, a);
    double fa = parser.evaluate(expression);

    parser.setVariable(DetectVariable, b);
    double fb = parser.evaluate(expression);

    sum = 0.5 * (fa + fb);

    // Sum of f(a+ih) for i=1 to n-1
    for (int i = 1; i < n; i++) {
        parser.setVariable(DetectVariable, a + i * h);
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

    std::string DetectVariable = "x";

    if (expression.find('y') != std::string::npos) {
        DetectVariable = "y";
    }

    // f(a) and f(b)
    parser.setVariable(DetectVariable, a);
    double fa = parser.evaluate(expression);

    parser.setVariable(DetectVariable, b);
    double fb = parser.evaluate(expression);

    sum = fa + fb;

    // Sum of 4*f(a+ih) for i=1,3,5,...,n-1 (odd indices)
    for (int i = 1; i < n; i += 2) {
        parser.setVariable(DetectVariable, a + i * h);
        sum += 4 * parser.evaluate(expression);
    }

    // Sum of 2*f(a+ih) for i=2,4,6,...,n-2 (even indices)
    for (int i = 2; i < n; i += 2) {
        parser.setVariable(DetectVariable, a + i * h);
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

    std::string DetectVariable = "x";

    if (expression.find('y') != std::string::npos) {
        DetectVariable = "y";
    }

    // f(a) and f(b)
    parser.setVariable(DetectVariable, a);
    double fa = parser.evaluate(expression);

    parser.setVariable(DetectVariable, b);
    double fb = parser.evaluate(expression);

    sum = fa + fb;

    // Sum of 3*f(a+ih) for i=1,2,4,5,7,8,...,n-2,n-1 (indices not divisible by 3)
    for (int i = 1; i < n; i++) {
        if (i % 3 != 0) {
            parser.setVariable(DetectVariable, a + i * h);
            sum += 3 * parser.evaluate(expression);
        }
    }

    // Sum of 2*f(a+ih) for i=3,6,9,...,n-3 (indices divisible by 3)
    for (int i = 3; i < n; i += 3) {
        parser.setVariable(DetectVariable, a + i * h);
        sum += 2 * parser.evaluate(expression);
    }

    return (3 * h / 8) * sum;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// 
// 
// 
// 
// Trapezoidal method for numerical integration
double trapezoidalIntegration2D(ExpressionParser& parser, const std::string& expression,
    double a, double b, int n,
    const std::vector<double>& x, const std::vector<double>& y) {
    // Input validation
    if (n <= 0) {
        throw std::runtime_error("Number of intervals must be positive");
    }
    if (a >= b) {
        throw std::runtime_error("Lower bound must be less than upper bound");
    }
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("X and Y vectors must have the same non-zero size");
    }
    // Check if x is actually the range [a,b] divided into n intervals
    // If not, we'll use the provided x vector as is
    double h;
    if (x.size() - 1 == n && std::abs(x[0] - a) < 1e-10 && std::abs(x[n] - b) < 1e-10) {
        h = (b - a) / n; // Use analytical step size if x matches our expectations
    }
    else {
        // Use the vector directly - variable step size
        h = 1.0; // We'll calculate individual step sizes below
    }
   

    double sum = 0.0;

    // First point - weight 1/2
    parser.setVariable("x", x[0]);
    parser.setVariable("y", y[0]);
    double first = parser.evaluate(expression);

    // Last point - weight 1/2
    parser.setVariable("x", x.back());
    parser.setVariable("y", y.back());
    double last = parser.evaluate(expression);

    // For equal step sizes, use the standard trapezoidal formula
    if (h == (b - a) / n) {
        sum = 0.5 * (first + last);

        // Interior points have weight 1
        for (size_t i = 1; i < x.size() - 1; i++) {
            parser.setVariable("x", x[i]);
            parser.setVariable("y", y[i]);
            sum += parser.evaluate(expression);
        }

        return h * sum;
    }
    // For variable step sizes, use the generalized trapezoidal rule
    else {
        sum = 0.0;
        for (size_t i = 1; i < x.size(); i++) {
            double step = x[i] - x[i - 1];

            parser.setVariable("x", x[i - 1]);
            parser.setVariable("y", y[i - 1]);
            double f_prev = parser.evaluate(expression);

            parser.setVariable("x", x[i]);
            parser.setVariable("y", y[i]);
            double f_current = parser.evaluate(expression);

            sum += (f_prev + f_current) * step / 2.0;
        }

        return sum;
    }
}
// Requires n to be an even number
double simpsonOneThird2D(ExpressionParser& parser, const std::string& expression,
    double a, double b, int n,
    const std::vector<double>& x = {}, const std::vector<double>& y = {}) {
    if (n <= 0 || n % 2 != 0) {
        throw std::runtime_error("Number of intervals must be positive and even");
    }
    if (a >= b) {
        throw std::runtime_error("Lower bound must be less than upper bound");
    }
    double h;
    if (x.size() - 1 == n && std::abs(x[0] - a) < 1e-10 && std::abs(x[n] - b) < 1e-10) {
        h = (b - a) / n; // Use analytical step size if x matches our expectations
    }
    else {
        // Use the vector directly - variable step size
        h = 1.0; // We'll calculate individual step sizes below
    }

    double sum = 0.0;

    // f(a) and f(b)
    parser.setVariable("x", x.front());
    parser.setVariable("y", y.front());
    double fa = parser.evaluate(expression);
    parser.setVariable("x", x.back());
    parser.setVariable("y", y.back());
    double fb = parser.evaluate(expression);

    sum = fa + fb;

    // Sum of 4*f(a+ih) for odd i (i=1,3,5,...,n-1)
    for (int i = 1; i < n; i += 2) {
        parser.setVariable("x", x[i]);
        parser.setVariable("y", y[i]);
        sum += 4 * parser.evaluate(expression);
    }

    // Sum of 2*f(a+ih) for even i (i=2,4,6,...,n-2)
    for (int i = 2; i < n; i += 2) {
        parser.setVariable("x", x[i]);
        parser.setVariable("y", y[i]);
        sum += 2 * parser.evaluate(expression);
    }

    return (h / 3) * sum;
}
// Requires n to be divisible by 3
double simpsonThreeEighth2D(ExpressionParser& parser, const std::string& expression,
    double a, double b, int n,
    const std::vector<double>& x = {}, const std::vector<double>& y = {}) {
    if (n <= 0 || n % 3 != 0) {
        throw std::runtime_error("Number of intervals must be positive and divisible by 3");
    }

    if (a >= b) {
        throw std::runtime_error("Lower bound must be less than upper bound");
    }


    double h;
    if (x.size() - 1 == n && std::abs(x[0] - a) < 1e-10 && std::abs(x[n] - b) < 1e-10) {
        h = (b - a) / n; // Use analytical step size if x matches our expectations
    }
    else {
        // Use the vector directly - variable step size
        h = 1.0; // We'll calculate individual step sizes below
    }

    
    double sum = 0.0;

   

    // f(a) and f(b)
    parser.setVariable("x", x.front());
    parser.setVariable("y", y.front());
    double fa = parser.evaluate(expression);

    parser.setVariable("x", x.back());
    parser.setVariable("y", y.back());
    double fb = parser.evaluate(expression);
    //std::cout << n << std::endl;
    sum = fa+fb;

    // Sum of 3*f(a+ih) for i=1,2,4,5,7,8,...,n-2,n-1 (indices not divisible by 3)
    for (int i = 1; i < n; i++) {
            
        if (i % 3 != 0) {
            parser.setVariable("x", x[i]);
            parser.setVariable("y", y[i]);
            sum += 3 * parser.evaluate(expression);
        }
    }

    
    // Sum of 2*f(a+ih) for i=3,6,9,...,n-3 (indices divisible by 3)
    for (int i = 3; i < n; i += 3) {
        parser.setVariable("x", x[i]);
        parser.setVariable("y", y[i]);
        sum += 2 * parser.evaluate(expression);
    }

    return (3 * h / 8) * sum;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
double trapezoidalIntegrationTable(ExpressionParser& parser, const std::string& expression,
    double a, double b, int n, const std::vector<double>& x = {}) {
    // Input validation
    if (n <= 0) {
        throw std::runtime_error("Number of intervals must be positive");
    }


    if (a >= b) {
        throw std::runtime_error("Lower bound must be less than upper bound");
    }

    // Calculate step size for uniform grid
    double h = (b - a) / n;
    
    // Compute trapezoidal sum
    double sum = 0.0;

    std::string DetectVariable="x";

    if (expression.find('y') != std::string::npos) {
        DetectVariable = "y";
    }
  

    //// First and last points have weight 1/2
    parser.setVariable(DetectVariable, x.front());
    double fa = parser.evaluate(expression);
    parser.setVariable(DetectVariable, x.back());
    double fb = parser.evaluate(expression);
   

    sum = 0.5 * (fa + fb);
    //// Interior points have weight 1
    for (size_t i = 1; i < x.size() - 1; i++) {
        parser.setVariable(DetectVariable, x[i]);
        
       
        sum += parser.evaluate(expression);
 
    }

    return h * sum;
}

double simpsonOneThirdTable(ExpressionParser& parser, const std::string& expression,
    double a, double b, int n, const std::vector<double>& x = {}) {
    // Input validation
    if (n <= 0 || n % 2 != 0) {
        throw std::runtime_error("Number of intervals must be positive and even");
    }
    if (a >= b) {
        throw std::runtime_error("Lower bound must be less than upper bound");
    }
    if (x.size() <= 2) {
        throw std::runtime_error("X vector must contain at least 3 points");
    }

    // Calculate step size for uniform grid
    double h = (b - a) / n;

    // Determine if expression uses x or y as the variable
    std::string DetectVariable = "x";
    if (expression.find('y') != std::string::npos) {
        DetectVariable = "y";
    }

    // First and last points
    parser.setVariable(DetectVariable, x.front());
    double fa = parser.evaluate(expression);
    parser.setVariable(DetectVariable, x.back());
    double fb = parser.evaluate(expression);

    double sum = fa + fb;

    // Sum of 4*f(a+ih) for i=1,3,5,...,n-1 (odd indices)
    for (size_t i = 1; i < x.size() - 1; i += 2) {
        parser.setVariable(DetectVariable, x[i]);
        sum += 4 * parser.evaluate(expression);
    }

    // Sum of 2*f(a+ih) for i=2,4,6,...,n-2 (even indices)
    for (size_t i = 2; i < x.size() - 1; i += 2) {
        parser.setVariable(DetectVariable, x[i]);
        sum += 2 * parser.evaluate(expression);
    }

    return (h / 3) * sum;
}

double simpsonThreeEighthTable(ExpressionParser& parser, const std::string& expression,
    double a, double b, int n, const std::vector<double>& x = {}) {
    // Input validation
    if (n <= 0 || n % 3 != 0) {
        throw std::runtime_error("Number of intervals must be positive and divisible by 3");
    }
    if (a >= b) {
        throw std::runtime_error("Lower bound must be less than upper bound");
    }
    if (x.size() <= 3) {
        throw std::runtime_error("X vector must contain at least 4 points");
    }

    // Calculate step size for uniform grid
    double h = (b - a) / n;

    // Determine if expression uses x or y as the variable
    std::string DetectVariable = "x";
    if (expression.find('y') != std::string::npos) {
        DetectVariable = "y";
    }

    // First and last points
    parser.setVariable(DetectVariable, x.front());
    double fa = parser.evaluate(expression);
    parser.setVariable(DetectVariable, x.back());
    double fb = parser.evaluate(expression);

    double sum = fa + fb;

    // Sum of 3*f(a+ih) for i=1,2,4,5,7,8,...,n-2,n-1 (indices not divisible by 3)
    for (size_t i = 1; i < x.size() - 1; i++) {
        if (i % 3 != 0) {
            parser.setVariable(DetectVariable, x[i]);
            sum += 3 * parser.evaluate(expression);
        }
    }

    // Sum of 2*f(a+ih) for i=3,6,9,...,n-3 (indices divisible by 3)
    for (size_t i = 3; i < x.size() - 1; i += 3) {
        parser.setVariable(DetectVariable, x[i]);
        sum += 2 * parser.evaluate(expression);
    }

    return (3 * h / 8) * sum;
}