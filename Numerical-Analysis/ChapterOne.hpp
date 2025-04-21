#include "FunctionHandler.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <limits>

// Forward declaration (assuming ExpressionParser is defined elsewhere)
class ExpressionParser;

// Class containing various numerical root-finding methods
class FunctionHandler {
private:
    // Internal function to check if a value is close to zero
    static bool isNearZero(double val, double epsilon = 1e-10) {
        return std::abs(val) < epsilon;
    }

public:
    // Bisection Method for finding roots of a function
    static double bisectionMethod(ExpressionParser& parser, const std::string& expression,
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
    static double secantMethod(ExpressionParser& parser, const std::string& expression,
        double x0, double x1, double tol = 1e-6, int maxIter = 100) {

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
    static std::vector<std::tuple<int, double, double, double>> secantMethodDetailed(
        ExpressionParser& parser, const std::string& expression,
        double x0, double x1, double tol = 1e-6, int maxIter = 100) {

        std::vector<std::tuple<int, double, double, double>> iterations;

        parser.setVariable("x", x0);
        double f0 = parser.evaluate(expression);
        parser.setVariable("x", x1);
        double f1 = parser.evaluate(expression);

        double x2;

        for (int iter = 0; iter < maxIter; iter++) {
            // Store iteration details
            iterations.push_back(std::make_tuple(iter, x1, f1, x1 - x0));

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

    // Newton-Raphson Method
    static double newtonMethod(ExpressionParser& parser,
        const std::string& expression, const std::string& derivativeExpression,
        double x0, double tol = 1e-6, int maxIter = 100) {

        double x = x0;

        for (int iter = 0; iter < maxIter; iter++) {
            // Evaluate function and its derivative at current point
            parser.setVariable("x", x);
            double fx = parser.evaluate(expression);
            double dfx = parser.evaluate(derivativeExpression);

            // Check if we're close enough to the root
            if (std::abs(fx) < tol) {
                return x;
            }

            // Check if derivative is too close to zero
            if (isNearZero(dfx)) {
                throw std::runtime_error("Derivative too close to zero in Newton's method");
            }

            // Calculate next approximation
            double xNext = x - fx / dfx;

            // Check if we've converged based on change in x
            if (std::abs(xNext - x) < tol) {
                return xNext;
            }

            x = xNext;
        }

        return x;  // Return best approximation after maxIter iterations
    }

    // False Position (Regula Falsi) Method
    static double falsePositionMethod(ExpressionParser& parser, const std::string& expression,
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

        double c, fc;

        for (int iter = 0; iter < maxIter; iter++) {
            // Calculate the false position point
            c = (a * fb - b * fa) / (fb - fa);

            // Evaluate function at the new point
            parser.setVariable("x", c);
            fc = parser.evaluate(expression);

            // Check if root is found
            if (std::abs(fc) < tol || std::abs(b - a) < tol) {
                break;
            }

            // Update interval boundaries
            if (fa * fc < 0) {
                b = c;
                fb = fc;
            }
            else {
                a = c;
                fa = fc;
            }
        }

        return c;
    }

    // Brent's Method - a combination of bisection, secant, and inverse quadratic interpolation
    static double brentMethod(ExpressionParser& parser, const std::string& expression,
        double a, double b, double tol = 1e-6, int maxIter = 100) {

        // Evaluate endpoints
        parser.setVariable("x", a);
        double fa = parser.evaluate(expression);
        parser.setVariable("x", b);
        double fb = parser.evaluate(expression);

        // Check if root is at endpoints
        if (std::abs(fa) < tol) return a;
        if (std::abs(fb) < tol) return b;

        // Ensure f(a) and f(b) have opposite signs
        if (fa * fb > 0) {
            throw std::runtime_error("Function must have opposite signs at interval endpoints");
        }

        // Ensure |f(b)| < |f(a)|
        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }

        double c = a;
        double fc = fa;
        double d = 0;  // Not used in first iteration
        double s = 0;
        bool mflag = true;

        for (int iter = 0; iter < maxIter; iter++) {
            // Check if f(a) and f(c) have opposite signs, and |f(c)| < |f(b)|
            if (std::abs(fc) < std::abs(fb)) {
                a = b; fa = fb;
                b = c; fb = fc;
                c = a; fc = fa;
            }

            double tol1 = 2 * std::numeric_limits<double>::epsilon() * std::abs(b) + 0.5 * tol;
            double xm = 0.5 * (c - b);

            // Check for convergence
            if (std::abs(xm) <= tol1 || isNearZero(fb)) {
                return b;  // Found the root
            }

            // Decide which method to use
            if (mflag && std::abs(s) > tol1 && std::abs(fb) < std::abs(fa)) {
                // Try inverse quadratic interpolation
                double r = fb / fa;
                double s2 = fc / fa;
                double p, q;

                if (isNearZero(a - c)) {
                    // Use secant method if a and c are too close
                    p = 2 * xm * r;
                    q = 1 - r;
                }
                else {
                    // Use inverse quadratic interpolation
                    q = fa / fc;
                    r = fb / fc;
                    s = fa / fb;
                    p = s * (2 * xm * q * (q - r) - (b - a) * (r - 1));
                    q = (q - 1) * (r - 1) * (s - 1);
                }

                // Adjust signs
                if (p > 0) q = -q;
                else p = -p;

                // Check if interpolation is acceptable
                if (2 * p < 3 * xm * q - std::abs(tol1 * q) && p < std::abs(0.5 * s * q)) {
                    s = p / q;
                    mflag = false;
                }
                else {
                    // Use bisection
                    s = xm;
                    mflag = true;
                }
            }
            else {
                // Use bisection
                s = xm;
                mflag = true;
            }

            // Update a, b, c
            a = b;
            fa = fb;

            if (std::abs(s) > tol1) {
                b += s;
            }
            else {
                b += (xm > 0 ? tol1 : -tol1);
            }

            // Evaluate function at new point
            parser.setVariable("x", b);
            fb = parser.evaluate(expression);

            // Adjust the interval
            if (fb * fc > 0) {
                c = a;
                fc = fa;
            }
        }

        return b;  // Return best approximation after maxIter
    }

    // Helper function to numerically calculate the derivative at a point
    static double numericalDerivative(ExpressionParser& parser, const std::string& expression,
        double x, double h = 1e-6) {

        parser.setVariable("x", x + h);
        double fxph = parser.evaluate(expression);

        parser.setVariable("x", x - h);
        double fxmh = parser.evaluate(expression);

        return (fxph - fxmh) / (2 * h);
    }

    // Newton method using numerical differentiation (no need for derivative expression)
    static double newtonMethodNumerical(ExpressionParser& parser, const std::string& expression,
        double x0, double tol = 1e-6, int maxIter = 100, double h = 1e-6) {

        double x = x0;

        for (int iter = 0; iter < maxIter; iter++) {
            // Evaluate function at current point
            parser.setVariable("x", x);
            double fx = parser.evaluate(expression);

            // Check if we're close enough to the root
            if (std::abs(fx) < tol) {
                return x;
            }

            // Calculate derivative numerically
            double dfx = numericalDerivative(parser, expression, x, h);

            // Check if derivative is too close to zero
            if (isNearZero(dfx)) {
                throw std::runtime_error("Derivative too close to zero in Newton's method");
            }

            // Calculate next approximation
            double xNext = x - fx / dfx;

            // Check if we've converged based on change in x
            if (std::abs(xNext - x) < tol) {
                return xNext;
            }

            x = xNext;
        }

        return x;  // Return best approximation after maxIter iterations
    }

    // Function to analyze the function behavior in an interval
    // Returns min value, max value, and count of sign changes
    static std::tuple<double, double, int> analyzeFunction(
        ExpressionParser& parser, const std::string& expression,
        double a, double b, int numPoints = 100) {

        double min_val = std::numeric_limits<double>::max();
        double max_val = std::numeric_limits<double>::lowest();
        int sign_changes = 0;
        double prev_val = 0;
        bool first_point = true;

        double step = (b - a) / numPoints;

        for (int i = 0; i <= numPoints; i++) {
            double x = a + i * step;
            parser.setVariable("x", x);
            double fx = parser.evaluate(expression);

            // Update min and max
            min_val = std::min(min_val, fx);
            max_val = std::max(max_val, fx);

            // Count sign changes
            if (!first_point) {
                if ((prev_val < 0 && fx >= 0) || (prev_val >= 0 && fx < 0)) {
                    sign_changes++;
                }
            }

            prev_val = fx;
            first_point = false;
        }

        return std::make_tuple(min_val, max_val, sign_changes);
    }

    // Function to find all roots in an interval using subdivision and root-finding
    static std::vector<double> findAllRoots(
        ExpressionParser& parser, const std::string& expression,
        double a, double b, int numSubintervals = 20, double tol = 1e-6) {

        std::vector<double> roots;
        double step = (b - a) / numSubintervals;

        for (int i = 0; i < numSubintervals; i++) {
            double subA = a + i * step;
            double subB = a + (i + 1) * step;

            parser.setVariable("x", subA);
            double fa = parser.evaluate(expression);
            parser.setVariable("x", subB);
            double fb = parser.evaluate(expression);

            // Check for root in this subinterval (sign change)
            if (fa * fb <= 0) {
                try {
                    // Use Brent's method for reliable root finding
                    double root = brentMethod(parser, expression, subA, subB, tol);

                    // Verify it's actually a root
                    parser.setVariable("x", root);
                    double froot = parser.evaluate(expression);

                    if (std::abs(froot) < tol) {
                        // Check if this root is distinct from previously found roots
                        bool distinct = true;
                        for (double r : roots) {
                            if (std::abs(r - root) < tol) {
                                distinct = false;
                                break;
                            }
                        }

                        if (distinct) {
                            roots.push_back(root);
                        }
                    }
                }
                catch (const std::exception& e) {
                    // Skip this subinterval if there's an error
                    continue;
                }
            }
        }

        return roots;
    }
};



