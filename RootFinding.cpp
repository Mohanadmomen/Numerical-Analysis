#include "rootfinding.h"
#include <stdexcept>

std::pair<double, double> findRootInterval(ExpressionParser& parser, const std::string& expression) {

    auto evaluate = [&](double x) {
        parser.setVariable("x", x);

        return parser.evaluate(expression);
    };

    double a = 0, b = 1, step = 1;

    double prevvalue, nextvalue;

    bool found = false;

    for (int i = 0; i < 10000; i++)
    {
        prevvalue = evaluate(a);
        nextvalue = evaluate(b);

        if (parser.isError(prevvalue) || parser.isError(nextvalue)) {
            a += step;
            b += step;
            continue;
        }
        if (nextvalue >= 0 && prevvalue < 0 || prevvalue > 0 && nextvalue <= 0) {

            found = true;
            break;
        }
        a += step;
        b += step;

    }

    if (found)
        return { floor(a), ceil(b) };
    else {
        //std::cout << "No positive interval is found in this range/n"<<std::endl;
        return { std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN() };
    }
}

// Main bisection method that takes external function to find interval
double bisectionMethod(ExpressionParser& parser, const std::string& expression,
                       double tol, int maxIter) {

    // Find suitable interval containing a root
    auto [a, b] = findRootInterval(parser, expression);
    if (std::isnan(a) || std::isnan(b))
        return std::numeric_limits<double>::quiet_NaN();
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
                    double tol , int maxIter ) {
    // Find suitable interval containing a root
    auto [a, b] = findRootInterval(parser, expression);

    if (std::isnan(a) || std::isnan(b))
        return std::numeric_limits<double>::quiet_NaN();

    // Initialize with the interval boundaries
    double x_prev = a;  // x₀
    double x_curr = b;  // x₁

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
        // Prevent division by zero
        if (std::abs(f_curr - f_prev) < 1e-15) {
            // If function values are too close, slightly perturb one of the points
            x_curr += 0.01;
            parser.setVariable("x", x_curr);
            f_curr = parser.evaluate(expression);
            continue;
        }

        // Calculate next approximation using secant formula
        x_next = x_curr - (f_curr * (x_curr - x_prev)) / (f_curr - f_prev);

        // If x_next is outside our initial interval, reset to midpoint
        if (x_next < a - 100 || x_next > b + 100) {
            x_next = (a + b) / 2;
        }

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


//############################### Newton Raphson ###########################################
double differentiate(ExpressionParser& parser, const std::string& expression, double x) {
    const double dx = 1e-6;
    // Evaluate f(x + dx)
    parser.setVariable("x", x + dx);
    double FxPlusdx = parser.evaluate(expression);
    // Evaluate f(x)
    parser.setVariable("x", x);
    double Fx = parser.evaluate(expression);
    // Return the derivative
    return (FxPlusdx - Fx) / dx;
}

double newton_raphson(ExpressionParser& parser, const std::string& expression,
                      double tol, int maxIter) {
    // Find suitable interval containing a root
    auto [a, b] = findRootInterval(parser, expression);
    if (std::isnan(a) || std::isnan(b))
        return std::numeric_limits<double>::quiet_NaN();

    // Evaluate at endpoints
    parser.setVariable("x", a);
    double fa = parser.evaluate(expression);
    parser.setVariable("x", b);
    double fb = parser.evaluate(expression);

    // Check if root is at endpoints
    if (std::abs(fa) < tol) return a;
    if (std::abs(fb) < tol) return b;

    // Initial guess (using midpoint as a start point)
    double x0 = (a + b) / 2.0;
    parser.setVariable("x", x0);  // Set x value before evaluating
    double f0 = parser.evaluate(expression);

    // Standard Newton-Raphson method
    for (int iter = 0; iter < maxIter; iter++) {
        // Compute the derivative using differentiate method
        double derivative = differentiate(parser, expression, x0);

        if (std::abs(derivative) <= 1e-9) {
            //std::cerr << "Division by zero detected in derivative" << std::endl;
            return std::numeric_limits<double>::quiet_NaN();
        }

        // Update the current guess
        double x1 = x0 - f0 / derivative;

        // Early exit if the result is good enough
        if (std::abs(x1 - x0) < tol) {
            return x1;
        }

        // Update for the next iteration
        x0 = x1;
        parser.setVariable("x", x0);  // Set x value before evaluating
        f0 = parser.evaluate(expression);
    }

    // If the root wasn't found within the max iterations, return the last approximation
    return x0;

}
