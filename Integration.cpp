#include "integration.h"
#include <cmath>
#include <stdexcept>

double formatResult(double result) {
    if (std::fabs(result) < EPSILON) {
        result = 0.0;
    }
    result = std::round(result * 1000000.0) / 1000000.0;
    return result;
}

// 1D Integration Methods
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

    parser.setVariable(DetectVariable, a);
    double fa = parser.evaluate(expression);
    parser.setVariable(DetectVariable, b);
    double fb = parser.evaluate(expression);

    sum = 0.5 * (fa + fb);
    for (int i = 1; i < n; i++) {
        parser.setVariable(DetectVariable, a + i * h);
        sum += parser.evaluate(expression);
    }

    return formatResult(h * sum);
}

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

    parser.setVariable(DetectVariable, a);
    double fa = parser.evaluate(expression);
    parser.setVariable(DetectVariable, b);
    double fb = parser.evaluate(expression);

    sum = fa + fb;
    for (int i = 1; i < n; i += 2) {
        parser.setVariable(DetectVariable, a + i * h);
        sum += 4 * parser.evaluate(expression);
    }
    for (int i = 2; i < n; i += 2) {
        parser.setVariable(DetectVariable, a + i * h);
        sum += 2 * parser.evaluate(expression);
    }

    return formatResult((h / 3) * sum);
}

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

    parser.setVariable(DetectVariable, a);
    double fa = parser.evaluate(expression);
    parser.setVariable(DetectVariable, b);
    double fb = parser.evaluate(expression);

    sum = fa + fb;
    for (int i = 1; i < n; i++) {
        if (i % 3 != 0) {
            parser.setVariable(DetectVariable, a + i * h);
            sum += 3 * parser.evaluate(expression);
        }
    }
    for (int i = 3; i < n; i += 3) {
        parser.setVariable(DetectVariable, a + i * h);
        sum += 2 * parser.evaluate(expression);
    }

    return formatResult((3 * h / 8) * sum);
}

// 2D Integration Methods
double trapezoidalIntegration2D(ExpressionParser& parser, const std::string& expression,
                                double a, double b, int n,
                                const std::vector<double>& x, const std::vector<double>& y) {
    if (n <= 0) {
        throw std::runtime_error("Number of intervals must be positive");
    }
    if (a >= b) {
        throw std::runtime_error("Lower bound must be less than upper bound");
    }
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("X and Y vectors must have the same non-zero size");
    }

    double h;

    h = (b - a) / n;

    double sum = 0.0;
    parser.setVariable("x", x[0]);
    parser.setVariable("y", y[0]);
    double first = parser.evaluate(expression);
    parser.setVariable("x", x.back());
    parser.setVariable("y", y.back());
    double last = parser.evaluate(expression);

    double result;
    if (h == (b - a) / n) {
        sum = 0.5 * (first + last);
        for (size_t i = 1; i < x.size() - 1; i++) {
            parser.setVariable("x", x[i]);
            parser.setVariable("y", y[i]);
            sum += parser.evaluate(expression);
        }
        result = h * sum;
    }
    else {
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
        result = sum;
    }

    return formatResult(result);
}

double simpsonOneThird2D(ExpressionParser& parser, const std::string& expression,
                         double a, double b, int n,
                         const std::vector<double>& x, const std::vector<double>& y) {
    if (n <= 0 || n % 2 != 0) {
        throw std::runtime_error("Number of intervals must be positive and even");
    }
    if (a >= b) {
        throw std::runtime_error("Lower bound must be less than upper bound");
    }

    double h;

    h = (b - a) / n;

    parser.setVariable("x", x.front());
    parser.setVariable("y", y.front());
    double fa = parser.evaluate(expression);
    parser.setVariable("x", x.back());
    parser.setVariable("y", y.back());
    double fb = parser.evaluate(expression);

    double sum = fa + fb;
    for (int i = 1; i < n; i += 2) {
        parser.setVariable("x", x[i]);
        parser.setVariable("y", y[i]);
        sum += 4 * parser.evaluate(expression);
    }
    for (int i = 2; i < n; i += 2) {
        parser.setVariable("x", x[i]);
        parser.setVariable("y", y[i]);
        sum += 2 * parser.evaluate(expression);
    }

    return formatResult((h / 3) * sum);
}

double simpsonThreeEighth2D(ExpressionParser& parser, const std::string& expression,
                            double a, double b, int n,
                            const std::vector<double>& x, const std::vector<double>& y) {
    if (n <= 0 || n % 3 != 0) {
        throw std::runtime_error("Number of intervals must be positive and divisible by 3");
    }
    if (a >= b) {
        throw std::runtime_error("Lower bound must be less than upper bound");
    }

    double h;

    h = (b - a) / n;

    parser.setVariable("x", x.front());
    parser.setVariable("y", y.front());
    double fa = parser.evaluate(expression);
    parser.setVariable("x", x.back());
    parser.setVariable("y", y.back());
    double fb = parser.evaluate(expression);

    double sum = fa + fb;
    for (int i = 1; i < n; i++) {
        if (i % 3 != 0) {
            parser.setVariable("x", x[i]);
            parser.setVariable("y", y[i]);
            sum += 3 * parser.evaluate(expression);
        }
    }
    for (int i = 3; i < n; i += 3) {
        parser.setVariable("x", x[i]);
        parser.setVariable("y", y[i]);
        sum += 2 * parser.evaluate(expression);
    }

    return formatResult((3 * h / 8) * sum);
}

// Table-based Integration Methods
double trapezoidalIntegrationTable(ExpressionParser& parser, const std::string& expression,
                                   double a, double b, int n,
                                   const std::vector<double>& x) {
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

    parser.setVariable(DetectVariable, x.front());
    double fa = parser.evaluate(expression);
    parser.setVariable(DetectVariable, x.back());
    double fb = parser.evaluate(expression);

    sum = 0.5 * (fa + fb);
    for (size_t i = 1; i < x.size() - 1; i++) {
        parser.setVariable(DetectVariable, x[i]);
        sum += parser.evaluate(expression);
    }

    return formatResult(h * sum);
}

double simpsonOneThirdTable(ExpressionParser& parser, const std::string& expression,
                            double a, double b, int n,
                            const std::vector<double>& x) {
    if (n <= 0 || n % 2 != 0) {
        throw std::runtime_error("Number of intervals must be positive and even");
    }
    if (a >= b) {
        throw std::runtime_error("Lower bound must be less than upper bound");
    }
    if (x.size() <= 2) {
        throw std::runtime_error("X vector must contain at least 3 points");
    }

    double h = (b - a) / n;
    std::string DetectVariable = "x";
    if (expression.find('y') != std::string::npos) {
        DetectVariable = "y";
    }

    parser.setVariable(DetectVariable, x.front());
    double fa = parser.evaluate(expression);
    parser.setVariable(DetectVariable, x.back());
    double fb = parser.evaluate(expression);

    double sum = fa + fb;
    for (size_t i = 1; i < x.size() - 1; i += 2) {
        parser.setVariable(DetectVariable, x[i]);
        sum += 4 * parser.evaluate(expression);
    }
    for (size_t i = 2; i < x.size() - 1; i += 2) {
        parser.setVariable(DetectVariable, x[i]);
        sum += 2 * parser.evaluate(expression);
    }

    return formatResult((h / 3) * sum);
}

double simpsonThreeEighthTable(ExpressionParser& parser, const std::string& expression,
                               double a, double b, int n,
                               const std::vector<double>& x) {
    if (n <= 0 || n % 3 != 0) {
        throw std::runtime_error("Number of intervals must be positive and divisible by 3");
    }
    if (a >= b) {
        throw std::runtime_error("Lower bound must be less than upper bound");
    }
    if (x.size() <= 3) {
        throw std::runtime_error("X vector must contain at least 4 points");
    }

    double h = (b - a) / n;
    std::string DetectVariable = "x";
    if (expression.find('y') != std::string::npos) {
        DetectVariable = "y";
    }

    parser.setVariable(DetectVariable, x.front());
    double fa = parser.evaluate(expression);
    parser.setVariable(DetectVariable, x.back());
    double fb = parser.evaluate(expression);

    double sum = fa + fb;
    for (size_t i = 1; i < x.size() - 1; i++) {
        if (i % 3 != 0) {
            parser.setVariable(DetectVariable, x[i]);
            sum += 3 * parser.evaluate(expression);
        }
    }
    for (size_t i = 3; i < x.size() - 1; i += 3) {
        parser.setVariable(DetectVariable, x[i]);
        sum += 2 * parser.evaluate(expression);
    }

    return formatResult((3 * h / 8) * sum);
}
