#include "Euler.h"

// Euler's Method - returns final y value
double euler(ExpressionParser& parser, const std::string& odeExpression, double x0, double y0, double h, double xEnd) {
    double x = x0;
    double y = y0;

    while (x < xEnd - 1e-10) {
        parser.setVariable("x", x);
        parser.setVariable("y", y);

        double slope = parser.evaluate(odeExpression);

        y = y + h * slope;
        x = x + h;
    }

    return y; // only final result
}

// Modified Euler's Method (Heun's Method) - returns final y value
double modifiedEuler(ExpressionParser& parser, const std::string& odeExpression, double x0, double y0, double h, double xEnd) {
    double x = x0;
    double y = y0;

    while (x < xEnd - 1e-10) {
        parser.setVariable("x", x);
        parser.setVariable("y", y);
        double k1 = parser.evaluate(odeExpression);

        double y_pred = y + h * k1;

        parser.setVariable("x", x + h);
        parser.setVariable("y", y_pred);
        double k2 = parser.evaluate(odeExpression);

        y = y + h * (k1 + k2) / 2;
        x = x + h;
    }

    return y; // only final result
}
