#include <iostream>
#include <cmath>
#include "FunctionHandler.hpp"
#define P 50

using namespace std;

ExpressionParser parser;

pair<double, double> CurveFittingLinear(int n);
pair<double, double> CurveFittingExponential(int n);
pair<double, double> CurveFittingLogarithmic(int n);
pair<double, double> CurveFittingPower(int n);
tuple<double, double, double> CurveFittingQuadratic(int n);

/**
 * @brief Prompts the user to choose a curve fitting model and fits the model to user-provided data.
 * 
 * Supports the following models:
 *   [1] Linear
 *   [2] Exponential
 *   [3] Logarithmic
 *   [4] Power
 *   [5] Exit
 * 
 * It then prints the best-fit equation for the selected model.
 */
void CurveFitting() {
    int ModelNum, n;
    double a, b, c;
    do {
        cout << "\n================================\n";
        cout << "[1] Linear       (a * x + b)\n";
        cout << "[2] Exponential  (a * exp(b * x))\n";
        cout << "[3] Logarithmic  (a * log(x) + b)\n";
        cout << "[4] Power        (a * x^b)\n";
        cout << "[5] Quadratic    (a * x^2 + b * x + c)\n";
        cout << "[6] Exit\n";
        cout << "================================\n\n";
        cout << "> Choose the model number: ";
        cin >> ModelNum;
    } while (ModelNum < 0 || ModelNum > 6);
    
    if (ModelNum != 6) {
        cout << "> How many data points? ";
        cin >> n;
    }

    cin.ignore(); // ignore newline after reading n

    switch(ModelNum) {
        case 1:
            tie(a, b) = CurveFittingLinear(n);
            cout << "Fitted Linear: y = " << a << "x + " << b << "\n";
            break;
        case 2:
            tie(a, b) = CurveFittingExponential(n);
            cout << "Fitted Exponential: y = " << a << " * exp(" << b << "x)\n";
            break;
        case 3:
            tie(a, b) = CurveFittingLogarithmic(n);
            cout << "Fitted Logarithmic: y = " << a << " * ln(x) + " << b << "\n";
            break;
        case 4:
            tie(a, b) = CurveFittingPower(n);
            cout << "Fitted Power: y = " << a << " * x^" << b << "\n";
            break;
        case 5:
            tie(a, b, c) = CurveFittingQuadratic(n);
            cout << "Fitted Quadratic: y = " << a << " * x^2 + " << b << " * x + " << c << "\n";
            break;
        case 6:
            return;
        default:
            cout << "Invalid model number.\n";
    }

}

/**
 * @brief Performs linear curve fitting on a set of (x, y) data points.
 * 
 * Uses the least squares method to compute the best-fit line in the form:
 *     y = a * x + b
 * 
 * @param n The number of data points to be entered.
 * @return pair<double, double> The coefficients (a, b):
 *         - a: slope of the line
 *         - b: y-intercept
 */
pair<double, double> CurveFittingLinear(int n) {
    double sumX = 0, sumX2 = 0, sumY = 0, sumXY = 0, a, b, x[P], y[P];

    for (int i = 0; i < n; i++) {
        string x_expr, y_expr;

        cout << "x[" << i << "] = ";
        getline(cin, x_expr);
        x[i] = parser.evaluate(x_expr);

        cout << "y[" << i << "] = ";
        getline(cin, y_expr);
        y[i] = parser.evaluate(y_expr);    
    }

    // Calculating summation of X, X^2, Y, XY
    for (int i = 0; i < n; i++) {
        sumX += x[i];
        sumY += y[i];
        sumX2 += x[i] * x[i];
        sumXY += x[i] * y[i];
    }

    // Calculating a and b
    a = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    b = (sumY - a * sumX) / n;

    return {a, b};
}

/**
 * @brief Performs exponential curve fitting on a set of (x, y) data points.
 * 
 * Fits the data to the model:
 *     y = a * exp(b * x)
 * 
 * The data is linearized by taking the natural log of y before applying least squares.
 * 
 * @param n The number of data points to be entered.
 * @return pair<double, double> The coefficients (a, b):
 *         - a: scaling factor
 *         - b: exponential rate
 */
pair<double, double> CurveFittingExponential(int n) {
    double sumX = 0, sumX2 = 0, sumY = 0, sumXY = 0, a, b, ln_a, x[P], y[P];

    for (int i = 0; i < n; i++) {
        string x_expr, y_expr;

        cout << "x[" << i << "] = ";
        getline(cin, x_expr);
        x[i] = parser.evaluate(x_expr);
    
        cout << "y[" << i << "] = ";
        getline(cin, y_expr);
        y[i] = parser.evaluate(y_expr);
        y[i] = log(y[i]); // transform y to ln(y)
    }

    for (int i = 0; i < n; i++) {
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumX2 += x[i] * x[i];
    }

    b = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    ln_a = (sumY - b * sumX) / n;
    a = exp(ln_a);

    return {a, b};
}

/**
 * @brief Performs logarithmic curve fitting on a set of (x, y) data points.
 * 
 * Fits the data to the model:
 *     y = a * ln(x) + b
 * 
 * Assumes x > 0 for all data points.
 * 
 * @param n The number of data points to be entered.
 * @return pair<double, double> The coefficients (a, b):
 *         - a: logarithmic scaling factor
 *         - b: vertical shift
 */
pair<double, double> CurveFittingLogarithmic(int n) {
    double sumX = 0, sumX2 = 0, sumY = 0, sumXY = 0, a, b, x[P], y[P];

    for (int i = 0; i < n; i++) {
        string x_expr, y_expr;

        cout << "x[" << i << "] = ";
        getline(cin, x_expr);
        x[i] = parser.evaluate(x_expr);
    
        cout << "y[" << i << "] = ";
        getline(cin, y_expr);
        y[i] = parser.evaluate(y_expr);    
        x[i] = log(x[i]); // transform x to ln(x)
    }

    for (int i = 0; i < n; i++) {
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumX2 += x[i] * x[i];
    }

    a = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    b = (sumY - a * sumX) / n;

    return {a, b};
}

/**
 * @brief Performs power curve fitting on a set of (x, y) data points.
 * 
 * Fits the data to the model:
 *     y = a * x^b
 * 
 * Linearizes the data using logarithms:
 *     ln(y) = ln(a) + b * ln(x)
 * 
 * Assumes x > 0 and y > 0 for all data points.
 * 
 * @param n The number of data points to be entered.
 * @return pair<double, double> The coefficients (a, b):
 *         - a: scaling constant
 *         - b: power exponent
 */
pair<double, double> CurveFittingPower(int n) {
    double sumX = 0, sumX2 = 0, sumY = 0, sumXY = 0, a, b, ln_a, x[P], y[P];

    for (int i = 0; i < n; i++) {
        string x_expr, y_expr;

        cout << "x[" << i << "] = ";
        getline(cin, x_expr);
        x[i] = parser.evaluate(x_expr);
    
        cout << "y[" << i << "] = ";
        getline(cin, y_expr);
        y[i] = parser.evaluate(y_expr);    
        x[i] = log(x[i]); // ln(x)
        y[i] = log(y[i]); // ln(y)
    }

    for (int i = 0; i < n; i++) {
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumX2 += x[i] * x[i];
    }

    b = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    ln_a = (sumY - b * sumX) / n;
    a = exp(ln_a);

    return {a, b};
}

/**
 * @brief Performs 2nd order polynomial curve fitting on a set of (x, y) data points.
 * 
 * Fits the data to a quadratic model using the least squares method:
 *     y = a * x^2 + b * x + c
 * 
 * Solves the normal equations using matrix methods.
 * 
 * @param n The number of data points to be entered.
 * @return tuple<double, double, double> The coefficients (a, b, c):
 *         - a: coefficient of x^2
 *         - b: coefficient of x
 *         - c: constant term (y-intercept)
 */
tuple<double, double, double> CurveFittingQuadratic(int n) {
    vector<double> x(n), y(n);
    double sumX = 0, sumX2 = 0, sumX3 = 0, sumX4 = 0;
    double sumY = 0, sumXY = 0, sumX2Y = 0;

    for (int i = 0; i < n; ++i) {
        string x_expr, y_expr;

        cout << "x[" << i << "] = ";
        getline(cin, x_expr);
        x[i] = parser.evaluate(x_expr);
    
        cout << "y[" << i << "] = ";
        getline(cin, y_expr);
        y[i] = parser.evaluate(y_expr);

        double xi = x[i], yi = y[i];
        sumX += xi;
        sumX2 += xi * xi;
        sumX3 += xi * xi * xi;
        sumX4 += xi * xi * xi * xi;
        sumY += yi;
        sumXY += xi * yi;
        sumX2Y += xi * xi * yi;
    }

    // Solving normal equations using Cramer's Rule
    double D = n * (sumX2 * sumX4 - sumX3 * sumX3) - sumX * (sumX * sumX4 - sumX2 * sumX3)
             + sumX2 * (sumX * sumX3 - sumX2 * sumX2);
             
    double Da = n * (sumX2 * sumX2Y - sumXY * sumX3) - sumX * (sumX * sumX2Y - sumXY * sumX2)
            + sumY * (sumX * sumX3 - sumX2 * sumX2);
    
    double Db = n * (sumXY * sumX4 - sumX2Y * sumX3) - sumY * (sumX * sumX4 - sumX2 * sumX3)
            + sumX2 * (sumX * sumX2Y - sumXY * sumX2);

    double Dc = sumY * (sumX2 * sumX4 - sumX3 * sumX3) - sumX * (sumXY * sumX4 - sumX2Y * sumX3)
            + sumX2 * (sumXY * sumX3 - sumX2 * sumX2Y);

    double a = Da / D;
    double b = Db / D;
    double c = Dc / D;

    return {a, b, c};
}
