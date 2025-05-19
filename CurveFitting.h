#ifndef CURVEFITTING_H
#define CURVEFITTING_H

#include <string>
#include <vector>
#include <tuple>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <limits>
#include <algorithm>
#include "FunctionHandler.h"

using namespace std;

class CurveFitter {
private:
    ExpressionParser parser;

    string formatDouble(double value, int precision = 4);
    string prettyCheckSign(double number);

public:
    enum Model {
        LINEAR,
        EXPONENTIAL,
        LOGARITHMIC,
        POWER,
        QUADRATIC,
        CUBIC
    };

    CurveFitter();

    string fit(const vector<double>& x, const vector<double>& y, Model model);

    pair<double, double> fitLinear(const vector<double>& x, const vector<double>& y);
    pair<double, double> fitExponential(const vector<double>& x, const vector<double>& y);
    pair<double, double> fitLogarithmic(const vector<double>& x, const vector<double>& y);
    pair<double, double> fitPower(const vector<double>& x, const vector<double>& y);
    tuple<double, double, double> fitQuadratic(const vector<double>& x, const vector<double>& y);
    tuple<double, double, double, double> fitCubic(const vector<double>& x, const vector<double>& y);

    pair<double, double> getLinearCoefficients(const vector<double>& x, const vector<double>& y);
    pair<double, double> getExponentialCoefficients(const vector<double>& x, const vector<double>& y);
    pair<double, double> getLogarithmicCoefficients(const vector<double>& x, const vector<double>& y);
    pair<double, double> getPowerCoefficients(const vector<double>& x, const vector<double>& y);
    tuple<double, double, double> getQuadraticCoefficients(const vector<double>& x, const vector<double>& y);
    tuple<double, double, double, double> getCubicCoefficients(const vector<double>& x, const vector<double>& y);
};

#endif // CURVEFITTING_H
