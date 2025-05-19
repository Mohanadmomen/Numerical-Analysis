#include "CurveFitting.h"

using namespace std;

CurveFitter::CurveFitter() {}

string CurveFitter::formatDouble(double value, int precision) {
    ostringstream oss;
    oss << fixed << setprecision(precision) << value;
    return oss.str();
}

string CurveFitter::prettyCheckSign(double number) {
    if (number < 0)
        return " - " + formatDouble(abs(number));
    return " + " + formatDouble(number);
}

string CurveFitter::fit(const vector<double>& x, const vector<double>& y, Model model) {
    int n = x.size();

    switch (model) {
    case LINEAR: {
        auto [a, b] = fitLinear(x, y);
        if (a == 0)
            return "Fitted Linear: y = " + formatDouble(b);
        else
            return "Fitted Linear: y = " + formatDouble(a) + "x" + prettyCheckSign(b);
    }
    case EXPONENTIAL: {
        auto [a, b] = fitExponential(x, y);
        return "Fitted Exponential: y = " + formatDouble(a) + " * exp(" + formatDouble(b) + "x)";
    }
    case LOGARITHMIC: {
        auto [a, b] = fitLogarithmic(x, y);
        return "Fitted Logarithmic: y = " + formatDouble(a) + " * ln(x)" + prettyCheckSign(b);
    }
    case POWER: {
        auto [a, b] = fitPower(x, y);
        return "Fitted Power: y = " + formatDouble(a) + " * x^" + formatDouble(b);
    }
    case QUADRATIC: {
        auto [a, b, c] = fitQuadratic(x, y);
        return "Fitted Quadratic: y = " + formatDouble(a) + " * x^2" + prettyCheckSign(b) + " * x" + prettyCheckSign(c);
    }
    case CUBIC: {
        auto [a, b, c, d] = fitCubic(x, y);
        return "Fitted Cubic: y = " + formatDouble(a) + " * x^3" + prettyCheckSign(b) + " * x^2" + prettyCheckSign(c) + " * x" + prettyCheckSign(d);
    }
    default:
        return "Invalid model";
    }
}

pair<double, double> CurveFitter::fitLinear(const vector<double>& x, const vector<double>& y) {
    int n = x.size();
    double sumX = 0, sumX2 = 0, sumY = 0, sumXY = 0;

    for (int i = 0; i < n; i++) {
        sumX += x[i];
        sumY += y[i];
        sumX2 += x[i] * x[i];
        sumXY += x[i] * y[i];
    }

    double denominator = (n * sumX2 - sumX * sumX);
    if (abs(denominator) < 1e-10) {
        return { 0, sumY / n };
    }

    double a = (n * sumXY - sumX * sumY) / denominator;
    double b = (sumY - a * sumX) / n;

    return { a, b };
}

pair<double, double> CurveFitter::fitExponential(const vector<double>& x, const vector<double>& y) {
    int n = x.size();
    vector<double> valid_x, valid_y, ln_y;

    for (int i = 0; i < n; i++) {
        if (y[i] > 0) {
            valid_x.push_back(x[i]);
            valid_y.push_back(y[i]);
            ln_y.push_back(log(y[i]));
        }
    }

    int valid_points = valid_x.size();
    if (valid_points == 0) return { 0, 0 };

    double sumX = 0, sumX2 = 0, sumY = 0, sumXY = 0;

    for (int i = 0; i < valid_points; i++) {
        sumX += valid_x[i];
        sumY += ln_y[i];
        sumXY += valid_x[i] * ln_y[i];
        sumX2 += valid_x[i] * valid_x[i];
    }

    double denominator = (valid_points * sumX2 - sumX * sumX);
    if (abs(denominator) < 1e-10) {
        return { exp(sumY / valid_points), 0 };
    }

    double b = (valid_points * sumXY - sumX * sumY) / denominator;
    double ln_a = (sumY - b * sumX) / valid_points;
    double a = exp(ln_a);

    return { a, b };
}

pair<double, double> CurveFitter::fitLogarithmic(const vector<double>& x, const vector<double>& y) {
    int n = x.size();
    vector<double> valid_x, valid_y, ln_x;

    for (int i = 0; i < n; i++) {
        if (x[i] > 0) {
            valid_x.push_back(x[i]);
            valid_y.push_back(y[i]);
            ln_x.push_back(log(x[i]));
        }
    }

    int valid_points = valid_x.size();
    if (valid_points == 0) return { 0, 0 };

    double sumX = 0, sumX2 = 0, sumY = 0, sumXY = 0;

    for (int i = 0; i < valid_points; i++) {
        sumX += ln_x[i];
        sumY += valid_y[i];
        sumXY += ln_x[i] * valid_y[i];
        sumX2 += ln_x[i] * ln_x[i];
    }

    double denominator = (valid_points * sumX2 - sumX * sumX);
    if (abs(denominator) < 1e-10) {
        return { 0, sumY / valid_points };
    }

    double a = (valid_points * sumXY - sumX * sumY) / denominator;
    double b = (sumY - a * sumX) / valid_points;

    return { a, b };
}

pair<double, double> CurveFitter::fitPower(const vector<double>& x, const vector<double>& y) {
    int n = x.size();
    vector<double> valid_x, valid_y, ln_x, ln_y;

    for (int i = 0; i < n; i++) {
        if (x[i] > 0 && y[i] > 0) {
            valid_x.push_back(x[i]);
            valid_y.push_back(y[i]);
            ln_x.push_back(log(x[i]));
            ln_y.push_back(log(y[i]));
        }
    }

    int valid_points = valid_x.size();
    if (valid_points == 0) return { 0, 0 };

    double sumX = 0, sumX2 = 0, sumY = 0, sumXY = 0;

    for (int i = 0; i < valid_points; i++) {
        sumX += ln_x[i];
        sumY += ln_y[i];
        sumXY += ln_x[i] * ln_y[i];
        sumX2 += ln_x[i] * ln_x[i];
    }

    double denominator = (valid_points * sumX2 - sumX * sumX);
    if (abs(denominator) < 1e-10) {
        return { exp(sumY / valid_points), 0 };
    }

    double b = (valid_points * sumXY - sumX * sumY) / denominator;
    double ln_a = (sumY - b * sumX) / valid_points;
    double a = exp(ln_a);

    return { a, b };
}

tuple<double, double, double> CurveFitter::fitQuadratic(const vector<double>& x, const vector<double>& y) {
    int n = x.size();
    double sumX = 0, sumX2 = 0, sumX3 = 0, sumX4 = 0;
    double sumY = 0, sumXY = 0, sumX2Y = 0;

    for (int i = 0; i < n; ++i) {
        double xi = x[i], yi = y[i];
        sumX += xi;
        sumX2 += xi * xi;
        sumX3 += xi * xi * xi;
        sumX4 += xi * xi * xi * xi;
        sumY += yi;
        sumXY += xi * yi;
        sumX2Y += xi * xi * yi;
    }

    double det = n * (sumX2 * sumX4 - sumX3 * sumX3) -
                 sumX * (sumX * sumX4 - sumX2 * sumX3) +
                 sumX2 * (sumX * sumX3 - sumX2 * sumX2);

    if (abs(det) < 1e-10) {
        auto [a_linear, b_linear] = fitLinear(x, y);
        return { 0, a_linear, b_linear };
    }

    double Da = n * (sumX2 * sumX2Y - sumXY * sumX3) -
                sumX * (sumX * sumX2Y - sumXY * sumX2) +
                sumY * (sumX * sumX3 - sumX2 * sumX2);

    double Db = n * (sumXY * sumX4 - sumX2Y * sumX3) -
                sumY * (sumX * sumX4 - sumX2 * sumX3) +
                sumX2 * (sumX * sumX2Y - sumXY * sumX2);

    double Dc = sumY * (sumX2 * sumX4 - sumX3 * sumX3) -
                sumX * (sumXY * sumX4 - sumX2Y * sumX3) +
                sumX2 * (sumXY * sumX3 - sumX2 * sumX2Y);

    double a = Da / det;
    double b = Db / det;
    double c = Dc / det;

    return { a, b, c };
}

tuple<double, double, double, double> CurveFitter::fitCubic(const vector<double>& x, const vector<double>& y) {
    int n = x.size();
    if (n < 4) {
        auto [a, b, c] = fitQuadratic(x, y);
        return { 0, a, b, c };
    }

    double sumX = 0, sumX2 = 0, sumX3 = 0, sumX4 = 0, sumX5 = 0, sumX6 = 0;
    double sumY = 0, sumXY = 0, sumX2Y = 0, sumX3Y = 0;

    for (int i = 0; i < n; ++i) {
        double xi = x[i], yi = y[i];
        double xi2 = xi * xi;
        double xi3 = xi2 * xi;

        sumX += xi;
        sumX2 += xi2;
        sumX3 += xi3;
        sumX4 += xi3 * xi;
        sumX5 += xi3 * xi2;
        sumX6 += xi3 * xi3;

        sumY += yi;
        sumXY += xi * yi;
        sumX2Y += xi2 * yi;
        sumX3Y += xi3 * yi;
    }

    double A[4][4] = {
        {n,    sumX,   sumX2,  sumX3},
        {sumX, sumX2,  sumX3,  sumX4},
        {sumX2,sumX3,  sumX4,  sumX5},
        {sumX3,sumX4,  sumX5,  sumX6}
    };

    double B[4] = { sumY, sumXY, sumX2Y, sumX3Y };
    double M[4][5];

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) M[i][j] = A[i][j];
        M[i][4] = B[i];
    }

    for (int i = 0; i < 4; i++) {
        int maxRow = i;
        double maxVal = abs(M[i][i]);
        for (int j = i + 1; j < 4; j++) {
            if (abs(M[j][i]) > maxVal) {
                maxVal = abs(M[j][i]);
                maxRow = j;
            }
        }

        if (maxVal < 1e-10) {
            auto [a_quad, b_quad, c_quad] = fitQuadratic(x, y);
            return { 0, a_quad, b_quad, c_quad };
        }

        if (maxRow != i) {
            for (int j = i; j < 5; j++) {
                double temp = M[i][j];
                M[i][j] = M[maxRow][j];
                M[maxRow][j] = temp;
            }
        }

        double pivot = M[i][i];
        for (int j = i; j < 5; j++) M[i][j] /= pivot;

        for (int j = 0; j < 4; j++) {
            if (j != i) {
                double factor = M[j][i];
                for (int k = i; k < 5; k++) M[j][k] -= factor * M[i][k];
            }
        }
    }

    double d = M[0][4];
    double c = M[1][4];
    double b = M[2][4];
    double a = M[3][4];

    return { a, b, c, d };
}

pair<double, double> CurveFitter::getLinearCoefficients(const vector<double>& x, const vector<double>& y) {
    return fitLinear(x, y);
}
pair<double, double> CurveFitter::getExponentialCoefficients(const vector<double>& x, const vector<double>& y) {
    return fitExponential(x, y);
}
pair<double, double> CurveFitter::getLogarithmicCoefficients(const vector<double>& x, const vector<double>& y) {
    return fitLogarithmic(x, y);
}
pair<double, double> CurveFitter::getPowerCoefficients(const vector<double>& x, const vector<double>& y) {
    return fitPower(x, y);
}
tuple<double, double, double> CurveFitter::getQuadraticCoefficients(const vector<double>& x, const vector<double>& y) {
    return fitQuadratic(x, y);
}
tuple<double, double, double, double> CurveFitter::getCubicCoefficients(const vector<double>& x, const vector<double>& y) {
    return fitCubic(x, y);
}
