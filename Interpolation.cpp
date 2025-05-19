#include "interpolation.h"
#include <stdexcept>
#include <cmath>

double lagrangeInterpolation(const std::vector<double>& x, const std::vector<double>& y, double xi) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("Input vectors must be non-empty and of the same size");
    }

    double result = 0.0;
    int n = x.size();

    for (int i = 0; i < n; i++) {
        double term = y[i];
        for (int j = 0; j < n; j++) {
            if (j != i) {
                if (x[i] == x[j]) {
                    throw std::runtime_error("Input x values must be distinct");
                }
                term *= (xi - x[j]) / (x[i] - x[j]);
            }
        }
        result += term;
    }

    return result;
}

std::vector<double> delta_finder(const std::vector<double>& x, const std::vector<double>& y, int n) {
    std::vector<double> delta;
    for (int i = 1; i < y.size(); i++) {
        if (fabs(((x.at(i + n - 1) - x.at(i - 1)))) < 1e-9) {
            throw std::runtime_error("Division by zero detected");
        }

        delta.push_back((y.at(i) - y.at(i - 1)) / (x.at(i + n - 1) - x.at(i - 1)));
    }
    return delta;
}

double helper_of_newton(std::vector<double> x, double xi) {
    double result = 1;
    for (int i = 0; i < x.size(); i++) {
        result *= (xi - x.at(i));
    }
    return result;
}

double newton_forward(std::vector<double> first_deltas, std::vector<double> x, double y0, double xi) {
    double result = 0.0;
    first_deltas.insert(first_deltas.begin(), y0);
    for (int i = 0; i < first_deltas.size(); i++) {
        std::vector<double> sliced(x.begin(), x.begin() + i);
        result += first_deltas[i] * helper_of_newton(sliced, xi);
    }
    return result;
}

double newton_backward(std::vector<double> last_deltas, std::vector<double> x, double yn, double xi) {
    double result = yn;
    for (int i = 0; i < last_deltas.size(); i++) {
        std::vector<double> sliced(x.end() - i - 1, x.end());
        result += last_deltas.at(i) * helper_of_newton(sliced, xi);
    }
    return result;
}

double NewtonInterpolation(const std::vector<double>& x, const std::vector<double>& y, double xi) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("Input vectors must be non-empty and of the same size");
    }
    double result = 0;
    std::vector<std::vector<double>> delta_vector;

    delta_vector.push_back(delta_finder(x, y, 1));
    for (int i = 2; i < y.size(); i++) {
        delta_vector.push_back(delta_finder(x, delta_vector.back(), i));
    }

    if (fabs(xi - x.at(0)) < fabs(xi - x.at(x.size() - 1))) {
        std::vector<double> first_deltas;
        for (int i = 0; i < delta_vector.at(0).size(); i++) {
            first_deltas.push_back(delta_vector.at(i).at(0));
        }
        result = newton_forward(first_deltas, x, y.at(0), xi);
    }
    else {
        std::vector<double> last_deltas;
        for (int i = 0; i < delta_vector.at(0).size(); i++) {
            last_deltas.push_back(delta_vector.at(i).back());
        }
        result = newton_backward(last_deltas, x, y.at(y.size() - 1), xi);
    }

    return result;
}
