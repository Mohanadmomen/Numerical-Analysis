#include <iostream>
#include <string>
#include <cmath>
#include <vector>







// Lagrange Interpolation Method
// Interpolates a function based on a set of known points
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



std::vector<double> delta_finder(const std::vector<double>& x, const std::vector<double>& y, int n) { /// nmax = size(y)
    std::vector<double> delta; /////// ysize = 2
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


// takes a vector that contains the first element of each
// delta vector , the values of x as a vector , first y and the required value to be found
double newton_forward(std::vector<double> first_deltas, std::vector<double> x, double y0, double xi) {
    double result = 0.0;
    first_deltas.insert(first_deltas.begin(), y0); // add y0 in the begaining
    // now multiply by (x - x0) (x-x1) ... (x-xi)
    for (int i = 0; i < first_deltas.size(); i++) {
        std::vector<double> sliced(x.begin(), x.begin() + i);
        result += first_deltas[i] * helper_of_newton(sliced, xi); // {i=0 , 1 , ...} ==> {y0 * 1 ,delta0 * (x - x0) ,delta1 * (x - x0) * (x - x1 ) , ...}
    }
    return result;
}

// takes a vector that contains the last element of each
// delta vector , the values of x as a vector , last y and the required value to be found
double newton_backward(std::vector<double> last_deltas, std::vector<double> x, double yn, double xi) {
    double result = yn; //will sum the rest in a loop

    // now multiply by (x - xn) (x - xn-1) ... (x-xi)
    for (int i = 0; i < last_deltas.size(); i++) {
        std::vector<double> sliced(x.end() - i - 1, x.end()); // i = {0 , 1 , 2 , ...} ====> {[xn] , [xn-1 , xn] , }
        result += last_deltas.at(i) * helper_of_newton(sliced, xi); // {i=0 , 1 , ...} ==> {delta0 * (x - xn) ,delta1 * (x - xn) * (x - xn-1 ) , ...}
    }
    return result;
}

double NewtonInterpolation(const std::vector<double>& x, const std::vector<double>& y, double xi) {
    if (x.size() != y.size() || x.empty()) {
        throw std::runtime_error("Input vectors must be non-empty and of the same size");
    }
    double result = 0;
    std::vector< std::vector<double> > delta_vector;
    std::vector< std::vector<double> > delta_n;

    // the number of deltas is less than the number of size of Y by 1
    // each index contains a vector
    delta_vector.push_back(delta_finder(x, y, 1)); /////##########
    for (int i = 2; i < y.size(); i++) {
        delta_vector.push_back(delta_finder(x, delta_vector.back(), i));
    }

    // Now choosing the proper method
    // Checking if the required value is closer to first value
    if (fabs(xi - x.at(0)) < fabs(xi - x.at(x.size() - 1))) {
        // applying newton forward
        std::vector<double> first_deltas; // making a vector of first deltas to pass it as an argument
        for (int i = 0; i < delta_vector.at(0).size(); i++) {
            first_deltas.push_back(delta_vector.at(i).at(0)); // Pushing delta[i][0]
        }
        result = newton_forward(first_deltas, x, y.at(0), xi);
    }
    else // Appling newton backward
    {
        std::vector<double> last_deltas; // making a vector of first last to pass it as an argument
        for (int i = 0; i < delta_vector.at(0).size(); i++) {
            last_deltas.push_back(delta_vector.at(i).back()); // Pushing delta[i][last element]
        }
        result = newton_backward(last_deltas, x, y.at(y.size() - 1), xi);
    }

    return result;
}