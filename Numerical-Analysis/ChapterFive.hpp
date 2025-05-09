#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <tuple>
#include <limits>
#include <algorithm>
#include <iomanip>

#define P 50

using namespace std;

ExpressionParser parser;

pair<double, double> CurveFittingLinear(int n);
pair<double, double> CurveFittingExponential(int n);
pair<double, double> CurveFittingLogarithmic(int n);
pair<double, double> CurveFittingPower(int n);
tuple<double, double, double> CurveFittingQuadratic(int n);
tuple<double, double, double, double> CurveFittingCubic(int n);
string CheckSign(double number);


string CheckSign(double number) {
    if (number < 0)
        return " - " + to_string(abs(number));
    return " + " + to_string(number);
}

// Improved number formatting 
string formatDouble(double value, int precision = 4) {
    ostringstream oss;
    oss << fixed << setprecision(precision) << value;
    return oss.str();
}

string prettyCheckSign(double number) {
    if (number < 0)
        return " - " + formatDouble(abs(number));
    return " + " + formatDouble(number);
}

void CurveFitting() {
    int ModelNum, n;
    double a, b, c, d;
    
    do {
        cout << "\n================================\n";
        cout << "[1] Linear       (a * x + b)\n";
        cout << "[2] Exponential  (a * exp(b * x))\n";
        cout << "[3] Logarithmic  (a * ln(x) + b)\n";
        cout << "[4] Power        (a * x^b)\n";
        cout << "[5] Quadratic    (a * x^2 + b * x + c)\n";
        cout << "[6] Cubic        (a * x^3 + b * x^2 + c * x + d)\n";
        cout << "[7] Exit\n";
        cout << "================================\n\n";
        cout << "> Choose the model number: ";
        cin >> ModelNum;
    } while (ModelNum < 1 || ModelNum > 7);
    
    if (ModelNum != 7) {
        cout << "> How many data points? ";
        cin >> n;
        
        // Input validation
        while (n <= 0 || n > P) {
            cout << "Invalid number of data points. Please enter a value between 1 and " << P << ": ";
            cin >> n;
        }
    }

    cin.ignore(); // ignore newline after reading n

    switch(ModelNum) {
        case 1: {
            tie(a, b) = CurveFittingLinear(n);
            if (a == 0)
                cout << "Fitted Linear: y = " << formatDouble(b) << "\n";
            else
                cout << "Fitted Linear: y = " << formatDouble(a) << "x" << prettyCheckSign(b) << "\n";
            break;
        }
        case 2: {
            tie(a, b) = CurveFittingExponential(n);
            cout << "Fitted Exponential: y = " << formatDouble(a) << " * exp(" << formatDouble(b) << "x)\n";
            break;
        }
        case 3: {
            tie(a, b) = CurveFittingLogarithmic(n);
            cout << "Fitted Logarithmic: y = " << formatDouble(a) << " * ln(x)" << prettyCheckSign(b) << "\n";
            break;
        }
        case 4: {
            tie(a, b) = CurveFittingPower(n);
            cout << "Fitted Power: y = " << formatDouble(a) << " * x^" << formatDouble(b) << "\n";
            break;
        }
        case 5: {
            tie(a, b, c) = CurveFittingQuadratic(n);
            cout << "Fitted Quadratic: y = " << formatDouble(a) << " * x^2" << prettyCheckSign(b) << " * x" << prettyCheckSign(c) << "\n";
            break;
        }
        case 6: {
            tie(a, b, c, d) = CurveFittingCubic(n);
            cout << "Fitted Cubic: y = " << formatDouble(a) << " * x^3" << prettyCheckSign(b) << " * x^2" 
                 << prettyCheckSign(c) << " * x" << prettyCheckSign(d) << "\n";
            break;
        }
        case 7:
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

    // Check for potential division by zero
    double denominator = (n * sumX2 - sumX * sumX);
    if (abs(denominator) < 1e-10) {
        cout << "Warning: Division by near-zero value. Results may be inaccurate." << endl;
        return {0, sumY / n}; // Return average y as constant model
    }

    // Calculating a and b
    a = (n * sumXY - sumX * sumY) / denominator;
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
    double sumX = 0, sumX2 = 0, sumY = 0, sumXY = 0, a, b, ln_a;
    double x[P], y[P], ln_y[P];
    int valid_points = 0;  // count valid data points

    for (int i = 0; i < n; i++) {
        string x_expr, y_expr;
        double x_val, y_val;

        cout << "x[" << i << "] = ";
        getline(cin, x_expr);
        x_val = parser.evaluate(x_expr);

        while (true) {
            cout << "y[" << i << "] = ";
            getline(cin, y_expr);
            y_val = parser.evaluate(y_expr);

            if (y_val <= 0) {
                cout << "Warning: Non-positive y value detected at point " << i << " (y = " << y_val << ").\n";
                cout << "Choose an option:\n";
                cout << "1. Exclude this point\n";
                cout << "2. Enter the point again\n";
                cout << "Your choice: ";

                int choice;
                cin >> choice;
                cin.ignore(); // clear newline from input buffer

                if (choice == 1) {
                    cout << "Point " << i << " excluded from fitting.\n";
                    y_val = NAN; // mark as excluded
                    break; // exit input loop
                } else if (choice == 2) {
                    cout << "Please re-enter y[" << i << "] again.\n";
                    continue; // loop again for re-entry
                } else {
                    cout << "Invalid choice. Please choose again.\n";
                    continue;
                }
            } else {
                // valid y
                break;
            }
        }

        // Only store valid points
        if (!isnan(y_val)) {
            x[valid_points] = x_val;
            y[valid_points] = y_val;
            ln_y[valid_points] = log(y_val);
            valid_points++;
        }
    }

    // Now compute sums using only valid points
    for (int i = 0; i < valid_points; i++) {
        sumX += x[i];
        sumY += ln_y[i];
        sumXY += x[i] * ln_y[i];
        sumX2 += x[i] * x[i];
    }

    if (valid_points == 0) {
        cout << "Error: No valid data points for fitting.\n";
        return {0, 0};
    }

    double denominator = (valid_points * sumX2 - sumX * sumX);
    if (abs(denominator) < 1e-10) {
        cout << "Warning: Division by near-zero value. Results may be inaccurate.\n";
        return {exp(sumY / valid_points), 0};
    }

    b = (valid_points * sumXY - sumX * sumY) / denominator;
    ln_a = (sumY - b * sumX) / valid_points;
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
    double sumX = 0, sumX2 = 0, sumY = 0, sumXY = 0, a, b;
    double x[P], y[P], ln_x[P];
    int valid_points = 0;  // count of valid data points

    for (int i = 0; i < n; i++) {
        string x_expr, y_expr;
        double x_val, y_val;

        while (true) {
            cout << "x[" << i << "] = ";
            getline(cin, x_expr);
            x_val = parser.evaluate(x_expr);

            if (x_val <= 0) {
                cout << "Warning: Non-positive x value detected at point " << i << " (x = " << x_val << ").\n";
                cout << "Choose an option:\n";
                cout << "1. Exclude this point\n";
                cout << "2. Enter the point again\n";
                cout << "Your choice: ";

                int choice;
                cin >> choice;
                cin.ignore();  // clear input buffer

                if (choice == 1) {
                    cout << "Point " << i << " excluded from fitting.\n";
                    x_val = NAN; // mark as excluded
                    break;
                } else if (choice == 2) {
                    cout << "Please re-enter x[" << i << "] again.\n";
                    continue; // ask again
                } else {
                    cout << "Invalid choice. Please choose again.\n";
                    continue;
                }
            } else {
                // valid x
                break;
            }
        }

        if (isnan(x_val)) {
            continue; // skip this point
        }

        ln_x[valid_points] = log(x_val);
        x[valid_points] = x_val;

        cout << "y[" << i << "] = ";
        getline(cin, y_expr);
        y_val = parser.evaluate(y_expr);

        y[valid_points] = y_val;
        valid_points++;
    }

    if (valid_points == 0) {
        cout << "Error: No valid data points for fitting.\n";
        return {0, 0};
    }

    for (int i = 0; i < valid_points; i++) {
        sumX += ln_x[i];
        sumY += y[i];
        sumXY += ln_x[i] * y[i];
        sumX2 += ln_x[i] * ln_x[i];
    }

    double denominator = (valid_points * sumX2 - sumX * sumX);
    if (abs(denominator) < 1e-10) {
        cout << "Warning: Division by near-zero value. Results may be inaccurate.\n";
        return {0, sumY / valid_points}; // Return constant model
    }

    a = (valid_points * sumXY - sumX * sumY) / denominator;
    b = (sumY - a * sumX) / valid_points;

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
    double sumX = 0, sumX2 = 0, sumY = 0, sumXY = 0, a, b, ln_a;
    double x[P], y[P], ln_x[P], ln_y[P];
    int valid_points = 0;  // count valid data points

    for (int i = 0; i < n; i++) {
        string x_expr, y_expr;
        double x_val, y_val;

        // Input x[i] with validation
        while (true) {
            cout << "x[" << i << "] = ";
            getline(cin, x_expr);
            x_val = parser.evaluate(x_expr);

            if (x_val <= 0) {
                cout << "Warning: Non-positive x value detected at point " << i << " (x = " << x_val << ").\n";
                cout << "Choose an option:\n";
                cout << "1. Exclude this point\n";
                cout << "2. Enter the point again\n";
                cout << "Your choice: ";

                int choice;
                cin >> choice;
                cin.ignore();

                if (choice == 1) {
                    cout << "Point " << i << " excluded from fitting.\n";
                    x_val = NAN;  // mark excluded
                    break;
                } else if (choice == 2) {
                    cout << "Please re-enter x[" << i << "] again.\n";
                    continue;
                } else {
                    cout << "Invalid choice. Please choose again.\n";
                    continue;
                }
            } else {
                break;
            }
        }

        if (isnan(x_val)) continue;  // skip this point

        // Input y[i] with validation
        while (true) {
            cout << "y[" << i << "] = ";
            getline(cin, y_expr);
            y_val = parser.evaluate(y_expr);

            if (y_val <= 0) {
                cout << "Warning: Non-positive y value detected at point " << i << " (y = " << y_val << ").\n";
                cout << "Choose an option:\n";
                cout << "1. Exclude this point\n";
                cout << "2. Enter the point again\n";
                cout << "Your choice: ";

                int choice;
                cin >> choice;
                cin.ignore();

                if (choice == 1) {
                    cout << "Point " << i << " excluded from fitting.\n";
                    y_val = NAN;  // mark excluded
                    break;
                } else if (choice == 2) {
                    cout << "Please re-enter y[" << i << "] again.\n";
                    continue;
                } else {
                    cout << "Invalid choice. Please choose again.\n";
                    continue;
                }
            } else {
                break;
            }
        }

        if (isnan(y_val)) continue;  // skip this point

        // Store valid values
        x[valid_points] = x_val;
        y[valid_points] = y_val;
        ln_x[valid_points] = log(x_val);
        ln_y[valid_points] = log(y_val);
        valid_points++;
    }

    if (valid_points == 0) {
        cout << "Error: No valid data points for fitting.\n";
        return {0, 0};
    }

    for (int i = 0; i < valid_points; i++) {
        sumX += ln_x[i];
        sumY += ln_y[i];
        sumXY += ln_x[i] * ln_y[i];
        sumX2 += ln_x[i] * ln_x[i];
    }

    double denominator = (valid_points * sumX2 - sumX * sumX);
    if (abs(denominator) < 1e-10) {
        cout << "Warning: Division by near-zero value. Results may be inaccurate.\n";
        return {exp(sumY / valid_points), 0};
    }

    b = (valid_points * sumXY - sumX * sumY) / denominator;
    ln_a = (sumY - b * sumX) / valid_points;
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

    // Checking for numerical stability issues
    double det = n * (sumX2 * sumX4 - sumX3 * sumX3) - 
                sumX * (sumX * sumX4 - sumX2 * sumX3) + 
                sumX2 * (sumX * sumX3 - sumX2 * sumX2);
                
    if (abs(det) < 1e-10) {
        cout << "Warning: Matrix is nearly singular. Results may be inaccurate." << endl;
        // Fall back to linear or constant fit
        auto [a_linear, b_linear] = CurveFittingLinear(n);
        return {0, a_linear, b_linear};
    }

    // Solving normal equations using Cramer's Rule
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

    return {a, b, c};
}

tuple<double, double, double, double> CurveFittingCubic(int n) {
    if (n < 4) {
        cout << "Warning: Need at least 4 points for cubic fitting. Falling back to quadratic." << endl;
        auto [a, b, c] = CurveFittingQuadratic(n);
        return {0, a, b, c};
    }

    vector<double> x(n), y(n);
    double sumX = 0, sumX2 = 0, sumX3 = 0, sumX4 = 0, sumX5 = 0, sumX6 = 0;
    double sumY = 0, sumXY = 0, sumX2Y = 0, sumX3Y = 0;

    for (int i = 0; i < n; ++i) {
        string x_expr, y_expr;

        cout << "x[" << i << "] = ";
        getline(cin, x_expr);
        x[i] = parser.evaluate(x_expr);

        cout << "y[" << i << "] = ";
        getline(cin, y_expr);
        y[i] = parser.evaluate(y_expr);

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

    // Coefficient matrix for the normal equations
    double A[4][4] = {
        {n,    sumX,   sumX2,  sumX3},
        {sumX, sumX2,  sumX3,  sumX4},
        {sumX2,sumX3,  sumX4,  sumX5},
        {sumX3,sumX4,  sumX5,  sumX6}
    };

    // Right-hand side
    double B[4] = {sumY, sumXY, sumX2Y, sumX3Y};

    // Implement Gaussian elimination with pivoting for better numerical stability
    double M[4][5]; // Augmented matrix
    
    // Initialize augmented matrix
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            M[i][j] = A[i][j];
        }
        M[i][4] = B[i];
    }
    
    // Perform row operations
    for (int i = 0; i < 4; i++) {
        // Find pivot
        int maxRow = i;
        double maxVal = abs(M[i][i]);
        for (int j = i + 1; j < 4; j++) {
            if (abs(M[j][i]) > maxVal) {
                maxVal = abs(M[j][i]);
                maxRow = j;
            }
        }
        
        // Check for singularity
        if (maxVal < 1e-10) {
            cout << "Warning: Matrix is singular. Falling back to quadratic fit." << endl;
            auto [a_quad, b_quad, c_quad] = CurveFittingQuadratic(n);
            return {0, a_quad, b_quad, c_quad};
        }
        
        // Swap rows if needed
        if (maxRow != i) {
            for (int j = i; j < 5; j++) {
                double temp = M[i][j];
                M[i][j] = M[maxRow][j];
                M[maxRow][j] = temp;
            }
        }
        
        // Normalize pivot row
        double pivot = M[i][i];
        for (int j = i; j < 5; j++) {
            M[i][j] /= pivot;
        }
        
        // Eliminate column
        for (int j = 0; j < 4; j++) {
            if (j != i) {
                double factor = M[j][i];
                for (int k = i; k < 5; k++) {
                    M[j][k] -= factor * M[i][k];
                }
            }
        }
    }
    
    // Extract solutions
    double d = M[0][4]; // constant term
    double c = M[1][4]; // coefficient of x
    double b = M[2][4]; // coefficient of x^2
    double a = M[3][4]; // coefficient of x^3

    return {a, b, c, d};  // Corresponds to ax^3 + bx^2 + cx + d
}
