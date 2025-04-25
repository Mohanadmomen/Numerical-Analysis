#pragma once

#include <iostream>
#include <vector>
#include <functional>
#include <string>
#include <utility>

class ChapterFour {
private:
    ExpressionParser parser;

    // Helper function to evaluate the expression at a specific point (x, y)
    double evaluateFunction(const std::string& expression, double x, double y) {
        parser.setVariable("x", x);
        parser.setVariable("y", y);
        return parser.evaluate(expression);
    }

public:
    ChapterFour() {}

    // Structure to store the solution points
    struct Solution {
        std::vector<double> x_values;
        std::vector<double> y_values;
    };

    // Euler's Method for solving first-order differential equations: dy/dx = f(x,y)
    Solution euler(const std::string& odeExpression, double x0, double y0, double h, double xEnd) {
        Solution result;

        // Initialize solution vectors with initial conditions
        result.x_values.push_back(x0);
        result.y_values.push_back(y0);

        double x = x0;
        double y = y0;

        // Continue until we reach or exceed the final x value
        while (x < xEnd - 1e-10) {
            // Set variables in the parser
            parser.setVariable("x", x);
            parser.setVariable("y", y);

            // Compute the slope using f(x,y)
            double slope = parser.evaluate(odeExpression);

            // Euler update: y_next = y + h * f(x,y)
            y = y + h * slope;
            x = x + h;

            // Store solution points
            result.x_values.push_back(x);
            result.y_values.push_back(y);
        }

        return result;
    }

    // Modified Euler Method (Heun's Method) for solving first-order differential equations: dy/dx = f(x,y)
    Solution modifiedEuler(const std::string& odeExpression, double x0, double y0, double h, double xEnd) {
        Solution result;

        // Initialize solution vectors with initial conditions
        result.x_values.push_back(x0);
        result.y_values.push_back(y0);

        double x = x0;
        double y = y0;

        // Continue until we reach or exceed the final x value
        while (x < xEnd - 1e-10) {
            // Step 1: Compute k1 = f(x, y)
            parser.setVariable("x", x);
            parser.setVariable("y", y);
            double k1 = parser.evaluate(odeExpression);

            // Step 2: Predictor: y* = y + h * k1
            double y_pred = y + h * k1;

            // Step 3: Compute k2 = f(x + h, y*)
            parser.setVariable("x", x + h);
            parser.setVariable("y", y_pred);
            double k2 = parser.evaluate(odeExpression);

            // Step 4: Modified Euler update: y_next = y + h * (k1 + k2) / 2
            y = y + h * (k1 + k2) / 2;
            x = x + h;

            // Store solution points
            result.x_values.push_back(x);
            result.y_values.push_back(y);
        }

        return result;
    }

    // Print the solution to the console
    void printSolution(const Solution& solution, const std::string& methodName) {
        std::cout << "\n" << methodName << " Solution:\n";
        std::cout << "-------------------\n";
        std::cout << "   x   |   y   \n";
        std::cout << "-------------------\n";

        for (size_t i = 0; i < solution.x_values.size(); ++i) {
            printf(" %.4f | %.6f\n", solution.x_values[i], solution.y_values[i]);
        }
        std::cout << "-------------------\n";
    }

    // Get the solution as a formatted string
    std::string getSolutionString(const Solution& solution, const std::string& methodName) {
        std::string result = methodName + " Solution:\n\n";
        result += "   x   |   y   \n";
        result += "-------------------\n";

        for (size_t i = 0; i < solution.x_values.size(); ++i) {
            char buffer[50];
            sprintf_s(buffer, " %.4f | %.6f\n", solution.x_values[i], solution.y_values[i]);
            result += buffer;
        }

        return result;
    }

    // Compare the results between Euler and Modified Euler methods
    void compareResults(const Solution& eulerSolution, const Solution& modifiedEulerSolution) {
        std::cout << "\nComparison between Euler and Modified Euler methods:\n";
        std::cout << "-------------------------------------------------------\n";
        std::cout << "   x   |  Euler y  | Modified y | Difference\n";
        std::cout << "-------------------------------------------------------\n";

        size_t n = std::min(eulerSolution.x_values.size(), modifiedEulerSolution.x_values.size());

        for (size_t i = 0; i < n; ++i) {
            double diff = std::abs(eulerSolution.y_values[i] - modifiedEulerSolution.y_values[i]);
            printf(" %.4f | %.6f | %.6f | %.6f\n",
                eulerSolution.x_values[i],
                eulerSolution.y_values[i],
                modifiedEulerSolution.y_values[i],
                diff);
        }
        std::cout << "-------------------------------------------------------\n";
    }
};