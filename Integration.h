#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <vector>
#include <string>
#include "FunctionHandler.h"

// Small epsilon value for comparing doubles and handling near-zero results
const double EPSILON = 1e-12;

// Helper function to round and handle near-zero results
double formatResult(double result);

// 1D Integration Methods
double trapezoidalIntegration(ExpressionParser& parser, const std::string& expression,
                              double a, double b, int n);
double simpsonOneThird(ExpressionParser& parser, const std::string& expression,
                       double a, double b, int n);
double simpsonThreeEighth(ExpressionParser& parser, const std::string& expression,
                          double a, double b, int n);

// 2D Integration Methods
double trapezoidalIntegration2D(ExpressionParser& parser, const std::string& expression,
                                double a, double b, int n,
                                const std::vector<double>& x, const std::vector<double>& y);
double simpsonOneThird2D(ExpressionParser& parser, const std::string& expression,
                         double a, double b, int n,
                         const std::vector<double>& x = {}, const std::vector<double>& y = {});
double simpsonThreeEighth2D(ExpressionParser& parser, const std::string& expression,
                            double a, double b, int n,
                            const std::vector<double>& x = {}, const std::vector<double>& y = {});

// Table-based Integration Methods
double trapezoidalIntegrationTable(ExpressionParser& parser, const std::string& expression,
                                   double a, double b, int n,
                                   const std::vector<double>& x = {});
double simpsonOneThirdTable(ExpressionParser& parser, const std::string& expression,
                            double a, double b, int n,
                            const std::vector<double>& x = {});
double simpsonThreeEighthTable(ExpressionParser& parser, const std::string& expression,
                               double a, double b, int n,
                               const std::vector<double>& x = {});

#endif // INTEGRATION_H
