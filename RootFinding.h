#ifndef ROOTFINDING_H
#define ROOTFINDING_H

#include <string>
#include <utility>
#include <cmath>
#include <iostream>
#include "FunctionHandler.h"

// Root interval finder
std::pair<double, double> findRootInterval(ExpressionParser& parser, const std::string& expression);

// Root-finding methods
double bisectionMethod(ExpressionParser& parser, const std::string& expression,
                       double tol = 1e-6, int maxIter = 100);

double secantMethod(ExpressionParser& parser, const std::string& expression,
                    double tol = 1e-6, int maxIter = 100);

double newton_raphson(ExpressionParser& parser, const std::string& expression,
                      double tol = 1e-6, int maxIter = 100);

// Utility
double differntiate(ExpressionParser& parser, const std::string& expression, double x);

#endif // ROOTFINDING_H
