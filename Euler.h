#ifndef EULER_H
#define EULER_H

#include <string>
#include <cmath>
#include "functionHandler.h"

// Euler's Method - returns final y value
double euler(ExpressionParser& parser, const std::string& odeExpression, double x0, double y0, double h, double xEnd);

// Modified Euler's Method (Heun's Method) - returns final y value
double modifiedEuler(ExpressionParser& parser, const std::string& odeExpression, double x0, double y0, double h, double xEnd);

#endif // EULER_H
