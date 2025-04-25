#include <iostream>
#include <iomanip>
#include "FunctionHandler.hpp"
#include "ChapterThree.hpp"
#include "ChapterTwo.hpp"
#include "ChapterOne.hpp"
#include "ChapterFive.hpp"

int main() {
    ExpressionParser parser;

    // Example expression: 1/(x+1)
    // Result should be: 1.405357143 by using Trapezoidal method
    
    // std::cout << "Write the equation: ";
    // std::string expression;
    // std::getline(std::cin, expression);

    // Generate C++ function
    // std::string cppFunction = parser.generateCppFunction(expression);
    // std::cout << "Generated C++ function:\n" << cppFunction << std::endl;

    // Evaluate for x = 1    
    // parser.setVariable("x", 1);
    // std::cout << parser.evaluate(expression);

    // Evaluate for different x values
    // std::cout << "Evaluating the expression for different x values:" << std::endl;
    // for (double x = 0; x <= 5; x += 1) {
    //    parser.setVariable("x", x);
    //    try {
    //        double result = parser.evaluate(expression);
    //        std::cout << "f(" << x << ") = " << result << std::endl;
    //    }
    //    catch (const std::exception& e) {
    //        std::cout << "Error: " << e.what() << std::endl;
    //    }
    // }

    CurveFitting();
    // double result = trapezoidalIntegration(parser, "exp(0-x^2)", 0, 1, 10);
    // std::cout << std::setprecision(10) << "Root = " << result << '\n';

    return 0;
}
