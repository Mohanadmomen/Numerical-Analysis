#include <iostream>
#include <iomanip>
#include "FunctionHandler.hpp"
#include "ChapterThree.hpp"
#include "ChapterTwo.hpp"
#include "ChapterOne.hpp"
int main() {
    ExpressionParser parser;

    // Example expression: 1/(x+1)
    

    std::string expression;
    std::cin >> expression;

    //// Generate C++ function
    //std::string cppFunction = parser.generateCppFunction(expression);
    //std::cout << "Generated C++ function:\n" << cppFunction << std::endl;

    

    double result = trapezoidalIntegration(parser, expression, 0,3,6);
    std::cout << std::setprecision(10) << result << '\n';
    //std::cout << n << "\t" << result << std::endl;







    //// Evaluate for different x values
    //std::cout << "Evaluating the expression for different x values:" << std::endl;
    //for (double x = 0; x <= 5; x += 1) {
    //    parser.setVariable("x", x);
    //    try {
    //        double result = parser.evaluate(expression);
    //        std::cout << "f(" << x << ") = " << result << std::endl;
    //    }
    //    catch (const std::exception& e) {
    //        std::cout << "Error: " << e.what() << std::endl;
    //    }
    //}

    return 0;
}