#ifndef FUNCTIONHANDLER_H
#define FUNCTIONHANDLER_H

#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <cmath>
#include <cctype>
#include <stdexcept>
#include <map>
#include <functional>
#include <limits>

class ExpressionParser {
private:
    enum TokenType {
        NUMBER,
        VARIABLE,
        OPERATOR,
        UNARY_OPERATOR, // Add unary operator type
        LEFT_PAREN,
        RIGHT_PAREN,
        FUNCTION,
        COMMA
    };

    struct Token {
        TokenType type;
        std::string value;
    };

    std::map<std::string, double> variables;
    std::map<std::string, std::function<double(const std::vector<double>&)>> functions;
    std::map<std::string, int> functionArgCount;

    // Define constants that might not be available
    const double PI = 3.14159265358979323846;
    const double E = 2.71828182845904523536;

    // Value to return when errors occur
    const double ERROR_VALUE = std::numeric_limits<double>::quiet_NaN();

    // Small epsilon value for comparing doubles
    const double EPSILON = 1e-12;

    // Operator precedence
    std::map<char, int> precedence = {
        {'+', 1}, {'-', 1},
        {'*', 2}, {'/', 2}, {'%', 2},
        {'^', 3}
    };

    // Higher precedence for unary operators
    const int UNARY_PRECEDENCE = 4;

    // Helper function to check if a value is approximately equal to a multiple of PI
    bool isMultipleOfPi(double value) {
        double ratio = value / PI;
        double nearestInt = round(ratio);
        return fabs(ratio - nearestInt) < EPSILON;
    }

    // Helper function to check if a value is approximately equal to 0
    bool isApproxZero(double value) {
        return fabs(value) < EPSILON;
    }

    // Helper function to check if a value is approximately equal to 1
    bool isApproxOne(double value) {
        return fabs(value - 1.0) < EPSILON;
    }

    void setupFunctions() {
        // Modified sin function to handle multiples of PI exactly
        functions["sin"] = [this](const std::vector<double>& args) {
            try {
                // Check if the argument is close to a multiple of PI
                if (isMultipleOfPi(args[0])) {
                    // sin(n*PI) is exactly 0 for integer n
                    return 0.0;
                }
                return sin(args[0]);
            }
            catch (...) {
                return ERROR_VALUE;
            }
            };
        functionArgCount["sin"] = 1;

        // Modified cos function to handle multiples of PI exactly
        functions["cos"] = [this](const std::vector<double>& args) {
            try {
                // Check if the argument is close to a multiple of PI
                if (isMultipleOfPi(args[0])) {
                    // cos(n*PI) is 1 for even n and -1 for odd n
                    return (fmod(round(args[0] / PI), 2.0) == 0) ? 1.0 : -1.0;
                }
                return cos(args[0]);
            }
            catch (...) {
                return ERROR_VALUE;
            }
            };
        functionArgCount["cos"] = 1;

        functions["tan"] = [this](const std::vector<double>& args) {
            try {
                // Check if the argument is close to a multiple of PI
                if (isMultipleOfPi(args[0])) {
                    // tan(n*PI) is exactly 0 for integer n
                    return 0.0;
                }
                // Check if arg is close to PI/2 + n*PI
                double halfPiMultiple = fmod(args[0] - PI / 2, PI);
                if (isApproxZero(halfPiMultiple)) {
                    return ERROR_VALUE; // tan is undefined at these points
                }
                return tan(args[0]);
            }
            catch (...) {
                return ERROR_VALUE;
            }
            };
        functionArgCount["tan"] = 1;

        // Add cotangent function
        functions["cot"] = [this](const std::vector<double>& args) {
            try {
                // Check if the argument is close to a multiple of PI
                if (isMultipleOfPi(args[0])) {
                    // cot(n*PI) would be 1/0, which is undefined
                    // Return a special value instead of throwing
                    return ERROR_VALUE;
                }
                double tanVal = tan(args[0]);
                if (isApproxZero(tanVal)) {
                    return ERROR_VALUE; // Division by zero
                }
                return 1.0 / tanVal;
            }
            catch (...) {
                return ERROR_VALUE;
            }
            };
        functionArgCount["cot"] = 1;

        // Add cosecant function
        functions["cosec"] = [this](const std::vector<double>& args) {
            // Alternative name
            return functions["csc"](args);
            };
        functionArgCount["cosec"] = 1;

        // Add cosecant function (csc)
        functions["csc"] = [this](const std::vector<double>& args) {
            try {
                // Check if the argument is close to a multiple of PI
                if (isMultipleOfPi(args[0])) {
                    // csc(n*PI) would be 1/0, which is undefined
                    return ERROR_VALUE;
                }
                double sinVal = sin(args[0]);
                if (isApproxZero(sinVal)) {
                    return ERROR_VALUE; // Division by zero
                }
                return 1.0 / sinVal;
            }
            catch (...) {
                return ERROR_VALUE;
            }
            };
        functionArgCount["csc"] = 1;

        // Add secant function (sec)
        functions["sec"] = [this](const std::vector<double>& args) {
            try {
                // Check if the argument is close to a multiple of PI/2 + PI*n
                double halfPiMultiple = fmod(args[0], PI);
                if (isApproxZero(fmod(halfPiMultiple + PI / 2, PI))) {
                    return ERROR_VALUE; // Undefined value
                }
                double cosVal = cos(args[0]);
                if (isApproxZero(cosVal)) {
                    return ERROR_VALUE; // Division by zero
                }
                return 1.0 / cosVal;
            }
            catch (...) {
                return ERROR_VALUE;
            }
            };
        functionArgCount["sec"] = 1;

        // Add inverse trigonometric functions
        functions["asin"] = [this](const std::vector<double>& args) {
            try {
                if (args[0] < -1.0 || args[0] > 1.0) {
                    return ERROR_VALUE; // Out of domain
                }
                // Handle special cases for exact values
                if (isApproxZero(args[0])) return 0.0;
                if (isApproxOne(args[0])) return PI / 2;
                if (isApproxOne(-args[0])) return -PI / 2;
                if (isApproxZero(args[0] - 0.5)) return PI / 6;
                if (isApproxZero(args[0] + 0.5)) return -PI / 6;
                if (isApproxZero(args[0] - (sqrt(2) / 2))) return PI / 4;
                if (isApproxZero(args[0] + (sqrt(2) / 2))) return -PI / 4;
                if (isApproxZero(args[0] - (sqrt(3) / 2))) return PI / 3;
                if (isApproxZero(args[0] + (sqrt(3) / 2))) return -PI / 3;

                return asin(args[0]);
            }
            catch (...) {
                return ERROR_VALUE;
            }
            };
        functionArgCount["asin"] = 1;
        functions["arcsin"] = functions["asin"]; // Alternative name
        functionArgCount["arcsin"] = 1;

        functions["acos"] = [this](const std::vector<double>& args) {
            try {
                if (args[0] < -1.0 || args[0] > 1.0) {
                    return ERROR_VALUE; // Out of domain
                }
                // Handle special cases for exact values
                if (isApproxZero(args[0])) return PI / 2;
                if (isApproxOne(args[0])) return 0.0;
                if (isApproxOne(-args[0])) return PI;
                if (isApproxZero(args[0] - 0.5)) return PI / 3;
                if (isApproxZero(args[0] + 0.5)) return 2 * PI / 3;
                if (isApproxZero(args[0] - (sqrt(2) / 2))) return PI / 4;
                if (isApproxZero(args[0] + (sqrt(2) / 2))) return 3 * PI / 4;
                if (isApproxZero(args[0] - (sqrt(3) / 2))) return PI / 6;
                if (isApproxZero(args[0] + (sqrt(3) / 2))) return 5 * PI / 6;

                return acos(args[0]);
            }
            catch (...) {
                return ERROR_VALUE;
            }
            };
        functionArgCount["acos"] = 1;
        functions["arccos"] = functions["acos"]; // Alternative name
        functionArgCount["arccos"] = 1;

        functions["atan"] = [this](const std::vector<double>& args) {
            try {
                // Handle special cases for exact values
                if (isApproxZero(args[0])) return 0.0;
                if (isApproxZero(args[0] - 1.0)) return PI / 4;
                if (isApproxZero(args[0] + 1.0)) return -PI / 4;
                if (isApproxZero(args[0] - sqrt(3))) return PI / 3;
                if (isApproxZero(args[0] + sqrt(3))) return -PI / 3;
                if (isApproxZero(args[0] - (1 / sqrt(3)))) return PI / 6;
                if (isApproxZero(args[0] + (1 / sqrt(3)))) return -PI / 6;

                return atan(args[0]);
            }
            catch (...) {
                return ERROR_VALUE;
            }
            };
        functionArgCount["atan"] = 1;
        functions["arctan"] = functions["atan"]; // Alternative name
        functionArgCount["arctan"] = 1;

        // Two-argument arctangent
        functions["atan2"] = [this](const std::vector<double>& args) {
            try {
                // Handle the case where both arguments are zero
                if (isApproxZero(args[0]) && isApproxZero(args[1])) {
                    return ERROR_VALUE; // atan2(0,0) is undefined
                }
                return atan2(args[0], args[1]);
            }
            catch (...) {
                return ERROR_VALUE;
            }
            };
        functionArgCount["atan2"] = 2;

        // Add arccotangent function
        functions["acot"] = [this](const std::vector<double>& args) {
            try {
                // arccot(x) = PI/2 - arctan(x)
                if (isApproxZero(args[0])) return PI / 2;
                return PI / 2 - atan(args[0]);
            }
            catch (...) {
                return ERROR_VALUE;
            }
            };
        functionArgCount["acot"] = 1;
        functions["arccot"] = functions["acot"]; // Alternative name
        functionArgCount["arccot"] = 1;

        // Add arcsecant function
        functions["asec"] = [this](const std::vector<double>& args) {
            try {
                if (args[0] > -1.0 && args[0] < 1.0) {
                    return ERROR_VALUE; // Out of domain
                }
                // arcsec(x) = arccos(1/x)
                return acos(1.0 / args[0]);
            }
            catch (...) {
                return ERROR_VALUE;
            }
            };
        functionArgCount["asec"] = 1;
        functions["arcsec"] = functions["asec"]; // Alternative name
        functionArgCount["arcsec"] = 1;

        // Add arccosecant function
        functions["acsc"] = [this](const std::vector<double>& args) {
            try {
                if (args[0] > -1.0 && args[0] < 1.0) {
                    return ERROR_VALUE; // Out of domain
                }
                // arccsc(x) = arcsin(1/x)
                return asin(1.0 / args[0]);
            }
            catch (...) {
                return ERROR_VALUE;
            }
            };
        functionArgCount["acsc"] = 1;
        functions["arccsc"] = functions["acsc"]; // Alternative name
        functionArgCount["arccsc"] = 1;
        functions["acosec"] = functions["acsc"]; // Another alternative name
        functionArgCount["acosec"] = 1;
        functions["arccosec"] = functions["acsc"]; // Another alternative name
        functionArgCount["arccosec"] = 1;

        // Modified mathematical functions with error handling
        functions["sqrt"] = [this](const std::vector<double>& args) {
            try {
                if (args[0] < 0) {
                    return ERROR_VALUE; // Negative input
                }
                return sqrt(args[0]);
            }
            catch (...) {
                return ERROR_VALUE;
            }
            };
        functionArgCount["sqrt"] = 1;

        functions["log"] = [this](const std::vector<double>& args) {
            try {
                if (args[0] <= 0) {
                    return ERROR_VALUE; // Non-positive input
                }
                return log10(args[0]);
            }
            catch (...) {
                return ERROR_VALUE;
            }
            };
        functionArgCount["log"] = 1;

        functions["ln"] = [this](const std::vector<double>& args) {
            try {
                if (args[0] <= 0) {
                    return ERROR_VALUE; // Non-positive input
                }
                return log(args[0]);
            }
            catch (...) {
                return ERROR_VALUE;
            }
            };
        functionArgCount["ln"] = 1;

        functions["exp"] = [this](const std::vector<double>& args) {
            try {
                return exp(args[0]);
            }
            catch (...) {
                return ERROR_VALUE;
            }
            };
        functionArgCount["exp"] = 1;

        functions["abs"] = [this](const std::vector<double>& args) {
            try {
                return fabs(args[0]);
            }
            catch (...) {
                return ERROR_VALUE;
            }
            };
        functionArgCount["abs"] = 1;

        functions["pow"] = [this](const std::vector<double>& args) {
            try {
                // Handle special cases
                if (isApproxZero(args[0]) && args[1] < 0) {
                    return ERROR_VALUE; // 0^negative is undefined
                }
                if (args[0] < 0 && fmod(args[1], 1.0) != 0) {
                    return ERROR_VALUE; // Negative base with non-integer exponent
                }
                return pow(args[0], args[1]);
            }
            catch (...) {
                return ERROR_VALUE;
            }
            };
        functionArgCount["pow"] = 2;
    }


    std::vector<Token> tokenize(const std::string& expression) {
        try {
            std::vector<Token> tokens;
            std::string current;

            for (size_t i = 0; i < expression.length(); ++i) {
                char c = expression[i];

                if (std::isspace(c)) {
                    continue;
                }
                else if (std::isdigit(c) || c == '.') {
                    current.clear();
                    while (i < expression.length() && (std::isdigit(expression[i]) || expression[i] == '.')) {
                        current += expression[i++];
                    }
                    --i;
                    tokens.push_back({ NUMBER, current });
                }
                else if (std::isalpha(c)) {
                    current.clear();
                    while (i < expression.length() && (std::isalpha(expression[i]) || std::isdigit(expression[i]))) {
                        current += expression[i++];
                    }
                    --i;

                    if (functions.find(current) != functions.end()) {
                        tokens.push_back({ FUNCTION, current });

                        // Check if the next character is ^ and handle "exp^x" case
                        size_t j = i + 1;
                        while (j < expression.length() && std::isspace(expression[j])) j++;

                        if (j < expression.length() && expression[j] == '^') {
                            // If it's exp^x, convert to exp(x)
                            if (current == "exp") {
                                // Skip the ^ character
                                i = j;
                                tokens.pop_back(); // Remove the FUNCTION token

                                // Convert to pow(e, x)
                                tokens.push_back({ FUNCTION, "pow" });
                                tokens.push_back({ LEFT_PAREN, "(" });
                                tokens.push_back({ NUMBER, std::to_string(E) }); // e constant
                                tokens.push_back({ COMMA, "," });
                            }
                        }
                    }
                    else {
                        tokens.push_back({ VARIABLE, current });
                    }
                }
                else if (c == '(') {
                    tokens.push_back({ LEFT_PAREN, "(" });
                }
                else if (c == ')') {
                    tokens.push_back({ RIGHT_PAREN, ")" });
                }
                else if (c == ',') {
                    tokens.push_back({ COMMA, "," });
                }
                else if (precedence.find(c) != precedence.end() || c == '^') {
                    std::string op(1, c);

                    // Handle unary minus (and potentially unary plus)
                    if ((c == '-' || c == '+') &&
                        (tokens.empty() ||
                            tokens.back().type == LEFT_PAREN ||
                            tokens.back().type == COMMA ||
                            tokens.back().type == OPERATOR ||
                            tokens.back().type == UNARY_OPERATOR)) {
                        // This is a unary operator
                        tokens.push_back({ UNARY_OPERATOR, op });
                    }
                    else {
                        tokens.push_back({ OPERATOR, op });
                    }
                }
                else {
                    // Instead of throwing, return a special token that will indicate error
                    //std::cerr << "Warning: Invalid character in expression: " << c << std::endl;
                    tokens.push_back({ NUMBER, "NaN" });
                    return tokens;
                }
            }

            return tokens;
        }
        catch (...) {
            // Return a token list with just NaN if any exception occurs
            std::vector<Token> errorTokens;
            errorTokens.push_back({ NUMBER, "NaN" });
            return errorTokens;
        }
    }

    std::vector<Token> shuntingYard(const std::vector<Token>& tokens) {
        try {
            std::vector<Token> output;
            std::stack<Token> operators;
            std::stack<int> argCounts; // To keep track of function argument counts

            for (size_t i = 0; i < tokens.size(); ++i) {
                const auto& token = tokens[i];

                // Check if we already have an error token
                if (token.type == NUMBER && token.value == "NaN") {
                    output.push_back(token);
                    return output;
                }

                switch (token.type) {
                case NUMBER:
                case VARIABLE:
                    output.push_back(token);
                    break;
                case FUNCTION:
                    operators.push(token);
                    argCounts.push(0); // Initialize argument count
                    break;
                case COMMA:
                    while (!operators.empty() && operators.top().type != LEFT_PAREN) {
                        output.push_back(operators.top());
                        operators.pop();
                    }

                    // Increment argument count for the current function
                    if (!argCounts.empty()) {
                        int count = argCounts.top();
                        argCounts.pop();
                        argCounts.push(count + 1);
                    }
                    break;
                case UNARY_OPERATOR:
                    // Unary operators have the highest precedence
                    operators.push(token);
                    break;
                case OPERATOR: {
                    char op = token.value[0];
                    while (!operators.empty() &&
                        (operators.top().type == OPERATOR || operators.top().type == UNARY_OPERATOR) &&
                        ((op != '^' &&
                            (operators.top().type == UNARY_OPERATOR ?
                                UNARY_PRECEDENCE : precedence[operators.top().value[0]]) >= precedence[op]) ||
                            (op == '^' &&
                                (operators.top().type == UNARY_OPERATOR ?
                                    UNARY_PRECEDENCE : precedence[operators.top().value[0]]) > precedence[op]))) {
                        output.push_back(operators.top());
                        operators.pop();
                    }
                    operators.push(token);
                    break;
                }
                case LEFT_PAREN:
                    operators.push(token);
                    break;
                case RIGHT_PAREN:
                    while (!operators.empty() && operators.top().type != LEFT_PAREN) {
                        output.push_back(operators.top());
                        operators.pop();
                    }

                    if (!operators.empty() && operators.top().type == LEFT_PAREN) {
                        operators.pop();
                    }
                    else {
                        // Mismatched parentheses
                        //std::cerr << "Warning: Mismatched parentheses" << std::endl;
                        output.clear();
                        output.push_back({ NUMBER, "NaN" });
                        return output;
                    }

                    // If we have a function call
                    if (!operators.empty() && operators.top().type == FUNCTION) {
                        // Get the argument count (plus one for the last argument)
                        int count = 1;
                        if (!argCounts.empty()) {
                            count += argCounts.top();
                            argCounts.pop();
                        }

                        // Store the function and its argument count
                        Token funcToken = operators.top();
                        funcToken.value += "," + std::to_string(count); // Append arg count to function name
                        output.push_back(funcToken);
                        operators.pop();
                    }
                    break;
                }
            }

            while (!operators.empty()) {
                if (operators.top().type == LEFT_PAREN || operators.top().type == RIGHT_PAREN) {
                    // Mismatched parentheses
                    //std::cerr << "Warning: Mismatched parentheses" << std::endl;
                    output.clear();
                    output.push_back({ NUMBER, "NaN" });
                    return output;
                }
                output.push_back(operators.top());
                operators.pop();
            }

            return output;
        }
        catch (...) {
            // Return a token list with just NaN if any exception occurs
            std::vector<Token> errorTokens;
            errorTokens.push_back({ NUMBER, "NaN" });
            return errorTokens;
        }
    }

    double evaluateRPN(const std::vector<Token>& rpn) {
        try {
            std::stack<double> values;

            for (const auto& token : rpn) {
                // Check if we already have an error token
                if (token.type == NUMBER && token.value == "NaN") {
                    return ERROR_VALUE;
                }

                switch (token.type) {
                case NUMBER:
                    try {
                        values.push(std::stod(token.value));
                    }
                    catch (...) {
                        return ERROR_VALUE;
                    }
                    break;
                case VARIABLE:
                    if (variables.find(token.value) != variables.end()) {
                        values.push(variables[token.value]);
                    }
                    else {
                        //std::cerr << "Warning: Undefined variable: " << token.value << std::endl;
                        return ERROR_VALUE;
                    }
                    break;
                case UNARY_OPERATOR: {
                    if (values.empty()) {
                        //std::cerr << "Warning: Invalid expression (not enough operands for unary operator)" << std::endl;
                        return ERROR_VALUE;
                    }

                    double a = values.top(); values.pop();

                    switch (token.value[0]) {
                    case '+': values.push(a); break;      // Unary plus (no effect)
                    case '-': values.push(-a); break;     // Unary minus (negation)
                    default:
                        //std::cerr << "Warning: Unknown unary operator: " << token.value << std::endl;
                        return ERROR_VALUE;
                    }
                    break;
                }
                case OPERATOR: {
                    if (values.size() < 2) {
                        //std::cerr << "Warning: Invalid expression (not enough operands)" << std::endl;
                        return ERROR_VALUE;
                    }

                    double b = values.top(); values.pop();
                    double a = values.top(); values.pop();

                    switch (token.value[0]) {
                    case '+': values.push(a + b); break;
                    case '-': values.push(a - b); break;
                    case '*': values.push(a * b); break;
                    case '/':
                        if (isApproxZero(b)) {
                            //std::cerr << "Warning: Division by zero" << std::endl;
                            return ERROR_VALUE;
                        }
                        values.push(a / b);
                        break;
                    case '%':
                        if (isApproxZero(b)) {
                            //std::cerr << "Warning: Modulo by zero" << std::endl;
                            return ERROR_VALUE;
                        }
                        values.push(fmod(a, b));
                        break;
                    case '^':
                        try {
                            // Special cases for power
                            if (isApproxZero(a) && b < 0) {
                                //std::cerr << "Warning: Zero raised to negative power" << std::endl;
                                return ERROR_VALUE;
                            }
                            if (a < 0 && fmod(b, 1.0) != 0) {
                                //std::cerr << "Warning: Negative base with non-integer exponent" << std::endl;
                                return ERROR_VALUE;
                            }
                            values.push(pow(a, b));
                        }
                        catch (...) {
                            return ERROR_VALUE;
                        }
                        break;
                    }
                    break;
                }
                case FUNCTION: {
                    // Parse function name and argument count
                    size_t commaPos = token.value.find(',');
                    std::string funcName = token.value;
                    int argCount = 1; // Default to 1 argument

                    if (commaPos != std::string::npos) {
                        funcName = token.value.substr(0, commaPos);
                        argCount = std::stoi(token.value.substr(commaPos + 1));
                    }
                    else if (functionArgCount.find(funcName) != functionArgCount.end()) {
                        argCount = functionArgCount[funcName];
                    }

                    // Check if we have enough values on the stack
                    if (values.size() < argCount) {
                        //std::cerr << "Warning: Not enough arguments for function: " << funcName << std::endl;
                        return ERROR_VALUE;
                    }

                    // Collect arguments
                    std::vector<double> args(argCount);
                    for (int i = argCount - 1; i >= 0; --i) {
                        args[i] = values.top();
                        values.pop();
                    }

                    // Apply function
                    if (functions.find(funcName) != functions.end()) {
                        values.push(functions[funcName](args));
                    }
                    else {
                        //std::cerr << "Warning: Unknown function: " << funcName << std::endl;
                        return ERROR_VALUE;
                    }
                    break;
                }
                default:
                    //std::cerr << "Warning: Invalid token in RPN" << std::endl;
                    return ERROR_VALUE;
                }
            }

            if (values.size() != 1) {
                //std::cerr << "Warning: Invalid expression (incorrect number of values at end)" << std::endl;
                return ERROR_VALUE;
            }

            return values.top();
        }
        catch (...) {
            return ERROR_VALUE;
        }
    }

public:
    ExpressionParser() {
        setupFunctions();

        // Set mathematical constants
        variables["pi"] = PI;
        variables["e"] = E;
    }

    void setVariable(const std::string& name, double value) {
        variables[name] = value;
    }

    double evaluate(const std::string& expression) {
        try {
            // Replace any standalone negative symbols at the beginning of the expression
            std::string processedExpr = expression;

            // For debugging
            // std::cout << "Evaluating: " << processedExpr << std::endl;

            auto tokens = tokenize(processedExpr);
            auto rpn = shuntingYard(tokens);
            return evaluateRPN(rpn);
        }
        catch (...) {
            return ERROR_VALUE;
        }
    }

    // Evaluate expression with specific x and y values
    double evaluate(const std::string& expression, double x_val, double y_val) {
        // Set x and y variables before evaluation
        setVariable("x", x_val);
        setVariable("y", y_val);
        return evaluate(expression);
    }

    // Check if a value is NaN (error value)
    bool isError(double value) {
        return std::isnan(value);
    }
};

#endif // FUNCTIONHANDLER_H
