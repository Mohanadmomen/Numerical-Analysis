#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <cmath>
#include <cctype>
#include <stdexcept>
#include <map>
#include <functional>
#pragma once 

class ExpressionParser {
private:
    enum TokenType {
        NUMBER,
        VARIABLE,
        OPERATOR,
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

    // Small epsilon value for comparing doubles
    const double EPSILON = 1e-12;

    // Operator precedence
    std::map<char, int> precedence = {
        {'+', 1}, {'-', 1},
        {'*', 2}, {'/', 2}, {'%', 2},
        {'^', 3}
    };

    // Helper function to check if a value is approximately equal to a multiple of PI
    bool isMultipleOfPi(double value) {
        double ratio = value / PI;
        double nearestInt = round(ratio);
        return fabs(ratio - nearestInt) < EPSILON;
    }

    void setupFunctions() {
        // Modified sin function to handle multiples of PI exactly
        functions["sin"] = [this](const std::vector<double>& args) {
            // Check if the argument is close to a multiple of PI
            if (isMultipleOfPi(args[0])) {
                // sin(n*PI) is exactly 0 for integer n
                return 0.0;
            }
            return sin(args[0]);
            };
        functionArgCount["sin"] = 1;

        // Modified cos function to handle multiples of PI exactly
        functions["cos"] = [this](const std::vector<double>& args) {
            // Check if the argument is close to a multiple of PI
            if (isMultipleOfPi(args[0])) {
                // cos(n*PI) is 1 for even n and -1 for odd n
                return (fmod(round(args[0] / PI), 2.0) == 0) ? 1.0 : -1.0;
            }
            return cos(args[0]);
            };
        functionArgCount["cos"] = 1;

        functions["tan"] = [this](const std::vector<double>& args) {
            // Check if the argument is close to a multiple of PI
            if (isMultipleOfPi(args[0])) {
                // tan(n*PI) is exactly 0 for integer n
                return 0.0;
            }
            return tan(args[0]);
            };
        functionArgCount["tan"] = 1;

        functions["sqrt"] = [](const std::vector<double>& args) { return sqrt(args[0]); };
        functionArgCount["sqrt"] = 1;

        functions["log"] = [](const std::vector<double>& args) { return log10(args[0]); };
        functionArgCount["log"] = 1;

        functions["ln"] = [](const std::vector<double>& args) { return log(args[0]); };
        functionArgCount["ln"] = 1;

        functions["exp"] = [](const std::vector<double>& args) { return exp(args[0]); };
        functionArgCount["exp"] = 1;

        functions["abs"] = [](const std::vector<double>& args) { return fabs(args[0]); };
        functionArgCount["abs"] = 1;

        functions["pow"] = [](const std::vector<double>& args) { return pow(args[0], args[1]); };
        functionArgCount["pow"] = 2;
    }


    std::vector<Token> tokenize(const std::string& expression) {
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
                tokens.push_back({ OPERATOR, op });
            }
            else {
                throw std::runtime_error("Invalid character in expression: " + std::string(1, c));
            }
        }

        return tokens;
    }

    std::vector<Token> shuntingYard(const std::vector<Token>& tokens) {
        std::vector<Token> output;
        std::stack<Token> operators;
        std::stack<int> argCounts; // To keep track of function argument counts

        for (size_t i = 0; i < tokens.size(); ++i) {
            const auto& token = tokens[i];

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
            case OPERATOR: {
                char op = token.value[0];
                while (!operators.empty() && operators.top().type == OPERATOR &&
                    ((op != '^' && precedence[op] <= precedence[operators.top().value[0]]) ||
                        (op == '^' && precedence[op] < precedence[operators.top().value[0]]))) {
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
                    throw std::runtime_error("Mismatched parentheses");
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
                throw std::runtime_error("Mismatched parentheses");
            }
            output.push_back(operators.top());
            operators.pop();
        }

        return output;
    }

    double evaluateRPN(const std::vector<Token>& rpn) {
        std::stack<double> values;

        for (const auto& token : rpn) {
            switch (token.type) {
            case NUMBER:
                values.push(std::stod(token.value));
                break;
            case VARIABLE:
                if (variables.find(token.value) != variables.end()) {
                    values.push(variables[token.value]);
                }
                else {
                    throw std::runtime_error("Undefined variable: " + token.value);
                }
                break;
            case OPERATOR: {
                if (values.size() < 2) {
                    throw std::runtime_error("Invalid expression");
                }

                double b = values.top(); values.pop();
                double a = values.top(); values.pop();

                switch (token.value[0]) {
                case '+': values.push(a + b); break;
                case '-': values.push(a - b); break;
                case '*': values.push(a * b); break;
                case '/':
                    if (b == 0) throw std::runtime_error("Division by zero");
                    values.push(a / b);
                    break;
                case '%': values.push(fmod(a, b)); break;
                case '^': values.push(pow(a, b)); break;
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
                    throw std::runtime_error("Not enough arguments for function: " + funcName);
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
                    throw std::runtime_error("Unknown function: " + funcName);
                }
                break;
            }
            default:
                throw std::runtime_error("Invalid token in RPN");
            }
        }

        if (values.size() != 1) {
            throw std::runtime_error("Invalid expression");
        }

        return values.top();
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
        auto tokens = tokenize(expression);
        auto rpn = shuntingYard(tokens);
        return evaluateRPN(rpn);
    }

    // Evaluate expression with specific x and y values
    double evaluate(const std::string& expression, double x_val, double y_val) {
        // Set x and y variables before evaluation
        setVariable("x", x_val);
        setVariable("y", y_val);
        return evaluate(expression);
    }
};