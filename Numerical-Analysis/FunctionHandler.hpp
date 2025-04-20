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

    // Operator precedence
    std::map<char, int> precedence = {
        {'+', 1}, {'-', 1},
        {'*', 2}, {'/', 2}, {'%', 2},
        {'^', 3}
    };

    void setupFunctions() {
        functions["sin"] = [](const std::vector<double>& args) { return sin(args[0]); };
        functions["cos"] = [](const std::vector<double>& args) { return cos(args[0]); };
        functions["tan"] = [](const std::vector<double>& args) { return tan(args[0]); };
        functions["sqrt"] = [](const std::vector<double>& args) { return sqrt(args[0]); };
        functions["log"] = [](const std::vector<double>& args) { return log(args[0]); };
        functions["exp"] = [](const std::vector<double>& args) { return exp(args[0]); };
        functions["abs"] = [](const std::vector<double>& args) { return fabs(args[0]); };
        functions["pow"] = [](const std::vector<double>& args) { return pow(args[0], args[1]); };
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
                while (i < expression.length() && std::isalpha(expression[i])) {
                    current += expression[i++];
                }
                --i;

                if (functions.find(current) != functions.end()) {
                    tokens.push_back({ FUNCTION, current });
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

        for (const auto& token : tokens) {
            switch (token.type) {
            case NUMBER:
            case VARIABLE:
                output.push_back(token);
                break;
            case FUNCTION:
                operators.push(token);
                break;
            case COMMA:
                while (!operators.empty() && operators.top().type != LEFT_PAREN) {
                    output.push_back(operators.top());
                    operators.pop();
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

                if (!operators.empty() && operators.top().type == FUNCTION) {
                    output.push_back(operators.top());
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
                auto func = functions[token.value];
                std::vector<double> args;
                args.push_back(values.top());
                values.pop();
                values.push(func(args));
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
    }

    void setVariable(const std::string& name, double value) {
        variables[name] = value;
    }

    double evaluate(const std::string& expression) {
        auto tokens = tokenize(expression);
        auto rpn = shuntingYard(tokens);
        return evaluateRPN(rpn);
    }

    // Generate C++ function from expression
    std::string generateCppFunction(const std::string& expression, const std::string& funcName = "mathFunction") {
        std::string code = "double " + funcName + "(double x) {\n";
        code += "    return " + expression + ";\n";
        code += "}\n";
        return code;
    }
};