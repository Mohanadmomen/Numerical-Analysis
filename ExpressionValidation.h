#ifndef EXPRESSIONVALIDATION_H
#define EXPRESSIONVALIDATION_H

#include <iostream>
#include <string>
#include <stack>
#include <set>
#include <cmath>

bool isValidFunction(const std::string& name) {
    // List of valid function names - only these exact function names are allowed
    static const std::set<std::string> validFunctions = {
        // Basic trigonometric functions
        "sin", "cos", "tan", "cot", "cosec", "csc", "sec",

        // Inverse trigonometric functions
        "asin", "arcsin", "acos", "arccos", "atan", "arctan", "atan2",
        "acot", "arccot", "asec", "arcsec", "acsc", "arccsc", "acosec", "arccosec",

        // Mathematical functions
        "sqrt", "log", "ln", "exp", "abs", "pow"
    };

    return validFunctions.find(name) != validFunctions.end();
}

bool isValidVariable(const std::string& name) {
    // Only allow 'x' and 'y' and 'Pi' as variables
    return name == "x" || name == "y" || name == "pi";
}

bool validateExpression(const std::string& expression) {
    if (expression.empty()) {
        return false;
    }

    // Check for capital X or Y and return 0
    for (size_t i = 0; i < expression.length(); ++i) {
        if (expression[i] == 'X' || expression[i] == 'Y') {
            return false;
        }
    }

    // Check for balanced parentheses
    std::stack<char> parentheses;
    bool expectOperand = true;       // Flag to track if we expect an operand next
    bool expectOperator = false;     // Flag to track if we expect an operator next
    bool lastWasOperand = false;     // Track if last token was an operand (number/variable)
    bool lastWasOperator = false;    // Track if last token was an operator
    bool lastWasFunction = false;    // Track if last token was a function name
    bool lastWasCloseParen = false;  // Track if last token was a closing parenthesis
    std::stack<int> functionCommas;  // Track comma count for nested functions
    bool emptyFunction = false;      // Track if we have an empty function call

    for (size_t i = 0; i < expression.length(); ++i) {
        char c = expression[i];

        // Skip whitespace
        if (std::isspace(c)) {
            continue;
        }

        // Handle numbers
        if (std::isdigit(c) || c == '.') {
            if (!expectOperand && !lastWasFunction && !lastWasOperator) {
                // Implicit multiplication like 2x is not allowed
                if (lastWasOperand || lastWasCloseParen) {
                    return false;
                }
            }

            // Start of number
            // Consume the entire number
            i++;
            while (i < expression.length() && (std::isdigit(expression[i]) || expression[i] == '.')) {
                i++;
            }
            i--; // Adjust for the loop increment

            lastWasOperand = true;
            lastWasOperator = false;
            lastWasFunction = false;
            lastWasCloseParen = false;
            expectOperand = false;
            expectOperator = true;
            continue;
        }

        // Handle variables and functions
        if (std::isalpha(c)) {
            if (!expectOperand && !lastWasFunction && !lastWasOperator) {
                // Implicit multiplication like 2x is not allowed
                if (lastWasOperand || lastWasCloseParen) {
                    return false;
                }
            }

            std::string identifier;
            // Consume the entire identifier
            while (i < expression.length() && (std::isalpha(expression[i]) || std::isdigit(expression[i]))) {
                identifier += expression[i++];
            }
            i--; // Adjust for the loop increment

            // Check if next non-whitespace character is '(' to determine if this is a function
            size_t j = i + 1;
            while (j < expression.length() && std::isspace(expression[j])) j++;

            if (j < expression.length() && expression[j] == '(') {
                // This is a function call
                // Check if it's a valid function name
                if (!isValidFunction(identifier)) {
                    return false;
                }

                lastWasFunction = true;
                lastWasOperand = false;
                lastWasOperator = false;
                lastWasCloseParen = false;
                functionCommas.push(0); // Initialize comma count for this function
            }
            else {
                // This is a variable - check if it's a valid variable name (x or y)
                if (!isValidVariable(identifier)) {
                    return false;
                }

                lastWasFunction = false;
                lastWasOperand = true;
                lastWasOperator = false;
                lastWasCloseParen = false;
                expectOperand = false;
                expectOperator = true;
            }
            continue;
        }

        // Handle open parenthesis
        if (c == '(') {
            if (!expectOperand && !lastWasFunction && !lastWasOperator) {
                // Implicit multiplication like x(3) is not allowed
                if (lastWasOperand || lastWasCloseParen) {
                    return false;
                }
            }

            parentheses.push(c);

            // Check for empty function immediately
            size_t j = i + 1;
            while (j < expression.length() && std::isspace(expression[j])) j++;

            if (j < expression.length() && expression[j] == ')') {
                // Found an empty function like tan()
                emptyFunction = true;
            }

            // Reset expectations within the new scope
            expectOperand = true;
            expectOperator = false;
            lastWasOperator = false;
            lastWasOperand = false;
            lastWasFunction = false;
            lastWasCloseParen = false;

            continue;
        }

        // Handle close parenthesis
        if (c == ')') {
            if (parentheses.empty()) {
                return false;
            }

            // Cannot have a closing parenthesis right after an operator or comma
            if (lastWasOperator) {
                return false;
            }

            parentheses.pop();

            lastWasOperator = false;
            lastWasOperand = false;
            lastWasFunction = false;
            lastWasCloseParen = true;
            expectOperand = false;
            expectOperator = true;

            // Update function comma tracking
            if (!functionCommas.empty()) {
                functionCommas.pop();
            }

            continue;
        }

        // Handle commas (for function arguments)
        if (c == ',') {
            if (parentheses.empty()) {
                return false;
            }

            if (lastWasOperator) {
                return false;
            }

            // Track comma count for the current function call
            if (!functionCommas.empty()) {
                int count = functionCommas.top();
                functionCommas.pop();
                functionCommas.push(count + 1);
            }

            expectOperand = true;
            expectOperator = false;
            lastWasOperator = true; // Treat comma similar to operator for syntax checking
            lastWasOperand = false;
            lastWasFunction = false;
            lastWasCloseParen = false;
            continue;
        }

        // Handle operators
        if (c == '+' || c == '-' || c == '*' || c == '/' || c == '%' || c == '^') {
            // Special case for unary +/- (can appear after other operators or at start)
            if ((c == '+' || c == '-') && (i == 0 || expectOperand || lastWasOperator || expression[i - 1] == '(' || expression[i - 1] == ',')) {
                // This is a unary operator
                expectOperand = true;
                lastWasOperator = true;
                lastWasOperand = false;
                lastWasFunction = false;
                lastWasCloseParen = false;
                continue;
            }

            // Binary operators must be preceded by an operand
            if (expectOperand) {
                return false;
            }

            expectOperand = true;
            expectOperator = false;
            lastWasOperator = true;
            lastWasOperand = false;
            lastWasFunction = false;
            lastWasCloseParen = false;
            continue;
        }

        // If we get here, we encountered an invalid character
        return false;
    }

    // Check for unclosed parentheses
    if (!parentheses.empty()) {
        return false;
    }

    // Check that the expression doesn't end with an operator
    if (lastWasOperator && !expression.empty()) {
        return false;
    }

    // Check for empty functions
    if (emptyFunction) {
        return false;
    }

    return true;
}

bool validateInterpolationInput(const std::vector<double>& x_values, const std::vector<double>& y_values) {
    // Check if vectors are empty and have the same size
    if (x_values.empty() || y_values.empty() || x_values.size() != y_values.size()) {
        return false;
    }

    // All x's must be distinct
    for (size_t i = 0; i < x_values.size(); ++i) {
        for (size_t j = i + 1; j < x_values.size(); ++j) {
            if (std::fabs(x_values[i] - x_values[j]) < 1e-9) {
                return false;
            }
        }
    }

    // All y's must be distinct
    for (size_t i = 0; i < y_values.size(); ++i) {
        for (size_t j = i + 1; j < y_values.size(); ++j) {
            if (std::fabs(y_values[i] - y_values[j]) < 1e-9) {
                return false;
            }
        }
    }

    // All checks passed
    return true;
}

#endif // EXPRESSIONVALIDATION_H
