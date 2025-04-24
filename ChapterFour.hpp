#pragma once

#include <iostream>
#include <vector>
#include <functional>
#include <string>
#include <utility>

class ChapterFour {
private:
    ExpressionParser parser;

    // دالة مساعدة لتقييم الدالة عند نقطة محددة
    double evaluateFunction(const std::string& expression, double x, double y) {
        parser.setVariable("x", x);
        parser.setVariable("y", y);
        return parser.evaluate(expression);
    }

public:
    ChapterFour() {}

    // هيكل لتخزين نقاط الحل
    struct Solution {
        std::vector<double> x_values;
        std::vector<double> y_values;
    };

    // طريقة أويلر لحل معادلة تفاضلية من الدرجة الأولى: dy/dx = f(x,y)
    Solution euler(const std::string& odeExpression, double x0, double y0, double h, double xEnd) {
        Solution result;

        // تهيئة متجهات الحل بالشروط الابتدائية
        result.x_values.push_back(x0);
        result.y_values.push_back(y0);

        double x = x0;
        double y = y0;

        // الاستمرار حتى نصل أو نتجاوز نقطة النهاية
        while (x < xEnd - 1e-10) {
            // تعيين المتغيرات في المحلل
            parser.setVariable("x", x);
            parser.setVariable("y", y);

            // حساب الميل باستخدام f(x,y)
            double slope = parser.evaluate(odeExpression);

            // تحديث أويلر: y_next = y + h * f(x,y)
            y = y + h * slope;
            x = x + h;

            // تخزين نقاط الحل
            result.x_values.push_back(x);
            result.y_values.push_back(y);
        }

        return result;
    }

    // طريقة أويلر المعدلة (طريقة هيون) لحل معادلة تفاضلية من الدرجة الأولى: dy/dx = f(x,y)
    Solution modifiedEuler(const std::string& odeExpression, double x0, double y0, double h, double xEnd) {
        Solution result;

        // تهيئة متجهات الحل بالشروط الابتدائية
        result.x_values.push_back(x0);
        result.y_values.push_back(y0);

        double x = x0;
        double y = y0;

        // الاستمرار حتى نصل أو نتجاوز نقطة النهاية
        while (x < xEnd - 1e-10) {
            // تعيين المتغيرات في المحلل للنقطة الحالية
            parser.setVariable("x", x);
            parser.setVariable("y", y);

            // الخطوة 1: حساب k1 = f(x,y)
            double k1 = parser.evaluate(odeExpression);

            // الخطوة 2: حساب المتنبئ: y* = y + h * k1
            double y_pred = y + h * k1;

            // الخطوة 3: حساب k2 = f(x+h, y*)
            parser.setVariable("x", x + h);
            parser.setVariable("y", y_pred);
            double k2 = parser.evaluate(odeExpression);

            // الخطوة 4: تحديث أويلر المعدل: y_next = y + h * (k1 + k2) / 2
            y = y + h * (k1 + k2) / 2;
            x = x + h;

            // تخزين نقاط الحل
            result.x_values.push_back(x);
            result.y_values.push_back(y);
        }

        return result;
    }

    // طباعة الحل على وحدة التحكم
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

    // الحصول على سلسلة الحل المنسقة
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

    // مقارنة النتائج بين طريقتي أويلر وأويلر المعدلة
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