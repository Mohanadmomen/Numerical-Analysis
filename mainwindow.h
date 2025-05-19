#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "FunctionHandler.h"

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void plotScatter(const std::vector<double>& x, const std::vector<double>& y);
    void plotCurve(double a, double b, const QString& methodName);
    void plotQuadraticCurve(double a, double b, double c, const QString& methodName);


    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////// \brief Root_Finding Tab Signals Declarations ////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    void on_functionLineEdit_textChanged(const QString &arg1);

    void on_bisectionRadioButton_toggled(bool checked);

    void on_secantRadioButton_toggled(bool checked);

    void on_newtonRadioButton_2_toggled(bool checked);

    void on_calculateButton_clicked();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////// \brief Interpolation Tab Signals Declarations ////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    void on_xValuesLineEdit_textChanged(const QString &arg1);

    void on_yValuesLineEdit_textChanged(const QString &arg1);

    void on_interpolateLineEdit_textChanged(const QString &arg1);

    void on_lagrangeRadioButton_toggled(bool checked);

    void on_newtonRadioButton_toggled(bool checked);

    void on_calculateButton_2_clicked();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////// \brief Integration Tab Signals Declarations /////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    void on_functionLineEdit_2_textChanged(const QString &arg1);

    void on_lowerLimitLineEdit_textChanged(const QString &arg1);

    void on_upperLimitLineEdit_textChanged(const QString &arg1);

    void on_nValuesLineEdit_3_textChanged(const QString &arg1);

    void on_xValuesLineEdit_5_textChanged(const QString &arg1);

    void on_yValuesLimitLineEdit_2_textChanged(const QString &arg1);

    void on_nValuesLineEdit_2_textChanged(const QString &arg1);

    void on_xValuesLineEdit_6_textChanged(const QString &arg1);

    void on_nValuesLineEdit_4_textChanged(const QString &arg1);

    void on_tabWidget_2_tabBarClicked(int index);

    void on_trapezoidRadioButton_toggled(bool checked);

    void on_simpson_1RadioButton_toggled(bool checked);

    void on_simpson_2RadioButton_toggled(bool checked);

    void on_calculateButton_3_clicked();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////// \brief Euler Tab Signals Declarations ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    void on_eulerRadioButton_clicked(bool checked);

    void on_calculateButton_4_clicked();

    void on_diffEqLineEdit_editingFinished();

    void on_initialXLineEdit_editingFinished();

    void on_initialYLineEdit_editingFinished();

    void on_finalXLineEdit_editingFinished();

    void on_stepsLineEdit_editingFinished();

    void on_ModifiedEulerRadioButton_clicked(bool checked);

    void on_xValuesLineEdit_2_editingFinished();

    void on_yValuesLineEdit_2_editingFinished();

    void on_nValuesLineEdit_5_editingFinished();

    void on_LinearRadioButton_toggled(bool checked);

    void on_ExponentialRadioButton_toggled(bool checked);

    void on_LogarithmicRadioButton_toggled(bool checked);

    void on_PowerRadioButton_toggled(bool checked);

    void on_QuadraticRradioButton_toggled(bool checked);

    void on_calculateButton_5_clicked();

private:
    Ui::MainWindow *ui;

    ExpressionParser parser;

    bool functionValid = false;
    bool methodValid = false;

    QString functionExpression = ""; // Stores the entered function expression
    QString selectedMethod = "";     // Stores the selected root-finding method (Bisection or Secant or newton)

    // Vectors for x and y values
    std::vector<double> xValuesOne;
    std::vector<double> yValuesOne;

    // Interpolation point
    double xiValue = 0.0;
    // To store selected interpolation method (either "lagrange" or "newton")
    QString interpolationMethod = "";

    bool xValuesValid = false;
    bool yValuesValid = false;
    bool xiValueValid = false;
    bool methodInterpolationValid = false;

    QString functionExpression_2 = ""; // Stores the entered function expression
    void handleFunctionInput();
    void handle2DTableInput();
    void handleXOnlyTableInput();

    double lowerLimit = 0.0;
    double upperLimit = 0.0;
    double intervals = 0.0;
    int currentTab = 0;

    bool integrationFunctionValid = false;
    bool integrationMethodValid = false;
    bool lowerLimitValid = false;
    bool upperLimitValid = false;
    bool intervalsValid = false;
    bool xValuesOneValid = false;
    bool yValuesOneValid = false;

    void updateIntegrationCalculateButton();

    // Vectors for x and y values
    std::vector<double> xValues;
    std::vector<double> yValues;

    // To store selected Euler method (either "Normal Euler" or "Modified Euler")
    QString EulerMethod = "";
    QString FunctionExpressionForEuler = "";
    double XnodeValue = 0.0;
    double XfinalValue = 0.0;
    double YnodeValue = 0.0;
    double StepValue = 0.0;

    std::vector<double> xValuesTwo;
    std::vector<double> yValuesTwo;
    QString FittingMethod = "";  // Like "Linear", "Exponential", etc.
    QString FittingResult = "";

    int nValues = 0;

};
#endif // MAINWINDOW_H
