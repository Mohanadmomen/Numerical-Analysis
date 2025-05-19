#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "RootFinding.h"
#include "Interpolation.h"
#include "Integration.h"
#include "Euler.h"
#include "CurveFitting.h"
#include "qcustomplot.h"
#include "ExpressionValidation.h"
#include <QMessageBox>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
    , parser()
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// \brief Root Finding Tab Signals Implementation ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

void MainWindow::on_functionLineEdit_textChanged(const QString &arg1)
{
    // Validate expression
    bool valid = validateExpression(arg1.toStdString());
    functionValid = valid;

    if (valid) {
        functionExpression = arg1;
        ui->functionLineEdit->setStyleSheet("");
    } else {
        ui->functionLineEdit->setStyleSheet("QLineEdit { border: 2px solid red; }");
    }

    // Only enable once both the text is valid AND a method is selected
    ui->calculateButton->setEnabled(functionValid && methodValid);
}


void MainWindow::on_bisectionRadioButton_toggled(bool checked)
{
    if (checked) {
        // Set the selected method to Bisection
        selectedMethod = "bisection";
        ui->methodResultLabel->setText("Bisection Method");
        ui->rootResultLabel->setText("--");
    }
    methodValid = checked;
    ui->calculateButton->setEnabled(functionValid && methodValid);
}


void MainWindow::on_secantRadioButton_toggled(bool checked)
{
    if (checked) {
        // Set the selected method to Secant
        selectedMethod = "secant";
        ui->methodResultLabel->setText("Secant Method");
        ui->rootResultLabel->setText("--");
    }
    methodValid = checked;
    ui->calculateButton->setEnabled(functionValid && methodValid);
}


void MainWindow::on_newtonRadioButton_2_toggled(bool checked)
{
    if (checked) {
        // Set the selected method to Newton
        selectedMethod = "newton";
        ui->methodResultLabel->setText("Newton Method");
        ui->rootResultLabel->setText("--");
    }
    methodValid = checked;
    ui->calculateButton->setEnabled(functionValid && methodValid);
}


void MainWindow::on_calculateButton_clicked()
{
    try {
        // Ensure the function expression is valid
        if (functionExpression.isEmpty()) {
            QMessageBox::warning(this, "Input Error", "Please enter a valid function expression.");
            return;
        }

        double result = 0.0;

        // Based on the selected method, call the appropriate function
        if (selectedMethod == "bisection") {
            result = bisectionMethod(parser, functionExpression.toStdString());
            // Display the result
            ui->rootResultLabel->setText(QString::number(result));
        }
        else if (selectedMethod == "secant") {
            result = secantMethod(parser, functionExpression.toStdString());
            // Display the result
            ui->rootResultLabel->setText(QString::number(result));
        }
        else if (selectedMethod == "newton") {
            result = newton_raphson(parser, functionExpression.toStdString());
            // Display the result
            ui->rootResultLabel->setText(QString::number(result));
        }
        else{
            QMessageBox::warning(this, "Input Error", "Please select a valid method.");
        }

    } catch (const std::exception& e) {
        QMessageBox::critical(this, "Error", QString::fromStdString(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// \brief Interpolation Tab Signals Implementation ///////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

void MainWindow::on_xValuesLineEdit_textChanged(const QString &arg1)
{
    // try to parse a comma‑list of doubles
    QStringList parts = arg1.split(",", Qt::SkipEmptyParts);
    std::vector<double> temp;
    bool ok = !parts.isEmpty();
    for (const QString &p : parts) {
        bool conv = false;
        double v = p.trimmed().toDouble(&conv);
        if (!conv) { ok = false; break; }
        temp.push_back(v);
    }

    if (ok) {
        xValuesOne = std::move(temp);
        ui->xValuesLineEdit->setStyleSheet("");
        xValuesValid = true;
    } else {
        ui->xValuesLineEdit->setStyleSheet("QLineEdit { border: 2px solid red; }");
        xValuesValid = false;
    }

    ui->calculateButton_2->setEnabled(xValuesValid && yValuesValid && xiValueValid && methodInterpolationValid);
}


void MainWindow::on_yValuesLineEdit_textChanged(const QString &arg1)
{
    QStringList parts = arg1.split(",", Qt::SkipEmptyParts);
    std::vector<double> temp;
    bool ok = !parts.isEmpty();
    for (const QString &p : parts) {
        bool conv = false;
        double v = p.trimmed().toDouble(&conv);
        if (!conv) { ok = false; break; }
        temp.push_back(v);
    }

    if (ok) {
        yValuesOne = std::move(temp);
        ui->yValuesLineEdit->setStyleSheet("");
        yValuesValid = true;
    } else {
        ui->yValuesLineEdit->setStyleSheet("QLineEdit { border: 2px solid red; }");
        yValuesValid = false;
    }

    ui->calculateButton_2->setEnabled(xValuesValid && yValuesValid && xiValueValid && methodInterpolationValid);
}


void MainWindow::on_interpolateLineEdit_textChanged(const QString &arg1)
{
    bool conv = false;
    double v = arg1.trimmed().toDouble(&conv);

    if (conv) {
        xiValue = v;
        ui->interpolateLineEdit->setStyleSheet("");
        xiValueValid = true;
    } else {
        ui->interpolateLineEdit->setStyleSheet("QLineEdit { border: 2px solid red; }");
        xiValueValid = false;
    }

    ui->calculateButton_2->setEnabled(xValuesValid && yValuesValid && xiValueValid && methodInterpolationValid);
}


void MainWindow::on_lagrangeRadioButton_toggled(bool checked)
{
    // Set the interpolation method to lagrange
    if (checked) {
        interpolationMethod = "lagrange";
        ui->methodResultLabel_2->setText("Lagrange Method");
        ui->valueResultLabel->setText("--");
    }
    methodInterpolationValid = checked;
    ui->calculateButton_2->setEnabled(xValuesValid && yValuesValid && xiValueValid && methodInterpolationValid);
}


void MainWindow::on_newtonRadioButton_toggled(bool checked)
{
    // Set the interpolation method to newton
    if (checked) {
        interpolationMethod = "newton";
        ui->methodResultLabel_2->setText("Newton Method");
        ui->valueResultLabel->setText("--");
    }
    methodInterpolationValid = checked;
    ui->calculateButton_2->setEnabled(xValuesValid && yValuesValid && xiValueValid && methodInterpolationValid);
}


void MainWindow::on_calculateButton_2_clicked()
{
    try {
        // Add validation check on x and y values
        if (!validateInterpolationInput(xValuesOne, yValuesOne)) {
            QMessageBox::warning(this, "Input Error", "Please enter valid x, y values");
            return;
        }

        double result = 0.0;

        // Call the appropriate interpolation function based on the selected method
        if (interpolationMethod == "lagrange") {
            result = lagrangeInterpolation(xValuesOne, yValuesOne, xiValue);
            // Display the result
            ui->valueResultLabel->setText(QString::number(result));
        } else if (interpolationMethod == "newton") {
            result = NewtonInterpolation(xValuesOne, yValuesOne, xiValue);
            // Display the result
            ui->valueResultLabel->setText(QString::number(result));
        } else{
            QMessageBox::warning(this, "Input Error", "Please select a valid method.");
        }

    } catch (const std::exception& e) {
        QMessageBox::critical(this, "Error", QString::fromStdString(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// \brief Integration Tab Signals Implementation ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

// --- Function Expression Input ---
void MainWindow::on_functionLineEdit_2_textChanged(const QString &arg1)
{
    integrationFunctionValid = validateExpression(arg1.toStdString());
    if (integrationFunctionValid) {
        functionExpression_2 = arg1;
        ui->functionLineEdit_2->setStyleSheet("");
    } else {
        ui->functionLineEdit_2->setStyleSheet("QLineEdit { border: 2px solid red; }");
    }
    updateIntegrationCalculateButton();
}

// --- Lower Limit Input ---
void MainWindow::on_lowerLimitLineEdit_textChanged(const QString &arg1)
{
    bool ok;
    double val = arg1.trimmed().toDouble(&ok);
    lowerLimitValid = ok;
    if (ok) {
        lowerLimit = val;
        ui->lowerLimitLineEdit->setStyleSheet("");
    } else {
        ui->lowerLimitLineEdit->setStyleSheet("QLineEdit { border: 2px solid red; }");
    }
    updateIntegrationCalculateButton();
}

// --- Upper Limit Input ---
void MainWindow::on_upperLimitLineEdit_textChanged(const QString &arg1)
{
    bool ok;
    double val = arg1.trimmed().toDouble(&ok);
    upperLimitValid = ok;
    if (ok) {
        upperLimit = val;
        ui->upperLimitLineEdit->setStyleSheet("");
    } else {
        ui->upperLimitLineEdit->setStyleSheet("QLineEdit { border: 2px solid red; }");
    }
    updateIntegrationCalculateButton();
}

// --- Intervals Input for Function Tab ---
void MainWindow::on_nValuesLineEdit_3_textChanged(const QString &arg1)
{
    bool ok;
    int n = arg1.trimmed().toInt(&ok);
    intervalsValid = ok && (n > 0);
    if (intervalsValid) {
        intervals = n;
        ui->nValuesLineEdit_3->setStyleSheet("");
    } else {
        ui->nValuesLineEdit_3->setStyleSheet("QLineEdit { border: 2px solid red; }");
    }
    updateIntegrationCalculateButton();
}

// --- X Values Table (2D) ---
void MainWindow::on_xValuesLineEdit_5_textChanged(const QString &arg1)
{
    QStringList parts = arg1.split(",", Qt::SkipEmptyParts);
    std::vector<double> temp;
    bool ok = !parts.isEmpty();
    for (const QString &p : parts) {
        bool conv;
        double v = p.trimmed().toDouble(&conv);
        if (!conv) { ok = false; break; }
        temp.push_back(v);
    }
    xValuesOneValid = ok;
    if (ok) {
        xValues = temp;
        ui->xValuesLineEdit_5->setStyleSheet("");
    } else {
        ui->xValuesLineEdit_5->setStyleSheet("QLineEdit { border: 2px solid red; }");
    }
    updateIntegrationCalculateButton();
}

// --- Y Values Table (2D) ---
void MainWindow::on_yValuesLimitLineEdit_2_textChanged(const QString &arg1)
{
    QStringList parts = arg1.split(",", Qt::SkipEmptyParts);
    std::vector<double> temp;
    bool ok = !parts.isEmpty();
    for (const QString &p : parts) {
        bool conv;
        double v = p.trimmed().toDouble(&conv);
        if (!conv) { ok = false; break; }
        temp.push_back(v);
    }
    yValuesOneValid = ok;
    if (ok) {
        yValues = temp;
        ui->yValuesLimitLineEdit_2->setStyleSheet("");
    } else {
        ui->yValuesLimitLineEdit_2->setStyleSheet("QLineEdit { border: 2px solid red; }");
    }
    updateIntegrationCalculateButton();
}

// --- Intervals Input for 2D Table Tab ---
void MainWindow::on_nValuesLineEdit_2_textChanged(const QString &arg1)
{
    bool ok;
    int n = arg1.trimmed().toInt(&ok);
    intervalsValid = ok && (n > 0);
    if (intervalsValid) {
        intervals = n;
        ui->nValuesLineEdit_2->setStyleSheet("");
    } else {
        ui->nValuesLineEdit_2->setStyleSheet("QLineEdit { border: 2px solid red; }");
    }
    updateIntegrationCalculateButton();
}

// --- X Values Table (X‑only) ---
void MainWindow::on_xValuesLineEdit_6_textChanged(const QString &arg1)
{
    QStringList parts = arg1.split(",", Qt::SkipEmptyParts);
    std::vector<double> temp;
    bool ok = !parts.isEmpty();
    for (const QString &p : parts) {
        bool conv;
        double v = p.trimmed().toDouble(&conv);
        if (!conv) { ok = false; break; }
        temp.push_back(v);
    }
    xValuesOneValid = ok;
    if (ok) {
        xValues = temp;
        ui->xValuesLineEdit_6->setStyleSheet("");
    } else {
        ui->xValuesLineEdit_6->setStyleSheet("QLineEdit { border: 2px solid red; }");
    }
    updateIntegrationCalculateButton();
}

// --- Intervals Input for X‑only Table Tab ---
void MainWindow::on_nValuesLineEdit_4_textChanged(const QString &arg1)
{
    bool ok;
    int n = arg1.trimmed().toInt(&ok);
    intervalsValid = ok && (n > 0);
    if (intervalsValid) {
        intervals = n;
        ui->nValuesLineEdit_4->setStyleSheet("");
    } else {
        ui->nValuesLineEdit_4->setStyleSheet("QLineEdit { border: 2px solid red; }");
    }
    updateIntegrationCalculateButton();
}

// --- Method Selection ---
void MainWindow::on_trapezoidRadioButton_toggled(bool checked)
{
    integrationMethodValid = checked;
    if (checked) {
        selectedMethod = "trapezoid";
        ui->methodResultLabel_3->setText("Trapezoid Rule");
        ui->integralResultLabel->setText("--");
    }
    updateIntegrationCalculateButton();
}

void MainWindow::on_simpson_1RadioButton_toggled(bool checked)
{
    integrationMethodValid = checked;
    if (checked) {
        selectedMethod = "simpson1";
        ui->methodResultLabel_3->setText("Simpson(1/3) Rule");
        ui->integralResultLabel->setText("--");
    }
    updateIntegrationCalculateButton();
}

void MainWindow::on_simpson_2RadioButton_toggled(bool checked)
{
    integrationMethodValid = checked;
    if (checked) {
        selectedMethod = "simpson2";
        ui->methodResultLabel_3->setText("Simpson(3/8) Rule");
        ui->integralResultLabel->setText("--");
    }
    updateIntegrationCalculateButton();
}

// --- Tab Change ---
void MainWindow::on_tabWidget_2_tabBarClicked(int index)
{
    currentTab = index;
    updateIntegrationCalculateButton();
}

// --- Calculate Button ---
void MainWindow::on_calculateButton_3_clicked()
{
    switch (currentTab) {
    case 0: handleFunctionInput();     break;
    case 1: handle2DTableInput();      break;
    case 2: handleXOnlyTableInput();   break;
    }
}

// --- Update Calculate Button State ---
void MainWindow::updateIntegrationCalculateButton()
{
    bool enable = false;
    switch (currentTab) {
    case 0:
        enable = integrationFunctionValid && lowerLimitValid && upperLimitValid
                 && intervalsValid && integrationMethodValid
                 && (lowerLimit < upperLimit);
        break;
    case 1:
        enable = xValuesOneValid && yValuesOneValid && intervalsValid && integrationMethodValid;
        break;
    case 2:
        enable = xValuesOneValid && intervalsValid && integrationFunctionValid && integrationMethodValid;
        break;
    }
    ui->calculateButton_3->setEnabled(enable);
}

////////////////////////////////////////////////////////////////////////////////////////////

void MainWindow::handleFunctionInput()
{
    if (functionExpression_2.isEmpty()) {
        QMessageBox::warning(this, "Input Error", "Please enter a valid function expression.");
        return;
    }
    if (lowerLimit >= upperLimit) {
        QMessageBox::warning(this, "Input Error", "Lower limit must be less than upper limit.");
        return;
    }
    if (selectedMethod.isEmpty()) {
        QMessageBox::warning(this, "Input Error", "Please select an integration method.");
        return;
    }

    try {
        double result = 0.0;
        if (selectedMethod == "trapezoid") {
            result = trapezoidalIntegration(parser, functionExpression_2.toStdString(), lowerLimit, upperLimit, intervals);
            ui->integralResultLabel->setText(QString::number(result));
        }
        else if (selectedMethod == "simpson1") {
            result = simpsonOneThird(parser, functionExpression_2.toStdString(), lowerLimit, upperLimit, intervals);
            ui->integralResultLabel->setText(QString::number(result));
        }
        else if (selectedMethod == "simpson2") {
            result = simpsonThreeEighth(parser, functionExpression_2.toStdString(), lowerLimit, upperLimit, intervals);
            ui->integralResultLabel->setText(QString::number(result));
        }
    }
    catch (const std::runtime_error& e) {
        QMessageBox::critical(this, "Calculation Error", QString::fromStdString(e.what()));
    }
}


void MainWindow::handle2DTableInput()
{
    //QVector<double> xValues = getXValuesFromTable();
    //QVector<double> yValues = getYValuesFromTable();

    // Ensure that x and y values are provided and are not empty
    if (xValues.size() != yValues.size() || xValues.empty()) {
        QMessageBox::warning(this, "Input Error", "X and Y values are invalid or empty.");
        return;
    }
    if (selectedMethod.isEmpty()) {
        QMessageBox::warning(this, "Input Error", "Please select an integration method.");
        return;
    }

    try {
        double result = 0.0;
        if (selectedMethod == "trapezoid") {
            result = trapezoidalIntegration2D(parser, functionExpression_2.toStdString(), lowerLimit, upperLimit, intervals,  xValues, yValues);
            ui->integralResultLabel->setText(QString::number(result));
        }
        else if (selectedMethod == "simpson1") {
            result = simpsonOneThird2D(parser, functionExpression_2.toStdString(), lowerLimit, upperLimit, intervals,  xValues, yValues);
            ui->integralResultLabel->setText(QString::number(result));
        }
        else if (selectedMethod == "simpson2") {
            result = simpsonThreeEighth2D(parser, functionExpression_2.toStdString(), lowerLimit, upperLimit, intervals,  xValues, yValues);
            ui->integralResultLabel->setText(QString::number(result));
        }
    }
    catch (const std::runtime_error& e) {
        QMessageBox::critical(this, "Calculation Error", QString::fromStdString(e.what()));
    }
}


void MainWindow::handleXOnlyTableInput()
{
    if (xValues.empty()) {
        QMessageBox::warning(this, "Input Error", "X values are invalid or empty.");
        return;
    }
    if (functionExpression.isEmpty()) {
        QMessageBox::warning(this, "Input Error", "Please enter a valid function to evaluate.");
        return;
    }
    if (selectedMethod.isEmpty()) {
        QMessageBox::warning(this, "Input Error", "Please select an integration method.");
        return;
    }

    try {
        double result = 0.0;
        if (selectedMethod == "trapezoid") {
            result = trapezoidalIntegrationTable(parser, functionExpression_2.toStdString(), lowerLimit, upperLimit, intervals,  xValues);
            ui->integralResultLabel->setText(QString::number(result));
        }
        else if (selectedMethod == "simpson1") {
            result = simpsonOneThirdTable(parser, functionExpression_2.toStdString(), lowerLimit, upperLimit, intervals,  xValues);
            ui->integralResultLabel->setText(QString::number(result));
        }
        else if (selectedMethod == "simpson2") {
            result = simpsonThreeEighthTable(parser, functionExpression_2.toStdString(), lowerLimit, upperLimit, intervals,  xValues);
            ui->integralResultLabel->setText(QString::number(result));
        }
    }
    catch (const std::runtime_error& e) {
        QMessageBox::critical(this, "Calculation Error", QString::fromStdString(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// \brief Euler Tab Signals Implementation ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////


void MainWindow::on_eulerRadioButton_clicked(bool checked)
{
    if (checked) {
        // Set the selected method to Euler
        EulerMethod = "Euler";
        ui->methodResultLabel_4->setText("Euler Method");
        ui->finalYResultLabel->setText("--");

    }
}


void MainWindow::on_ModifiedEulerRadioButton_clicked(bool checked)
{
    if (checked) {
        // Set the selected method to Modified Euler
        EulerMethod = "Modified Euler";
        ui->methodResultLabel_4->setText("Modified Euler Method");
        ui->finalYResultLabel->setText("--");

    }
}


void MainWindow::on_diffEqLineEdit_editingFinished()
{
    FunctionExpressionForEuler = ui->diffEqLineEdit->text();  // Get input expression

    // Add validation check
    if (!validateExpression(FunctionExpressionForEuler.toStdString())) {
        QMessageBox::warning(this, "Error", "Invalid expression");
        FunctionExpressionForEuler.clear();
    }
}


void MainWindow::on_initialXLineEdit_editingFinished()
{
    XnodeValue = ui->initialXLineEdit->text().trimmed().toDouble();
}


void MainWindow::on_initialYLineEdit_editingFinished()
{
    YnodeValue = ui->initialYLineEdit->text().trimmed().toDouble();
}


void MainWindow::on_finalXLineEdit_editingFinished()
{
    XfinalValue = ui->finalXLineEdit->text().trimmed().toDouble();
}


void MainWindow::on_stepsLineEdit_editingFinished()
{
    StepValue = ui->stepsLineEdit->text().trimmed().toDouble();
}


void MainWindow::on_calculateButton_4_clicked()
{
    // Ensure the function expression is valid
    if (FunctionExpressionForEuler.isEmpty()) {
        QMessageBox::warning(this, "Input Error", "Please fill all the information needed");
        return;
    }

    double result = 0.0;

    // Based on the selected method, call the appropriate function
    if (EulerMethod == "Euler") {
        result =  euler(parser, FunctionExpressionForEuler.toStdString(), XnodeValue, YnodeValue, StepValue, XfinalValue);
        // Display the result
        ui->finalYResultLabel->setText(QString::number(result));
    }
    else if (EulerMethod == "Modified Euler") {
        result = modifiedEuler(parser, FunctionExpressionForEuler.toStdString(), XnodeValue, YnodeValue, StepValue, XfinalValue);
        // Display the result
        ui->finalYResultLabel->setText(QString::number(result));
    }else{
        QMessageBox::warning(this, "Input Error", "Please Select a method");
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// \brief CurveFitting Tab Signals Implementation ///////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////


void MainWindow::on_xValuesLineEdit_2_editingFinished()
{
    QString xText = ui->xValuesLineEdit_2->text();
    xValuesTwo.clear();
    QStringList xStringList = xText.split(",", Qt::SkipEmptyParts);
    for (const QString& x : xStringList) {
        xValuesTwo.push_back(x.trimmed().toDouble());
    }
}


void MainWindow::on_yValuesLineEdit_2_editingFinished()
{
    QString yText = ui->yValuesLineEdit_2->text();
    yValuesTwo.clear();
    QStringList yStringList = yText.split(",", Qt::SkipEmptyParts);
    for (const QString& y : yStringList) {
        yValuesTwo.push_back(y.trimmed().toDouble());
    }
}


void MainWindow::on_nValuesLineEdit_5_editingFinished()
{
    nValues = ui->nValuesLineEdit_5->text().toInt();
}


void MainWindow::on_LinearRadioButton_toggled(bool checked)
{
    if (checked) {
        // Set the selected method to Linear
        FittingMethod = "Linear";
        ui->methodResultLabel_8->setText("Linear");
        ui->finalYResultLabel_5->setText("--");
    }
}


void MainWindow::on_ExponentialRadioButton_toggled(bool checked)
{
    // Set the selected method to Exponential
    FittingMethod = "Exponential";
    ui->methodResultLabel_8->setText("Exponential");
    ui->finalYResultLabel_5->setText("--");
}


void MainWindow::on_LogarithmicRadioButton_toggled(bool checked)
{
    // Set the selected method to Logarithmic
    FittingMethod = "Logarithmic";
    ui->methodResultLabel_8->setText("Logarithmic");
    ui->finalYResultLabel_5->setText("--");
}


void MainWindow::on_PowerRadioButton_toggled(bool checked)
{
    // Set the selected method to Power
    FittingMethod = "Power";
    ui->methodResultLabel_8->setText("Power");
    ui->finalYResultLabel_5->setText("--");
}


void MainWindow::on_QuadraticRradioButton_toggled(bool checked)
{
    // Set the selected method to Quadratic
    FittingMethod = "Quadratic";
    ui->methodResultLabel_8->setText("Quadratic");
    ui->finalYResultLabel_5->setText("--");
}


void MainWindow::on_calculateButton_5_clicked()
{
    try {
        if (xValuesTwo.size() != yValuesTwo.size() || xValuesTwo.empty()) {
            QMessageBox::warning(this, "Input Error", "X and Y must have the same non-zero size.");
            return;
        }

        // 1. Clear previous plots
        ui->plotWidget->clearGraphs();

        // 2. Plot original scatter points
        plotScatter(xValuesTwo, yValuesTwo);

        CurveFitter fitter;

        // 3. Based on selected method, fit and plot curve
        if (FittingMethod == "Linear") {
            //FittingResult = fitter.fit(xValuesTwo, yValuesTwo, CurveFitter::LINEAR);
            auto [a, b] = fitter.getLinearCoefficients(xValuesTwo, yValuesTwo);
            // Format equation string: y = a·e^(bx)
            FittingResult = QString("y = %1·x + %2").arg(a, 0, 'g', 4).arg(b, 0, 'g', 4);

            // Display the result
            ui->finalYResultLabel_5->setText(FittingResult);
            plotCurve(a, b, "Linear");
        }
        else if (FittingMethod == "Exponential") {
            //FittingResult = fitter.fit(xValuesTwo, yValuesTwo, CurveFitter::EXPONENTIAL);
            auto [a, b] = fitter.getExponentialCoefficients(xValuesTwo, yValuesTwo);
            // Format equation string: y = ax + b
            FittingResult = QString("y = %1·e^(%2x)").arg(a, 0, 'g', 4).arg(b, 0, 'g', 4);
            // Display the result
            ui->finalYResultLabel_5->setText(FittingResult);
            plotCurve(a, b, "Exponential");
        }
        else if (FittingMethod == "Logarithmic") {
            //FittingResult = fitter.fit(xValuesTwo, yValuesTwo, CurveFitter::LOGARITHMIC);
            auto [a, b] = fitter.getLogarithmicCoefficients(xValuesTwo, yValuesTwo);
            // Format equation string: y = a + b·ln(x)
            FittingResult = QString("y = %1 + %2·ln(x)").arg(a, 0, 'g', 4).arg(b, 0, 'g', 4);
            // Display the result
            ui->finalYResultLabel_5->setText(FittingResult);
            plotCurve(a, b, "Logarithmic");
        }
        else if (FittingMethod == "Power") {
            //FittingResult = fitter.fit(xValuesTwo, yValuesTwo, CurveFitter::POWER);
            auto [a, b] = fitter.getPowerCoefficients(xValuesTwo, yValuesTwo);
            // Format equation string: y = a·x^b
            FittingResult = QString("y = %1·x^%2").arg(a, 0, 'g', 4).arg(b, 0, 'g', 4);
            // Display the result
            ui->finalYResultLabel_5->setText(FittingResult);
            plotCurve(a, b, "Power");
        }
        else if (FittingMethod == "Quadratic") {
            //FittingResult = fitter.fit(xValuesTwo, yValuesTwo, CurveFitter::QUADRATIC);
            auto [a, b, c] = fitter.getQuadraticCoefficients(xValuesTwo, yValuesTwo);
            // Format equation string: y = ax² + bx + c
            FittingResult = QString("y = %1x² + %2x + %3").arg(a, 0, 'g', 4).arg(b, 0, 'g', 4).arg(c, 0, 'g', 4);
            // Display the result
            ui->finalYResultLabel_5->setText(FittingResult);
            // Special case: Quadratic needs 3 parameters
            // So we need a special plotCurve for Quadratic
            plotQuadraticCurve(a, b, c, "Quadratic");
        }

        // Redraw
        ui->plotWidget->replot();
    }
    catch (const std::runtime_error& e) {
        QMessageBox::critical(this, "Calculation Error", QString::fromStdString(e.what()));
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MainWindow::plotScatter(const std::vector<double>& x, const std::vector<double>& y) {
    // Manually convert std::vector to QVector
    QVector<double> xQVec(x.begin(), x.end());
    QVector<double> yQVec(y.begin(), y.end());

    // Clear any existing graphs on the plot widget
    ui->plotWidget->clearGraphs();

    // Add a new graph to the plot
    ui->plotWidget->addGraph();
    ui->plotWidget->graph(0)->setData(xQVec, yQVec);  // Use QVector here

    // Set appearance of the scatter plot
    ui->plotWidget->graph(0)->setPen(QPen(Qt::black));  // Set color of points
    ui->plotWidget->graph(0)->setLineStyle(QCPGraph::lsNone);  // No lines between points
    ui->plotWidget->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 5));  // Circle marker, size 5

    // Rescale axes to fit the data
    ui->plotWidget->rescaleAxes();

    // Replot to update the graph
    ui->plotWidget->replot();
}


void MainWindow::plotCurve(double a, double b, const QString& methodName) {
    // Number of points to plot the curve
    const int points = 100;

    QVector<double> x(points), y(points);

    // Define your x range (for example from minX to maxX)
    double minX = -10; // Change as needed
    double maxX =  10;

    double step = (maxX - minX) / (points - 1);

    for (int i = 0; i < points; ++i) {
        x[i] = minX + i * step;
        if (methodName == "Linear") {
            y[i] = a * x[i] + b;
        }
        else if (methodName == "Exponential") {
            y[i] = a * std::exp(b * x[i]);
        }
        else if (methodName == "Logarithmic") {
            if (x[i] > 0) // log is undefined for x <= 0
                y[i] = a + b * std::log(x[i]);
            else
                y[i] = 0; // or you can skip these points if needed
        }
        else if (methodName == "Power") {
            if (x[i] > 0) // Power model assumes positive x
                y[i] = a * std::pow(x[i], b);
            else
                y[i] = 0;
        }
        else {
            // Default: straight line
            y[i] = a * x[i] + b;
        }
    }

    // Now plot
    ui->plotWidget->addGraph();
    ui->plotWidget->graph()->setData(x, y);
    ui->plotWidget->graph()->setPen(QPen(Qt::red)); // You can change the color
    ui->plotWidget->graph()->setName(methodName); // Name it
    ui->plotWidget->rescaleAxes();
    ui->plotWidget->replot();
}

void MainWindow::plotQuadraticCurve(double a, double b, double c, const QString& methodName) {
    const int points = 100;
    QVector<double> x(points), y(points);

    double minX = -10; // or based on your xValuesTwo range
    double maxX = 10;
    double step = (maxX - minX) / (points - 1);

    for (int i = 0; i < points; ++i) {
        x[i] = minX + i * step;
        y[i] = a * x[i] * x[i] + b * x[i] + c;
    }

    ui->plotWidget->addGraph();
    ui->plotWidget->graph()->setData(x, y);
    ui->plotWidget->graph()->setPen(QPen(Qt::blue)); // different color
    ui->plotWidget->graph()->setName(methodName);
    ui->plotWidget->rescaleAxes();
}
