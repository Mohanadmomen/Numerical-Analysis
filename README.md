@@ -0,0 +1,57 @@
# Numerical-Analysis

Welcome to the Numerical Analysis repository! This project covers a range of numerical methods and techniques, categorized into five main chapters. These methods are essential for solving various mathematical problems that arise in real-world applications.

## 📚 Chapters Overview

### Chapter 1: Root Finding Methods
This chapter includes methods for finding the roots of functions. The key techniques covered are:

- **Bisection Method**: A reliable method for finding the root of a function within a specified interval.
- **Newton-Raphson Method**: A fast convergence method for finding roots using derivatives.
- **Secant Method**: A method similar to Newton's, but does not require the calculation of derivatives.

### Chapter 2: Interpolation
In this chapter, we explore various interpolation methods that estimate values between known data points:
- **Lagrange Interpolation**: A polynomial interpolation method using Lagrange polynomials.
- **Newton Interpolation (Forward and Backward)**: Techniques to build interpolation polynomials based on Newton’s formula, both for forward and backward differences.

### Chapter 3: Numerical Integration
This chapter introduces techniques for approximating definite integrals:
- **Trapezoidal Rule**: A simple method for approximating the integral using trapezoids.
- **Simpson's 1/3 Rule**: A more accurate method for integration using quadratic approximations.
- **Simpson's 3/8 Rule**: Another integration rule that improves accuracy by using cubic approximations.

### Chapter 4: Numerical Solution of Ordinary Differential Equations (ODEs)
This chapter focuses on solving ODEs using:
- **Euler's Method**: A straightforward approach to solving initial value problems.
- **Modified Euler's Method**: An improved version of Euler’s method that enhances accuracy.

### Chapter 5: Curve Fitting Methods
In this chapter, we cover various methods for fitting a curve to a set of data points:
- **Linear Regression**: Fitting a straight line to the data.
- **Exponential Fit**: Fitting an exponential curve to the data.
- **Logarithmic Fit**: Fitting a logarithmic curve to the data.
- **Power Fit**: Fitting a power-law curve to the data.
- **Quadratic Fit**: Fitting a parabolic curve to the data.

## 🖥️ GUI with Qt

This project also features a **Qt-based GUI** that allows users to interact with the numerical methods more easily. The GUI provides the following features:

- **Input Fields**: Enter the necessary parameters for the selected method, such as intervals for root-finding methods, data points for interpolation, and limits for integration methods.
- **Result Display**: View the results of the numerical methods, including graphical representations of the methods (e.g., root-finding results, interpolation curves, integration areas).
- **Plotting**: Visualization of interpolation curves, integration regions, and curve fitting models.
- **Interactive Controls**: Modify parameters dynamically and observe how the results change in real-time.

## 🛠️ Technologies Used

- **C++**: The primary programming language for implementing numerical methods and the GUI.
- **Qt5/Qt6**: Framework for building the GUI.

## 📚 Learning Objectives
The project aims to:

- Understand and implement fundamental numerical methods for solving mathematical problems.

- Develop and use a GUI to interact with numerical methods and visualize results.