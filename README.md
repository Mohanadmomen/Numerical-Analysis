# Numerical Analysis Tool

**Numerical Analysis Tool** is a cross‑platform Qt/C++ application that provides interactive, visual implementations of fundamental numerical methods. Whether you need to find roots of nonlinear equations, solve ODEs, interpolate data, approximate integrals, or fit curves, this tool has you covered—all in an easy‑to‑use tabbed interface.

## Features

- **Root Finding**  
  - Bisection Method  
  - Secant Method  
  - Newton–Raphson Method  

- **Ordinary Differential Equations (ODE) Solver**  
  - Euler’s Method  
  - Modified Euler’s Method  

- **Polynomial Interpolation**  
  - Lagrange Interpolation  
  - Newton’s Divided Difference  

- **Numerical Integration**  
  - Trapezoidal Rule  
  - Simpson’s 1/3 Rule  
  - Simpson’s 3/8 Rule  
  - Support for function input, 2D data tables, or X‑only tables  

- **Curve Fitting** (Least Squares)  
  - Linear  
  - Logarithmic  
  - Exponential  
  - Power  
  - Quadratic  
  - Cubic  

- **Validation & Ease‑of‑Use**  
  - Real‑time input validation with red‑border error feedback  
  - Dynamic enabling/disabling of Calculate buttons  
  - Clear result display and error messages  

## Getting Started

### Prerequisites

- Qt 5 or 6  
- C++17 compiler (GCC, Clang, MSVC)  
- CMake (or Qt qmake)  

### Building

```bash
git clone https://github.com/<your‑username>/numerical-analysis-tool.git
cd numerical-analysis-tool
mkdir build && cd build
cmake ..           # or qmake ../Numerical_Analysis_Tool.pro
make
