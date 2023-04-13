# Nelder-Mead Simplex Method for Optimization

This is an iterative numerical method to find the (local) minimum of an objective function.

For an optimization problem with $n$ variables, a simplex of $n+1$ points is created which gradually contracts around the local optimal value. This method does not involve computing the derivatives/gradient of the objective function.

A "simplex" is basically a generalization of what a triangle is in 2-dimensional space. For example, a simplex in 3D would be a tetrahedron, and a simplex in 1D would be a line segment.

## Theory

For detailed information on the theory, check the [DOCUMENTATION](NM_Simplex.pdf) PDF.

The algorithm is very simple to understand. I have adapted my explanation from the original paper by JA Nelder and R Mead, which is available for free online:

````
Nelder, J. A., & Mead, R. (1965). A Simplex Method for Function Minimization. The Computer Journal, 7(4), 308-313.
````

## Usage

See [example.c](example.c) for a sample usage of the procedures. The example is the "Motivating Problem" in the [DOCUMENTATION](NM_Simplex.pdf) PDF.

You can run `./test.sh` (after giving execute permissions) to see the example in action.

