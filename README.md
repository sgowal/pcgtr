# Preconditioned Conjugate Gradient and Trust-Region Solver (PCG-TR)
### or Efficient Numerical Unconstrained Optimization for Embedded Systems

Especially indented for embedded systems, this C library solves any unconstrained minimization problem for which the Jacobian and Hessian of the cost function can be computed. The library weighs 24kB on a 32bits machine.

This library is based on a standard optimization strategy that combines Preconditioned Conjugate Gradient and Trust Region methods.
It is very efficient when the Hessian matrix is sparse (see examples are provided in the examples folder).

Its API is simple: it offers 4 functions:
* `int pcgtr_solver_init(pcgtr_solver_t *solver, int n, int sparse);`
* `void pcgtr_solver_halloc(pcgtr_solver_t *solver, int row, int col);`
* `void pcgtr_solver_solve(pcgtr_solver_t *solver, matrix_t *x);`
* `void pcgtr_solver_destroy(pcgtr_solver_t *solver);`
All functions are heavily documented in the code.

### Example : Rosenbrock function

Solving the Rosenbrock function, given an initial starting point uniformly sampled between 0 and 2 for all n dimensions, yields the correct final solution in an average time of (using a single core of an Intel i7 2.7GHz processor with 100 runs): 

 n   |   Using sparseness   |   Not using sparseness
-----+----------------------+----------------------------------
   2 |        8.61 us       |           6.32 us
   4 |       18.39 us       |          13.77 us
   8 |       40.33 us       |          31.38 us
  16 |       82.32 us       |          88.44 us
  32 |      138.69 us       |         236.20 us
  64 |      262.38 us       |         863.00 us
 128 |      522.55 us       |       3'343.59 us 
 256 |    1'035.23 us       |      12'935.90 us
 512 |    1'970.24 us       |      49'911.36 us
1024 |    4'025.77 us       |     219'008.00 us
2048 |    9'179.77 us       |   1'055'101.06 us
4096 |   19'243.42 us       |   4'708'783.17 us (10 runs only)

### Compilation

You can type `make` to compile the library, which will be located in the build folder.
You can type `make examples` to compile the examples.
You can execute the first example with `examples/example1 -h` to get help or to try the Rosenbrock example `examples/example2 2` with 2 variables.

### Enjoy!

### License

Licensed under the MIT license.

Copyright (c) 2013 Sven Gowal

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.