# N-body
We simulate the classical n-body problem, i.e. how do n massive object evolve in time, given some initial state, when the only acting force on them is gravity. This is a standard exercise for a Numerical Methods class and some interesting initial condtions can be found in `setup.pdf`. 

### Goals of this project:
 - get familiar with `boost/odeint`, professional ODE solver. By this I mean: get familiar with syntax, bsc. classes and functions.  For this, implement some simple ODE solvers (Explicit Euler, Explicit Midpoint, Velocity Verlet), and compare it to them.
 - test out Python's animation library. Seems interesting to visualize data that evolves in time. 
 - learn how to use HDF5 in C++.

### Currently working on:
 - `ode/int` solution of the problem. The solution for Pluto comes out completely wrong, while the rest works just fine.. Needs more investigation. 

### Things that still need to be done:
 - Tests! Currently there are no tests for the implementation... The results of the simulation seem fine, but we should probably make some basic analysis of the precision and how well it handels what kind of problems. A reasonable reference would probably be an analytical solution or otherwise `ode/int`, which seems to be well tested.
 - The implementation assumes in general that we are dealing with a 3D problem, i.e. that the phase space coordinates have the form (x,y,z,v_x,v_y,v_z). This is fine, but it needs to be documented to make sure that one can modify 2D problems to still be solvable with this implementation. One might solve this problem by reimplementing the naive solver using `point_type.hpp`.
