# Q.M-Nuclear-Physics-Phase-shifts-from-R-matrix-theory

## The stationary Schrödinger equation

If the Coulomb and Nuclear potential are considered, the radial stationary Schrödinger equation can be written as

![SchrodingerEq](src/SchrodingerEq.png)

In order to find the solutions of the equation above, it is useful to divide the configuration space into two regions: i) The internal region with r<a, where the nuclear and Coulomb potential interaction dominate; ii) The external region with r>a where the Nuclear potential vanishes. 

The solution for the external region can be written in terms of the Hankel functions 
![Hankel](src/Hankel.png)

which are constructed by Coulomb functions and the cosine and sine of a phase shift.

The phase shift can be obtained using the method of R-matrix.

Firstly C matrix is created as

![Cmatrix](src/Cmatrix.png)

where the kinetic elements are calculated using the Gauss-Legendre quadrature plus a Blonch operator L<sub>i,j</sub>

For i=j

![KineticMatrix1](src/KineticMatrix1.png)

For i≠j

![KineticMAtrix2](src/KineticMAtrix2.png)

The potential matrix elements take the simple form

![Potential](src/Potential.png)
