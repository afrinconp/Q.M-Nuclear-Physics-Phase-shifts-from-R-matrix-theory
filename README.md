# Q.M-Nuclear Physics: Collision Phase shifts of two nuclei using the R-matrix Theory

## The stationary Schrödinger equation

If the Coulomb and Nuclear potential are considered, the radial stationary Schrödinger equation can be written as

![SchrodingerEq](src/SchrodingerEq.png)

In order to find the solutions of the equation above, it is useful to divide the configuration space into two regions: i) The internal region with r<a, where the nuclear and Coulomb potential interaction dominate; ii) The external region with r>a where the Nuclear potential vanishes. 

The solution for the external region can be written in terms of the Hankel functions 
![Hankel](src/Hankel.png)

which are constructed by Coulomb functions and the cosine and sine of a phase shift.

Here:

![k](src/k.png)

![eta](src/eta.png)

Z<sub>1</sub> and Z<sub>2</sub> are the atomic number of projectile and target respectively

The phase shift can be obtained using the method of R-matrix.

## R-Matrix 

Firstly C matrix is created as

![Cmatrix](src/Cmatrix.png)

where the kinetic elements are calculated using the Gauss-Legendre quadrature plus a Blonch operator L<sub>i,j</sub>

For i=j

![KineticMatrix1](src/KineticMatrix1.png)

For i≠j

![KineticMAtrix2](src/KineticMAtrix2.png)

The potential matrix elements take the simple form

![Potential](src/Potential.png)

φ<sub>i</sub> are the Lagrange-Legendre functions evaluated in a (f<sub>i</sub>(a)), these fucntions are defined as

![Lagrange-Legendre-functions](src/Lagrange-Legendre-functions.png)

Here x<sub>i</sub> are the roots of Legendre polynomial of order N  (P<sub>N</sub>), which satisfies  

![Roots-Legendre](src/Roots-Legendre.png)

Then, the R-matrix is made following the equation below

![Rmatrix](src/Rmatrix.png)

## Scattering matrix

Once we know the R-matrix, it is possible to find the scaterring matrix appearing in the asymptotic scattering wave function.

![Scattering-matrix](src/Scattering-matrix.png)

Then, the phase shift is easy to find due to the fact that

![Scatter2](src/Scatter2.png)

## Python Program

I constructed a program in python that calculates the shift phase in function of the energy (MeV) (MatrixPhaseShiftsInDegrees.py), for the 12C and a proton. The potential which is used in this program is:

V<sub>N</sub>(r)=-73.8exp(-(r/2.70)<sup>2</sup>)

V<sub>C</sub>(r)=6e<sup>2</sup>/r

The results are showed below for l=0, N=7 and N=10

![12Cp](12Cp.png)

Besides, instead of having 12C and a proton, I  obtained the result for other  for two α particles and let be l=4

The potential, which must be used is:

V<sub>N</sub>(r)=-73.8exp(-(r/2.70)<sup>2</sup>)

V<sub>C</sub>(r)=4e<sup>2</sup>erf(r/1.33)/r.

![Alpha shift α-α](AlphaShiftα-α.png)

You can get more information if you click the paper below:

[The R-matrix theory](https://iopscience.iop.org/article/10.1088/0034-4885/73/3/036301/meta)



