In this repository you can find a Molecular Dynamics simulation for liquid Argon and Monte Carlo simulation for the same species. A brief overview of both methods is given in the paragraphs that follow.

Note that the codes and files required for the MC simulation are all inside the folder named "MC".

# MD-for-LJ-liquid-Argon

This program simulates liquid Argon using a Lennard Jones (LJ) potential and allows to calculate the radial distribution function (RDF), the pressure and the energy.

Moreover it incorporates a refinement of the employed Boundary Conditions which allow to calculate the velocity autocorrelation function (VACF) and the mean square displancement (MSD), providing insights into the system dynamics.

## Example of the Energy conversion over the iterations

The following plot shows the convergence of the different energy contributions as a function of the simulation length, as an example of some of the results that can be obtained by performing such a MD simulation.
![Energy convergence for LJ fluid (Argon) as a function of the simulation length](figures/energy_convergence.jpg)


## Comments on use

To run the code make sure you download the the two files "leap-conf.data" & "leap-lj.data" and save them on the same folder as the ".f" file. The leap-lj.data file contains physical values for the simulation of liquid Argon and technical parameters that are useful to run the code, namely Nconf which can be adjusted to lower values when modyfying the code to see whether the new implementations are correct. (Don't forget to set it back to a big enough value to get illustrative results).

Moreover, the file leap-conf.data contains the initial conditions for the positions and velocities of the n atoms.


# MC-for-LJ-liquid-Argon

This program simulates liquid Argon on the (N,V,T) ensemble using a Lennard Jones (LJ) potential, allowing to calculate the radial distribution function (RDF), the pressure and the energy of the system. 

Moreover, the code can be used to compute the RDF, the energy and the pressure of the system for different densities, leading to the equation of state of our LJ fluid, which relates pressure as a function of density.

## Example of the equation of state

The following plot shows the numerically-obtained equation of state for liquid Argon at T=2 (in reduced units), as an example of some of the results that can be obtained by performing such a MC simulation.
![Equation of state for LJ fluid (Argon) using MC methods](figures/energy_convergence.jpg)
