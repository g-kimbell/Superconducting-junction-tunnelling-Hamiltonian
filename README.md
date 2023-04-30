This project is a model for superconducting junctions based on tunnelling Hamiltonian Green functions.
The code has been written by Graham Kimbell with support from Xavier Montiel on the theory aspects.

This model considers structures with periodic boundary conditions in x-y and broken translational symmetry along z. In other words,
we consider 2-D layers of materials defined in momentum space in x-y stacked along z. The Hamiltonian at every k-point is defined in
a spin x nambu space basis for each layer, giving a 4x4 matrix for a single layer, and a 4n x 4n matrix for n layers. The Hamiltonian
can consist of the electronic disperson, s-wave or d-wave BCS coupling, applied field at any angle, internal exchange field at any 
angle, and Rashba spin-orbit coupling. An interlayer coupling allows tunnelling between layers.

The model is mostly of interest for ferromagnetic superconducting junctions, where we can show properties such as singlet-triplet 
conversion at the S/F interfaces, as well as more complex behaviour such as competing opposite sign singlet and triplet currents in 
ferromagnetic junctions and anomalous Josephson currents in junctions with finite spin chirality.

Different junctions can be quickly defined, and then the superconducting properties calculated. The superconducting gap can be solved
self-consistently. From this the critical temperatre and critical field (both first order and second order) can be found for junctions.

Additionally, the critical current can be calculated through the structure as a function of the phase difference of the superconductors,
and there are additional tools for quickly finding 0-pi transitions as a function of ferromagnet angles.

If you are interested in utilising this code for yourself I am happy to help - you can contact me at grahamkimbell [at] outlook.com
