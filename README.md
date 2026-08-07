<div align="center">

# Model Hamiltonian


<img width="843" height="423" alt="image" src="https://github.com/user-attachments/assets/45307ed2-9d8c-4884-93c6-b11ff923d885" />
</div>




# Summary
A Model Hamiltonian for Hydrogen Bonded Complexes comprising of donor, acceptor, and hydrogen-bond coordinate. This model is capable of capturing the coupling patterns resulting different energy flow dynamics from excited donor or acceptor modes. 
# Theory

The model describes vibrational energy transfer in hydrogen-bonded complexes using a reduced three-mode Hamiltonian comprising:

- **H-bond Donor (D)**
- **H-bond Acceptor (A)**
- **Hydrogen-Bond Coordinate (HB)**

The Hamiltonian is expressed as

**H = H<sub>D</sub> + H<sub>A</sub> + H<sub>HB</sub> + H<sub>coup</sub>**

where each mode is represented by a Morse oscillator, and the coupling Hamiltonian includes:

- Bilinear donor–acceptor interaction
- Cubic nonlinear donor–acceptor interactions
- Linear donor/acceptor–HB couplings
- Quadratic donor/acceptor–HB couplings

These interaction terms capture the anharmonic vibrational couplings responsible for:

- Intramolecular vibrational redistribution (IVR)
- Vibrational predissociation (VP)
- Energy transfer between donor, acceptor, and hydrogen-bond modes

The mode-specific energy currents are obtained from the coupling Hamiltonian as


<p align="center">
<b>
J<sub>i</sub> = dE<sub>i</sub>/dt = -q̇<sub>i</sub> ∂H<sub>coup</sub>/∂q<sub>i</sub>
</b>
</p>
which allows direct analysis of the energy flow during classical trajectory simulations.

## Dynamics

The equations of motion are propagated using classical molecular dynamics.

- **Trajectory length:** 3 ps
- **Time step:** 0.2 fs
- **Number of trajectories:** 200
- **Initial excitations:** OH donor or OH acceptor stretch

The coupling parameters are obtained from the molecular geometry through numerical differentiation of the potential energy surface.

> **Note:** A complete derivation of the Hamiltonian, energy-current expressions, and parameterization is available in the accompanying publication.
