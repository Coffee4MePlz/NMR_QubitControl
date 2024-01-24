# NMR Qubit Control MATLAB Code

## Overview
This repository contains MATLAB scripts for controlling qubits using Nuclear Magnetic Resonance (NMR) techniques, focusing on modulated radiofrequency pulses (mostly for Fourier Series modulation and optimization of Fourier coefficients), Quantum Process Tomography, and noise analysis in quantum systems. The code is based on extensive research conducted in the field of Quantum Information Processing using NMR, as detailed in the accompanying Master's thesis. The codes in this repository where developed at UFABC (Federal University of The São Paulo ABC) at the [Quantum information Lab](https://www.quantumufabc.org/) and financed by CAPEs .

### Files of interest in the Repository

1. `TUTORIAL_ModulatedPulses.m`
2. `modulados_HomQbits_Par_withNoise.m`
3. `genChoi_map_exp.m`
4. `MasterThesis_TheoryForModulatedPulses_GustavoCafe.pdf`


## `TUTORIAL_ModulatedPulses.m`
This tutorial script introduces modulated pulses in qubit control and serves as a tutorial for the `modulados_HomQbits_Par_withNoise.m` file. 

### Key Concepts:
- Fundamentals of pulse modulation for qubit control.
- Practical MATLAB implementation in `modulados_HomQbits_Par_withNoise.m` and result analysis.

## `modulados_HomQbits_Par_withNoise.m`
This advanced script applies modulated pulses to homogenous qubits with added noise, it is the main file of interest, and is used by the tutorial file. In it there is an embedded optimization routine that finds modulated pulses by Fourier Series optimization. It demonstrates noise simulation and its impact on qubit control for various systems and molecules (chloroform, trifluor, etc.). Care should be taken for the use of diferent systems, i.e., heteronuclear and homonuclear. 

### Features:
- Noise simulation in quantum systems.
- Techniques for controlling heteronuclear and homonuclear qubits in noisy environments.

## `genChoi_map_exp.m`
This script is for generating the Choi matrix from experimental data, i.e. from experimental Quantum Process Tomography (QPT). It directly links to the concepts discussed in the master's thesis, particularly the challenges and methodologies in QPT and noise analysis in NMR quantum information protocols.

### Focus Areas:
- Generating Choi matrices from experimental quantum process tomographies.
- Understanding the implications of noise and imperfections in quantum systems.


## How to Use
1. Clone or download this repository to your local machine.
2. Open MATLAB and navigate to the downloaded scripts.
3. Start with `TUTORIAL_ModulatedPulses.m` for basic concepts.
4. Proceed to `modulados_HomQbits_Par_withNoise.m` and `genChoi_map_exp.m` for advanced topics and practical applications.

## Contributing
Contributions to enhance, optimize, or expand the capabilities of this repository are welcome. Please ensure that your contributions are well-documented and include tests or examples.

## License
This project is open-sourced, developed at UFABC (Federal University of The São Paulo ABC) at the [Quantum information Lab](https://www.quantumufabc.org/) and financed by CAPEs .
