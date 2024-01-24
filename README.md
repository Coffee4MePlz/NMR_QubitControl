# Matlab codes for NMR Quantum Computing
# NMR Qubit Control MATLAB Code

## Overview
This repository contains MATLAB scripts for controlling qubits using Nuclear Magnetic Resonance (NMR) techniques, focusing on modulated pulses, process tomography, and noise analysis in quantum systems. The code is based on extensive research conducted in the field of quantum information processing using NMR, as detailed in the accompanying master's thesis.

### Files in the Repository

1. `TUTORIAL_ModulatedPulses.m`
2. `modulados_HomQbits_Par_withNoise.m`
3. `genChoi_map_exp.m`

## `TUTORIAL_ModulatedPulses.m`
This tutorial script introduces modulated pulses in qubit control. It demonstrates the implementation and analysis of these techniques in MATLAB, providing a foundation for understanding pulse modulation in quantum systems.

### Key Concepts:
- Fundamentals of pulse modulation for qubit control.
- Practical MATLAB implementation and result analysis.

## `modulados_HomQbits_Par_withNoise.m`
This advanced script applies modulated pulses to homogenous qubits with added noise, building on the concepts from the tutorial. It demonstrates noise simulation and its impact on qubit control.

### Features:
- Noise simulation in quantum systems.
- Techniques for controlling homogenous qubits in noisy environments.

## `genChoi_map_exp.m`
This script is for generating the Choi matrix from experimental data, a crucial component in Quantum Process Tomography (QPT). It directly links to the concepts discussed in the master's thesis, particularly the challenges and methodologies in QPT and noise analysis in NMR quantum systems.

### Focus Areas:
- Generating Choi matrices from experimental quantum process tomographies.
- Understanding the implications of noise and imperfections in quantum systems.

## Prerequisites
- MATLAB (preferably recent versions for compatibility)
- Basic knowledge of NMR and quantum computing concepts
- Familiarity with MATLAB programming and matrix operations

## How to Use
1. Clone or download this repository to your local machine.
2. Open MATLAB and navigate to the downloaded scripts.
3. Start with `TUTORIAL_ModulatedPulses.m` for basic concepts.
4. Proceed to `modulados_HomQbits_Par_withNoise.m` and `genChoi_map_exp.m` for advanced topics and practical applications.

## Contributing
Contributions to enhance, optimize, or expand the capabilities of this repository are welcome. Please ensure that your contributions are well-documented and include tests or examples.

## License
This project is open-sourced, developed at UFABC (Federal University of The São Paulo ABC) at the [Quantum information Lab](https://www.quantumufabc.org/) and financed by CAPEs .
