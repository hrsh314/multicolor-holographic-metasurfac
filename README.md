# Multicolour Holographic Metasurfaces

This repository contains the design, phase library, and Lumerical FDTD simulation pipeline for a multicolour hologram utilizing geometric phase metasurfaces. The workflow applies the Gerchberg-Saxton (GS) algorithm for Computer-Generated Holography (CGH) to extract target phase masks for individual colours.  

## Device Architecture
The metasurface operates with Right-Circularly Polarized (RCP) incident light, evaluating the Left-Circularly Polarized (LCP) light on the transmitted side.  

* **Material Stack:** The device utilizes a straightforward two-layer architecture consisting of a Silicon Dioxide (SiO2) substrate with Polycrystalline Silicon (Poly Si) meta-atoms patterned on top.
* **Metasurface Construction:** The design scales from independent single-colour metasurfaces up to a fully integrated multicolour architecture.
* **Spatial Multiplexing:** To achieve multicolour projection, both the red-diffracting and green-diffracting meta-atoms are placed diagonally within a single unit cell and spatially multiplexed.
* **Unit Cell Periodicity:** The unit cells were rigorously evaluated at periodicities of both **P = 400 nm** and **P = 350 nm**.  

## Simulation Workflow & Usage

The simulation pipeline is executed in four sequential steps to ensure pure colour projection and eliminate crosstalk:

### 1. Meta-Atom Simulation & Optimization
Simulate the unit meta-atoms for both red and green colours to generate Phase Conversion Ratio (PCR) and Diffraction Efficiency (DE) heatmaps. Using these graphs, select specific length and width dimensions ensuring that the red meta-atom diffracts exclusively red light, and the green meta-atom diffracts exclusively green light.

### 2. Geometric Phase Shift Extraction
Using the optimal length and width dimensions identified in the previous step, systematically rotate the meta-atoms to extract the geometric phase response and achieve a full 2π phase shift.

### 3. Single-Colour Validation
Create individual, single-colour metasurfaces utilizing the extracted phase data. Run independent simulations to project the red flower and the green leaves separately, verifying that each target image reconstructs successfully with its respective phase mask.

### 4. Multicolour Integration & Final Testing
Once the single-colour projections are validated, place both the red and green meta-atoms diagonally into a single unit cell to generate the final multicolour metasurface. Test the combined structure by injecting the colours individually: verify that injecting only red light projects the red flower, and injecting only green light projects the green leaves. Confirming these isolated projections guarantees that injecting both wavelengths simultaneously will successfully project the complete, multicolour holographic image.
