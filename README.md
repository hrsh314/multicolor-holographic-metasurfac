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
<img width="1365" height="732" alt="Picture2" src="https://github.com/user-attachments/assets/b84a20dd-21da-4ffd-a751-ed86983533b3" />



### 2. Geometric Phase Shift Extraction
Using the optimal length and width dimensions identified in the previous step, systematically rotate the meta-atoms to extract the geometric phase response and achieve a full 2π phase shift.

<img width="528" height="258" alt="Picture4" src="https://github.com/user-attachments/assets/dbba53a4-96da-4e43-95a7-acb8bc1aa14f" />
<img width="492" height="242" alt="Picture3" src="https://github.com/user-attachments/assets/7fef245e-4416-4361-aee6-d9be269d4786" />


### 3. Single-Colour Validation
Create individual, single-colour metasurfaces utilizing the extracted phase data. Run independent simulations to project the red flower and the green leaves separately, verifying that each target image reconstructs successfully with its respective phase mask.
<img width="349" height="305" alt="Picture5" src="https://github.com/user-attachments/assets/6d2265bb-976c-456b-9118-26f3a124c193" />
<img width="345" height="307" alt="Picture6" src="https://github.com/user-attachments/assets/5f391b20-94be-4835-a942-3eaf63ecfa68" />


### 4. Multicolour Integration & Final Testing
Once the single-colour projections are validated, place both the red and green meta-atoms diagonally into a single unit cell to generate the final multicolour metasurface. Test the combined structure by injecting the colours individually: verify that injecting only red light projects the red flower, and injecting only green light projects the green leaves. Confirming these isolated projections guarantees that injecting both wavelengths simultaneously will successfully project the complete, multicolour holographic image.
<img width="430" height="439" alt="Picture7" src="https://github.com/user-attachments/assets/c746c115-2db4-4057-acf5-d27fc4bf3fe4" />

