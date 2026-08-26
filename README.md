# Multicolour Holographic Metasurfaces

This repository contains the design, phase library, and Lumerical FDTD simulation pipeline for a multicolour hologram utilizing geometric phase metasurfaces. The workflow applies the Gerchberg-Saxton (GS) algorithm for Computer-Generated Holography (CGH) to extract target phase masks for individual colours.  

## Device Architecture
The metasurface operates with Right-Circularly Polarized (RCP) incident light, evaluating the Left-Circularly Polarized (LCP) light on the transmitted side.  

* **Spatial Multiplexing:** The meta-atoms are arranged in a diagonal segmented quadrant layout, where two diagonally opposite areas are mapped for red light and the remaining two are mapped for green.
* **Stack Integration:** An active liquid crystal layer is integrated directly below the metasurface layer to modulate the optical response.
* **Unit Cell Periodicity:** The unit cells were rigorously evaluated at periodicities of both **P = 400 nm** and **P = 350 nm**.  

## Holographic Reconstruction
The phase masks map target images with a CGH pixel size of 20. To achieve a high-quality floating holographic effect, the white background is computationally removed from the target projection images before applying the phase retrieval algorithms.  

* **Red Reconstruction:** Projects a floating flower hologram at **λ = 633 nm**, achieving a Root Mean Square Error (RMSE) of 0.1018.  
* **Green Reconstruction:** Projects a leaf structure at **λ = 532 nm**, achieving an RMSE of 0.0731.  

## Simulation Workflow & Usage

The simulation pipeline is executed in three primary sequential steps:

### 1. Diffraction Efficiency (DE) & PCR Optimization
Run the PCR simulation script to determine the ideal meta-atom dimensions for minimizing cross-talk between colours.
* The script creates a unit cell meta-atom and sweeps its length and width.
* It outputs heatmaps of DE and PCR for both 633 nm and 532 nm wavelengths.
* **Selection Criteria:** The optimal length and width for the red meta-atom are selected at the exact point where the DE for red light is maximized while the DE for green light is minimized. The exact inverse logic is used to select the dimensions for the green meta-atom. 

### 2. Geometric Phase Library Generation
Once the precise length and width dimensions are established, generate the lookup tables:

# Multicolour Holographic Metasurfaces

This repository contains the design, phase library, and Lumerical FDTD simulation pipeline for a multicolour hologram utilizing geometric phase metasurfaces. The workflow applies the Gerchberg-Saxton (GS) algorithm for Computer-Generated Holography (CGH) to extract target phase masks for individual colours.  

## Device Architecture
The metasurface operates with Right-Circularly Polarized (RCP) incident light, evaluating the Left-Circularly Polarized (LCP) light on the transmitted side[cite: 2].  

* **Spatial Multiplexing:** The meta-atoms are arranged in a diagonal segmented quadrant arrangement, where two diagonally opposite areas are mapped for red light and two for green light.
* **Stack Integration:** An active liquid crystal layer is integrated directly below the metasurface layer to modulate the optical response.
* **Unit Cell Periodicity:** The unit cells were rigorously evaluated at periodicities of both **P = 400 nm** and **P = 350 nm**[cite: 2].  

## Holographic Reconstruction
The phase masks map target images with a CGH pixel size of 20[cite: 2]. To achieve a high-quality floating holographic effect, the white background is computationally removed from the target projection images before applying the phase retrieval algorithms.  

* **Red Reconstruction:** Projects a floating flower hologram at **λ = 633 nm**, achieving a Root Mean Square Error (RMSE) of 0.1018[cite: 2].  
* **Green Reconstruction:** Projects a leaf structure at **λ = 532 nm**, achieving an RMSE of 0.0731[cite: 2].  

## Simulation Workflow & Usage

The simulation pipeline is executed in three primary sequential steps:

### 1. Diffraction Efficiency (DE) & PCR Optimization
Run the PCR simulation script to determine the ideal meta-atom dimensions for minimizing cross-talk between colours.
* The script creates a unit cell meta-atom and sweeps its length and width.
* It outputs heatmaps of DE and PCR for both 633 nm and 532 nm wavelengths.
* **Selection Criteria:** The optimal length and width for the red meta-atom are selected at the exact point where the DE for red light is maximized while the DE for green light is minimized. The exact inverse logic is used to select the dimensions for the green meta-atom. 

### 2. Geometric Phase Library Generation
Once the precise length and width dimensions are established, generate the lookup tables:

```bash
python Geometric_phase_lib.py
```bash
python Geometric_phase_lib.py
