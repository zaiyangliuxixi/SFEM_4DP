# SFEM-4DP

This directory contains essential files and folders for conducting finite element analysis (FEA) for 4D printing. The project supports both two-dimensional and three-dimensional simulations. 

## Lib Folders
- **+three_dim_fea/**  
  Contains MATLAB code and functions specifically developed for conducting three-dimensional FEA simulations for 4D printing applications.

- **+two_dim_fea/**  
  Holds MATLAB scripts and classes designed for two-dimensional FEA simulations.

## MATLAB Files
- **FEA_4D_Printing.mlapp**  
  The main MATLAB App Designer file that provides a graphical user interface (GUI) for managing simulations, selecting material properties, and running FEA analyses.

- **meshObject.m**  
  Handles the meshing process for objects in preparation for finite element analysis.

- **applyConstrain.m**  
  Script for applying boundary conditions and constraints.

- **setSubregion.m**  
  Divides the meshed object into subregions, allowing the assignment of different material properties (e.g., Young’s modulus, Poisson’s ratio, density, measured strain) to each subregion.

- **runSimulation.m**  
  Executes the simulation based on the model's geometry, boundary conditions, and assigned material properties.

- **generateVideo.m**  
  A script to create video visualizations of the simulation results over time.


## Executable File
- **SFEM_4DP_installer.exe**  
  If you are not involved in the underlying development, you can simply download and use this executable file to run the application. This is a standalone installer for the SFEM-4DP (Strain-based Finite Element Model for 4D Printing) project.

---

### How to Use
1. If you're developing or modifying the code, open **FEA_4D_Printing.mlapp** in MATLAB R2024a to start the graphical interface.
2. For end-users, download and run **SFEM_4DP_installer.exe** to access all the functionalities.

---

This project was developed using MATLAB version R2024a.
