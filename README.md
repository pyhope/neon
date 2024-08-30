# Analyze which sites neon occupies in MgO

## Possible sites: oxygen sites, cation sites, and interstitial sites
### Package requirements: numpy, matplotlib, MDAnalysis

**Important note:** There is a bug in the MDAnalysis module where it cannot correctly handle periodic boundary conditions when the left box boundary (xlo, ylo, zlo) has a non-zero value (e.g., the first values of the BOX BOUNDS are 1.0897046722 for the .dump files in the directory). I converted the file to the LAMMPS data format using atomsk, which adjusted xlo, ylo, and zlo to 0 (and xhi, yhi, zhi, along with the coordinates of all atoms, were also modified accordingly), allowing the code to run correctly. Therefore, when running this code, please ensure that the left box boundary is set to 0.

### Useage:
* **Step 1:** Smooth the trajectory using OVITO (remove thermal vibrational displacements of atoms)
* **Step 2:** Extract one frame from the smoothed trajectory
* **Step 3:** Change the filename in `neighbor_analysis.py` and run `python3 neighbor_analysis.py`
