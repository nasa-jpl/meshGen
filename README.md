## meshGen
- Python mesh generator to create mesh files for AERO FSI (Fluid-Structure Interaction) suite running Embedded framework
- Primary input is STL file(s) representing the surface mesh of desired bodies
- Also includes the ability to run user-defined modules to create parametric meshes
  - meshGen_parachuteDGB.py is an example where gmsh is used to create parametrically a DGB parachute mesh
    - Note this has examples of writing more output files for AERO-S, including frames, decomposition, etc.
  - Users can create additional modules (and update the calling in the runMeshGen.py script) for other problems of interest
- Note that some functions (like converting .top files to .exo format for viewing in Paraview) use a shell script to run AERO-S or XP2EXO utilities in the AERO suite. The assumption is that the user has built these tools, and that sourcing the config.sh file will enable them to be used via the environment variables as noted in the functions calling them. If the user has not built AERO and/or adjusted the config.sh file appropriately for their setup, then it is likely these will not run (but should be caught by the try/except block to fail gracefully). 

## Usage
- Use "pip install -r requirements.txt" to make sure all required packages are installed
- To run meshGen, update settings in runMeshGen.py script as appropriate, and run "python runMeshGen.py"
- Output will be written to "output" folder (and warning that previous content here will be deleted)

## References
- AERO code suite (Stanford University)
  - https://bitbucket.org/frg/workspace/projects/AERO
  - https://frg.bitbucket.io/
- Select publications using AERO, FSI, supersonic parachute inflation
  - https://doi.org/10.2514/6.2026-3822
  - https://doi.org/10.2514/1.J064791
  - https://doi.org/10.2514/6.2022-2746
  - https://doi.org/10.2514/6.2022-0351
