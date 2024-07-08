Hayward Fault Rupture Simulation using AxiSEM3D

This repository contains scripts and instructions for performing a finite rupture simulation on the Hayward Fault using AxiSEM3D.
Steps to Prepare and Run the Simulation

1- Prepare Simulation Folder:
 
   Before running the simulation, prepare the necessary input files by executing the following command:

        bash

        python prepare_input.py --eventfilename input_files/hayward_rupture.pkl

        This command will generate the required input files for the Hayward Fault rupture simulation.

    Run the Simulation:
        Ensure that AxiSEM3D is installed and available in the current working directory.
        Run the simulation using your preferred method or script.

Additional Notes

    AxiSEM3D Installation: Make sure that AxiSEM3D is correctly installed and configured in your environment before running the simulation.
    Input Files: The input file hayward_rupture.pkl should be located in the input_files directory. Adjust paths accordingly if your setup differs.

