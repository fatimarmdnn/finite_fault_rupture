# Hayward Fault Rupture Simulation using AxiSEM3D

This repository contains scripts and instructions for performing a finite rupture simulation on the Hayward Fault using AxiSEM3D.
Steps to Prepare and Run the Simulation

## 1- Prepare Simulation Folder:
 
   Before running the simulation, prepare the necessary input files by executing the following command:
 ```
        python prepare_input.py --eventfilename input_files/hayward_rupture.pkl
```

## 2- Run the Simulation: (e.g., on Archer2)

 ```
        sbatch submit.slurm
```
## Notes

- If you're using a different rupture model, update the `--eventfilename` argument in the `prepare_input.py` command accordingly.
  You can generate a `.pkl` file similar to the one above using your own rupture parameters.
- Check AxiSEM3D's documentation for output interpretation and post-processing tools.




