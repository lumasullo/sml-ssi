# sml-ssi
Repository with code to study Single-molecule Localization with Sequential Structured Illumination. It contains code to numerically calculate CRB of different configurations and simulate (Monte Carlo) experiments and measurements.

- How to reproduce the figures contained in Masullo et al (2021):

1) Clone the repository
2) Install the necessary packages (matplotlib, numpy, scipy, configparser)
3) Run the  [method]_plots.py   

where [method] = ot, otmin, minflux, rastmin, rastmax, camera

This should produce the different figures in the manuscript.

- How to explore new configurations and produce new simulated data or reproduce the data in Masullo et al (2021):

1) Go the the folder sml-ssi/crb/[method] and open the desired script. They work as follows:
   - crb.py produces a 2D CRB map
   - sigma_vs_fov.py calculates the average CRB as a function of the FOV size
   - sigma_vs_N.py calculates the average CRB as a function of N (detected photons)
   - sigma_vs_sbr.py calculates the average CRB as a function of SBR (signal-to-background ratio)
   
2) The scripts will generate preview plots and also save the results in the folder "results" together with the metadata of the simulation parameters.
You will have to move these results to the desired folder. For example I moved them into the sml_ssi/figure_[method] folders.

3) If you want to explore other parameters than the ones analyzed in Masullo et al (2021) you have to change the simulation parameters in the scripts.
All parameters are defined at the beginning of the script and saved into a config file.

4) These scripts are just short files that define arrays of variables and call and use the functions from sml-ssi/tools/tools_simulations.py
That is where you can find all the functions that are actually used to perform all the calculations. 
These functions are implementations of the theory described in Masullo et al (2021).

5) The folder sml-ssi/simulations contains scripts to simulate single-molecule localization experiments for each method.

NOTE: The first lines of the scripts automatically fix the working directory, which should be [path]/sml-ssi where [path] is the absolute path to the folder where you cloned the repo. This should work fine for any Python environment, but has only been checked in Anaconda + Spyder.

If you need help, or have comments or suggestions, please write me to: lu.masullo@gmail.com
