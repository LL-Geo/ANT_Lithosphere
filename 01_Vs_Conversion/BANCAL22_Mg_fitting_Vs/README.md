This directory is the Bayesian invesrion procedure required to calibrate experimental parameterisations of anelasticity. The version provided here is as-used for the manuscript this repository pertains to. For the full, up-to-date, version of the code which is periodically updated, see BANCAL22 here: https://github.com/JamesHazzard/BANCAL22.

The structure of this directory differs from the other directories because this piece of code was written as separate product. 

algorithm: 
- contains GASWAM.py algorithm procedure
- allows flexibility to code in your own algorithm procedure

anelasticity_parameterisation:
- contains YT16 anelasticity parameterisation (Yamauchi and Takei, 2016, JGR)
- allows flexibility to use any chosen parameterisation

data:
- contains plate, adiabat, attenuation, viscosity, and potential temperature information
- import this data from the "inversiondata" directory

options:
- choose your algorithm by adding the name of the algorithm module you wish to use e.g. "GASWAM" to the file algorithm_selection.txt
- choose your anelasticity parameterisation by adding the name of the parameterisation you wish to use e.g. "YT16" to the file parameterisation_selection.txt
- choose your data constraints to use in the inversion procedure by adding a "1" beneath the name of each data type you wish to use, and a "0" beneath the name of each data type you wish not to use  in the file data_selection.txt

output:
- stores output of inversion runs including full history of samples
- archives the data used in the given inversion run

how to use:
- import chosen data into relevant data directories from "inversiondata" directory
- use options directory to choose your inversion set up
- once options/data imported, run "python3 setup.py" to set up the code
- after setting up, use "python3 main.py" to begin the inversion procedure

References:
See Richards, F.D. et al. (2020, JGR), DOI: https://doi.org/10.1029/2019JB019062, for a deterministic version of this inversion procedure