# Genetic algorithm Assigment

This program is thought as application to try to set the parametres of a model that predicts the evoluion of the COVID-19 pandemia. 


## Files

The folder sent contains several files used in the execution of the assigment: `GA_covid.c` program that consists on the main function to execute and that contains the principal genetic algorithm functions; `Fit_documentation.c` program that contains the functions needed to set the data for each individual and that calcuates the fitness function; the scripts `RKF78` and `RKF78`  that are responsable for solving the Runge-Kutta system in order to get the ODE results needed. and the `Report` explaining the results accuired.  
`solution` executables are also included for if needed.

## Usage

`GA_covid.c` must be first compiled. gcc compiler is recomended with -Ofast or -O3 optimizer for both files. The form followed can be:

    gcc -lm -O3 GA_covid.c -o solution

Then, `solution`  executable can be executed to obtain the results.

## Authors
@author: Adrià Colás Gallardo.        NIU: 1419956
@author: Xavier Olivé Fernández.    NIU: 1394407

