# SIO_MonteCarlo

A small project done in the Simulation and Optimization course at HEIG-VD. 
The idea is to use three algorithms in order to compute statistically the integral of a function over a given interval:
* Uniform sampling
* Importance sampling
* Uniform sampling with the control variable method

Three methods can be used in order to create the samplings (that contain the estimated area, 
its confidence interval at 95% and more):

* Specify the size of the sample,
* Continue sampling until a given width of the CI is reached,
* Continue sampling until a given elapsed time limit is reached.

The results are printed in the console and can be exported as CSV if the option is enabled (EXPORT_CSV).
