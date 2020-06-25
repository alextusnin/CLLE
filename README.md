# CLLE
Numerical solver of LLE equation for modeling micro-resonator frequency combs

The first version of the solver is based on the conventional Split-step Fourier method. User should initialize all the simulation parameters in the main function. 
The step adaptative method uses Numerical recipes library. It is way faster than Split-step because it does not use expensive evaluation of exponents. 
The sript "read_data.py" provides with a small routine on python to facilitate data analysis in python
