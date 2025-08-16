This code is based on Xlife++'s infrastructure.
It has a heavy dependency on it. 

The code is scinded into two main parts : 
-the resolution of wave equation with Xlife++
-the display of information with Matlab

Solving the wave equation
the /Sim/ directory and simulation files contains the algorithm and they handle the input and output of the solution

the /wave_opt/ directory solves one problem you input on the "data.txt" file and output a /trail/ directory that contains all the info you could want to restart the simulation and see the solution. 

the /batt_wave/ directory solves multiple problems. you input key parameters in the "data.txt" file and you can access 9 other parameters in the "data_batterie.txt" file. The code then outputs multiple /trail/ directories

Displaying the solution
in the /matlab/ directory by adding the /trail/ to the root, you can read the solution (see the (x,t) graph, the error on the (x,t) graph, the error and energy rate over time along with multiples frames of the simulation and/or a video of the simulation). 
