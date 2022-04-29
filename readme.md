# PIC/Flip Fluid Simulator


# Introduction

This program simulate the velocity and pressure field of an incomplesible fluid in 2D using the PIC/Flip or the PIC algorithm. The boundaries conditions are homogeneous, i.e., v=0 in the boundaries.

The simulation is rendered in OpenGL and visualize the value of a parameter of the fluid, ink color is the standard.

You have to press 's' to start the simulation and you can toggle the visualization of the velocity field and/or the particles position pressing 'g'.



# Options

- You can change the simulation algorithm turning `false` the parameter `flipEnabled` in [Fluid2.h](src/Fluid2.h). If it is false you are using PIC and it is true you are using PIC/Flip.
- You can modify the parameter of the fluid like viscosity and density and the value of gravity or time step in [Scene.cpp](src/Scene.cpp). If you increase the time step or the viscosity, the algorithm could become unstable.
- You can change the number of cell in each axis in [Scene.cpp](src/Scene.cpp). If you inclease these number, the computation time will increase too and if you decrease these number, the algorithm could become unstable.

# Questions

Write to the mail: r.lozanoc93@gmail.com.

# Miscellaneous

Freeglut libraries are included in external but you can check [the freeglut project](http://freeglut.sourceforge.net/) for more information and credits.