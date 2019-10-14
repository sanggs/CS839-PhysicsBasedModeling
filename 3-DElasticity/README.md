# Path to the code
https://github.com/sanggs/CS839-PhysicsBasedModeling/tree/master/3-DElasticity/tests/assignment2

# Implementation
The simulation is made on a 3D "X". The ends of the arms of this 3D "X" are used as handles to induce a stretch and the body reacts to this based on Linear Elasticity.
If the macro in main.cpp and AnimatedMesh.cpp is changed to TWO_DIMENSION, the code compiles Linear Elasticity on a 2D memberane stretched in a similar fashion. 

# Platform
This code has been tested on OSX.

# How to run
Run the following commands from within this directory

cd ../../
cmake -G Xcode

1. Open the Xcode project created in this directory in Xcode.
2. Select the assignment2 build option and run. 
3. This generates the 'assignment2.usda' file in the CS839-PhysicsBasedModeling/3-DElasticity/tests/assignment2/Debug folder.

# Visualise the usda file
Use the usdview API provided by PIXAR USD to visualise this file.

