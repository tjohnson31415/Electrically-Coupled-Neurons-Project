This is the result of the final project for a course on the Biophsyics of
Neurons. My goal was to reproduce and validate the results of the following
paper:

Loewenstein, Y., Yarom, Y., & Sompolinsky, H. (2001). The generation of
oscillations in networks of electrically coupled cells. Proceedings of the
National Academy of Sciences, 98(14), 8095-8100.
http://www.pnas.org/content/98/14/8095.full

The authors present a number of dynamic models to describe the spiking
oscillations of neurons in electrically coupled networks. Such cells are
coupled by gap junctions that allow ions in the cytoplasm to equilibrate
between cells. It is logical to think that such coupling would simply damp the
spiking behavior, but simulations of the presented models reveal synchronous
oscillations in networks of two or more cells with very interesting dynamics.

The models are presented as systems of coupled differential equations.
Simulations were carried out in MATLAB using the appropriate ode* solver (ode45
was used for the simple models and ode15s for the stiffer problems of coupled
networks). For a summary of the models and results refer to either
Presentation.pdf or the Report.pdf which describe my simulations. For detailed
background on the models and their relevance to studying neurons refer to the
original paper.

Matlab Code
=========== 
Herein is contained the Matlab scripts and functions that I wrote to produce
the simulations. 

parameters.m -- defines the variables and functions used in the models 

*Model*.m -- these files define the systems of equations for each of the models
using the variables defined in the parameters file. 

Run.m and RunCellNetwork.m -- scripts with the commands to simulate the models
and produce the graphics used in Report.pdf

StackedPlot.m -- utility function used to produce the plots for the network models

License-ish
===========
Feel free to take and use the code without limitations.
