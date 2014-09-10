# The Generation of Oscillations in Networks of Electrically Coupled Cells

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

Copyright 2014 Travis Johnson

## Matlab Code

Herein is contained the Matlab scripts and functions that I wrote to produce
the simulations. 

- parameters.m -- defines the variables and functions used in the models 

- *Model*.m -- these files define the systems of equations for each of the models
using the variables defined in the parameters file. 

- Run.m and RunCellNetwork.m -- scripts with the commands to simulate the models 
  and produce the graphics used in Report.pdf

- StackedPlot.m -- utility function used to produce the plots for the network models

## License

The MIT License (MIT)

Copyright (c) 2014 Travis Johnson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
