Program for Simulating the Phase Frustrated Kuramoto Oscillator
=================================================================

Background Information
----------------------
The Kuramoto Oscillator models synchonisaiton in complex networks. When a parameter 
called the phase frustration is introduced, it forces adjacent oscillators out of 
synchronisation. For small values of this parameter, a phenomena called 'cluster 
synchonisation' is observed, whereby clusters belonging to the same symmetry group 
become synchronised.

Program Functionality 
---------------------
The Phase Frustrated Kuramoto model is implemented using numerical integration
(Runge-Kutta algo). Both single and double layer networks can be implemented, 
with the latter case allowing  for different inter- and intra-layer coupling.

Two examples have been implemtnted in the 'Examples' folder.

Visualisation of Graph Symmetries
----------------------------------
You might wonder what I mean in 'Background Information' by 'belonging to the 
same symmetry group'. My friend Bricker has written some code which calculates 
these group decompositions for arbitrary graphs, and so I refer you to his 
GitHub repo: 
[ostlerb/Graph-Geometric-Decomposition](https://github.com/ostlerb/Graph-Geometric-Decomposition).


Dependencies
-------------
1. Matplotlib
2. Sympy
3. Numpy
