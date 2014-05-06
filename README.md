hydro-gravity
=============


code workflow
-----
* setup: add spectral-tools to matlab's path
* run hydro simulation
  * design initial conditions
  * run stage1 - a coarse-grain hydro simulation on a small grid to get an idea how the solution looks like.
  * run stage2 - a high-res simulation on the entire time interval, writes a high-res output.
  * run stage3 - specify time points you plan to evaluate metric at, evaluate high-res hydro at these and neighbouring points for use with metric evaluation. 
* run gravity simulation

references
----------
* A\. Adams, N. Benjamin, and W. Musial - "Dynamical spacetimes from Numerical Hydrodynamics"
* W\. Musial - "The Trial of the Holographic Principle"