# FBPIC-Electric Work Profiler (FBPIC-EWP)
A simple code to estimate the work performed by the Laser and Plasma fields on FBPIC (Fourier-Bessel Particle-In-Cell) tracked electrons.

## Overview and Motivation
FBPIC employs a set of 2D cylindrical grids to represent the fields. Each grid corresponds to a distinct azimuthal mode, denoted with an integer m. More specifically, m=0 represents all fields that are independent of the azimuthal angle $\theta$,
while m>1 refers to the fields that vary proportionally to $\cos(\theta)$ and $\sin(\theta)$. Ideally, m=1 represents the laser.

When studying regimes such as Direct Laser Acceleration, it can be relevant to estimate the contribution of both the Laser and the Plasma fields on the accelerated electrons. This code simply exploits FBPIC mode decomposition to calculate the work of both fields. This is done as follows:
* FBPIC-EWP retrieves the total field ($E_T$) on all selected particles, along with their coordinates $(x, y, z)$ and normalized momenta ($u$) for every FBPIC iteration. The total electric field is defined as the sum of the Laser and Plasma fields (i.e., $E_L$ and $E_W$ respectively)
* Disposing of the 3D electric field map, a 2D linear interpolation is performed for m=0 in order to calculate the Plasma electric field on each electron. Subsequently, the Laser electric field is retrieved as $E_L = E_T - E_W$

