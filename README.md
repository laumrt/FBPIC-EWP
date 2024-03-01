# FBPIC-Electric Work Profiler
A simple code to estimate the work performed by the Laser and Plasma field on (FBPIC) Fourier-Bessel Particle-In-Cell code tracked electrons.

## Overview
FBPIC employs a set of 2D cylindrical grids to represent the fields. Each grid corresponds to a distinct azimuthal mode, denoted with an integer m. More specifically, m=0 represents all fields that are independent of the azimuthal angle $\theta$,
while m>1 refers to the fields that vary proportionally to $\cos(\theta)$ and $\sin(\theta)$. Ideally, m=1 represents the laser.
