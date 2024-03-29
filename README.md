[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10776689.svg)](https://doi.org/10.5281/zenodo.10776689)

# FBPIC-Electric Work Profiler (FBPIC-EWP)
FBPIC-EWP, short for FBPIC-Electric Work Profiler, is a tool designed to estimate the work performed by the Laser and Plasma fields on tracked electrons in FBPIC ([Fourier-Bessel Particle-In-Cell](https://github.com/fbpic/fbpic)) simulations.

## Overview and Motivation
FBPIC-EWP operates within the context of FBPIC simulations, which employs a set of 2D cylindrical grids to represent fields with distinct azimuthal modes. While m=0 corresponds to fields independent of the azimuthal angle $\theta$, m>1 represents fields varying proportionally to $\cos(\theta)$ and $\sin(\theta)$. In an idealized scenario, m=1 represents the laser field.

The motivation behind FBPIC-EWP stems from the need to understand the contribution of both the Laser and Plasma fields to regimes like Direct Laser Acceleration. By estimating the work performed by these fields, researchers can gain insights into the acceleration mechanisms of tracked electrons.

## Key Features
FBPIC-EWP workflow is the following:
* __Field Retrieval Process__: FBPIC-EWP retrieves the total field ($\mathbf{E}_T$) on all selected particles, along with their coordinates $(x, y, z)$ and normalized momenta ($u$) for every FBPIC iteration. The total electric field is defined as the sum of the Laser and Plasma fields (i.e., $\mathbf{E}_L$ and $\mathbf{E}_W$ respectively)
* __Laser and Plasma field Calculation__: Disposing of the 3D electric field map, a 2D linear interpolation is performed for m=0 in order to calculate the Plasma electric field on each electron. Subsequently, the Laser electric field is retrieved as $\mathbf{E}_L = \mathbf{E}_T - \mathbf{E}_W$
* __Work Calculation__: Lastly, the work of each electric field component is calculated as $W_{W,L} = -e \int_{0}^{t} \mathbf{E}_{W,L} \cdot \mathbf{v}\,dt'$ in the time interval $[0, t]$. The results are stored in an HDF5 file.

FBPIC-EWP also allows to calculate the Laser and Plasma magnetic fields on each particle, through the `magnetic_field` parameter in the `main.py` module. `magnetic_field` is set to `False` by default.

## Requirements
* Python 3.x
* NumPy
* SciPy
* numba
* h5py
* [openPMD-viewer](https://github.com/openPMD/openPMD-viewer)
* tqdm

## Files
* `func.py`: Python module containing all the functions necessary to perform the work calculation (`interp_2D`, `particles_loop`, `spatial_lowpass_filter` and `shift_half_step`)
* `main.py`: Python module that allows to select the simulation and the electrons to study. It performs the actual work calcuation using the `func.py` module

## Installation
You can install FBPIC-EWP cloning the source
```python
git clone https://github.com/laumrt/FBPIC-EWP.git
cd FBPIC-EWP
pip install .
```
or directly with PiPy
```
pip install git+https://github.com/laumrt/FBPIC-EWP.git
```

## Contributions
Your contributions and suggestions for improvements are welcome.
If you encounter any issues or have ideas for enhancements, please feel free to open an issue or submit a pull request. 
