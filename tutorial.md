# Code Explanation (main.py)
## 1. Imports and Initial Setup
First, we import the necessary libraries and define the output file name (`file_out`) and simulation path (`hdf5_path`). Subsequently, the code imports the simulation (`sim = LpaDiagnostics(...)`) and it is possible to set `magnetic_field = True`, in case to calculate the laser and plasma magnetic field on each tracked particle. 

```python
import numpy as np
from scipy.constants import c, e, m_e, pi
import time
import h5py
from tqdm import tqdm
from openpmd_viewer import OpenPMDTimeSeries, ParticleTracker
from openpmd_viewer.addons import LpaDiagnostics
from func import particles_loop, spatial_lowpass_filter, shift_half_step

# Defining the output file name
file_out = 'work_result_test.h5'

# Selecting the simulation
hdf5_path = ".../hdf5"
sim = LpaDiagnostics(hdf5_path, backend='h5py', check_all_files=False)

# Set to True if you want to decompose the magnetic field in m=0,1 modes (False by default)
magnetic_field = False
```
## 2. Particle Selection
Select the species of particles to study that has to correspond to a species defined in the original FBPIC simulation (in this example `species='elec_L'`). Here, we import the ids of the tracked particles at the last iteration with `id_ = sim.get_particle(['id'], t = sim.t[-1], species=species)[0]`. We then randomly choose a subset of ten particles to study (i.e., `id_array = np.random.choice(id_, size=10, replace=False)`).
```python
# Selecting the particles to study
species = 'elec_L'
id_ = sim.get_particle(['id'], t = sim.t[-1], species=species)[0]
id_array = np.random.choice(id_, size=10, replace=False)
pt_select = ParticleTracker(sim, species=species, preserve_particle_index=True, select=id_array)
N_p = pt_select.N_selected  # number of selected electrons
```
## 3. Timestep Parameters
FBPIC-EWP automatically defines the timestep (`dt`) from the original FBPIC simulation timestep.
```python
# Define the time range
t = sim.t  # (s)
N_t = t.size
dt = t[0:2].ptp()  # (s), time step
```
## 4. Data Import
Import the coordinates (`x, y, z`), momenta (`ux, uy, uz`), and total electric (`Ex_t, Ey_t, Ez_t`) and magnetic (`Bx_t, By_t, Bz_t`) fields of each particle for every iteration.
```python
# Importing the coordinates, momenta and total EM-field of each particle for every iteration
if magnetic_field:
    print('(1/3) Importing the coordinates, momenta and total EM-field of each particle for every iteration\n')
    x, y, z, ux, uy, uz, Ex_t, Ey_t, Ez_t, Bx_t, By_t, Bz_t, w = sim.iterate(sim.get_particle,
       ['x', 'y', 'z', 'ux', 'uy', 'uz', 'E/x', 'E/y', 'E/z', 'B/x', 'B/y', 'B/z', 'w'],
       spec=species, select=pt_select)
    Bx_t, By_t, Bz_t = np.nan_to_num(Bx_t), np.nan_to_num(By_t), np.nan_to_num(Bz_t)
else:
    print('(1/3) Importing the coordinates, momenta and total E-field of each particle for every iteration\n')
    x, y, z, ux, uy, uz, Ex_t, Ey_t, Ez_t, w = sim.iterate(sim.get_particle,
       ['x', 'y', 'z', 'ux', 'uy', 'uz', 'E/x', 'E/y', 'E/z', 'w'],
       spec=species, select=pt_select)
```
## 5. Field Calculation Loop
For each time step, extract and filter the electric and magnetic field components. Use a JIT-compiled loop to interpolate the fields onto the particle positions and save the electric field components for each particle. Here, ('Ex_las, Ey_las, Ez_las') and ('Ex_wake, Ey_wake, Ez_wake') are arrays containing the laser and plasma electric field at each iteration for each selected particle. Analogously, ('Bx_las, By_las, Bz_las') and ('Bx_wake, By_wake, Bz_wake') are arrays containing the laser and plasma magnetic field at each iteration for each selected particle. 
```python
print('(2/3) Calculating the electric fields on each particle... \n')

# Electric field on particle calculation (iterate over each time step in the specified time range)
for it in tqdm(range(N_t)):
    
    ################################ ELECTRIC FIELD ##########################################
    
    # Extract electric field components in cylindrical coordinates (V/m)
    Er_w_map = sim.get_field(t=t[it], field='E', coord='r', m=0, theta=pi/2)[0]
    Ez_w_map, info = sim.get_field(t=t[it], field='E', coord='z', m=0, theta=pi/2)

    # Lowpass filter on the radial field
    Er_w_map = spatial_lowpass_filter(Er_w_map, info.dz)

    # Extract radial and axial coordinates for field interpolation
    r_int, dr_int = info.r, info.dr  # m
    z_int, dz_int = info.z, info.dz  # m
    
    # JIT particle loop
    Ex_l, Ey_l, Ez_l, Ex_w, Ey_w, Ez_w = particles_loop(
        x[it, :], y[it, :], z[it, :], Ex_t[it, :], Ey_t[it, :], Ez_t[it, :], 
        r_int.min(), z_int.min(), dr_int, dz_int,
        Er_w_map, Ez_w_map
    )

    # Saving all the electric fields (V/m)
    Ex_las[it, :] = Ex_l
    Ey_las[it, :] = Ey_l
    Ez_las[it, :] = Ez_l
    Ex_wake[it, :] = Ex_w
    Ey_wake[it, :] = Ey_w
    Ez_wake[it, :] = Ez_w
    
    ################################# MAGNETIC FIELD ##########################################
    
    if magnetic_field:
        # Extract magnetic field components in cylindrical coordinates (T)
        Br_w_map = sim.get_field(t=t[it], field='B', coord='r', m=0, theta=pi/2)[0]
        Bz_w_map, info = sim.get_field(t=t[it], field='B', coord='z', m=0, theta=pi/2)

        # Lowpass filter on the radial field
        Br_w_map = spatial_lowpass_filter(Br_w_map, info.dz)

        # Extract radial and axial coordinates for field interpolation
        r_int, dr_int = info.r, info.dr  # m
        z_int, dz_int = info.z, info.dz  # m

        # JIT particle loop
        Bx_l, By_l, Bz_l, Bx_w, By_w, Bz_w = particles_loop(
            x[it, :], y[it, :], z[it, :], Bx_t[it, :], By_t[it, :], Bz_t[it, :],
            r_int.min(), z_int.min(), dr_int, dz_int,
            Br_w_map, Bz_w_map
        )

        # Saving all the magnetic fields (T)
        Bx_las[it, :] = Bx_l
        By_las[it, :] = By_l
        Bz_las[it, :] = Bz_l
        Bx_wake[it, :] = Bx_w
        By_wake[it, :] = By_w
        Bz_wake[it, :] = Bz_w

```
